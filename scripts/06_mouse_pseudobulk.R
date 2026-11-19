#Baier et al.
#06: Mouse pseudobulk, age/cognition differential expression, pathway enrichment, 
# three-group patterns, 126-label reassignment, and mouse-level program scores

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "06_mouse_bulk")
ckdir <- file.path(out, "checkpoints")
dir.create(ckdir, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(Seurat)
library(Matrix)
library(data.table)
library(edgeR)
library(fgsea)
library(msigdbr)

set.seed(20261002)
NULL_EPS <- 1e-10
#Aggregate SoupX-corrected atlas counts by mouse and Allen subclass.
atlas_file <- hca_atlas_input
atlas <- readRDS(atlas_file)
eligibility <- data.table::fread(file.path(root, "05_differential_abundance",
                                           "05_baiermouse_subclasses_allgroups_statistical_eligibility.csv"))
primary_types <- eligibility[statistical_tier == "Primary", allen_subclass]
stopifnot(length(primary_types) == 11L, !anyDuplicated(primary_types))

counts <- SeuratObject::LayerData(atlas, assay = "RNA", layer = "counts")
cell_meta <- atlas[[]]
stopifnot(identical(colnames(counts), rownames(cell_meta)),
          all(c("mouse_id", "allen_subclass", "age_group", "primary_group",
                "sequence_batch") %in% names(cell_meta)))

sample_info <- unique(as.data.table(cell_meta)[, .(
  mouse_id, sequence_batch, age_group, primary_group
)])
setorder(sample_info, mouse_id)
stopifnot(nrow(sample_info) == 13L, !anyDuplicated(sample_info$mouse_id))
ids <- sample_info$mouse_id

make_pseudobulk <- function(L) {
  idx <- which(cell_meta$allen_subclass == L & cell_meta$mouse_id %in% ids)
  stopifnot(length(idx) > 0L)
  mouse <- factor(cell_meta$mouse_id[idx], levels = ids)
  design <- Matrix::sparse.model.matrix(~ 0 + mouse)
  colnames(design) <- ids
  pb <- counts[, idx, drop = FALSE] %*% design
  colnames(pb) <- ids
  rownames(pb) <- rownames(counts)
  pb
}

pb <- setNames(lapply(primary_types, make_pseudobulk), primary_types)
stopifnot(all(vapply(pb, function(x) identical(colnames(x), ids), logical(1))))

n_cells <- as.data.table(cell_meta)[allen_subclass %in% primary_types,
                                    .(n_cells = .N), by = .(allen_subclass, mouse_id)]
qc <- rbindlist(lapply(names(pb), function(L) {
  data.table(lineage = L, mouse_id = ids,
             total_UMIs = as.numeric(Matrix::colSums(pb[[L]])))
}))
qc <- merge(qc, n_cells, by.x = c("lineage", "mouse_id"),
            by.y = c("allen_subclass", "mouse_id"), all.x = TRUE)
qc[is.na(n_cells), n_cells := 0L]
qc[, UMIs_per_cell := total_UMIs / pmax(n_cells, 1L)]
fwrite(qc, file.path(out, "06_baiermouse_alllineages_pseudobulk_library_qc.csv"))
fwrite(sample_info, file.path(out, "06_baiermouse_allmice_pseudobulk_sample_metadata.csv"))
saveRDS(pb, file.path(out, "06_baiermouse_alllineages_pseudobulk_raw_counts.rds"), compress = FALSE)
rm(atlas, counts, cell_meta)
gc()

# Use mouse metadata carried by the de novo annotated atlas.
meta <- as.data.frame(sample_info)
meta$group <- factor(
  ifelse(grepl("young", tolower(meta$age_group)), "Young",
         ifelse(grepl("intact", tolower(meta$primary_group)), "Intact", "Impaired")),
  levels = c("Young", "Intact", "Impaired")
)
meta$batch <- factor(meta$sequence_batch)
stopifnot(sum(meta$group == "Young") == 4L,
          sum(meta$group == "Intact") == 4L,
          sum(meta$group == "Impaired") == 5L, !anyNA(meta$batch))
print(table(meta$group, meta$batch))

#Retain Hallmark and GOBP memberships from the original msigdbr specification.
hd <- msigdbr::msigdbr(species = "Mus musculus", collection = "H")
hs <- split(hd$gene_symbol, hd$gs_name)
stopifnot(length(hs) == 50L)
saveRDS(hs, file.path(out, "06_baiermouse_alllineages_hallmark_gene_sets.rds"))
gd <- msigdbr::msigdbr(species = "Mus musculus", collection = "C5",
                       subcollection = "GO:BP")
gs <- split(gd$gene_symbol, gd$gs_name)
saveRDS(gs, file.path(out, "06_baiermouse_alllineages_gobp_gene_sets.rds"))

#Primary all-13-mouse edgeR models
as_num <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }
fit_one <- function(x, group, batch) {
  design <- model.matrix(~ 0 + group + batch)
  colnames(design) <- sub("^group", "", colnames(design))
  stopifnot(qr(design)$rank == ncol(design))
  y <- edgeR::DGEList(x)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  list(fit = edgeR::glmQLFit(y, design, robust = TRUE), design = design)
}

CON <- list(
  Cog = c(Intact = 1, Impaired = -1),
  AgeP = c(Young = 1, Intact = -1),
  YvImp = c(Young = 1, Impaired = -1),
  AgePool = c(Young = 1, Intact = -4/9, Impaired = -5/9)
)

de <- function(f, weights) {
  v <- setNames(rep(0, ncol(f$design)), colnames(f$design))
  stopifnot(all(names(weights) %in% names(v)))
  v[names(weights)] <- weights
  tab <- edgeR::topTags(edgeR::glmQLFTest(f$fit, contrast = v),
                        n = Inf, sort.by = "none")$table
  data.table(gene = rownames(tab), logFC = tab$logFC,
             logCPM = tab$logCPM, F = tab$F, PValue = tab$PValue,
             FDR = tab$FDR,
             s = sign(tab$logFC) * sqrt(pmax(tab$F, 0)))
}

# Independent, stable random streams for each enrichment family. This avoids
# changes in Hallmark FDR when an unrelated analysis or permutation is skipped.
fgsea_seed <- function(lineage, contrast, analysis = "observed", split = 0L) {
  i <- match(lineage, primary_types)
  k <- match(contrast, c(names(CON), "GOBP"))
  stopifnot(!is.na(i), !is.na(k))
  section <- switch(analysis, observed = 0L, reassignment = 100000L,
                    two_group_GO = 300000L)
  as.integer(20261002L + section + i * 1000L + k + 10L * as.integer(split))
}

gsea <- function(x, paths = hs, min_size = 10, eps = 0, seed) {
  ranks <- sort(setNames(x$s, x$gene)[is.finite(x$s)], decreasing = TRUE)
  set.seed(seed)
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    paths, ranks, minSize = min_size, maxSize = 500, eps = eps,
    BPPARAM = BiocParallel::SerialParam()))
  as.data.table(fg)[, .(pathway, NES, pval, padj)]
}

# Prespecified fixed lineage/contrast order for reproducible observed FGSEA.
# Reset the random seed once before the entire sequence of comparisons.
OBSERVED_FGSEA_LINEAGE_ORDER <- c(
  "Microglia NN",
  "Oligo NN",
  "Endo NN",
  "CHOR NN",
  "OPC NN",
  "Astro-TE NN",
  "Lymphoid NN",
  "Ependymal NN",
  "DG Glut",
  "Peri NN",
  "BAM NN"
)
stopifnot(
  setequal(OBSERVED_FGSEA_LINEAGE_ORDER, names(pb)),
  identical(names(CON), c("Cog", "AgeP", "YvImp", "AgePool"))
)

# First fit all observed edgeR models; these do not alter the FGSEA RNG stream.
obs_file <- file.path(out, "06_baiermouse_alllineages_fourcontrasts_edger_checkpoint.rds")
obs <- lapply(names(pb), function(L) {
  message("Primary mouse models: ", L)
  f <- fit_one(as_num(pb[[L]]), meta$group, meta$batch)
  list(DE = lapply(CON, function(w) de(f, w)))
})
names(obs) <- names(pb)

# Start the deterministic observed FGSEA RNG stream before the first call.
# Do not reseed between comparisons or alter the invocation parameters.
set.seed(20261002)
for (L in OBSERVED_FGSEA_LINEAGE_ORDER) {
  message("Observed Hallmark FGSEA: ", L)
  hall <- rbindlist(lapply(names(CON), function(k) {
    x <- obs[[L]]$DE[[k]]
    ranks <- sort(setNames(x$s, x$gene)[is.finite(x$s)],
                  decreasing = TRUE)
    fg <- as.data.table(suppressWarnings(
      fgsea::fgseaMultilevel(
        pathways = hs, stats = ranks,
        minSize = 10, maxSize = 500, eps = 0
      )
    ))
    fg[, .(pathway, NES, pval, padj)][, contrast := k]
  }))
  obs[[L]]$H <- hall[, lineage := L]
}

saveRDS(obs, obs_file)

H <- rbindlist(lapply(obs, `[[`, "H"))
fwrite(H, file.path(out, "06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv"))
fwrite(rbindlist(lapply(names(obs), function(L) {
  rbindlist(lapply(names(obs[[L]]$DE), function(k) {
    copy(obs[[L]]$DE[[k]])[, `:=`(lineage = L, contrast = k)]
  }))
})), file.path(out, "06_baiermouse_alllineages_fourcontrasts_gene_effects.csv"))

#Classify Hallmarks by the complete young/intact/impaired expression pattern
W <- dcast(H, lineage + pathway ~ contrast, value.var = c("NES", "padj"))
W[, pattern := fcase(
  padj_Cog < .05 & padj_AgeP < .05 & sign(NES_Cog) == sign(NES_AgeP), "Progressive",
  padj_Cog < .05 & padj_AgeP < .05, "Intact-specific",
  padj_Cog < .05 & padj_YvImp < .05 & sign(NES_YvImp) == sign(NES_Cog), "Impairment-specific",
  padj_Cog < .05, "Cognition-only (other)",
  padj_AgeP < .05 & padj_YvImp < .05 & sign(NES_AgeP) == sign(NES_YvImp), "Aging-shared",
  default = "ns"
)]
fwrite(W, file.path(out, "06_baiermouse_alllineages_threegroups_pathway_patterns.csv"))

#126-assignment permutation null
aged <- ids[meta$group != "Young"]
combos <- combn(aged, 4)
stopifnot(ncol(combos) == 126L)
lins <- names(pb)
stopifnot(all(lins %in% names(pb)))
checkpoint_file <- function(L) file.path(
  ckdir, paste0("06_baiermouse_", tolower(gsub("[^A-Za-z0-9]+", "_", L)), "_cognition_label_reassignment_checkpoint.rds"))
rho <- function(x, y) suppressWarnings(
  cor(x, y, method = "spearman", use = "complete.obs"))

for (L in lins) {
  ck <- checkpoint_file(L)
  message("126 assignments: ", L)
  cnt <- as_num(pb[[L]])
  res <- lapply(seq_len(ncol(combos)), function(i) {
    g <- as.character(meta$group)
    g[g != "Young"] <- "Impaired"
    g[ids %in% combos[, i]] <- "Intact"
    f <- fit_one(cnt, factor(g, levels = c("Young", "Intact", "Impaired")),
                 meta$batch)
    dc <- de(f, CON$Cog)
    da <- de(f, CON$AgeP)
    hc <- gsea(dc, eps = NULL_EPS,
               seed = fgsea_seed(L, "Cog", "reassignment", i))
    ha <- gsea(da, eps = NULL_EPS,
               seed = fgsea_seed(L, "AgeP", "reassignment", i))
    m <- merge(hc[, .(pathway, c = NES)], ha[, .(pathway, a = NES)],
               by = "pathway")
    list(split = i,
         is_true = setequal(combos[, i], ids[meta$group == "Intact"]),
         gene_rho = rho(dc$logFC, da$logFC),
         hall_rho = rho(m$c, m$a),
         cogNES = hc[, .(pathway, NES)])
  })
  saveRDS(res, ck)
}

#Shared-reference correlation null (125 alternative labels for each lineage)
summ <- rbindlist(lapply(lins, function(L) {
  r <- readRDS(checkpoint_file(L))
  is_true <- vapply(r, `[[`, logical(1), "is_true")
  stopifnot(sum(is_true) == 1L)
  t <- r[[which(is_true)]]
  n <- r[!is_true]
  g <- vapply(n, `[[`, numeric(1), "gene_rho")
  h <- vapply(n, `[[`, numeric(1), "hall_rho")
  data.table(lineage = L,
             obs_gene_rho = t$gene_rho,
             null_gene_mean = mean(g),
             null_gene_q025 = as.numeric(quantile(g, .025)),
             null_gene_q975 = as.numeric(quantile(g, .975)),
             obs_hall_rho = t$hall_rho,
             null_hall_mean = mean(h),
             null_hall_q025 = as.numeric(quantile(h, .025)),
             null_hall_q975 = as.numeric(quantile(h, .975)),
             hall_rho_minus_null = t$hall_rho - mean(h),
             p_hall_upper = (1 + sum(h >= t$hall_rho)) / (1 + length(h)),
             p_hall_lower = (1 + sum(h <= t$hall_rho)) / (1 + length(h)))
}))
fwrite(summ, file.path(out, "06_baiermouse_alllineages_cognition_label_reassignment_null_summary.csv"))

#Nominal two-sided empirical Hallmark P values
perm <- rbindlist(lapply(lins, function(L) {
  r <- readRDS(checkpoint_file(L))
  z <- rbindlist(lapply(r, function(x) {
    copy(x$cogNES)[, `:=`(split = x$split, is_true = x$is_true)]
  }))
  z[, .(obs_NES = NES[is_true],
        perm_p = (1 + sum(abs(NES[!is_true]) >= abs(NES[is_true]))) /
          (1 + sum(!is_true))), by = pathway][, lineage := L]
}))
fwrite(perm, file.path(out, "06_baiermouse_alllineages_cognition_label_reassignment_pathway_p.csv"))

#Mouse-level Hallmark scores for the five-lineage trajectories
LT <- c("Astro-TE NN", "Microglia NN", "OPC NN", "Oligo NN", "Peri NN")
PROG <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_E2F_TARGETS",
          "HALLMARK_DNA_REPAIR")
stopifnot(all(LT %in% names(pb)), all(PROG %in% names(hs)))
lcpm <- function(L) {
  m <- as_num(pb[[L]])
  edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(m)), log = TRUE)
}
zs <- function(m, genes) {
  genes <- intersect(genes, rownames(m))
  colMeans(t(scale(t(m[genes, , drop = FALSE]))))
}
scores <- rbindlist(lapply(LT, function(L) {
  m <- lcpm(L)
  rbindlist(lapply(PROG, function(p) {
    data.table(lineage = L, program = p, mouse_id = ids,
               group = meta$group, score = zs(m, hs[[p]]))
  }))
}))
scores[, z := as.numeric(scale(score)), by = .(lineage, program)]
score_means <- scores[, .(mean_z = mean(z), se = sd(z) / sqrt(.N), n = .N),
                      by = .(lineage, program, group)]
score_means <- merge(score_means,
                     W[, .(lineage, program = pathway, pattern)],
                     by = c("lineage", "program"), all.x = TRUE)
stopifnot(nrow(scores) == 5L * 8L * 13L,
          all(scores[, uniqueN(mouse_id), by = .(lineage, program)]$V1 == 13L))
fwrite(scores, file.path(out, "06_baiermouse_alllineages_threegroups_program_scores_by_mouse.csv"))
fwrite(score_means, file.path(out, "06_baiermouse_alllineages_threegroups_program_group_means.csv"))

#Two group cognition GO analysis
two_group_de <- rbindlist(lapply(names(pb), function(L) {
  keep_mouse <- meta$group != "Young"
  cnt <- as_num(pb[[L]][, keep_mouse, drop = FALSE])
  info <- meta[keep_mouse, , drop = FALSE]
  group <- factor(tolower(as.character(info$group)),
                  levels = c("intact", "impaired"))
  batch <- factor(info$sequence_batch)
  design <- model.matrix(~ batch + group)
  stopifnot(qr(design)$rank == ncol(design))
  y <- edgeR::DGEList(counts = cnt)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  tab <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = "groupimpaired"),
                        n = Inf, sort.by = "PValue")$table
  data.table(lineage = L, gene = rownames(tab), logFC = tab$logFC,
             logCPM = tab$logCPM, F = tab$F, PValue = tab$PValue,
             FDR = tab$FDR,
             s = sign(tab$logFC) * sqrt(pmax(tab$F, 0)))
}))
fwrite(two_group_de, file.path(out, "06_baiermouse_alllineages_cognition_twogroup_gene_effects.csv"))
two_group_go <- rbindlist(lapply(names(pb), function(L) {
  message("Original two-group GO: ", L)
  gsea(two_group_de[lineage == L], paths = gs, min_size = 15,
       eps = 0, seed = fgsea_seed(L, "GOBP", "two_group_GO"))[, lineage := L]
}))
fwrite(two_group_go, file.path(out, "06_baiermouse_alllineages_cognition_twogroup_gobp_enrichment.csv"))

#Continuous cognitive-flexibility sensitivity (astrocytes; aged mice only)
# Secondary association; the mouse is the replicate (n=4 Intact, n=5 Impaired).
# The standardized CF coefficient is positive for higher cf51to61, as observed,
# and does not duplicate the primary categorical intact-vs-impaired contrast.
# No resampling or FGSEA is performed; existing outputs above are unchanged.
flex_behavior_file <- file.path(root, "01_mouse_behavior",
                                "01_baiermouse_allmice_behavior_mouse_level.csv")
if (!file.exists(flex_behavior_file))
  stop("Missing upstream mouse-level behavioral data: ", flex_behavior_file)
flex_behavior <- data.table::fread(flex_behavior_file)
stopifnot(all(c("mouse_id", "cf51to61") %in% names(flex_behavior)),
          !anyDuplicated(flex_behavior$mouse_id))
flex_aged <- which(as.character(meta$group) %in% c("Intact", "Impaired"))
flex_info <- meta[flex_aged, , drop = FALSE]
flex_info$cf51to61 <- as.numeric(flex_behavior$cf51to61[
  match(flex_info$mouse_id, flex_behavior$mouse_id)])
stopifnot(nrow(flex_info) == 9L,
          sum(as.character(flex_info$group) == "Intact") == 4L,
          sum(as.character(flex_info$group) == "Impaired") == 5L,
          all(is.finite(flex_info$cf51to61)),
          stats::sd(flex_info$cf51to61) > 0)
flex_info$CF_z <- as.numeric(scale(flex_info$cf51to61))
flex_info$batch <- droplevels(factor(flex_info$sequence_batch))
flex_design <- stats::model.matrix(~ batch + CF_z, data = flex_info)
if (qr(flex_design)$rank != ncol(flex_design))
  stop("Continuous cognitive-flexibility design is not full rank")
flex_y <- edgeR::DGEList(counts = as_num(
  pb[["Astro-TE NN"]][, flex_aged, drop = FALSE]))
flex_keep <- edgeR::filterByExpr(flex_y, design = flex_design)
if (sum(flex_keep) < 100L)
  stop("Too few expressed astrocyte genes for continuous CF analysis")
flex_y <- edgeR::calcNormFactors(
  flex_y[flex_keep, , keep.lib.sizes = FALSE], method = "TMM")
flex_y <- edgeR::estimateDisp(flex_y, flex_design, robust = TRUE)
flex_fit <- edgeR::glmQLFit(flex_y, flex_design, robust = TRUE)
flex_tab <- edgeR::topTags(
  edgeR::glmQLFTest(flex_fit, coef = "CF_z"),
  n = Inf, sort.by = "PValue")$table
flex_de <- data.table(
  gene = rownames(flex_tab), logFC = flex_tab$logFC,
  logCPM = flex_tab$logCPM, F = flex_tab$F,
  PValue = flex_tab$PValue, FDR = flex_tab$FDR)
data.table::fwrite(
  flex_de, file.path(out,
    "06_baiermouse_astrocyte_cognitive_flexibility_continuous_gene_effects.csv"))
data.table::fwrite(
  data.table(mouse_id = flex_info$mouse_id,
             group = as.character(flex_info$group),
             sequence_batch = as.character(flex_info$sequence_batch),
             cf51to61 = flex_info$cf51to61, CF_z = flex_info$CF_z),
  file.path(out,
    "06_baiermouse_astrocyte_cognitive_flexibility_continuous_model_design.csv"))
message("Continuous astrocyte CF model: ", nrow(flex_de),
        " genes tested across nine aged mice (batch adjusted)")
