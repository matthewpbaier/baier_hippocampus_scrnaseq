#Baier et al.
#07: Mouse pseudobulk sensitivity analyses

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "07_mouse_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261003)

#Load libraries
library(data.table)
library(edgeR)
library(fgsea)
library(limma)

#Load outputs from scripts 01 and 06
pb <- readRDS(file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_pseudobulk_raw_counts.rds"))
obs <- readRDS(file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_fourcontrasts_edger_checkpoint.rds"))
hs <- readRDS(file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_hallmark_gene_sets.rds"))
meta <- fread(file.path(root, "06_mouse_bulk", "06_baiermouse_allmice_pseudobulk_sample_metadata.csv"))
behavior <- fread(file.path(root, "01_mouse_behavior", "01_baiermouse_allmice_behavior_mouse_level.csv"))

ids <- colnames(pb[[1]])
stopifnot(length(hs) == 50L, length(ids) == 13L,
          all(vapply(pb, function(x) identical(colnames(x), ids), logical(1))),
          all(names(pb) %in% names(obs)))
meta <- meta[match(ids, mouse_id)]
stopifnot(identical(as.character(meta$mouse_id), ids),
          all(c("mouse_id", "sequence_batch", "age_group", "primary_group") %in% names(meta)))
meta[, group := factor(ifelse(grepl("young", tolower(age_group)), "Young",
                       ifelse(grepl("intact", tolower(primary_group)), "Intact", "Impaired")),
                       levels = c("Young", "Intact", "Impaired"))]
meta[, batch := factor(sequence_batch)]
stopifnot(identical(as.integer(table(meta$group)), c(4L, 4L, 5L)))

behavior <- behavior[match(ids, mouse_id)]
stopifnot(identical(as.character(behavior$mouse_id), ids),
          all(behavior$sequenced),
          identical(as.character(behavior$group), as.character(meta$group)),
          all(is.finite(behavior$weight_loss_pct)))
meta[, weight_loss_abs_pct := abs(behavior$weight_loss_pct)]

AST <- "Astro-TE NN"
PER <- "Peri NN"
END <- "Endo NN"
OX <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
FA <- "HALLMARK_FATTY_ACID_METABOLISM"
stopifnot(all(c(AST, PER, END) %in% names(pb)),
          all(c(OX, FA) %in% names(hs)))


#edgeR, FGSEA and mouse-score helpers
#Positive denotes intact-favored (cognition) or younger age (age)
as_num <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

fit_mouse <- function(lineage, keep = rep(TRUE, length(ids)), covariate = NULL,
                      remove_genes = character()) {
  stopifnot(lineage %in% names(pb), is.logical(keep), length(keep) == length(ids),
            !anyNA(keep))
  cnt <- as_num(pb[[lineage]][, keep, drop = FALSE])
  if (length(remove_genes)) {
    cnt <- cnt[!rownames(cnt) %in% remove_genes, , drop = FALSE]
  }
  group <- droplevels(meta$group[keep])
  batch <- droplevels(meta$batch[keep])
  stopifnot(all(c("Intact", "Impaired") %in% levels(group)))
  design <- if (nlevels(batch) > 1L) model.matrix(~ 0 + group + batch) else model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  if (!is.null(covariate)) {
    stopifnot(length(covariate) == length(ids), all(is.finite(covariate[keep])),
              sd(covariate[keep]) > 0)
    design <- cbind(design, covariate_z = as.numeric(scale(covariate[keep])))
  }
  if (qr(design)$rank != ncol(design)) {
    stop("Model is rank-deficient: ", lineage, "; mice=", paste(ids[keep], collapse = ","))
  }
  y <- edgeR::DGEList(cnt)
  expressed <- edgeR::filterByExpr(y, design = design)
  if (!any(expressed)) stop("No expressed genes after filtering: ", lineage)
  y <- edgeR::calcNormFactors(y[expressed, , keep.lib.sizes = FALSE], method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  list(fit = edgeR::glmQLFit(y, design, robust = TRUE), design = design)
}

contrasts <- list(
  Cog = c(Intact = 1, Impaired = -1),
  AgeP = c(Young = 1, Intact = -1)
)

get_de <- function(fitted, contrast) {
  w <- contrasts[[contrast]]
  if (is.null(w) || !all(names(w) %in% colnames(fitted$design))) {
    stop("Unavailable contrast: ", contrast)
  }
  v <- setNames(numeric(ncol(fitted$design)), colnames(fitted$design))
  v[names(w)] <- w
  tab <- edgeR::topTags(edgeR::glmQLFTest(fitted$fit, contrast = v),
                        n = Inf, sort.by = "none")$table
  data.table(gene = rownames(tab), logFC = tab$logFC, logCPM = tab$logCPM,
             F = tab$F, PValue = tab$PValue, FDR = tab$FDR,
             s = sign(tab$logFC) * sqrt(pmax(tab$F, 0)))
}

run_hallmarks <- function(de_table, gene_sets = hs, min_size = 10L) {
  z <- de_table[is.finite(s) & !is.na(gene)]
  ranks <- sort(setNames(z$s, z$gene), decreasing = TRUE)
  fg <- as.data.table(suppressWarnings(fgsea::fgseaMultilevel(
    pathways = gene_sets, stats = ranks, minSize = min_size, maxSize = 500,
    eps = 0, BPPARAM = BiocParallel::SerialParam()
  )))
  fg[, .(pathway, NES, pval, padj, size,
         leading_edge = vapply(leadingEdge, paste, character(1), collapse = ";"))]
}

canonical_primary <- function(lineage, axis) {
  x <- as.data.table(obs[[lineage]]$H)
  ans <- x[contrast == axis, .(pathway, NES, pval, padj)]
  stopifnot(nrow(ans) == 50L, !anyDuplicated(ans$pathway))
  ans
}

use_canonical_primary <- function(estimated, lineage, axis) {
  baseline <- canonical_primary(lineage, axis)
  idx <- match(estimated$pathway, baseline$pathway)
  stopifnot(nrow(estimated) == 50L, !anyNA(idx))
  estimated[, `:=`(
    NES = baseline$NES[idx],
    pval = baseline$pval[idx],
    padj = baseline$padj[idx]
  )]
  estimated
}

log_cpm <- function(lineage) {
  x <- edgeR::DGEList(as_num(pb[[lineage]]))
  edgeR::cpm(edgeR::calcNormFactors(x, method = "TMM"), log = TRUE)
}

zscore_program <- function(logcpm, genes) {
  genes <- intersect(genes, rownames(logcpm))
  if (!length(genes)) stop("No signature genes present in pseudobulk matrix")
  z <- t(scale(t(logcpm[genes, , drop = FALSE])))
  z <- z[apply(z, 1L, function(x) all(is.finite(x))), , drop = FALSE]
  if (!nrow(z)) stop("No finite gene-standardized signature genes")
  colMeans(z)
}

#Behavioral weight sensitivity
aged <- meta$group != "Young"
stopifnot(sum(aged) == 9L)
nonventricular <- setdiff(names(pb), c("CHOR NN", "Ependymal NN"))
stopifnot(length(nonventricular) == 9L)
weight_pathways <- c(OX, FA, "HALLMARK_PEROXISOME", "HALLMARK_GLYCOLYSIS")
stopifnot(all(weight_pathways %in% names(hs)))

weight_results <- rbindlist(lapply(nonventricular, function(L) {
  message("Weight adjustment: ", L)
  f0 <- fit_mouse(L, keep = aged)
  f1 <- fit_mouse(L, keep = aged, covariate = meta$weight_loss_abs_pct)
  rbindlist(list(
    run_hallmarks(get_de(f0, "Cog"))[, model := "Batch adjusted"],
    run_hallmarks(get_de(f1, "Cog"))[, model := "Batch + weight adjusted"]
  ))[, lineage := L]
}))
fwrite(weight_results, file.path(out, "07_baiermouse_alllineages_weightadjusted_weight_hallmarks_all50.csv"))
fwrite(weight_results[pathway %in% weight_pathways,
                      .(lineage, model, pathway, NES, pval, padj, size)],
       file.path(out, "07_baiermouse_alllineages_weightadjusted_weight_sensitivity.csv"))

#Independent aged-only limma-voom/CAMERA sensitivity (inter-gene correlation 0.01)
run_camera <- function(lineage) {
  cnt <- as_num(pb[[lineage]][, aged, drop = FALSE])
  info <- copy(meta[aged])
  info[, camera_group := factor(group, levels = c("Impaired", "Intact"))]
  X <- model.matrix(~ batch + camera_group, data = info)
  stopifnot(qr(X)$rank == ncol(X), "camera_groupIntact" %in% colnames(X))
  y <- edgeR::DGEList(cnt)
  keep <- edgeR::filterByExpr(y, design = X)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  v <- limma::voom(y, design = X, plot = FALSE)
  index <- limma::ids2indices(hs, rownames(v))
  size <- vapply(index, length, integer(1))
  index <- index[size >= 10L & size <= 500L]
  contrast <- setNames(numeric(ncol(X)), colnames(X))
  contrast["camera_groupIntact"] <- 1
  cam <- as.data.table(limma::camera(
    v, index = index, design = X, contrast = contrast,
    inter.gene.cor = 0.01, sort = FALSE
  ), keep.rownames = "pathway")
  cam[, .(lineage = lineage, pathway, camera_direction = Direction,
          camera_P = PValue, camera_FDR = p.adjust(PValue, "BH"),
          n_genes = NGenes, camera_sign = fifelse(Direction == "Up", 1L, -1L))]
}

camera_all <- rbindlist(lapply(nonventricular, function(L) {
  message("CAMERA: ", L)
  run_camera(L)
}))
fwrite(camera_all, file.path(out, "07_baiermouse_alllineages_cognition_camera_hallmarks_aged.csv"))

fgsea_aged <- weight_results[model == "Batch adjusted",
                              .(lineage, pathway, aged_NES = NES,
                                aged_FGSEA_FDR = padj)]
comparison <- merge(fgsea_aged, camera_all, by = c("lineage", "pathway"))
comparison[, support := fcase(
  aged_FGSEA_FDR < .05 & camera_FDR < .05 & sign(aged_NES) == camera_sign,
  "FGSEA + CAMERA",
  aged_FGSEA_FDR < .05, "FGSEA only",
  camera_FDR < .05, "CAMERA only",
  default = "Neither"
)]
primary_hallmarks <- rbindlist(lapply(names(obs), function(L) {
  as.data.table(obs[[L]]$H)[contrast == "Cog",
    .(lineage = L, pathway, primary_13m_NES = NES, primary_13m_FDR = padj)]
}))
comparison <- merge(comparison, primary_hallmarks,
                    by = c("lineage", "pathway"), all.x = TRUE)
fwrite(comparison, file.path(out, "07_baiermouse_alllineages_cognition_fgsea_camera_comparison.csv"))

#Dissociation-response and immediate-early gene sensitivity
IEG <- c("Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3",
         "Atf3", "Ier2", "Ier3", "Ier5", "Nr4a1", "Dusp1", "Zfp36",
         "Klf2", "Klf4", "Klf6", "Hspa1a", "Hspa1b", "Hspa8", "Hsph1",
         "Dnajb1", "Socs3", "Btg2", "Cebpd", "Ppp1r15a", "Gadd45b", "Gadd45g")
stopifnot(length(IEG) == 29L, !anyDuplicated(IEG))

#Preserve original standardized mean across genes, rather than an AddModuleScore
diss_scores <- rbindlist(lapply(names(pb), function(L) {
  lc <- log_cpm(L)
  data.table(lineage = L, mouse_id = colnames(lc),
             diss_score = zscore_program(lc, IEG),
             n_genes = length(intersect(IEG, rownames(lc))))
}))
fwrite(diss_scores, file.path(out, "07_baiermouse_alllineages_dissociation_dissociation_scores_by_mouse_lineage.csv"))
hs_no_IEG <- lapply(hs, setdiff, y = IEG)

handling <- rbindlist(lapply(nonventricular, function(L) {
  message("Dissociation / IEG sensitivity: ", L)
  original <- obs[[L]]$DE$Cog
  d <- diss_scores[lineage == L][match(ids, mouse_id), diss_score]
  stopifnot(length(d) == length(ids), all(is.finite(d)))
  adjusted <- get_de(fit_mouse(L, covariate = d), "Cog")
  # Use the exact upstream primary FGSEA results instead of re-estimating them.
  primary <- as.data.table(obs[[L]]$H)[contrast == "Cog",
      .(pathway, NES, pval, padj)]
  primary[, `:=`(size = NA_integer_, leading_edge = NA_character_,
                 model = "Primary 13 mice")]
  rbindlist(list(
    primary,
    run_hallmarks(original, hs_no_IEG)[, model := "IEG removed"],
    run_hallmarks(adjusted)[, model := "Dissociation adjusted"]
  ))[, lineage := L]
}))
fwrite(handling, file.path(out, "07_baiermouse_alllineages_dissociation_dissociation_ieg_hallmarks_all50.csv"))
immune_pathways <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                     "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA",
                     "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS")
fwrite(handling[lineage == END & pathway %in% immune_pathways,
                .(lineage, pathway, model, NES, pval, padj)],
       file.path(out, "07_baiermouse_endothelial_cognition_endothelial_immune_sensitivity.csv"))


#Endothelial mouse-level program scores and immune-marker adjustment
end_lc <- log_cpm(END)
endo_programs <- immune_pathways[1:4]
endo_scores <- rbindlist(lapply(endo_programs, function(path) {
  data.table(mouse_id = ids, pathway = path,
             score = zscore_program(end_lc, hs[[path]]))
}))
immune_markers <- c("Ptprc", "Lst1", "Tyrobp", "Cd74", "C1qa")
marker_genes <- intersect(immune_markers, rownames(end_lc))
stopifnot(length(marker_genes) >= 3L)
marker_score <- colMeans(end_lc[marker_genes, , drop = FALSE])
endo_scores <- merge(endo_scores,
                     meta[, .(mouse_id, group, batch)], by = "mouse_id")
endo_scores[, immune_marker_logCPM := marker_score[match(mouse_id, ids)]]
fwrite(endo_scores, file.path(out, "07_baiermouse_endothelial_cognition_endothelial_program_scores_by_mouse.csv"))

fit_score <- function(x, with_marker) {
  x <- as.data.frame(x)
  x$group <- factor(as.character(x$group), levels = c("Intact", "Impaired"))
  x$batch <- droplevels(factor(x$batch))
  predictors <- c("batch", "group",
                  if (with_marker) "immune_marker_logCPM")
  form <- stats::reformulate(predictors, response = "score")
  fit <- stats::lm(form, data = x)
  stopifnot(qr(model.matrix(fit))$rank == ncol(model.matrix(fit)))
  coef <- summary(fit)$coefficients["groupImpaired", ]
  ci <- confint(fit, "groupImpaired")
  data.table(estimate_impaired_minus_intact = unname(coef["Estimate"]),
             estimate_intact_minus_impaired = -unname(coef["Estimate"]),
             ci_low_intact_minus_impaired = -unname(ci[2]),
             ci_high_intact_minus_impaired = -unname(ci[1]),
             p = unname(coef["Pr(>|t|)"]), n = nrow(x))
}

endo_models <- rbindlist(lapply(c(FALSE, TRUE), function(add_marker) {
  z <- endo_scores[group != "Young",
                   fit_score(.SD, add_marker), by = pathway]
  z[, model := if (add_marker) "Batch + immune marker" else "Batch adjusted"]
  z[, FDR_four := p.adjust(p, "BH")]
  z
}))
fwrite(endo_models, file.path(out, "07_baiermouse_endothelial_cognition_endothelial_score_models.csv"))

#Pericyte astrocyte transcript sensitivity
astro_marker_genes <- c("Aqp4", "Slc1a2", "Slc1a3", "Gja1", "Aldoc", "Gpc5",
                        "Atp1a2", "Sparcl1", "Ntsr2", "Mlc1", "Gjb6", "Slc4a4")
pericyte_marker_genes <- c("Pdgfrb", "Rgs5", "Kcnj8", "Abcc9", "Vtn",
                          "Higd1b", "Ndufa4l2", "Atp13a5")
astro_lc <- log_cpm(AST)
peri_lc <- log_cpm(PER)
astro_score_peri <- zscore_program(peri_lc, astro_marker_genes)
peri_score_astro <- zscore_program(astro_lc, pericyte_marker_genes)
marker_table <- data.table(mouse_id = ids, group = meta$group,
                           astro_marker_in_pericytes = astro_score_peri,
                           peri_marker_in_astrocytes = peri_score_astro)
fwrite(marker_table, file.path(out, "07_baiermouse_pericyte_cognition_pericyte_astrocyte_marker_scores.csv"))
common <- intersect(rownames(astro_lc), rownames(peri_lc))
astro_enriched <- common[rowMeans(astro_lc[common, , drop = FALSE]) -
                         rowMeans(peri_lc[common, , drop = FALSE]) > 2]

peri_de <- list(
  Primary = obs[[PER]]$DE$Cog,
  Astro_marker_adjusted = get_de(fit_mouse(PER, covariate = astro_score_peri), "Cog"),
  Astro_enriched_genes_removed = get_de(fit_mouse(PER, remove_genes = astro_enriched), "Cog")
)
peri_gsea <- rbindlist(lapply(names(peri_de), function(model) {
  result <- run_hallmarks(peri_de[[model]])
  if (model == "Primary")
    result <- use_canonical_primary(result, PER, "Cog")
  result[, model := model]
  result
}))
peri_gsea[, lineage := PER]
fwrite(peri_gsea, file.path(out, "07_baiermouse_pericyte_cognition_pericyte_astro_transcript_all50.csv"))
fwrite(peri_gsea[pathway %in% c(OX, FA),
                  .(lineage, model, pathway, NES, pval, padj, size)],
       file.path(out, "07_baiermouse_pericyte_cognition_pericyte_astro_transcript_sensitivity.csv"))

astro_gsea <- run_hallmarks(obs[[AST]]$DE$Cog)
overlap <- rbindlist(lapply(c(OX, FA), function(p) {
  pe <- strsplit(peri_gsea[model == "Primary" & pathway == p, leading_edge], ";", fixed = TRUE)[[1]]
  ae <- strsplit(astro_gsea[pathway == p, leading_edge], ";", fixed = TRUE)[[1]]
  data.table(pathway = p, n_pericyte_le = length(pe), n_astro_le = length(ae),
             n_shared = length(intersect(pe, ae)),
             jaccard = length(intersect(pe, ae)) / length(union(pe, ae)),
             frac_pericyte_le_astro_enriched = mean(pe %in% astro_enriched))
}))
fwrite(overlap, file.path(out, "07_baiermouse_pericyte_cognition_pericyte_astro_leading_edge_overlap.csv"))

#Astrocyte OXPHOS leave-one-mouse-out
ox_loo <- rbindlist(lapply(c("none", ids), function(drop_id) {
  message("Astrocyte OXPHOS omission: ", drop_id)
  keep <- if (drop_id == "none") rep(TRUE, length(ids)) else ids != drop_id
  fitted <- fit_mouse(AST, keep = keep)
  rbindlist(lapply(c("Cog", "AgeP"), function(axis) {
    res <- run_hallmarks(get_de(fitted, axis))
    if (drop_id == "none")
      res <- use_canonical_primary(res, AST, axis)
    setorder(res, -NES)
    res[, rank_positive_NES := seq_len(.N)]
    z <- res[pathway == OX]
    stopifnot(nrow(z) == 1L)
    z[, .(dropped_mouse = drop_id,
          dropped_group = if (drop_id == "none") "none" else as.character(meta$group[match(drop_id, ids)]),
          contrast = axis, OX_NES = NES, OX_padj = padj,
          OX_rank = rank_positive_NES)]
  }))
}))
fwrite(ox_loo, file.path(out, "07_baiermouse_astrocyte_cognition_astrocyte_oxphos_leave_one_mouse_out.csv"))
stopifnot(nrow(ox_loo) == (length(ids) + 1L) * 2L,
          all(ox_loo$OX_rank >= 1L & ox_loo$OX_rank <= 50L))

#Astrocyte mitochondrial-gene and neuronal ambient-RNA sensitivity
astro_primary_de <- as.data.table(obs[[AST]]$DE$Cog)
astro_primary_h <- as.data.table(obs[[AST]]$H)[contrast == "Cog",
    .(pathway, NES, pval, padj)]
astro_no_mt_h <- run_hallmarks(
  astro_primary_de[!grepl("^mt-", gene)])[, .(pathway, NES, pval, padj)]
astro_mt_compare <- rbindlist(list(
  copy(astro_primary_h)[, model := "Primary"],
  copy(astro_no_mt_h)[, model := "mt-encoded genes excluded"]
))
fwrite(astro_mt_compare, file.path(out, "07_baiermouse_astrocyte_cognition_astrocyte_mt_gene_exclusion_all50.csv"))

#Strongest DG-enriched genes with DG minus astrocyte mean logCPM >4
dg <- grep("^DG Glut", names(pb), value = TRUE)
stopifnot(length(dg) == 1L)
dg_lc <- log_cpm(dg)
shared <- intersect(rownames(dg_lc), rownames(astro_lc))
dg_mean <- rowMeans(dg_lc[shared, , drop = FALSE])
astro_mean <- rowMeans(astro_lc[shared, , drop = FALSE])
ambient_genes <- head(names(sort(dg_mean[dg_mean - astro_mean > 4],
                                 decreasing = TRUE)), 50L)

n_ambient_tested <- sum(ambient_genes %in% obs[[AST]]$DE$Cog$gene)
stopifnot(length(ambient_genes) == 50L, n_ambient_tested == 49L)
message("Neuronal ambient signature: ", length(ambient_genes),
        " genes scored; ", n_ambient_tested,
        " represented in the astrocyte cognition DE test universe")
ambient_score <- colMeans(astro_lc[ambient_genes, , drop = FALSE])
ambient_results <- run_hallmarks(
  get_de(fit_mouse(AST, covariate = ambient_score), "Cog"))
astro_ambient <- merge(astro_primary_h,
                       ambient_results[, .(pathway, ambient_NES = NES,
                                           ambient_padj = padj)], by = "pathway")
fwrite(astro_ambient, file.path(out, "07_baiermouse_astrocyte_cognition_astrocyte_ambient_adjusted_all50.csv"))
fwrite(data.table(gene = ambient_genes),
       file.path(out, "07_baiermouse_astrocyte_cognition_neuronal_ambient_signature_genes.csv"))
fwrite(data.table(mouse_id = ids, group = meta$group,
                  astro_ambient_score = ambient_score),
       file.path(out, "07_baiermouse_astrocyte_cognition_astrocyte_ambient_score_by_mouse.csv"))

prop_file <- file.path(root, "05_differential_abundance",
                       "05_baiermouse_subclasses_allgroups_proportions_by_mouse.csv")
prop <- fread(prop_file)
proxy_labels <- c("DG Glut", "DG-PIR Ex IMN", "CA1-ProS Glut", "CA3 Glut",
                  "CHOR NN", "Ependymal NN")
proxy_mouse <- dcast(prop[allen_subclass %in% proxy_labels],
                     mouse_id ~ allen_subclass, value.var = "proportion", fill = 0)
proxies <- intersect(proxy_labels, names(proxy_mouse))
stopifnot(length(proxies) > 0L)
proxy_mouse <- proxy_mouse[match(ids, mouse_id)]
stopifnot(identical(as.character(proxy_mouse$mouse_id), ids))

spearman_test <- function(a, b) {
  good <- is.finite(a) & is.finite(b)
  if (sum(good) < 4L || uniqueN(a[good]) < 2L || uniqueN(b[good]) < 2L) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  t <- suppressWarnings(cor.test(a[good], b[good], method = "spearman",
                                 exact = FALSE))
  c(rho = unname(t$estimate), p = t$p.value)
}
dissection <- rbindlist(lapply(proxies, function(p) {
  x <- proxy_mouse[[p]]
  cf <- spearman_test(x[aged], behavior$cf51to61[aged])
  ox <- spearman_test(x[aged], zscore_program(astro_lc, hs[[OX]])[aged])
  data.table(proxy = p, mean_young = mean(x[meta$group == "Young"]),
             mean_intact = mean(x[meta$group == "Intact"]),
             mean_impaired = mean(x[meta$group == "Impaired"]),
             rho_CF = cf[["rho"]], p_CF = cf[["p"]],
             rho_astro_OXPHOS = ox[["rho"]], p_astro_OXPHOS = ox[["p"]])
}))
fwrite(dissection, file.path(out, "07_baiermouse_alllineages_cognition_dissection_proxy_correlations.csv"))
