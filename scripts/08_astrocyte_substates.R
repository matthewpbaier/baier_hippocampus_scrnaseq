#Baier et al.
#08: Astrocyte transcriptional states, abundance, metabolic expression,
# state-identity pathways, and within-state analyses 

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "08_astrocyte_states")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(Seurat)
library(Matrix)
library(data.table)
library(edgeR)
library(fgsea)
library(limma)
library(UCell)
library(msigdbr)
library(readxl)

set.seed(120926)
STATES <- paste0("cAST", 0:5)
MAJOR <- STATES[1:2]
MINOR <- STATES[3:6]
MIN_SCORE_CELLS <- 10L   
MIN_PB_CELLS <- 25L      
MIN_REF_CELLS <- 5L      
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, paste0("08_", name, ".csv")))
require_file <- function(p) if (!file.exists(p)) stop("Required input not found: ", p)

atlas_file <- hca_atlas_input
meta_file <- file.path(root, "06_mouse_bulk", "06_baiermouse_allmice_pseudobulk_sample_metadata.csv")
beh_file <- file.path(root, "01_mouse_behavior", "01_baiermouse_allmice_behavior_mouse_level.csv")
hall_file <- file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_hallmark_gene_sets.rds")
go_file <- file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_gobp_gene_sets.rds")
whole_de_file <- file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_cognition_twogroup_gene_effects.csv")
astro_file <- file.path(out, "08_frozen_clean_astrocytes.rds")
for (p in c(meta_file, beh_file, hall_file, go_file, whole_de_file)) require_file(p)

mouse <- fread(meta_file)
need <- c("mouse_id", "age_group", "primary_group", "sequence_batch")
stopifnot(all(need %in% names(mouse)), nrow(mouse) == 13L,
          !anyDuplicated(mouse$mouse_id))
mouse[, mouse_id := as.character(mouse_id)]
mouse[, sequence_batch := as.character(sequence_batch)]
mouse[, group := fifelse(tolower(age_group) == "young", "Young",
                        fifelse(tolower(primary_group) == "intact", "Intact", "Impaired"))]
mouse[, group := factor(group, levels = c("Young", "Intact", "Impaired"))]
stopifnot(mouse[, sum(group == "Young")] == 4L,
          mouse[, sum(group == "Intact")] == 4L,
          mouse[, sum(group == "Impaired")] == 5L)
setorder(mouse, mouse_id)
beh <- fread(beh_file)
stopifnot(all(c("mouse_id", "cf51to61") %in% names(beh)))
mouse[, CF := as.numeric(beh$cf51to61[match(mouse_id, beh$mouse_id)])]
stopifnot(all(is.finite(mouse[group != "Young", CF])))
write_table(mouse[, .(mouse_id, group, sequence_batch, CF)], "mouse_manifest")
hs <- readRDS(hall_file)
gs <- readRDS(go_file)
stopifnot(length(hs) == 50L,
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION" %in% names(hs), length(gs) > 100L)

#Astrocyte subclustering
  require_file(atlas_file)
  obj <- readRDS(atlas_file)
  stopifnot(inherits(obj, "Seurat"), "allen_subclass" %in% names(obj[[]]),
            "mouse_id" %in% names(obj[[]]))
  a <- subset(obj, cells = colnames(obj)[as.character(obj$allen_subclass) == "Astro-TE NN"])
  rm(obj); gc()
  DefaultAssay(a) <- "RNA"
  stopifnot(ncol(a) > 10000L)
  a$mouse_id <- as.character(a$mouse_id)
  a$batch <- mouse$sequence_batch[match(a$mouse_id, mouse$mouse_id)]
  stopifnot(!anyNA(a$batch), setequal(unique(a$mouse_id), mouse$mouse_id))

  cluster_astro <- function(x) {
    x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize",
                               scale.factor = 10000, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst",
                                     nfeatures = 3000, verbose = FALSE)
    hvg <- Seurat::VariableFeatures(x)
    feat <- hvg[!toupper(hvg) %in% toupper(hs[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]) &
                  !grepl("^mt-|^MT-", hvg) & !grepl("^Rpl|^Rps|^RPL|^RPS", hvg)]
    if (length(feat) < 1000L) stop("Insufficient non-OXPHOS clustering features")
    x <- Seurat::ScaleData(x, features = feat, verbose = FALSE)
    x <- Seurat::RunPCA(x, features = feat, npcs = 30, verbose = FALSE)
    x <- Seurat::FindNeighbors(x, reduction = "pca", dims = 1:20, verbose = FALSE)
    x <- Seurat::RunUMAP(x, reduction = "pca", dims = 1:20,
                        seed.use = 120926, verbose = FALSE)
    x <- Seurat::FindClusters(x, resolution = 0.1, algorithm = 1,
                              random.seed = 120926, verbose = FALSE)
    list(object = x, features = feat)
  }

  initial <- cluster_astro(a)
  a <- initial$object
  raw <- as.character(Seurat::Idents(a))
  ranked <- names(sort(table(raw), decreasing = TRUE))
  lookup <- setNames(paste0("AST", seq_along(ranked) - 1L), ranked)
  a$initial_state <- unname(lookup[raw])
  write_table(data.table(raw_cluster = ranked, initial_state = unname(lookup),
                         n_cells = as.integer(table(factor(raw, levels = ranked)))),
              "initial_state_crosswalk")
  stopifnot("AST2" %in% a$initial_state)

  #Mixed lineage filter, applied to log-normalized RNA
  norm <- SeuratObject::LayerData(a, assay = "RNA", layer = "data")
  endo <- intersect(c("Pecam1", "Cldn5", "Kdr", "Tek", "Egfl7",
                      "Adgrl4", "Slco1a4", "Emcn"), rownames(norm))
  if (length(endo) < 5L) stop("Endothelial-marker coverage is unexpectedly low")
  endo_n <- as.numeric(Matrix::colSums(norm[endo, , drop = FALSE] > 0))
  audit <- data.table(cell_id = colnames(a),
                      initial_state = as.character(a$initial_state), endo_n = endo_n)
  audit[, remove := initial_state == "AST2" | endo_n >= 2L]
  audit[, reason := fifelse(initial_state == "AST2" & endo_n >= 2L,
                           "initial_AST2_and_vascular",
                           fifelse(initial_state == "AST2", "initial_AST2",
                                   fifelse(endo_n >= 2L, "vascular", "retained")))]
  write_table(audit, "cleaning_cell_ledger")
  a <- subset(a, cells = audit[remove == FALSE, cell_id])
  final <- cluster_astro(a)
  a <- final$object
  raw_clean <- as.character(Seurat::Idents(a))
  a$clean_state <- factor(paste0("cAST", raw_clean), levels = STATES)
  write_table(data.table(gene = final$features), "final_clustering_features")
  saveRDS(Seurat::DietSeurat(a, assays = "RNA", dimreducs = c("pca", "umap"),
                            graphs = NULL, misc = TRUE), astro_file)
  rm(a); gc()
astro <- readRDS(astro_file)
stopifnot(ncol(astro) == 10155L,
          setequal(as.character(unique(astro$clean_state)), STATES),
          setequal(as.character(unique(astro$mouse_id)), mouse$mouse_id))
write_table(data.table(cell_id = colnames(astro),
                       mouse_id = as.character(astro$mouse_id),
                       state = as.character(astro$clean_state)), "frozen_cell_assignments")

#Mouse x state ledger
md <- as.data.table(astro[[]], keep.rownames = "cell_id")
md[, `:=`(mouse_id = as.character(mouse_id), state = as.character(clean_state))]
md <- md[, intersect(c("cell_id", "mouse_id", "state", "nCount_RNA",
                        "nFeature_RNA", "percent.mt"), names(md)), with = FALSE]
md <- merge(md, mouse[, .(mouse_id, group, sequence_batch, CF)], by = "mouse_id",
            all.x = TRUE, sort = FALSE)
md <- md[match(colnames(astro), cell_id)]
stopifnot(identical(md$cell_id, colnames(astro)), !anyNA(md$group))
ledger <- CJ(mouse_id = mouse$mouse_id, state = STATES)
ledger <- merge(ledger, md[, .(n = .N), by = .(mouse_id, state)],
                by = c("mouse_id", "state"), all.x = TRUE, sort = FALSE)
ledger[is.na(n), n := 0L]
ledger <- merge(ledger, mouse[, .(mouse_id, group, sequence_batch, CF)], by = "mouse_id")
ledger[, N_astro := sum(n), by = mouse_id]
ledger[, fraction := n / N_astro]
stopifnot(all(ledger$N_astro > 0),
          all(abs(ledger[, sum(fraction), by = mouse_id]$V1 - 1) < 1e-10))
write_table(ledger, "state_counts_fractions_by_mouse")
pooled <- ledger[state %in% MINOR, .(n = sum(n), N_astro = unique(N_astro),
                                    group = unique(group), sequence_batch = unique(sequence_batch),
                                    CF = unique(CF)), by = mouse_id]
pooled[, `:=`(state = "cAST2_to_5", fraction = n / N_astro)]
write_table(pooled, "pooled_minor_fractions_by_mouse")
write_table(rbindlist(list(ledger[, .(mouse_id, group, state, fraction)],
                           pooled[, .(mouse_id, group, state, fraction)]))[
                             , .(mean_fraction = mean(fraction), sem_fraction = sd(fraction)/sqrt(.N),
                                 n_mice = .N), by = .(group, state)],
            "state_abundance_group_summary")
write_table(ledger[, .(n_mice_ge10 = sum(n >= 10L), n_mice_with_cells = sum(n > 0L)),
                   by = .(state, group)], "state_coverage_by_group")

tests <- rbindlist(list(ledger[, .(mouse_id, state, n, N_astro, group,
                                   sequence_batch)],
                        pooled[, .(mouse_id, state, n, N_astro, group,
                                    sequence_batch)]))
tests <- tests[as.character(group) %in% c("Intact", "Impaired")]
stopifnot(nrow(tests) == 63L, !anyDuplicated(tests, by = c("mouse_id", "state")))
tests[, `:=`(group = factor(as.character(group), levels = c("Intact", "Impaired")),
             batch = factor(sequence_batch),
             y = asin(sqrt((n + .5) / (N_astro + 1))))]
fit_abundance <- function(d) {
  d <- as.data.frame(d)
  d$group <- factor(
    as.character(d$group),
    levels = c("Intact", "Impaired")
  )
  d$batch <- droplevels(factor(d$batch))
  
  # Determine whether sequencing batch is estimable
  batch_adjusted <- FALSE
  
  if (nlevels(d$batch) > 1L) {
    dm <- model.matrix(~ batch + group, data = d)
    batch_adjusted <- qr(dm)$rank == ncol(dm)
  }
  
  # Preserve the original fallback to a group-only model
  predictors <- if (batch_adjusted) {
    c("batch", "group")
  } else {
    "group"
  }
  
  form <- stats::reformulate(predictors, response = "y")
  
  fit <- summary(stats::lm(form, data = d))$coefficients
  
  data.table(
    estimate_imp_minus_intact =
      unname(fit["groupImpaired", "Estimate"]),
    PValue =
      unname(fit["groupImpaired", "Pr(>|t|)"]),
    batch_adjusted = batch_adjusted
  )
}

da <- tests[, fit_abundance(.SD), by = state]
da[, FDR_seven_tests := p.adjust(PValue, "BH")]
da <- merge(da, tests[, .(mean_fraction = mean(n/N_astro)), by = .(state, group)],
            by = "state", allow.cartesian = TRUE)
da <- dcast(da, state + estimate_imp_minus_intact + PValue + FDR_seven_tests + batch_adjusted ~ group,
            value.var = "mean_fraction")
da[, estimate_intact_minus_imp := -estimate_imp_minus_intact]
write_table(da, "state_abundance_seven_test_models")
loo <- rbindlist(lapply(unique(tests$mouse_id), function(drop_id)
  tests[mouse_id != drop_id, fit_abundance(.SD), by = state][, dropped_mouse := drop_id]))
write_table(loo, "state_abundance_leave_one_mouse_out")

#AST1 fraction versus cognitive flexibility
cf <- ledger[state == "cAST1" & as.character(group) != "Young"]
stopifnot(nrow(cf) == 9L, all(is.finite(cf$CF)))
pooled_rho <- unname(suppressWarnings(cor(cf$fraction, cf$CF, method = "spearman")))
resid_rank_rho <- function(v, cog, grp, batch) {
  design <- model.matrix(~ factor(grp) + factor(batch))
  x <- lm.fit(design, rank(v, ties.method = "average"))$residuals
  y <- lm.fit(design, rank(cog, ties.method = "average"))$residuals
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  unname(cor(x, y))
}
partial_rho <- resid_rank_rho(cf$fraction, cf$CF, cf$group, cf$sequence_batch)
# The nine-mouse, within-group permutation space is small enough to enumerate
# exactly (4! x 5! = 2,880 assignments). No Monte Carlo RNG is required.
all_permutations <- function(values) {
  if (length(values) == 1L) return(matrix(values, nrow = 1L))
  do.call(rbind, lapply(seq_along(values), function(i) {
    cbind(values[i], all_permutations(values[-i]))
  }))
}
cf_groups <- unique(as.character(cf$group))
stopifnot(length(cf_groups) == 2L)
cf_indices <- lapply(cf_groups, function(g) which(as.character(cf$group) == g))
stopifnot(identical(sort(lengths(cf_indices)), c(4L, 5L)))
cf_perm_a <- all_permutations(cf$CF[cf_indices[[1L]]])
cf_perm_b <- all_permutations(cf$CF[cf_indices[[2L]]])
null_cf <- numeric(nrow(cf_perm_a) * nrow(cf_perm_b))
iter <- 0L
for (i in seq_len(nrow(cf_perm_a))) {
  for (j in seq_len(nrow(cf_perm_b))) {
    perm <- cf$CF
    perm[cf_indices[[1L]]] <- cf_perm_a[i, ]
    perm[cf_indices[[2L]]] <- cf_perm_b[j, ]
    iter <- iter + 1L
    null_cf[iter] <- resid_rank_rho(cf$fraction, perm,
                                    cf$group, cf$sequence_batch)
  }
}
stopifnot(length(null_cf) == 2880L, all(is.finite(null_cf)))
# Allow for floating-point differences at ties of identical rank statistics.
n_cf_extreme <- sum(abs(null_cf) >= abs(partial_rho) - 1e-10)
p_cf <- n_cf_extreme / length(null_cf)
write_table(data.table(n_mice = 9L, pooled_spearman_rho = pooled_rho,
                       adjusted_rank_residual_rho = partial_rho,
                       within_group_permutation_P = p_cf,
                       B = length(null_cf), n_extreme = n_cf_extreme,
                       permutation_method = "exact_enumeration"),
            "AST1_cognitive_flexibility")

#Mouse-paired state-versus-rest pseudobulk
get_counts <- function(x) {
  z <- tryCatch(SeuratObject::LayerData(x, assay = "RNA", layer = "counts"),
                error = function(e) NULL)
  if (is.null(z) || ncol(z) != ncol(x)) {
    x <- SeuratObject::JoinLayers(x, assay = "RNA")
    z <- SeuratObject::LayerData(x, assay = "RNA", layer = "counts")
  }
  if (!identical(colnames(z), colnames(x))) stop("RNA counts and cell order differ")
  if (!inherits(z, "sparseMatrix")) z <- Matrix::Matrix(z, sparse = TRUE)
  if (any(Matrix::colSums(z) <= 0)) stop("Zero-depth astrocyte found")
  z
}
cnt <- get_counts(astro)
keys <- as.vector(outer(mouse$mouse_id, STATES, paste, sep = "__"))
key <- factor(paste(md$mouse_id, md$state, sep = "__"), levels = keys)
stopifnot(!anyNA(key))
X <- Matrix::sparse.model.matrix(~ 0 + key)
colnames(X) <- keys
pb <- cnt %*% X

# Construct mouse-state metadata
keymeta <- data.table(
  mouse_id = rep(mouse$mouse_id, times = length(STATES)),
  state = rep(STATES, each = nrow(mouse)),
  n_cells = tabulate(as.integer(key), nbins = length(keys))
)

# Add the key column separately to avoid data.table's reserved argument
keymeta[, key := keys]

# Restore the intended column order
setcolorder(keymeta, c("key", "mouse_id", "state", "n_cells"))

# Verify alignment
stopifnot(
  nrow(keymeta) == length(keys),
  identical(colnames(pb), keymeta$key),
  sum(keymeta$n_cells) == ncol(astro)
)

cat("Mouse-state pseudobulk metadata verified:",
    nrow(keymeta), "libraries\n")


stopifnot(identical(colnames(pb), keymeta$key), sum(keymeta$n_cells) == ncol(astro))
write_table(keymeta, "mouse_state_pseudobulk_ledger")
all_pb <- do.call(cbind, lapply(mouse$mouse_id, function(m)
  Matrix::rowSums(pb[, keymeta$mouse_id == m, drop = FALSE])))
colnames(all_pb) <- mouse$mouse_id

state_vs_rest <- function(target, included) {
  n_target <- vapply(mouse$mouse_id, function(m)
    sum(keymeta$n_cells[keymeta$mouse_id == m & keymeta$state %in% included]), integer(1))
  n_total <- vapply(mouse$mouse_id, function(m)
    sum(keymeta$n_cells[keymeta$mouse_id == m]), integer(1))
  eligible <- n_target >= 10L & (n_total - n_target) >= 30L
  ids <- mouse$mouse_id[eligible]
  if (length(ids) < 4L)
    return(list(table = NULL, status = data.table(state = target, n_mice = length(ids),
                                                   status = "INSUFFICIENT_MOUSE_COVERAGE")))
  target_mat <- do.call(cbind, lapply(ids, function(m)
    Matrix::rowSums(pb[, keymeta$mouse_id == m & keymeta$state %in% included,
                       drop = FALSE])))
  colnames(target_mat) <- paste0(ids, "__Target")
  rest <- all_pb[, ids, drop = FALSE] - target_mat
  colnames(rest) <- paste0(ids, "__Rest")
  design <- model.matrix(~ mouse + type,
                         data.frame(mouse = factor(rep(ids, 2L)),
                                    type = factor(rep(c("Rest", "Target"),
                                                      each = length(ids)),
                                                  levels = c("Rest", "Target"))))
  y <- edgeR::DGEList(counts = cbind(rest, target_mat))
  keep <- edgeR::filterByExpr(y, design = design)
  if (sum(keep) < 100L)
    return(list(table = NULL, status = data.table(state = target, n_mice = length(ids),
                                                   status = "TOO_FEW_TESTED_GENES")))
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design, robust = FALSE)
  fit <- edgeR::glmQLFit(y, design, robust = FALSE)
  t <- edgeR::glmQLFTest(fit, coef = "typeTarget")$table
  d <- data.table(gene = rownames(t), state = target, logFC = t$logFC,
                  F = t$F, PValue = t$PValue, FDR = p.adjust(t$PValue, "BH"))
  d[, s := sign(logFC) * sqrt(pmax(F, 0))]
  list(table = d, status = data.table(state = target, n_mice = length(ids),
                                      status = "PAIRED_PSEUDOBULK"))
}
targets <- c(setNames(as.list(STATES), STATES), list(cAST2_to_5 = MINOR))
identity <- lapply(names(targets), function(nm) state_vs_rest(nm, targets[[nm]]))
names(identity) <- names(targets)
write_table(rbindlist(lapply(identity, `[[`, "status"), fill = TRUE), "state_identity_coverage")
identity_de <- rbindlist(lapply(identity, `[[`, "table"), fill = TRUE)
if (!nrow(identity_de)) stop("No state has eligible paired identity pseudobulk")
write_table(identity_de, "state_identity_gene_effects")
# Enrichment is run in a fixed state-major, collection-minor order.
# Seeding immediately before this stage isolates the FGSEA RNG stream from
# stochastic operations elsewhere in the pipeline (e.g. permutations).
reactome_identity_data <- msigdbr::msigdbr(
  species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
reactome_identity <- lapply(split(reactome_identity_data$gene_symbol,
                                  reactome_identity_data$gs_name), unique)
identity_collections <- list(
  HALLMARK = lapply(hs, unique),
  GO_BP = lapply(gs, unique),
  REACTOME = reactome_identity
)
stopifnot(length(identity_collections$HALLMARK) == 50L,
          length(identity_collections$GO_BP) > 100L,
          length(identity_collections$REACTOME) > 100L)

set.seed(3101L)
identity_enrichments <- rbindlist(lapply(unique(identity_de$state), function(st) {
  d <- identity_de[state == st]
  ranks <- sort(setNames(d$s, d$gene)[is.finite(d$s)], decreasing = TRUE)
  rbindlist(lapply(names(identity_collections), function(coll) {
    fg <- suppressWarnings(fgsea::fgseaMultilevel(
      identity_collections[[coll]], ranks, minSize = 10,
      maxSize = 450, eps = 0, nPermSimple = 1000))
    result <- as.data.table(fg)[, .(pathway, NES, pval, padj)]
    result[, `:=`(state = st, collection = coll)]
    result
  }))
}))

identity_go <- identity_enrichments[
  collection == "GO_BP", .(pathway, NES, pval, padj, state)]
identity_hallmark <- identity_enrichments[
  collection == "HALLMARK", .(pathway, NES, pval, padj, state)]
identity_reactome <- identity_enrichments[
  collection == "REACTOME", .(pathway, NES, pval, padj, state)]
write_table(identity_go, "state_identity_GO_BP_enrichment")
write_table(identity_hallmark, "state_identity_Hallmark_enrichment")
write_table(identity_reactome, "state_identity_Reactome_enrichment")

mc_file <- file.path(out, "Mouse.MitoCarta3.0.xls")
if (!file.exists(mc_file)) utils::download.file(
  "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls",
  destfile = mc_file, mode = "wb")
mc_sheet <- grep("^A Mouse", readxl::excel_sheets(mc_file), value = TRUE)[1]
mc <- as.data.table(readxl::read_excel(mc_file, sheet = mc_sheet))
mc_path <- grep("MitoPathways$", names(mc), value = TRUE)[1]
mc <- mc[!is.na(get(mc_path)) & get(mc_path) != "0",
         .(gene = as.character(Symbol), path = as.character(get(mc_path)))]
mc <- mc[, .(path = trimws(unlist(strsplit(path, "\\|")))), by = gene]

mito_patterns <- list(
  RESP_CHAIN = "^OXPHOS", TCA = "TCA cycle", PYRUVATE = "Pyruvate",
  BETA_OX = "Fatty acid oxidation|Carnitine", KETONE = "Ketone",
  MITO_TRANSL = "central dogma|ribosome|Translation",
  MITO_IMPORT = "Protein import")
mito_modules <- lapply(mito_patterns, function(p) unique(mc[grepl(p, path), gene]))
reactome_data <- msigdbr::msigdbr(species = "Mus musculus", collection = "C2",
                                   subcollection = "CP:REACTOME")
reactome <- split(reactome_data$gene_symbol, reactome_data$gs_name)
# The existing Script 06 GO:BP collection is the same version used by mouse DE.
kegg <- msigdbr::msigdbr(species = "Mus musculus", collection = "C2",
                         subcollection = "CP:KEGG_LEGACY")
kegg <- split(kegg$gene_symbol, kegg$gs_name)
pick <- function(collection, patterns) {
  hits <- unique(unlist(lapply(patterns, grep, x = names(collection), value = TRUE)))
  unique(unlist(collection[hits], use.names = FALSE))
}
other_modules <- list(
  PEROX_LIPID = unique(c(pick(reactome, "^REACTOME_PEROXISOMAL_LIPID_METABOLISM"),
                         pick(gs, "BETA_OXIDATION_USING_ACYL_COA_OXIDASE"))),
  PEROX_BIOGEN = pick(reactome, "^REACTOME_PEROXISOMAL_PROTEIN_IMPORT"),
  LIPID_STORE = unique(c(pick(reactome, "^REACTOME_TRIGLYCERIDE_BIOSYNTHESIS"),
                         pick(gs, c("^GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT$",
                                    "^GOBP_LIPID_DROPLET_ORGANIZATION$")))),
  GLYCOLYSIS = unique(c(pick(reactome, "^REACTOME_GLYCOLYSIS$"),
                        pick(kegg, "^KEGG_GLYCOLYSIS_GLUCONEOGENESIS$"))),
  CHOLESTEROL = pick(reactome, "^REACTOME_CHOLESTEROL_BIOSYNTHESIS$"),
  ROS_DETOX = pick(reactome, "^REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES$"),
  GLUTATHIONE = unique(c(pick(reactome, "^REACTOME_GLUTATHIONE_CONJUGATION$"),
                         pick(kegg, "^KEGG_GLUTATHIONE_METABOLISM$"))))

#Assign mitochondrial genes and exclude them from non-mitochondrial modules
mito_members <- rbindlist(lapply(names(mito_modules), function(k)
  data.table(module = k, gene = mito_modules[[k]])))
mito_members[, priority := match(module, names(mito_modules))]
setorder(mito_members, gene, priority)
mito_members <- unique(mito_members, by = "gene")[, .(module, gene)]
mito_genes <- unique(mito_members$gene)
other_members <- rbindlist(lapply(names(other_modules), function(k)
  data.table(module = k, gene = setdiff(other_modules[[k]], mito_genes))))
module <- unique(rbindlist(list(mito_members, other_members)))
write_table(module, "metabolic_module_membership_USED")
mods <- split(module$gene, module$module)
required_mods <- c("RESP_CHAIN", "BETA_OX", "ROS_DETOX")
if (!all(required_mods %in% names(mods)))
  stop("Missing original named metabolic programs: ",
       paste(setdiff(required_mods, names(mods)), collapse = ", "))
ANTIOX8 <- c("Sod2", "Sod3", "Prdx3", "Prdx6", "Gpx1", "Txn1", "Glrx5", "Park7")
sets <- c(hs[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME")],
          mods[intersect(c("RESP_CHAIN", "TCA", "BETA_OX", "PEROX_LIPID",
                           "ROS_DETOX", "GLUTATHIONE", "GLYCOLYSIS"), names(mods))],
          list(ANTIOX8 = ANTIOX8))
sets <- lapply(sets, function(x) intersect(unique(x), rownames(cnt)))
write_table(rbindlist(lapply(names(sets), function(k)
  data.table(program = k, n_detected = length(sets[[k]]),
             genes = paste(sets[[k]], collapse = ";")))), "program_coverage")
stopifnot(length(sets$ANTIOX8) == 8L,
          length(sets$HALLMARK_OXIDATIVE_PHOSPHORYLATION) > 100L,
          length(sets$RESP_CHAIN) >= 5L,
          length(sets$BETA_OX) >= 5L)

#UCell respiratory chain set
uc <- UCell::ScoreSignatures_UCell(matrix = cnt,
                                   features = list(Respiration = sets$RESP_CHAIN),
                                   maxRank = 1500, ncores = 1)
if (!"Respiration_UCell" %in% colnames(uc)) stop("Unexpected UCell column names")
um <- as.data.table(Seurat::Embeddings(astro, "umap"), keep.rownames = "cell_id")
stopifnot(ncol(Seurat::Embeddings(astro, "umap")) == 2L)
setnames(um, old = names(um)[2:3], new = c("U1", "U2"))
um[, state := as.character(astro$clean_state)[match(cell_id, colnames(astro))]]
um[, Respiration_UCell := as.numeric(uc[match(cell_id, rownames(uc)),
                                           "Respiration_UCell"])]
stopifnot(nrow(um) == 10155L, !anyNA(um$Respiration_UCell))
write_table(um, "UMAP_and_respiratory_UCell")

#log1p(CP10K) within each cell, mean of measured genes within program, 
#then mouse x state average
lib <- as.numeric(Matrix::colSums(cnt))
program_genes <- unique(unlist(sets, use.names = FALSE))
lognorm <- cnt[program_genes, , drop = FALSE] %*%
  Matrix::Diagonal(x = 1e4 / lib)
dimnames(lognorm) <- list(program_genes, colnames(cnt))
lognorm@x <- log1p(lognorm@x)
score_base <- md[, .(cell_id, mouse_id, group, sequence_batch, state)]
score_base <- score_base[match(colnames(cnt), cell_id)]
stopifnot(identical(score_base$cell_id, colnames(cnt)))
for (k in names(sets)) {
  g <- sets[[k]]
  if (length(g) >= 3L) score_base[[k]] <- as.numeric(Matrix::colMeans(lognorm[g, , drop = FALSE]))
}
scores <- melt(score_base, id.vars = c("cell_id", "mouse_id", "group", "sequence_batch", "state"),
               variable.name = "program", value.name = "score")
mouse_scores <- scores[, .(score = mean(score), n_cells = .N),
                       by = .(mouse_id, group, sequence_batch, state, program)]
write_table(mouse_scores, "mouse_state_program_scores")
fig4g <- mouse_scores[state %in% MAJOR & program %in%
  c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "BETA_OX", "ANTIOX8")]
stopifnot(nrow(fig4g) == 13L * 2L * 3L)
write_table(fig4g, "Figure4g_three_programs_mouse_state")

#Within-state program and paired AST1-minus-AST0 group models
#BH: 6 state x program comparisons, 3 paired contrasts
aged_scores <- fig4g[as.character(group) != "Young" & n_cells >= MIN_SCORE_CELLS]
fit_score <- function(d, value_col = "score") {
  d <- as.data.frame(d)
  d$group <- factor(as.character(d$group), levels = c("Intact", "Impaired"))
  d$batch <- factor(d$sequence_batch)
  design <- model.matrix(~ batch + group, d)
  if (qr(design)$rank < ncol(design)) stop("Within-state score model rank deficient")
  co <- summary(lm(reformulate(c("batch", "group"), response = value_col), d))$coefficients
  data.table(estimate_intact_minus_impaired = -unname(co["groupImpaired", "Estimate"]),
             SE = unname(co["groupImpaired", "Std. Error"]),
             PValue = unname(co["groupImpaired", "Pr(>|t|)"]), n_mice = nrow(d))
}
state_fit <- aged_scores[, fit_score(.SD), by = .(state, program)]
state_fit[, FDR_six := p.adjust(PValue, method = "BH")]
write_table(state_fit, "within_state_program_group_models")
paired <- dcast(fig4g[n_cells >= MIN_SCORE_CELLS],
                mouse_id + group + sequence_batch + program ~ state, value.var = "score")
paired[, AST1_minus_AST0 := cAST1 - cAST0]
write_table(paired, "mouse_paired_AST1_minus_AST0")
paired_fit <- paired[as.character(group) != "Young" & is.finite(AST1_minus_AST0),
                     fit_score(.SD, value_col = "AST1_minus_AST0"), by = program]
paired_fit[, FDR_three := p.adjust(PValue, method = "BH")]
write_table(paired_fit, "paired_state_group_models")

#Descriptive six-state/three-compartment mixture on CP10K
MIX <- sets[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME",
              "RESP_CHAIN", "BETA_OX", "ROS_DETOX")]
MIX <- MIX[vapply(MIX, length, integer(1)) >= 10L]
mg <- unique(unlist(MIX, use.names = FALSE))
linear <- cnt[mg, , drop = FALSE] %*% Matrix::Diagonal(x = 1e4 / lib)
state_sum <- as.matrix(linear %*% X)
state_avg <- sweep(state_sum, 2, pmax(keymeta$n_cells, 1L), "/")
colnames(state_avg) <- keymeta$key
M <- Matrix::sparse.model.matrix(~ 0 + factor(md$mouse_id, levels = mouse$mouse_id))
colnames(M) <- mouse$mouse_id
whole <- sweep(as.matrix(linear %*% M), 2, as.numeric(table(factor(md$mouse_id,
                                                              levels = mouse$mouse_id))), "/")
colnames(whole) <- mouse$mouse_id

max_recon_error <- max(vapply(mouse$mouse_id, function(m) {
  ix <- which(keymeta$mouse_id == m)
  nn <- sum(keymeta$n_cells[ix])
  pred <- as.numeric(state_sum[, ix, drop = FALSE] %*% rep(1 / nn, length(ix)))
  max(abs(pred - whole[, m]))
}, numeric(1)))
if (!is.finite(max_recon_error) || max_recon_error >= 1e-7)
  stop("Mouse-specific linear mixture identity failed: ", max_recon_error)
write_table(data.table(check = "mouse_specific_linear_CP10K_identity",
                       max_absolute_error = max_recon_error, tolerance = 1e-7),
            "composition_math_check")
ageI <- mouse[as.character(group) == "Intact", mouse_id]
ageM <- mouse[as.character(group) == "Impaired", mouse_id]
aged_ids <- c(ageI, ageM)
frac <- dcast(ledger, mouse_id ~ state, value.var = "fraction")
F <- as.matrix(frac[, ..STATES]); rownames(F) <- frac$mouse_id
ref6 <- sapply(STATES, function(st) {
  ix <- which(keymeta$mouse_id %in% aged_ids & keymeta$state == st &
                keymeta$n_cells >= MIN_REF_CELLS)
  if (length(ix) < 2L) return(rep(NA_real_, length(mg)))
  rowMeans(state_avg[, ix, drop = FALSE])
})
rownames(ref6) <- mg
minor_mean <- sapply(aged_ids, function(m) {
  ix <- which(keymeta$mouse_id == m & keymeta$state %in% MINOR)
  n <- sum(keymeta$n_cells[ix])
  if (n < MIN_REF_CELLS) return(rep(NA_real_, length(mg)))
  rowSums(state_sum[, ix, drop = FALSE]) / n
})
rownames(minor_mean) <- mg
ok_minor <- colSums(is.finite(minor_mean)) == nrow(minor_mean)
ref_minor <- if (sum(ok_minor) >= 2L)
  rowMeans(minor_mean[, ok_minor, drop = FALSE]) else rep(NA_real_, length(mg))
ref3 <- cbind(cAST0 = ref6[, "cAST0"], cAST1 = ref6[, "cAST1"],
              cAST2_to_5 = ref_minor)
P6 <- rbind(Intact = colMeans(F[ageI, , drop = FALSE]),
            Impaired = colMeans(F[ageM, , drop = FALSE]))
P3 <- cbind(cAST0 = P6[, "cAST0"], cAST1 = P6[, "cAST1"],
            cAST2_to_5 = rowSums(P6[, MINOR, drop = FALSE]))
observed <- cbind(Intact = rowMeans(whole[, ageI, drop = FALSE]),
                  Impaired = rowMeans(whole[, ageM, drop = FALSE]))
miss <- matrix(NA_real_, length(mg), 2, dimnames = list(mg, c("Intact", "Impaired")))
pred6 <- if (all(is.finite(ref6))) ref6 %*% t(P6) else miss
pred3 <- if (all(is.finite(ref3))) ref3 %*% t(P3) else miss
colnames(pred6) <- colnames(pred3) <- c("Intact", "Impaired")
log_effect <- function(m, pc = 0.1)
  if (anyNA(m)) rep(NA_real_, nrow(m)) else log2((m[, "Intact"] + pc) /
                                                  (m[, "Impaired"] + pc))
mix_genes <- data.table(gene = mg, observed = log_effect(observed),
                        six_state_pred = log_effect(pred6),
                        three_comp_pred = log_effect(pred3))
write_table(mix_genes, "composition_reconstruction_gene_effects")
mix_summary <- rbindlist(lapply(names(MIX), function(k) {
  d <- mix_genes[gene %in% MIX[[k]]]
  data.table(program = k, n_genes = nrow(d), observed = median(d$observed),
             six_state_pred = median(d$six_state_pred),
             three_comp_pred = median(d$three_comp_pred))
}))
write_table(mix_summary, "composition_reconstruction_program_summary")
write_table(data.table(state = STATES, n_aged_mice_ref = vapply(STATES, function(st)
  sum(keymeta$mouse_id %in% aged_ids & keymeta$state == st &
        keymeta$n_cells >= MIN_REF_CELLS), integer(1))),
  "composition_reference_coverage")

#Within-state edgeR cognition/aging gene effects for matched-gene sensitivity
fit_state_de <- function(st, comparison = c("cognition", "aging")) {
  comparison <- match.arg(comparison)
  eligible <- keymeta[state == st & n_cells >= MIN_PB_CELLS]
  info <- merge(eligible[, .(key, mouse_id)], mouse[, .(mouse_id, group,
                          sequence_batch)], by = "mouse_id")
  if (comparison == "cognition") {
    info <- info[as.character(group) %in% c("Intact", "Impaired")]
    info[, cmp := factor(as.character(group), levels = c("Intact", "Impaired"))]
  } else {
    info <- info[as.character(group) %in% c("Young", "Intact")]
    info[, cmp := factor(as.character(group), levels = c("Young", "Intact"))]
  }
  if (length(unique(info$cmp)) < 2L || min(table(info$cmp)) < 3L)
    stop("Insufficient within-state mouse coverage: ", st, " / ", comparison)
  setorder(info, mouse_id)
  dm <- model.matrix(~ factor(sequence_batch) + cmp, info)
  if (qr(dm)$rank != ncol(dm)) dm <- model.matrix(~ cmp, info)
  y <- edgeR::DGEList(counts = as.matrix(pb[, info$key, drop = FALSE]))
  keep <- edgeR::filterByExpr(y, design = dm)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- edgeR::estimateDisp(y, dm, robust = TRUE)
  fit <- edgeR::glmQLFit(y, dm, robust = TRUE)
  coef <- grep("^cmp", colnames(dm))
  stopifnot(length(coef) == 1L)
  tab <- edgeR::glmQLFTest(fit, coef = coef)$table
  #cognition test coef = Impaired-Intact; aging = Intact-Young;
  #negate both to orient toward preserved cognition/younger age.
  data.table(state = st, comparison = comparison, gene = rownames(tab),
             logFC_favorable = -tab$logFC, logCPM = tab$logCPM,
             F = tab$F, PValue = tab$PValue, FDR = p.adjust(tab$PValue, "BH"))
}
state_de <- rbindlist(lapply(MAJOR, function(st)
  rbindlist(lapply(c("cognition", "aging"), function(k) fit_state_de(st, k)))))
write_table(state_de, "within_state_cognition_aging_gene_effects")


#mm39 gene lengths, additional AST1 abundance QC covariates and
#length/expression-matched within-state programs 
glen_file <- file.path(root, "reference", "mouse_mm39_gene_lengths.csv")
if (file.exists(glen_file)) {
  gl <- fread(glen_file)
} else {
  if (requireNamespace("TxDb.Mmusculus.UCSC.mm39.knownGene", quietly = TRUE)) {
    txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene::TxDb.Mmusculus.UCSC.mm39.knownGene
  } else {
    txdb <- TxDb.Mmusculus.UCSC.mm39.refGene::TxDb.Mmusculus.UCSC.mm39.refGene
  }
  gr <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
  sym <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, names(gr),
                                "SYMBOL", "ENTREZID")
  gl <- data.table(gene = unname(as.character(sym)),
                   length = as.numeric(IRanges::width(gr)))
  gl <- gl[!is.na(gene), .(length = max(length)), by = gene]
  dir.create(dirname(glen_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(gl, glen_file)
}
stopifnot(all(c("gene", "length") %in% names(gl)), !anyDuplicated(gl$gene),
          all(gl$length > 0))
write_table(gl, "mouse_mm39_gene_lengths_USED")

long_genes <- intersect(gl[length > 100000, gene], rownames(cnt))
if (length(long_genes) < 10L) stop("Unexpectedly few long genes in mm39 reference")
md[, long_frac := as.numeric(Matrix::colSums(cnt[long_genes, , drop = FALSE]) / lib)]
if (!"nFeature_RNA" %in% names(md))
  md[, nFeature_RNA := as.numeric(Matrix::colSums(cnt > 0))]
if (!"percent.mt" %in% names(md)) {
  mt <- grep("^mt-", rownames(cnt), value = TRUE)
  md[, percent.mt := 100 * as.numeric(Matrix::colSums(cnt[mt, , drop = FALSE])) / lib]
}
qc_mouse <- md[, .(nFeature = median(as.numeric(nFeature_RNA)),
                   long_frac = median(long_frac),
                   pct_mt = median(as.numeric(`percent.mt`)),
                   cAST1 = mean(state == "cAST1")),
               by = .(mouse_id, group, sequence_batch)]
write_table(qc_mouse, "astrocyte_QC_by_mouse")


#Whole-lineage metabolic-program sensitivity 
lineage_pb <- readRDS(file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_pseudobulk_raw_counts.rds"))
s5_ids <- colnames(lineage_pb[["Astro-TE NN"]])
s5_info <- copy(mouse[match(s5_ids, mouse_id)])
stopifnot(identical(s5_info$mouse_id, s5_ids),
          identical(colnames(lineage_pb[["Peri NN"]]), s5_ids))
s5_info[, batch := factor(sequence_batch)]
s5_info[, group := factor(as.character(group), levels = c("Young", "Intact", "Impaired"))]

#Gene-standardized cytosolic-ribosomal logCPM
s5_astro <- as.matrix(lineage_pb[["Astro-TE NN"]])
s5_lcpm <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(s5_astro),
                                               method = "TMM"), log = TRUE)
s5_ribo <- grep("^R[Pp][SsLl][0-9]+[AaXxYy]?[0-9]*$",
                rownames(s5_lcpm), value = TRUE)
s5_ribo_z <- t(scale(t(s5_lcpm[s5_ribo, , drop = FALSE])))
s5_ribo_z <- s5_ribo_z[rowSums(!is.finite(s5_ribo_z)) == 0L, , drop = FALSE]
stopifnot(nrow(s5_ribo_z) > 10L)
s5_info[, ribo_score := colMeans(s5_ribo_z)]
s5_info[, ribo_frac := as.numeric(Matrix::colSums(s5_astro[s5_ribo, , drop = FALSE]) /
                                        Matrix::colSums(s5_astro))]
s5_info[, `:=`(pct_mt = qc_mouse$pct_mt[match(mouse_id, qc_mouse$mouse_id)],
               nFeature = qc_mouse$nFeature[match(mouse_id, qc_mouse$mouse_id)],
               long_frac = qc_mouse$long_frac[match(mouse_id, qc_mouse$mouse_id)],
               cAST1 = qc_mouse$cAST1[match(mouse_id, qc_mouse$mouse_id)])]

s5_ambient <- fread(file.path(root, "07_mouse_sensitivity",
                              "07_baiermouse_astrocyte_cognition_astrocyte_ambient_score_by_mouse.csv"))
s5_diss <- fread(file.path(root, "07_mouse_sensitivity",
                           "07_baiermouse_alllineages_dissociation_dissociation_scores_by_mouse_lineage.csv"))
s5_info[, ambient := s5_ambient$astro_ambient_score[
  match(mouse_id, s5_ambient$mouse_id)]]
s5_info[, weight_loss := abs(as.numeric(
  beh$weight_loss_pct[match(mouse_id, beh$mouse_id)]))]
write_table(s5_info[, .(mouse_id, group, batch, nFeature, long_frac, pct_mt,
                        cAST1, ribo_score, ribo_frac, ambient, weight_loss)],
            "S5_mouse_covariate_values")

ANTIOX30 <- c("Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx4", "Prdx1",
              "Prdx2", "Prdx3", "Prdx4", "Prdx5", "Prdx6", "Txn1", "Txn2",
              "Txnrd1", "Txnrd2", "Gsr", "Gclc", "Gclm", "Gss", "Mt1",
              "Mt2", "Mt3", "Glrx", "Glrx5", "Selenow", "Park7",
              "Srxn1", "Nqo1", "Hmox1")
s5_hallmark_names <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                       "HALLMARK_FATTY_ACID_METABOLISM",
                       "HALLMARK_PEROXISOME", "HALLMARK_ADIPOGENESIS",
                       "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
s5_module_names <- c("RESP_CHAIN", "TCA", "BETA_OX", "PEROX_LIPID",
                     "ROS_DETOX", "GLUTATHIONE", "GLYCOLYSIS", "CHOLESTEROL")
s5_sets <- c(hs[s5_hallmark_names], mods[s5_module_names],
             list(ANTIOXIDANT30 = ANTIOX30))
stopifnot(length(s5_sets) == 14L, !anyNA(names(s5_sets)),
          all(lengths(s5_sets) >= 5L))
write_table(rbindlist(lapply(names(s5_sets), function(k)
  data.table(pathway = k, gene = unique(s5_sets[[k]])))),
  "S5_program_membership")

s5_fit <- function(lineage, covariate = NULL, omit_mouse = NULL,
                   drop_mt_genes = FALSE) {
  keep <- if (is.null(omit_mouse)) rep(TRUE, nrow(s5_info)) else
    s5_ids != omit_mouse
  cnt <- as.matrix(lineage_pb[[lineage]][, keep, drop = FALSE])
  storage.mode(cnt) <- "double"
  if (drop_mt_genes)
    cnt <- cnt[!grepl("^mt-", rownames(cnt)), , drop = FALSE]
  g <- droplevels(s5_info$group[keep])
  b <- droplevels(s5_info$batch[keep])
  design <- if (nlevels(b) > 1L) model.matrix(~ 0 + g + b) else
    model.matrix(~ 0 + g)
  colnames(design) <- sub("^g", "", colnames(design))
  if (!is.null(covariate)) {
    v <- covariate[keep]
    design <- cbind(design, covariate_z = as.numeric(scale(v)))
  }
  if (qr(design)$rank != ncol(design))
    stop("Rank-deficient sensitivity model: ", lineage, " / ", omit_mouse)
  y <- edgeR::DGEList(counts = cnt)
  expressed <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[expressed, , keep.lib.sizes = FALSE],
                              method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  v <- setNames(numeric(ncol(design)), colnames(design))
  v[c("Intact", "Impaired")] <- c(1, -1)
  tab <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = v),
                        n = Inf, sort.by = "none")$table
  data.table(gene = rownames(tab), logFC = tab$logFC,
             s = sign(tab$logFC) * sqrt(pmax(tab$F, 0)))
}

s5_enrich <- function(de, pathway_sets = s5_sets) {
  ranks <- sort(setNames(de$s, de$gene)[is.finite(de$s)], decreasing = TRUE)
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = pathway_sets, stats = ranks, minSize = 5L,
    maxSize = 500L, eps = 0,
    BPPARAM = BiocParallel::SerialParam()))
  as.data.table(fg)[, .(pathway, NES, pval, padj, size)]
}

s5_diss_for <- function(lineage) {
  d <- s5_diss[s5_diss[["lineage"]] == lineage]
  d$diss_score[match(s5_ids, d$mouse_id)]
}

set.seed(20261004)
s5_results <- rbindlist(lapply(c("Astro-TE NN", "Peri NN"), function(lineage) {
  covars <- list(primary = NULL,
                pct_mt = s5_info$pct_mt,
                ribosomal = s5_info$ribo_score,
                genes_per_cell = s5_info$nFeature,
                ambient = s5_info$ambient,
                dissociation = s5_diss_for(lineage),
                weight_loss = s5_info$weight_loss)
  z <- rbindlist(lapply(names(covars), function(k) {
    message("25A S5: ", lineage, "; model = ", k)
    d <- s5_enrich(s5_fit(lineage, covariate = covars[[k]]))
    d[, model := k]
    d
  }))
  no_mt <- s5_enrich(s5_fit(lineage, drop_mt_genes = TRUE))
  no_mt[, model := "mt_genes_removed"]
  omitted <- rbindlist(lapply(s5_ids[as.character(s5_info$group) != "Young"], function(m) {
    d <- s5_enrich(s5_fit(lineage, omit_mouse = m))
    d[, model := paste0("drop_", m)]
    d
  }))
  r <- rbindlist(list(z, no_mt, omitted))
  r[, lineage := lineage]
  r
}))
write_table(s5_results, "S5_metabolic_module_sensitivity_all_models")
write_table(s5_results[!grepl("^drop_", model)],
            "S5_metabolic_module_sensitivity_covariates")
write_table(s5_results[grepl("^drop_", model),
                       .(min_NES = min(NES), max_NES = max(NES),
                         max_FDR = max(padj), n_omissions = .N),
                       by = .(lineage, pathway)],
            "S5_metabolic_module_leave_one_aged_mouse_out")

#Ribosomal-score diagnostic mouse-level QC relationships
s5_aged <- s5_info[as.character(group) != "Young"]
s5_ribosomal_cor <- rbindlist(lapply(c("nFeature", "pct_mt", "long_frac",
                                      "cAST1", "ribo_frac"), function(k) {
  ct <- suppressWarnings(cor.test(s5_aged$ribo_score, s5_aged[[k]],
                                  method = "spearman", exact = FALSE))
  data.table(variable = k, n_mice = nrow(s5_aged),
             rho = unname(ct$estimate), PValue = ct$p.value)
}))
write_table(s5_ribosomal_cor, "ribosomal_mouse_QC_correlations")
qa <- qc_mouse[as.character(group) != "Young"]
qa[, `:=`(group = factor(as.character(group), levels = c("Intact", "Impaired")),
           batch = factor(sequence_batch), lg = qlogis(cAST1))]
if (!all(is.finite(qa$lg))) stop("AST1 fraction is zero or one in QC logit model")
qc_formulas <- list(
  base = lg ~ group + batch,
  nFeature = lg ~ group + batch + scale(nFeature),
  long_frac = lg ~ group + batch + scale(long_frac),
  pct_mt = lg ~ group + batch + scale(pct_mt)
)
qc_fit <- rbindlist(lapply(names(qc_formulas), function(k) {
  mod <- lm(qc_formulas[[k]], data = qa)
  co <- summary(mod)$coefficients
  data.table(model = k,
             intact_minus_impaired_logit = -unname(co["groupImpaired", "Estimate"]),
             SE = unname(co["groupImpaired", "Std. Error"]),
             PValue = unname(co["groupImpaired", "Pr(>|t|)"]), n_mice = nrow(qa))
}))
write_table(qc_fit, "AST1_abundance_QC_covariate_models")

wd <- fread(whole_de_file)[lineage == "Astro-TE NN"]
stopifnot(all(c("gene", "logFC", "PValue", "FDR", "logCPM") %in% names(wd)))
wcov <- merge(gl, wd[, .(gene, logCPM)], by = "gene")
matched_test <- function(effect, members, nperm = 2000L) {
  d <- merge(data.table(gene = names(effect), lfc = unname(effect)),
             wcov, by = "gene")[is.finite(lfc) & length > 0 & is.finite(logCPM)]
  if (!nrow(d)) return(NULL)
  d[, `:=`(zl = as.numeric(scale(log10(length))),
             ze = as.numeric(scale(logCPM)))]
  S <- d[gene %in% members]
  C <- d[!gene %in% members]
  if (nrow(S) < 8L || nrow(C) < 100L) return(NULL)
  nearest <- lapply(seq_len(nrow(S)), function(i)
    head(order((C$zl - S$zl[i])^2 + (C$ze - S$ze[i])^2), 15L))
  null_med <- replicate(nperm, median(C$lfc[vapply(nearest, function(ix)
    ix[sample.int(length(ix), 1L)], integer(1))]))
  obs <- median(S$lfc)
  data.table(n_genes = nrow(S), obs_median = obs,
             null_median = median(null_med), excess = obs - median(null_med),
             p_emp_upper = (1 + sum(null_med >= obs))/(1 + nperm),
             p_emp_lower = (1 + sum(null_med <= obs))/(1 + nperm))
}
key_programs <- intersect(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "RESP_CHAIN",
                            "BETA_OX", "HALLMARK_FATTY_ACID_METABOLISM",
                            "ROS_DETOX", "GLYCOLYSIS"), names(sets))
matched <- rbindlist(lapply(split(state_de, by = c("state", "comparison")), function(d) {
  effect <- setNames(d$logFC_favorable, d$gene)
  rbindlist(lapply(key_programs, function(k) {
    ans <- matched_test(effect, sets[[k]])
    if (is.null(ans)) return(NULL)
    ans[, `:=`(state = d$state[1], comparison = d$comparison[1], program = k)]
    ans
  }), fill = TRUE)
}), fill = TRUE)
write_table(matched, "matched_gene_within_state_programs")

#Whole-astrocyte respiratory-complex gene directions
complex_genes <- list(
  I = c("Ndufs1","Ndufv2","Ndufa6","Ndufb8","Ndufs2","Ndufv1","Ndufs8","Ndufs7","Ndufs6","Ndufb9","Ndufb5","Ndufa9","Ndufa2","Ndufs4","Ndufb6","Ndufa5","Ndufs5","Ndufb10","Ndufa10","Ndufa7","Ndufv3","Ndufb2","Ndufc2","Ndufa8","Ndufb11","Ndufa12","Ndufs3","Ndufb7","Ndufc1","Ndufa3","Ndufab1","Ndufa13","Ndufb3","Ndufa11","mt-Nd5","mt-Nd4","mt-Nd2","Ndufa1","Ndufb4","mt-Nd3","mt-Nd1","mt-Nd4l","mt-Nd6"),
  II = c("Sdha","Sdhb","Sdhd","Sdhc"),
  III = c("Cyc1","Uqcrc1","Uqcrfs1","Uqcrc2","Uqcrq","Uqcr10","Uqcrh","Uqcrb","Uqcr11","mt-Cytb"),
  IV = c("Cox5a","Cox6a1","Cox7a2","Cox6c","Cox6b1","Cox4i1","Cox5b","Cox8b","Cox7c","Cox7a1","mt-Co2","Cox8a","Ndufa4","Cox4i2","Cox7b","Cox6a2","Cox6b2","mt-Co1","Cox8c","mt-Co3","Cox7a2l"),
  V = c("Atp5d","Atp5a1","Atp5e","Atp5b","Atp5o","Atp5k","Atp5j","Atp5c1","Atpif1","Atp5l","mt-Atp6","Atp5g3","Atp5h","Atp5pb","Atp5mpl","Atp5j2","Atp5g2","mt-Atp8","Atp5g1","Dmac2l","Atp5md")
)
complex_map <- rbindlist(lapply(names(complex_genes), function(k)
  data.table(gene = complex_genes[[k]], complex = k)))
stopifnot(!anyDuplicated(complex_map$gene))
complex_dt <- merge(complex_map, wd[, .(gene, logFC, PValue, FDR, logCPM)],
                    by = "gene")
complex_dt[, `:=`(intact_minus_impaired_log2FC = -logFC,
                   genome = fifelse(grepl("^mt-", gene), "mtDNA", "nuclear"),
                   positive_intact = -logFC > 0)]
write_table(complex_dt, "respiratory_complex_gene_effects")
write_table(complex_dt[, .(n_measured = .N, n_positive_intact = sum(positive_intact),
                            median_log2FC = median(intact_minus_impaired_log2FC)),
                         by = complex], "respiratory_complex_summary")