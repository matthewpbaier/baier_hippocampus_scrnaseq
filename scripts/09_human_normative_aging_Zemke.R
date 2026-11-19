#Baier et al.
#09: Normative human hippocampal aging (Zemke et al.; GEO GSE278576)

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "09_human_normative_aging")
for (p in c(out, file.path(out, "checkpoints")))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(Matrix)
library(Seurat)
library(edgeR)
library(fgsea)
library(msigdbr)
library(data.table)
library(dplyr)
library(tidyr)
library(readr)

set.seed(20261002)

#GEO files
metadata_name <- "GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz"
metadata_file <- file.path(hca_zemke_source_root, metadata_name)
if (!file.exists(metadata_file))
  stop("Missing Zemke metadata: ", metadata_file, "\nSee input/README.md (GSE278576).")
meta <- readr::read_tsv(metadata_file, show_col_types = FALSE)

#Original GEO metadata uses the literal column name 'bacrode'
needed <- c("orig.ident", "Age", "Gender", "age_group", "subclass", "bacrode")
if (!all(needed %in% names(meta)))
  stop("Missing original metadata columns: ",
       paste(setdiff(needed, names(meta)), collapse = ", "))

donors <- unique(as.character(meta[["orig.ident"]]))
stopifnot(length(donors) == 40L, !anyNA(donors), all(nzchar(donors)))
raw_names <- paste0("GSE278576_", donors, "_raw_feature_bc_matrix.h5")
raw_paths <- setNames(file.path(hca_zemke_source_root, raw_names), donors)
missing_h5 <- raw_paths[!file.exists(raw_paths)]
if (length(missing_h5))
  stop("Missing ", length(missing_h5), " Zemke donor H5 files in ",
       hca_zemke_source_root, ":\n", paste(basename(missing_h5), collapse = "\n"))

celltype_map <- c(Astro = "Astrocyte", Microglia = "Microglia",
                  OPC = "OPC", Oligo = "Oligodendrocyte",
                  Endo = "Endothelial", DG = "DG neuron")
aggregation_types <- c("Astrocyte", "Microglia", "OPC",
                       "Oligodendrocyte", "DG neuron", "Endothelial")
primary_types <- c("Astrocyte", "Microglia", "OPC", "Oligodendrocyte")

meta <- meta %>%
  mutate(donor_id = as.character(.data[["orig.ident"]]),
         Age = as.numeric(Age), Gender = as.character(Gender),
         barcode_raw = sub("^.*_", "", as.character(bacrode)),
         crossstudy_celltype = unname(celltype_map[as.character(subclass)]))

donor_info <- meta %>% distinct(donor_id, Age, Gender, age_group) %>%
  arrange(Age, donor_id)
donor_info <- donor_info %>%
  mutate(Gender = factor(Gender, levels = c("F", "M")),
         age_z = as.numeric(scale(Age)))
readr::write_csv(donor_info, file.path(out, "09_zemke_alllineages_aging_donor_design.csv"))

donor_order <- donor_info$donor_id

# Donor-by-cell-type coverage
coverage <- tidyr::crossing(
  donor_id = donor_order,
  crossstudy_celltype = aggregation_types
) %>%
  dplyr::left_join(
    meta %>%
      dplyr::filter(!is.na(crossstudy_celltype)) %>%
      dplyr::count(
        donor_id,
        crossstudy_celltype,
        name = "n_nuclei"
      ),
    by = c("donor_id", "crossstudy_celltype")
  ) %>%
  dplyr::mutate(
    n_nuclei = tidyr::replace_na(n_nuclei, 0L)
  ) %>%
  dplyr::left_join(
    donor_info,
    by = "donor_id"
  )

# Validate the complete 40-donor x 6-cell-type matrix
stopifnot(
  nrow(coverage) == 240L,
  !anyNA(coverage$n_nuclei),
  sum(coverage$n_nuclei) ==
    sum(!is.na(meta$crossstudy_celltype)),
  !anyDuplicated(coverage[, c("donor_id", "crossstudy_celltype")])
)

cat("Zemke donor coverage validated:",
    nrow(coverage), "donor-cell-type combinations\n")


readr::write_csv(coverage, file.path(out, "09_zemke_alllineages_aging_lineage_coverage.csv"))
coverage_primary <- coverage %>% filter(crossstudy_celltype %in% primary_types)

#Aggregate QC-filtered nuclei from donor-specific 10x H5 files
#Reuse the checkpoint location on later runs
new_ck <- file.path(out, "checkpoints")
dir.create(new_ck, recursive = TRUE, showWarnings = FALSE)

read_donor <- function(d, reference_genes = NULL) {
  ck <- file.path(new_ck, paste0("09_zemke_", d, "_alllineages_pseudobulk_counts.rds"))
  if (file.exists(ck)) {
    x <- readRDS(ck)
    if (!identical(colnames(x), aggregation_types))
      stop("Unexpected six-lineage pseudobulk columns for ", d)
    if (!is.null(reference_genes) && !identical(rownames(x), reference_genes))
      stop("Donor checkpoint gene universe mismatch for ", d)
    return(list(path = ck, genes = rownames(x)))
  }

  raw_file <- raw_paths[[d]]
  raw <- Seurat::Read10X_h5(raw_file, use.names = TRUE, unique.features = TRUE)
  if (is.list(raw)) {
    if (!"Gene Expression" %in% names(raw))
      stop("Missing Gene Expression modality: ", d)
    raw <- raw[["Gene Expression"]]
  }
  if (!is.null(reference_genes)) {
    if (!setequal(rownames(raw), reference_genes))
      stop("Gene universe differs among donor H5 files: ", d)
    raw <- raw[reference_genes, , drop = FALSE]
  } else {
    reference_genes <- rownames(raw)
  }

  md <- meta %>% filter(donor_id == d,
                        crossstudy_celltype %in% aggregation_types)
  idx <- match(md$barcode_raw, colnames(raw))
  if (anyNA(idx))
    stop("QC-filtered nuclei missing from original raw count matrix: ", d)
  design <- Matrix::sparse.model.matrix(~ 0 + factor(md$crossstudy_celltype,
                                         levels = aggregation_types))
  colnames(design) <- aggregation_types
  x <- raw[, idx, drop = FALSE] %*% design
  if (!identical(rownames(x), reference_genes)) stop("Gene order mismatch")
  stopifnot(identical(colnames(x), aggregation_types))
  saveRDS(x, file.path(new_ck, paste0("09_zemke_", d, "_alllineages_pseudobulk_counts.rds")))
  rm(raw, design, x)
  invisible(gc())
  list(path = file.path(new_ck, paste0("09_zemke_", d, "_alllineages_pseudobulk_counts.rds")),
       genes = reference_genes)
}

ck_paths <- setNames(character(length(donor_order)), donor_order)
reference_genes <- NULL
for (d in donor_order) {
  message("Zemke donor ", d)
  z <- read_donor(d, reference_genes)
  if (is.null(reference_genes)) reference_genes <- z$genes
  ck_paths[[d]] <- z$path
}

#Assemble one gene-by-donor pseudobulk matrix per lineage
first <- readRDS(ck_paths[[1]])
reference_genes <- rownames(first)
rm(first)
pb <- setNames(lapply(aggregation_types, function(ct) {
  mats <- lapply(donor_order, function(d) {
    z <- readRDS(ck_paths[[d]])
    if (!identical(rownames(z), reference_genes))
      stop("Gene mismatch during donor-pseudobulk assembly: ", d)
    x <- z[, ct, drop = FALSE]
    colnames(x) <- d
    x
  })
  do.call(cbind, mats)
}), aggregation_types)
if (!all(vapply(pb[primary_types], function(x)
  identical(colnames(x), donor_order) && all(Matrix::colSums(x) > 0), logical(1))))
  stop("Primary lineage donor coverage or ordering is invalid")
saveRDS(pb, file.path(out, "09_zemke_alllineages_aging_six_lineage_pseudobulk.rds"))

#Sex-adjusted continuous-age differential expression, donor as replicate
fit_aging <- function(counts, info) {
  stopifnot(identical(colnames(counts), info$donor_id))
  design <- model.matrix(~ Gender + age_z, data = info)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient age model")
  y <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = "age_z")
  res <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table
  res <- data.table(gene = rownames(res), as.data.table(res))
  # Original model: positive logFC = older. Favorable axis: younger = -logFC.
  res[, `:=`(logFC_younger = -logFC,
             younger_signed_sqrtF = -sign(logFC) * sqrt(pmax(F, 0)))]
  list(results = res, dge = y, design = design,
       n_genes_tested = sum(keep))
}

age_de <- setNames(lapply(primary_types, function(ct) {
  message("Sex-adjusted continuous-age DE: ", ct)
  fit_aging(pb[[ct]], donor_info)
}), primary_types)
for (ct in primary_types)
  data.table::fwrite(age_de[[ct]]$results,
    file.path(out, paste0("09_zemke_", tolower(gsub("[^A-Za-z0-9]+", "_", ct)), "_aging_differential_expression.csv")))

#Within-lineage Hallmark enrichment (oriented toward younger age)
hallmarks <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
sets <- split(hallmarks$gene_symbol, hallmarks$gs_name)
stopifnot(length(sets) == 50L)
saveRDS(sets, file.path(out, "09_zemke_alllineages_aging_human_hallmark_sets.rds"))

run_hallmark <- function(de, rank_column) {
  x <- de[, .(gene, statistic = get(rank_column))]
  x <- x[!is.na(gene) & is.finite(statistic)]
  setorder(x, -statistic)
  x <- x[!duplicated(gene)]
  ranks <- sort(setNames(x$statistic, x$gene), decreasing = TRUE)
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = sets, stats = ranks, minSize = 10, maxSize = 500,
    eps = 0, BPPARAM = BiocParallel::SerialParam()))
  as.data.table(fg)[, .(pathway, NES, pval, padj)]
}

hallmark_inference <- rbindlist(lapply(primary_types, function(ct) {
  z <- run_hallmark(age_de[[ct]]$results, "younger_signed_sqrtF")
  z[, crossstudy_celltype := ct]
  z
}))
hallmark_logfc <- rbindlist(lapply(primary_types, function(ct) {
  z <- run_hallmark(age_de[[ct]]$results, "logFC_younger")
  z[, crossstudy_celltype := ct]
  z
}))
fwrite(hallmark_inference,
       file.path(out, "09_zemke_alllineages_aging_hallmark_younger_age_inference.csv"))
fwrite(hallmark_logfc,
       file.path(out, "09_zemke_alllineages_aging_hallmark_younger_age_log_fc_rank.csv"))

#Checkpoint for cross-cohort integration
saveRDS(list(
  zemke_age_de = age_de,
  zemke_sample_info = donor_info,
  zemke_hallmark_inference = hallmark_inference,
  zemke_hallmark_logfc = hallmark_logfc,
  human_hallmark_sets = sets
), file.path(out, "09_zemke_alllineages_aging_normative_aging_checkpoint.rds"))
