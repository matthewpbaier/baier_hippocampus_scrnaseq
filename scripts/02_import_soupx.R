#Baier et. al
#02: Import of Cell Ranger processed data and SoupX

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
base_dir <- hca_output_root
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c("02_mouse_metadata", "02_qc", "objects")) dir.create(file.path(base_dir,d), recursive=TRUE, showWarnings=FALSE)

#Load libraries
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(scales)
library(SoupX)

# Read sequenced sample metadata and Cell Ranger matrices
raw_root <- hca_mouse_input_root
meta <- readr::read_csv(file.path(raw_root, "metadata", "sample_metadata.csv"))

meta %>% count(age_group, primary_group)
meta %>% count(sequence_batch, primary_group)

meta %>%
  select(
    mouse_id,
    sequence_id,
    sequence_batch,
    age_group,
    primary_group,
    cf51to61
  )

#Load Cell Ranger-processed matrices
matrix_files <- list.files(
  raw_root,
  pattern = "matrix\\.mtx(\\.gz)?$",
  recursive = TRUE,
  full.names = TRUE
)

matrix_files

matrix_inventory <- tibble(
  matrix_file = matrix_files
) %>%
  mutate(
    matrix_dir = dirname(matrix_file),

    matrix_type = case_when(
      str_detect(matrix_dir, "filtered_feature_bc_matrix") ~ "filtered",
      str_detect(matrix_dir, "raw_feature_bc_matrix") ~ "raw",
      TRUE ~ "unknown"
    )
  )

matrix_inventory

matrix_inventory %>%
  count(matrix_type)

seq_ids <- as.character(meta$sequence_id)

extract_seq_id <- function(path, seq_ids) {
  hits <- seq_ids[
    map_lgl(
      seq_ids,
      ~ str_detect(
        path,
        paste0(
          "(^|[/\\\\])",
          fixed(.x),
          "([/\\\\]|$)"
        )
      )
    )
  ]
  if (length(hits) == 1) {
    return(hits)
  }
  NA_character_
}

matrix_inventory <- matrix_inventory %>%
  mutate(
    sequence_id = map_chr(
      matrix_dir,
      extract_seq_id,
      seq_ids = seq_ids
    )
  )

matrix_inventory %>%
  select(
    sequence_id,
    matrix_type,
    matrix_dir
  ) %>%
  arrange(sequence_id, matrix_type)

matrix_inventory %>%
  count(sequence_id, matrix_type)

matrix_inventory %>%
  filter(is.na(sequence_id))

matrix_inventory %>%
  count(sequence_id) %>%
  filter(n != 2)

filtered_lookup <- matrix_inventory %>%
  filter(matrix_type == "filtered") %>%
  select(
    sequence_id,
    filtered_dir = matrix_dir
  )

raw_lookup <- matrix_inventory %>%
  filter(matrix_type == "raw") %>%
  select(
    sequence_id,
    raw_dir = matrix_dir
  )

meta <- meta %>%
  mutate(
    sequence_id = as.character(sequence_id)
  ) %>%
  left_join(
    filtered_lookup,
    by = "sequence_id"
  ) %>%
  left_join(
    raw_lookup,
    by = "sequence_id"
  )

meta %>%
  select(
    mouse_id,
    sequence_id,
    filtered_dir,
    raw_dir
  )

#Load filtered matrices into Seurat objects
read_10x_gene_expression <- function(path) {
  x <- Read10X(
    data.dir = path
  )

  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      x <- x[["Gene Expression"]]
    } else {

      stop(
        paste(
          "Multiple modalities found but no Gene Expression assay:",
          path
        )
      )
    }
  }
  x
}

seurat_list <- vector(
  "list",
  nrow(meta)
)

names(seurat_list) <- meta$mouse_id

for (i in seq_len(nrow(meta))) {
  counts <- read_10x_gene_expression(
    meta$filtered_dir[i]
  )

  obj <- CreateSeuratObject(
    counts = counts,
    project = meta$mouse_id[i],
    min.cells = 0,
    min.features = 0
  )

  obj <- RenameCells(
    obj,
    add.cell.id = meta$mouse_id[i]
  )

  obj$mouse_id <- meta$mouse_id[i]
  obj$sequence_id <- meta$sequence_id[i]
  obj$sequence_batch <- meta$sequence_batch[i]
  obj$age_group <- meta$age_group[i]
  obj$primary_group <- meta$primary_group[i]

  seurat_list[[i]] <- obj

  rm(counts, obj)
  gc()
}

#Verify number of cells from Cell Ranger
input_counts <- tibble(
  mouse_id = names(seurat_list),

  n_cells_cellranger =
    map_int(
      seurat_list,
      ncol
    ),

  n_features_matrix =
    map_int(
      seurat_list,
      nrow
    )
)

input_counts <- meta %>%
  select(
    mouse_id,
    sequence_id,
    sequence_batch,
    age_group,
    primary_group
  ) %>%
  left_join(
    input_counts,
    by = "mouse_id"
  )

input_counts

write_csv(
  input_counts,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_cellranger_filtered_counts.csv")
)

#Calculate QC metrics: % mitochondrial, % ribosomal, % hemoglobin
for (i in seq_along(seurat_list)) {

  obj <- seurat_list[[i]]

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[ab]")

  seurat_list[[i]] <- obj

  rm(obj)
}

#Pre-QC summary
qc_summary_list <- map(
  seurat_list,
  function(obj) {

    md <- obj[[]]

    tibble(
      mouse_id = unique(md$mouse_id),

      sequence_id = unique(md$sequence_id),
      sequence_batch = unique(md$sequence_batch),
      age_group = unique(md$age_group),
      primary_group = unique(md$primary_group),
      n_cells = ncol(obj),
      median_umi = median(md$nCount_RNA),
      q1_umi =quantile(md$nCount_RNA, 0.25),
      q3_umi = quantile(md$nCount_RNA, 0.75),
      median_features = median(md$nFeature_RNA),
      q1_features = quantile(md$nFeature_RNA, 0.25),
      q3_features = quantile(md$nFeature_RNA, 0.75),
      median_percent_mt = median(md$percent.mt),
      q95_percent_mt = quantile(md$percent.mt, 0.95),
      median_percent_ribo = median(md$percent.ribo),
      median_percent_hb =median(md$percent.hb)
    )
  }
)

qc_summary <- bind_rows(qc_summary_list)
qc_summary

write_csv(qc_summary,file.path(base_dir, "02_qc/02_baiermouse_allcells_initial_qc_sample_summary.csv"))

qc_cell <- map_dfr(
  seurat_list,
  function(obj) {

    obj[[]] %>%
      rownames_to_column(
        "cell_id"
      ) %>%
      select(
        cell_id,
        mouse_id,
        sequence_id,
        sequence_batch,
        age_group,
        primary_group,
        nCount_RNA,
        nFeature_RNA,
        percent.mt,
        percent.ribo,
        percent.hb
      )
  }
)

#Provisional QC prior to SoupX
preliminary_qc_summary <- qc_cell %>%
  group_by(mouse_id) %>%
  summarise(
    n_total = n(),

    n_feature_lt_200 =
      sum(nFeature_RNA < 200),

    pct_feature_lt_200 =
      100 * mean(nFeature_RNA < 200),

    n_count_lt_500 =
      sum(nCount_RNA < 500),

    pct_count_lt_500 =
      100 * mean(nCount_RNA < 500),

    n_mt_gt_20 =
      sum(percent.mt > 20),

    pct_mt_gt_20 =
      100 * mean(percent.mt > 20),

    n_mt_gt_30 =
      sum(percent.mt > 30),

    pct_mt_gt_30 =
      100 * mean(percent.mt > 30),

    n_fail_permissive =
      sum(
        nFeature_RNA < 200 |
        nCount_RNA < 500 |
        percent.mt > 30
      ),

    pct_fail_permissive =
      100 * mean(
        nFeature_RNA < 200 |
        nCount_RNA < 500 |
        percent.mt > 30
      ),

    .groups = "drop"
  )

preliminary_qc_summary

preSoup_list <- lapply(
  seurat_list,
  function(obj) {

    subset(
      obj,
      subset =
        nFeature_RNA >= 200 &
        nCount_RNA >= 500 &
        percent.mt <= 30
    )
  }
)

names(preSoup_list) <- names(seurat_list)

preSoup_retention <- tibble(
  mouse_id = names(seurat_list),

  before =
    map_int(seurat_list, ncol),

  after =
    map_int(preSoup_list, ncol)
) %>%

  mutate(
    removed = before - after,
    pct_retained =
      100 * after / before
  )

preSoup_retention

write_csv(
  preSoup_retention,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_presoupx_qc_retention.csv")
)

#Temporary clustering for SoupX
for (i in seq_along(preSoup_list)) {
  obj <- preSoup_list[[i]]

  obj <- NormalizeData(
    obj,
    verbose = FALSE
  )

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )

  obj <- ScaleData(
    obj,
    features = VariableFeatures(obj),
    verbose = FALSE
  )

  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = 30,
    verbose = FALSE
  )

  obj <- FindNeighbors(
    obj,
    dims = 1:20,
    verbose = FALSE
  )

  obj <- FindClusters(
    obj,
    resolution = 0.5,
    verbose = FALSE
  )

  preSoup_list[[i]] <- obj
}

temp_clusters <- map_dfr(
  preSoup_list,
  function(obj) {

    tibble(
      mouse_id =
        unique(obj$mouse_id),

      n_cells =
        ncol(obj),

      n_clusters =
        length(unique(Idents(obj)))
    )
  }
)

temp_clusters

#Alignment of cell barcodes for SoupX
for (i in seq_along(preSoup_list)) {
  obj <- preSoup_list[[i]]

  mouse <- unique(obj$mouse_id)

  obj$barcode_10x <- sub(
    paste0("^", mouse, "_"),
    "",
    colnames(obj)
  )

  preSoup_list[[i]] <- obj
}

head(
  preSoup_list[[1]]$barcode_10x
)

#SoupX
dir.create(
  file.path(base_dir, "02_qc/SoupX_corrected"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(base_dir, "02_qc/SoupX_diagnostics"),
  recursive = TRUE,
  showWarnings = FALSE
)

run_soupx_sample <- function(
    mouse_id,
    meta,
    preSoup_list,
    seed = 1234
) {

  i <- which(meta$mouse_id == mouse_id)

  if (length(i) != 1) {
    stop("Could not uniquely identify mouse: ", mouse_id)
  }

  obj <- preSoup_list[[mouse_id]]

  #Original 10X barcodes
  keep_barcodes <- unname(
    obj$barcode_10x
  )
  
  #Load raw droplets and Cell Ranger called cells
  tod <- Read10X(
    data.dir = meta$raw_dir[i]
  )

  toc_full <- Read10X(
    data.dir = meta$filtered_dir[i]
  )

  #Check features match
  stopifnot(
    identical(
      rownames(tod),
      rownames(toc_full)
    )
  )
  
  #Restrict filtered matrix to cells retained after permissive QC
  toc <- toc_full[
    ,
    keep_barcodes,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(toc),
      keep_barcodes
    )
  )

  #Temporary clustering
  clusters <- as.character(
    Idents(obj)
  )

  names(clusters) <- unname(
    obj$barcode_10x
  )

  clusters <- clusters[
    colnames(toc)
  ]

  stopifnot(
    !anyNA(clusters)
  )

  #Create SoupChannel
  sc <- SoupChannel(
    tod = tod,
    toc = toc
  )

  sc <- setClusters(
    sc,
    clusters
  )

  #Save ambient profile
  soup_profile <- sc$soupProfile %>%
    rownames_to_column("gene") %>%
    arrange(desc(est))

  set.seed(seed)

  #Contamination estimate
  sc <- autoEstCont(
    sc,
    doPlot = FALSE
  )

  rho <- median(
    sc$metaData$rho,
    na.rm = TRUE
  )

  #Adjust counts
  set.seed(seed)

  corrected <- adjustCounts(
    sc,
    roundToInt = TRUE
  )

  stopifnot(
    identical(
      rownames(corrected),
      rownames(toc)
    )
  )

  stopifnot(
    identical(
      colnames(corrected),
      colnames(toc)
    )
  )

  original_umis <- sum(toc)
  corrected_umis <- sum(corrected)

  pct_removed <-
    100 *
    (1 - corrected_umis / original_umis)

  gene_correction <- tibble(
    gene = rownames(toc),

    original_counts =
      Matrix::rowSums(toc),

    corrected_counts =
      Matrix::rowSums(corrected)
  ) %>%
    mutate(
      counts_removed =
        original_counts - corrected_counts,

      fraction_removed =
        if_else(
          original_counts > 0,
          counts_removed / original_counts,
          NA_real_
        )
    ) %>%
    arrange(desc(counts_removed))

  #Save corrected counts for downstream analysis
  saveRDS(
    corrected,
    file.path(
      file.path(base_dir, "02_qc/SoupX_corrected"),
      paste0("02_baiermouse_", mouse_id, "_soupx_corrected_counts.rds")
    )
  )

  write_csv(
    soup_profile,
    file.path(
      file.path(base_dir, "02_qc/SoupX_diagnostics"),
      paste0("02_baiermouse_", mouse_id, "_soupx_ambient_profile.csv")
    )
  )

  write_csv(
    gene_correction,
    file.path(
      file.path(base_dir, "02_qc/SoupX_diagnostics"),
      paste0("02_baiermouse_", mouse_id, "_soupx_gene_correction.csv")
    )
  )

  summary_row <- tibble(
    mouse_id = mouse_id,
    sequence_id = meta$sequence_id[i],
    sequence_batch = meta$sequence_batch[i],
    age_group = meta$age_group[i],
    primary_group = meta$primary_group[i],

    n_cells = ncol(toc),
    n_clusters = length(unique(clusters)),

    rho = rho,

    original_umis = original_umis,
    corrected_umis = corrected_umis,
    pct_umis_removed = pct_removed
  )

  rm(
    tod,
    toc_full,
    toc,
    sc,
    corrected
  )

  gc()

  list(
    summary = summary_row,

    top_soup =
      soup_profile %>%
      slice_head(n = 25) %>%
      mutate(mouse_id = mouse_id),

    top_correction =
      gene_correction %>%
      slice_head(n = 25) %>%
      mutate(mouse_id = mouse_id)
  )
}

soupx_results <- vector(
  "list",
  length = nrow(meta)
)

names(soupx_results) <- meta$mouse_id

for (mouse in meta$mouse_id) {

  soupx_results[[mouse]] <-
    run_soupx_sample(
      mouse_id = mouse,
      meta = meta,
      preSoup_list = preSoup_list,
      seed = 1234
    )
}

soupx_summary <- map_dfr(
  soupx_results,
  "summary"
)

soupx_summary

write_csv(
  soupx_summary,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_soupx_sample_summary.csv")
)

soupx_top_genes <- map_dfr(
  soupx_results,
  "top_soup"
)

write_csv(
  soupx_top_genes,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_soupx_ambient_top_genes.csv")
)

soupx_top_corrections <- map_dfr(
  soupx_results,
  "top_correction"
)

write_csv(
  soupx_top_corrections,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_soupx_corrected_top_genes.csv")
)

saveRDS(
  meta,
  file.path(base_dir, "02_mouse_metadata/02_baiermouse_allmice_sequencing_metadata.rds")
)

write_csv(
  meta,
  file.path(base_dir, "02_mouse_metadata/02_baiermouse_allmice_sequencing_metadata.csv")
)

s2_checkpoint <- list(

  meta = meta,

  input_counts = input_counts,

  qc_summary = qc_summary,

  preliminary_qc_summary =
    preliminary_qc_summary,

  preSoup_retention =
    preSoup_retention,

  temp_clusters =
    temp_clusters,

  soupx_summary =
    soupx_summary
)

saveRDS(
  s2_checkpoint,
  file.path(base_dir, "objects/02_baiermouse_allcells_soupx_checkpoint.rds")
)

soupx_cluster_assignments <- purrr::map_dfr(
  preSoup_list,
  function(obj) {

    tibble(
      mouse_id = unique(obj$mouse_id),
      barcode_10x = unname(obj$barcode_10x),
      soupx_temp_cluster =
        as.character(Idents(obj))
    )
  }
)

dim(soupx_cluster_assignments)
head(soupx_cluster_assignments)

write_csv(
  soupx_cluster_assignments,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_soupx_preliminary_clusters.csv")
)

saveRDS(
  soupx_cluster_assignments,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_soupx_preliminary_clusters.rds")
)

saveRDS(
  qc_cell,
  file.path(base_dir, "02_qc/02_baiermouse_allcells_presoupx_qc_metrics.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(base_dir, "02_qc/02_baiermouse_allcells_soupx_sessioninfo.txt")
)
