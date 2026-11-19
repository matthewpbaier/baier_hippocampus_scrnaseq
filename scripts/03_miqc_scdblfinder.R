source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))

#Baier et. al
#03: QC via miQC and scDblFinder on SoupX-processed objects

base_dir <- hca_output_root
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
miqc_dir <- file.path(base_dir, "03_final_qc", "miQC")
miqc_plot_dir <- file.path(miqc_dir, "plots")
miqc_obj_dir <- file.path(miqc_dir, "objects")
doublet_dir <- file.path(base_dir, "03_final_qc", "scDblFinder")
doublet_obj_dir <- file.path(doublet_dir, "objects")
for (d in c(miqc_dir, miqc_plot_dir, miqc_obj_dir, doublet_dir, doublet_obj_dir)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

#Load libraries
library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(miQC)
library(patchwork)
library(scDblFinder)

#Load checkpoint from 02_import_soupx.R
s2_checkpoint <- readRDS(file.path(base_dir,"objects","02_baiermouse_allcells_soupx_checkpoint.rds"))

meta <- s2_checkpoint$meta
soupx_summary <- s2_checkpoint$soupx_summary

meta <- meta %>% mutate(sequence_id = as.character(sequence_id))

soupx_dir <- file.path(base_dir, "02_qc","SoupX_corrected")

meta <- meta %>% mutate(soupx_file = file.path(soupx_dir, paste0("02_baiermouse_", mouse_id, "_soupx_corrected_counts.rds")))

meta %>% select(mouse_id,sequence_id, soupx_file)
stopifnot(all(file.exists(meta$soupx_file)))

test_counts <- readRDS(meta$soupx_file[1])
dim(test_counts)

grep("^mt-", rownames(test_counts), value = TRUE
)

rm(test_counts)
gc()

#Run miQC
run_miqc_sample <- function(
    mouse_id,
    meta,
    posterior_cutoff = 0.75
) {

  i <- which(meta$mouse_id == mouse_id)

  if (length(i) != 1) {
    stop(
      "Could not uniquely identify mouse: ",
      mouse_id
    )
  }

  #Load SoupX corrected counts
  counts <- readRDS(meta$soupx_file[i])

  #Create SingleCellExperiment
  sce <- SingleCellExperiment(assays = list(counts = counts))

  #Retain original 10X barcode
  colData(sce)$barcode_10x <- colnames(sce)

  #Retain sample-level metadata
  colData(sce)$mouse_id <- meta$mouse_id[i]
  colData(sce)$sequence_id <- meta$sequence_id[i]
  colData(sce)$sequence_batch <- meta$sequence_batch[i]
  colData(sce)$age_group <- meta$age_group[i]
  colData(sce)$primary_group <- meta$primary_group[i]

  #Retain QC metrics
  mt_genes <- grepl("^mt-", rownames(sce))
  feature_ctrls <- list(mito = rownames(sce)[mt_genes])
  sce <- addPerCellQC(sce, subsets = feature_ctrls)

  #Specify miQC model
  model <- mixtureModel(sce, model_type = "linear")

  #miQC's optional boundary corrections use data.frame subset assignments,
  #when no cells satisfy a correction, those assignments can fail on a
  #zero-row subset--apply each correction only when it can change calls,
  #otherwise retain the original posterior cutoff and filtering behavior.
  model_intercepts <- vapply(
    1:2,
    function(k) as.numeric(flexmix::parameters(model, component = k)[1]),
    numeric(1)
  )
  compromised_component <- if (model_intercepts[1] > model_intercepts[2]) 1L else 2L
  intact_component <- 3L - compromised_component
  mt_percent <- as.numeric(colData(sce)$subsets_mito_percent)
  posterior_compromised <- flexmix::posterior(model)[, compromised_component]
  below_intact_boundary <- mt_percent < flexmix::fitted(model)[, intact_component]
  apply_boundary <- any(below_intact_boundary)
  apply_left_cutoff <- any(
    posterior_compromised > posterior_cutoff & !below_intact_boundary
  )

  message(
    "miQC ", mouse_id, ": ",
    sum(posterior_compromised > posterior_cutoff & !below_intact_boundary),
    " cells marked for removal before optional left-cutoff correction"
  )

  #Diagnostic plots
  p_metrics <- plotMetrics(sce) + ggtitle(paste0(mouse_id, " — QC metrics"))
  p_model <- plotModel(sce, model) + ggtitle(paste0(mouse_id," — miQC model"))

  p_filter <- plotFiltering(
    sce, model,
    posterior_cutoff = posterior_cutoff,
    keep_all_below_boundary = apply_boundary,
    enforce_left_cutoff = apply_left_cutoff
  ) + ggtitle(paste0(mouse_id, " — filtering"))

  diagnostic_plot <- p_metrics | p_model | p_filter

  ggsave(file.path(miqc_plot_dir,paste0("03_baiermouse_", mouse_id, "_miqc_diagnostic.pdf")), diagnostic_plot, width = 15, height = 4.5, units = "in")

  #Apply filter
  cells_before <- colnames(sce)

  sce_filtered <- filterCells(
    sce,
    model,
    posterior_cutoff = posterior_cutoff,
    keep_all_below_boundary = apply_boundary,
    enforce_left_cutoff = apply_left_cutoff,
    verbose = TRUE
  )

  cells_after <- colnames(
    sce_filtered
  )

  retained <- cells_before %in%
    cells_after

  #Save barcode-level calls
  cell_calls <- tibble(
    mouse_id = mouse_id,
    barcode_10x = cells_before,
    miqc_keep = retained
  )

  write_csv(
    cell_calls,
    file.path(
      miqc_dir,
      paste0("03_baiermouse_", mouse_id, "_miqc_cell_calls.csv")
    )
  )

  #Sample summary
  summary_row <- tibble(
    mouse_id = mouse_id,
    sequence_id = meta$sequence_id[i],
    sequence_batch = meta$sequence_batch[i],
    age_group = meta$age_group[i],
    primary_group = meta$primary_group[i],

    n_before_miqc = ncol(sce),
    n_after_miqc = ncol(sce_filtered),

    n_removed_miqc =
      ncol(sce) -
      ncol(sce_filtered),

    pct_retained_miqc =
      100 *
      ncol(sce_filtered) /
      ncol(sce),

    median_genes_before =
      median(
        colData(sce)$detected
      ),

    median_mt_before =
      median(
        colData(sce)$subsets_mito_percent
      ),

    median_genes_after =
      median(
        colData(sce_filtered)$detected
      ),

    median_mt_after =
      median(
        colData(sce_filtered)$subsets_mito_percent
      )
  )

  #Save filtered object and model
  saveRDS(
    sce_filtered,
    file.path(
      miqc_obj_dir,
      paste0("03_baiermouse_", mouse_id, "_miqc_filtered_cells.rds")
    )
  )

  saveRDS(
    model,
    file.path(
      miqc_obj_dir,
      paste0("03_baiermouse_", mouse_id, "_miqc_fitted_model.rds")
    )
  )

  rm(
    counts,
    sce,
    sce_filtered,
    model
  )

  gc()

  return(
    summary_row
  )
}

miqc_results <- vector(
  "list",
  nrow(meta)
)

names(miqc_results) <-
  meta$mouse_id

for (mouse in meta$mouse_id) {

  miqc_results[[mouse]] <-
    run_miqc_sample(
      mouse_id = mouse,
      meta = meta,
      posterior_cutoff = 0.75
    )
}

miqc_summary <- bind_rows(
  miqc_results
)

miqc_summary

write_csv(
  miqc_summary,
  file.path(
    miqc_dir,
    "03_baiermouse_allcells_miqc_sample_summary.csv"
  )
)

miqc_overall <- miqc_summary %>%
  summarise(
    total_before =
      sum(n_before_miqc),

    total_after =
      sum(n_after_miqc),

    total_removed =
      sum(n_removed_miqc),

    pct_retained =
      100 *
      total_after /
      total_before
  )

miqc_overall

#Doublet removal via scDblFinder
run_scdblfinder_sample <- function(
    mouse_id,
    meta,
    miqc_obj_dir,
    seed = 1234
) {

  i <- which(meta$mouse_id == mouse_id)

  sce <- readRDS(
    file.path(
      miqc_obj_dir,
      paste0("03_baiermouse_", mouse_id, "_miqc_filtered_cells.rds")
    )
  )

  n_before <- ncol(sce)

  set.seed(seed)

  sce <- scDblFinder(
    sce
  )

  dbl_table <- table(
    colData(sce)$scDblFinder.class
  )

  print(dbl_table)

  calls <- tibble(
    mouse_id = mouse_id,
    barcode_10x = colnames(sce),
    scDblFinder_score =
      colData(sce)$scDblFinder.score,
    scDblFinder_class =
      colData(sce)$scDblFinder.class
  )

  write_csv(
    calls,
    file.path(
      doublet_dir,
      paste0("03_baiermouse_", mouse_id, "_doublet_calls.csv")
    )
  )

  #Retain singlets
  sce_singlet <- sce[
    ,
    colData(sce)$scDblFinder.class ==
      "singlet"
  ]

  n_after <- ncol(sce_singlet)

  summary_row <- tibble(
    mouse_id = mouse_id,
    sequence_id = meta$sequence_id[i],
    sequence_batch = meta$sequence_batch[i],
    age_group = meta$age_group[i],
    primary_group = meta$primary_group[i],

    n_before_doublet = n_before,
    n_singlets = n_after,
    n_doublets = n_before - n_after,

    pct_doublets =
      100 * (n_before - n_after) / n_before,

    pct_singlets =
      100 * n_after / n_before
  )

  #Save object with all calls
  saveRDS(
    sce,
    file.path(
      doublet_obj_dir,
      paste0("03_baiermouse_", mouse_id, "_doublet_all_cells.rds")
    )
  )

  #Save clean singlet object
  saveRDS(
    sce_singlet,
    file.path(
      doublet_obj_dir,
      paste0("03_baiermouse_", mouse_id, "_final_singlets.rds")
    )
  )

  rm(
    sce,
    sce_singlet
  )

  gc()

  summary_row
}

dbl_results <- vector(
  "list",
  nrow(meta)
)

names(dbl_results) <- meta$mouse_id

for (mouse in meta$mouse_id) {

  dbl_results[[mouse]] <-
    run_scdblfinder_sample(
      mouse_id = mouse,
      meta = meta,
      miqc_obj_dir = miqc_obj_dir,
      seed = 1234
    )
}

dbl_summary <- bind_rows(
  dbl_results
)

dbl_summary

write_csv(
  dbl_summary,
  file.path(
    doublet_dir,
    "03_baiermouse_allcells_doublets_sample_summary.csv"
  )
)

#Final QC summary
final_qc_summary <- miqc_summary %>%
  left_join(
    dbl_summary %>%
      select(
        mouse_id,
        n_before_doublet,
        n_singlets,
        n_doublets,
        pct_doublets
      ),
    by = "mouse_id"
  ) %>%
  mutate(
    pct_of_preSoup_remaining =
      100 * n_singlets /
      n_before_miqc
  )

final_qc_summary

final_overall <- final_qc_summary %>%
  summarise(
    post_SoupX_cells =
      sum(n_before_miqc),

    post_miQC_cells =
      sum(n_after_miqc),

    final_singlets =
      sum(n_singlets),

    miQC_removed =
      sum(n_removed_miqc),

    doublets_removed =
      sum(n_doublets),

    final_retention_pct =
      100 *
      final_singlets /
      post_SoupX_cells
  )

final_overall

final_qc_summary <- miqc_summary %>%
  left_join(
    dbl_summary %>%
      select(
        mouse_id,
        n_before_doublet,
        n_singlets,
        n_doublets,
        pct_doublets,
        pct_singlets
      ),
    by = "mouse_id"
  ) %>%
  mutate(
    pct_final_vs_preSoup =
      100 * n_singlets / n_before_miqc
  )

write_csv(
  final_qc_summary,
  file.path(
    base_dir,
    "03_final_qc",
    "03_baiermouse_allcells_final_qc_summary.csv"
  )
)

s3_final_checkpoint <- list(
  meta = meta,
  soupx_summary = soupx_summary,
  miqc_summary = miqc_summary,
  dbl_summary = dbl_summary,
  final_qc_summary = final_qc_summary,
  final_overall = final_overall
)

saveRDS(
  s3_final_checkpoint,
  file.path(
    base_dir,
    "objects",
    "03_baiermouse_allcells_final_qc_checkpoint.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    base_dir,
    "03_final_qc",
    "03_baiermouse_allcells_final_qc_sessioninfo.txt"
  )
)
