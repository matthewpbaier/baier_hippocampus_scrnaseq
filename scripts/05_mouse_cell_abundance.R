#Baier et. al
#05: Differential abundance analysis 

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
base_dir <- hca_output_root
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(Seurat)
library(tidyverse)
library(limma)
library(speckle)
library(patchwork)

#Set output directory
out_dir <- file.path(
  base_dir,
  "05_differential_abundance"
)

fig_dir <- file.path(
  out_dir,
  "figures"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

atlas_file <- hca_atlas_input
atlas <- readRDS(atlas_file)
ncol(atlas)

cell_meta <- atlas[[]] %>%
  tibble::rownames_to_column(
    "cell_id"
  )

stopifnot(
  "allen_subclass" %in%
    colnames(atlas[[]])
)

#Reconstruct annotation labels if needed
if (!"broad_lineage" %in% colnames(atlas[[]])) {

  atlas$broad_lineage <- dplyr::case_when(

    atlas$allen_subclass == "Microglia NN" ~
      "Microglia",

    atlas$allen_subclass %in% c(
      "BAM NN",
      "Monocytes NN",
      "DC NN",
      "Lymphoid NN"
    ) ~
      "Other immune",

    atlas$allen_subclass == "Oligo NN" ~
      "Oligodendrocyte",

    atlas$allen_subclass == "OPC NN" ~
      "OPC",

    atlas$allen_subclass %in% c(
      "Astro-TE NN",
      "Astro-NT NN"
    ) ~
      "Astrocyte",

    atlas$allen_subclass == "CHOR NN" ~
      "Choroid plexus",

    atlas$allen_subclass %in% c(
      "Ependymal NN",
      "Astroependymal NN"
    ) ~
      "Ependymal",

    atlas$allen_subclass == "Endo NN" ~
      "Endothelial",

    atlas$allen_subclass == "Peri NN" ~
      "Pericyte",

    atlas$allen_subclass == "SMC NN" ~
      "SMC",

    atlas$allen_subclass == "VLMC NN" ~
      "VLMC",

    stringr::str_detect(
      atlas$allen_subclass,
      "IMN"
    ) ~
      "Immature neuron",

    stringr::str_detect(
      atlas$allen_subclass,
      "Glut"
    ) ~
      "Excitatory neuron",

    stringr::str_detect(
      atlas$allen_subclass,
      "Gaba|GABA|Gly"
    ) ~
      "Inhibitory neuron",

    TRUE ~
      "Other/rare"
  )
}

if (!"anatomic_scope" %in% colnames(atlas[[]])) {

  atlas$anatomic_scope <- dplyr::case_when(

    atlas$allen_subclass %in% c(
      "CHOR NN",
      "Ependymal NN",
      "Astroependymal NN",
      "Hypendymal NN",
      "Tanycyte NN"
    ) ~
      "Adjacent ventricular",

    atlas$allen_subclass %in% c(
      "BAM NN",
      "Lymphoid NN",
      "DC NN",
      "Monocytes NN"
    ) ~
      "Immune/interface",

    atlas$allen_subclass %in% c(
      "Endo NN",
      "Peri NN",
      "SMC NN",
      "VLMC NN"
    ) ~
      "Vascular",

    atlas$allen_subclass %in% c(
      "Microglia NN",
      "Oligo NN",
      "OPC NN",
      "Astro-TE NN",
      "Astro-NT NN",
      "DG Glut",
      "HPF CR Glut",
      "CA1-ProS Glut",
      "CA3 Glut"
    ) ~
      "Hippocampal/parenchymal",

    TRUE ~
      "Other/reference-only"
  )
}

cell_meta <- atlas[[]] %>%
  tibble::rownames_to_column(
    "cell_id"
  )

required_cols <- c(
  "mouse_id",
  "sequence_batch",
  "age_group",
  "primary_group",
  "broad_lineage",
  "allen_subclass",
  "anatomic_scope"
)

stopifnot(all(required_cols %in% names(cell_meta)))

# Recompute the original 04b eligibility rule from recovered cells.
# Primary: >=100 cells in >=8 of the nine aged mice; secondary: >=50 in >=7.
ct <- data.table::as.data.table(cell_meta)
aged_ids <- unique(as.character(ct[
  tolower(as.character(primary_group)) %in% c("intact", "impaired"), mouse_id]))
counts_aged <- ct[mouse_id %in% aged_ids,
  .(n = .N), by = .(allen_subclass, mouse_id)]
all_counts <- merge(
  data.table::CJ(allen_subclass = sort(unique(as.character(ct$allen_subclass))),
                 mouse_id = sort(aged_ids)),
  counts_aged, by = c("allen_subclass", "mouse_id"), all.x = TRUE)
all_counts[is.na(n), n := 0L]
eligibility <- all_counts[, .(
  total_aged_cells = sum(n), n_aged_mice_present = sum(n > 0),
  n_aged_mice_ge100 = sum(n >= 100), n_aged_mice_ge50 = sum(n >= 50),
  min_cells = min(n), median_cells = median(n), max_cells = max(n)),
  by = allen_subclass]
eligibility[, statistical_tier := data.table::fcase(
  n_aged_mice_ge100 >= 8L, "Primary",
  n_aged_mice_ge50 >= 7L, "Secondary", default = "Descriptive")]
data.table::fwrite(eligibility, file.path(
  out_dir, "05_baiermouse_subclasses_allgroups_statistical_eligibility.csv"))

sample_info <- cell_meta %>%
  dplyr::distinct(
    mouse_id,
    sequence_batch,
    age_group,
    primary_group
  ) %>%
  dplyr::arrange(
    mouse_id
  )

stopifnot(
  nrow(sample_info) == 13
)

sample_info

readr::write_csv(
  sample_info,
  file.path(
    out_dir,
    "05_baiermouse_allmice_abundance_sample_metadata.csv"
  )
)

broad_counts <- cell_meta %>%
  dplyr::count(
    mouse_id,
    broad_lineage,
    name = "n_cells"
  ) %>%

  dplyr::group_by(
    mouse_id
  ) %>%

  dplyr::mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%

  dplyr::ungroup() %>%

  dplyr::left_join(
    sample_info,
    by = "mouse_id"
  )

readr::write_csv(
  broad_counts,
  file.path(
    out_dir,
    "05_baiermouse_alllineages_allgroups_proportions_by_mouse.csv"
  )
)

subclass_counts <- cell_meta %>%
  dplyr::count(
    mouse_id,
    allen_subclass,
    name = "n_cells"
  ) %>%

  dplyr::group_by(
    mouse_id
  ) %>%

  dplyr::mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%

  dplyr::ungroup() %>%

  dplyr::left_join(
    sample_info,
    by = "mouse_id"
  )

readr::write_csv(
  subclass_counts,
  file.path(
    out_dir,
    "05_baiermouse_subclasses_allgroups_proportions_by_mouse.csv"
  )
)

#Propeller helper
make_prop_list <- function(
    dat,
    celltype_col,
    sample_info,
    keep_types = NULL
) {

  dat <- dat %>%
    dplyr::filter(
      !is.na(.data[[celltype_col]]),
      !is.na(mouse_id)
    )

  sample_levels <-
    sample_info$mouse_id

  clusters <- factor(
    dat[[celltype_col]]
  )

  samples <- factor(
    dat$mouse_id,
    levels = sample_levels
  )

  prop_list <-
    speckle::getTransformedProps(
      clusters = clusters,
      sample = samples,
      transform = "logit"
    )

  for (
    nm in c(
      "Counts",
      "Proportions",
      "TransformedProps"
    )
  ) {

    prop_list[[nm]] <-
      prop_list[[nm]][
        ,
        sample_levels,
        drop = FALSE
      ]
  }

  if (!is.null(keep_types)) {

    keep_types <- intersect(
      keep_types,
      rownames(
        prop_list$Proportions
      )
    )

    for (
      nm in c(
        "Counts",
        "Proportions",
        "TransformedProps"
      )
    ) {

      prop_list[[nm]] <-
        prop_list[[nm]][
          keep_types,
          ,
          drop = FALSE
        ]
    }
  }

  prop_list
}

run_propeller_two_group <- function(
    prop_list,
    sample_info,
    group_col,
    reference_group,
    test_group,
    batch_col = "sequence_batch"
) {

  group <- factor(
    sample_info[[group_col]],
    levels = c(
      reference_group,
      test_group
    )
  )

  stopifnot(
    !anyNA(group)
  )

  group_design <- model.matrix(
    ~ 0 + group
  )

  colnames(
    group_design
  ) <- c(
    reference_group,
    test_group
  )

  batch <- factor(
    sample_info[[batch_col]]
  )

  if (nlevels(batch) > 1) {

    batch_design <- model.matrix(
      ~ batch
    )[
      ,
      -1,
      drop = FALSE
    ]

    colnames(
      batch_design
    ) <- paste0(
      "batch_",
      levels(batch)[-1]
    )

    design <- cbind(
      group_design,
      batch_design
    )

  } else {

    design <- group_design
  }

  rownames(design) <-
    sample_info$mouse_id

  stopifnot(
    identical(
      colnames(
        prop_list$Proportions
      ),
      rownames(design)
    )
  )

  contrast <- rep(
    0,
    ncol(design)
  )

  names(contrast) <-
    colnames(design)

  contrast[
    reference_group
  ] <- -1

  contrast[
    test_group
  ] <- 1

  res <- speckle::propeller.ttest(
    prop.list = prop_list,
    design = design,
    contrasts = contrast,
    robust = TRUE,
    trend = FALSE,
    sort = FALSE
  )

  res <- as.data.frame(
    res
  ) %>%
    tibble::rownames_to_column(
      "cell_type"
    )

  pmat <-
    prop_list$Proportions

  reference_mice <-
    sample_info$mouse_id[
      group == reference_group
    ]

  test_mice <-
    sample_info$mouse_id[
      group == test_group
    ]

  descriptive <- tibble::tibble(
    cell_type =
      rownames(pmat),

    mean_reference_prop =
      rowMeans(
        pmat[
          ,
          reference_mice,
          drop = FALSE
        ]
      ),

    mean_test_prop =
      rowMeans(
        pmat[
          ,
          test_mice,
          drop = FALSE
        ]
      ),

    median_reference_prop =
      apply(
        pmat[
          ,
          reference_mice,
          drop = FALSE
        ],
        1,
        median
      ),

    median_test_prop =
      apply(
        pmat[
          ,
          test_mice,
          drop = FALSE
        ],
        1,
        median
      )
  ) %>%

    dplyr::mutate(
      difference_pp =
        100 *
        (
          mean_test_prop -
          mean_reference_prop
        ),

      log2_ratio =
        log2(
          (
            mean_test_prop +
              1e-8
          ) /
          (
            mean_reference_prop +
              1e-8
          )
        ),

      reference_group =
        reference_group,

      test_group =
        test_group
    )

 res %>%
  dplyr::left_join(
    descriptive,
    by = "cell_type"
  ) %>%
  tibble::as_tibble()
}

#Intact vs. Impaired comparison
aged_info <- sample_info %>%
  dplyr::filter(
    primary_group %in%
      c(
        "intact",
        "impaired"
      )
  ) %>%

  dplyr::arrange(
    mouse_id
  )

table(
  aged_info$primary_group
)

table(
  aged_info$primary_group,
  aged_info$sequence_batch
)

aged_cells <- cell_meta %>%
  dplyr::filter(
    mouse_id %in%
      aged_info$mouse_id
  )

broad_prop_aged <- make_prop_list(
  dat = aged_cells,
  celltype_col = "broad_lineage",
  sample_info = aged_info
)

broad_da_aged <- run_propeller_two_group(
  prop_list = broad_prop_aged,
  sample_info = aged_info,
  group_col = "primary_group",
  reference_group = "intact",
  test_group = "impaired"
)

broad_da_aged %>%
  dplyr::arrange(FDR) %>%
  tibble::as_tibble() %>%
  print(n = Inf)

readr::write_csv(
  broad_da_aged,
  file.path(
    out_dir,
    "05_baiermouse_alllineages_cognition_differential_abundance.csv"
  )
)

# Exact original three-group abundance broad-lineage design: all 13 mice, two contrasts,
# logit-transformed proportions, with sequencing batch as covariate.
count_mat <- tidyr::expand_grid(
  broad_lineage = sort(unique(as.character(cell_meta$broad_lineage))),
  mouse_id = sample_info$mouse_id
) %>%
  dplyr::left_join(broad_counts %>% dplyr::select(broad_lineage, mouse_id, n_cells),
                   by = c("broad_lineage", "mouse_id")) %>%
  dplyr::mutate(n_cells = tidyr::replace_na(n_cells, 0L)) %>%
  tidyr::pivot_wider(names_from = mouse_id, values_from = n_cells)
raw_counts <- as.matrix(count_mat[, -1, drop = FALSE])
rownames(raw_counts) <- count_mat$broad_lineage
storage.mode(raw_counts) <- "numeric"
all_info <- sample_info[match(colnames(raw_counts), sample_info$mouse_id), ]
all_info$group <- factor(as.character(all_info$primary_group),
                         levels = c("young", "intact", "impaired"))
all_info$batch <- factor(all_info$sequence_batch)
all_design <- model.matrix(~ 0 + group + batch, data = all_info)
colnames(all_design) <- sub("^group", "", colnames(all_design))
all_props <- speckle::convertDataToList(raw_counts, data.type = "counts",
                                        transform = "logit")
contrast <- matrix(0, nrow = ncol(all_design), ncol = 2,
                   dimnames = list(colnames(all_design),
                                   c("Impaired_vs_Intact", "AgedIntact_vs_Young")))
contrast["impaired", "Impaired_vs_Intact"] <- 1
contrast["intact", "Impaired_vs_Intact"] <- -1
contrast["intact", "AgedIntact_vs_Young"] <- 1
contrast["young", "AgedIntact_vs_Young"] <- -1
broad_logit_three_group <- data.table::rbindlist(
  lapply(colnames(contrast), function(comp) {
    z <- speckle::propeller.ttest(all_props, all_design,
      contrasts = contrast[, comp], robust = TRUE, trend = FALSE, sort = FALSE)
    data.table::data.table(contrast = comp, lineage = rownames(z), z)
  }), fill = TRUE)
data.table::fwrite(broad_logit_three_group,
  file.path(out_dir, "05_baiermouse_alllineages_threegroups_logit_propeller.csv"))

primary_types <- eligibility %>%
  dplyr::filter(
    statistical_tier ==
      "Primary"
  ) %>%
  dplyr::pull(
    allen_subclass
  )

primary_types

length(
  primary_types
)

primary_sub_prop_aged <- make_prop_list(
  dat = aged_cells,
  celltype_col = "allen_subclass",
  sample_info = aged_info,
  keep_types = primary_types
)

primary_subclass_da <- run_propeller_two_group(
  prop_list = primary_sub_prop_aged,
  sample_info = aged_info,
  group_col = "primary_group",
  reference_group = "intact",
  test_group = "impaired"
)

primary_subclass_da %>%
  dplyr::arrange(
    FDR
  ) %>%
  print(
    n = Inf
  )

readr::write_csv(
  primary_subclass_da,
  file.path(
    out_dir,
    "05_baiermouse_subclasses_cognition_differential_abundance.csv"
  )
)

broad_da_aged %>%
  dplyr::select(
    cell_type,
    PropMean.intact,
    PropMean.impaired,
    PropRatio,
    difference_pp,
    Tstatistic,
    P.Value,
    FDR
  ) %>%
  dplyr::arrange(FDR, P.Value) %>%
  print(n = Inf)

primary_types <- eligibility %>%
  dplyr::filter(
    statistical_tier == "Primary"
  ) %>%
  dplyr::pull(
    allen_subclass
  )

primary_types
length(primary_types)

primary_sub_prop_aged <- make_prop_list(
  dat = aged_cells,
  celltype_col = "allen_subclass",
  sample_info = aged_info,
  keep_types = primary_types
)

dim(
  primary_sub_prop_aged$Proportions
)

rownames(
  primary_sub_prop_aged$Proportions
)

colnames(
  primary_sub_prop_aged$Proportions
)

primary_subclass_da <- run_propeller_two_group(
  prop_list = primary_sub_prop_aged,
  sample_info = aged_info,
  group_col = "primary_group",
  reference_group = "intact",
  test_group = "impaired"
)

primary_subclass_da %>%
  dplyr::select(
    cell_type,
    PropMean.intact,
    PropMean.impaired,
    PropRatio,
    difference_pp,
    Tstatistic,
    P.Value,
    FDR
  ) %>%
  dplyr::arrange(FDR, P.Value) %>%
  print(n = Inf)

secondary_types <- eligibility %>%
  dplyr::filter(
    statistical_tier == "Secondary"
  ) %>%
  dplyr::pull(
    allen_subclass
  )

secondary_types

secondary_sub_prop_aged <- make_prop_list(
  dat = aged_cells,
  celltype_col = "allen_subclass",
  sample_info = aged_info,
  keep_types = secondary_types
)

secondary_subclass_da <- run_propeller_two_group(
  prop_list = secondary_sub_prop_aged,
  sample_info = aged_info,
  group_col = "primary_group",
  reference_group = "intact",
  test_group = "impaired"
)

secondary_subclass_da %>%
  dplyr::select(
    cell_type,
    PropMean.intact,
    PropMean.impaired,
    PropRatio,
    difference_pp,
    Tstatistic,
    P.Value,
    FDR
  ) %>%
  dplyr::arrange(FDR, P.Value) %>%
  print(n = Inf, width = Inf)

#Sensitivity analysis with ventricular populations
aged_nonvent <- aged_cells %>%
  dplyr::filter(
    anatomic_scope != "Adjacent ventricular"
  )

ventricular_denominator_check <- aged_cells %>%
  dplyr::count(
    mouse_id,
    name = "all_cells"
  ) %>%

  dplyr::left_join(
    aged_nonvent %>%
      dplyr::count(
        mouse_id,
        name = "nonvent_cells"
      ),
    by = "mouse_id"
  ) %>%

  dplyr::mutate(
    ventricular_cells =
      all_cells - nonvent_cells,

    pct_ventricular =
      100 * ventricular_cells / all_cells
  ) %>%

  dplyr::left_join(
    aged_info,
    by = "mouse_id"
  ) %>%

  dplyr::arrange(
    primary_group,
    pct_ventricular
  ) %>%

  tibble::as_tibble()

ventricular_denominator_check <- ventricular_denominator_check %>%
  tibble::as_tibble()

ventricular_denominator_check %>%
  print(n = Inf, width = Inf)

primary_nonvent_types <- setdiff(
  primary_types,
  c(
    "CHOR NN",
    "Ependymal NN"
  )
)

primary_sub_prop_nonvent <- make_prop_list(
  dat = aged_nonvent,
  celltype_col = "allen_subclass",
  sample_info = aged_info,
  keep_types = primary_nonvent_types
)

primary_subclass_da_nonvent <- run_propeller_two_group(
  prop_list = primary_sub_prop_nonvent,
  sample_info = aged_info,
  group_col = "primary_group",
  reference_group = "intact",
  test_group = "impaired"
)

primary_subclass_da_nonvent %>%
  dplyr::select(
    cell_type,
    PropMean.intact,
    PropMean.impaired,
    PropRatio,
    difference_pp,
    Tstatistic,
    P.Value,
    FDR
  ) %>%
  dplyr::arrange(FDR, P.Value) %>%
  print(n = Inf, width = Inf)

denominator_sensitivity <- primary_subclass_da %>%
  dplyr::filter(
    !cell_type %in%
      c(
        "CHOR NN",
        "Ependymal NN"
      )
  ) %>%
  dplyr::select(
    cell_type,
    full_difference_pp = difference_pp,
    full_P = P.Value,
    full_FDR = FDR
  ) %>%
  dplyr::inner_join(
    primary_subclass_da_nonvent %>%
      dplyr::select(
        cell_type,
        nonvent_difference_pp = difference_pp,
        nonvent_P = P.Value,
        nonvent_FDR = FDR
      ),
    by = "cell_type"
  ) %>%
  dplyr::mutate(
    same_direction =
      sign(full_difference_pp) ==
      sign(nonvent_difference_pp),

    effect_ratio =
      nonvent_difference_pp /
      full_difference_pp
  )

denominator_sensitivity %>%
  print(n = Inf, width = Inf)

#Sensitivity analysis of impaired performer that reached criterion
extreme_info <- aged_info %>%
  dplyr::filter(
    mouse_id != "C291"
  )

extreme_cells <- aged_cells %>%
  dplyr::filter(
    mouse_id != "C291"
  )

table(
  extreme_info$primary_group
)

extreme_prop <- make_prop_list(
  dat = extreme_cells,
  celltype_col = "allen_subclass",
  sample_info = extreme_info,
  keep_types = primary_types
)

primary_subclass_da_extreme <- run_propeller_two_group(
  prop_list = extreme_prop,
  sample_info = extreme_info,
  group_col = "primary_group",
  reference_group = "intact",
  test_group = "impaired"
)

primary_subclass_da_extreme <- tibble::as_tibble(
  primary_subclass_da_extreme
)

primary_subclass_da_extreme %>%
  dplyr::select(
    cell_type,
    PropMean.intact,
    PropMean.impaired,
    PropRatio,
    difference_pp,
    Tstatistic,
    P.Value,
    FDR
  ) %>%
  dplyr::arrange(FDR, P.Value) %>%
  print(n = Inf, width = Inf)

extreme_sensitivity <- primary_subclass_da %>%
  dplyr::select(
    cell_type,
    primary_difference_pp = difference_pp,
    primary_P = P.Value,
    primary_FDR = FDR
  ) %>%
  dplyr::inner_join(
    primary_subclass_da_extreme %>%
      dplyr::select(
        cell_type,
        extreme_difference_pp = difference_pp,
        extreme_P = P.Value,
        extreme_FDR = FDR
      ),
    by = "cell_type"
  ) %>%
  dplyr::mutate(
    same_direction =
      sign(primary_difference_pp) ==
      sign(extreme_difference_pp)
  )

extreme_sensitivity %>%
  print(n = Inf, width = Inf)

#Cognitive flexibility sensitivity analysis
# Use the full-precision measurements generated by Script 01.
behavior <- readr::read_csv(
  file.path(base_dir, "01_mouse_behavior", "01_baiermouse_allmice_behavior_mouse_level.csv"),
  show_col_types = FALSE
)
cf_lookup <- behavior %>%
  dplyr::filter(sequenced, age_group == "Aged") %>%
  dplyr::transmute(mouse_id = as.character(mouse_id),
                   CF51to61 = as.numeric(cf51to61))

stopifnot(
  !anyDuplicated(cf_lookup$mouse_id),
  setequal(cf_lookup$mouse_id, aged_info$mouse_id),
  all(is.finite(cf_lookup$CF51to61))
)
aged_info_cf <- aged_info %>%
  dplyr::left_join(cf_lookup, by = "mouse_id")
readr::write_csv(
  aged_info_cf,
  file.path(out_dir, "05_baiermouse_agedmice_cognition_flexibility_inputs.csv")
)

run_continuous_prop <- function(
    prop_list,
    sample_info,
    predictor,
    batch_col = "sequence_batch"
) {

  predictor_z <- as.numeric(
    scale(sample_info[[predictor]])
  )

  batch <- factor(
    sample_info[[batch_col]]
  )

  design <- model.matrix(
    ~ predictor_z + batch
  )

  rownames(design) <-
    sample_info$mouse_id

  pmat <- prop_list$TransformedProps[
    ,
    sample_info$mouse_id,
    drop = FALSE
  ]

  fit <- limma::lmFit(
    pmat,
    design
  )

  fit <- limma::eBayes(
    fit,
    robust = TRUE,
    trend = FALSE
  )

  res <- limma::topTable(
    fit,
    coef = "predictor_z",
    number = Inf,
    sort.by = "none"
  ) %>%
    tibble::rownames_to_column(
      "cell_type"
    ) %>%
    dplyr::rename(
      FDR = adj.P.Val
    )

  raw_p <- prop_list$Proportions[
    ,
    sample_info$mouse_id,
    drop = FALSE
  ]

  rho <- apply(
    raw_p,
    1,
    function(x) {
      suppressWarnings(
        cor(
          x,
          sample_info[[predictor]],
          method = "spearman"
        )
      )
    }
  )

  res %>%
    dplyr::left_join(
      tibble::tibble(
        cell_type = names(rho),
        spearman_rho = as.numeric(rho)
      ),
      by = "cell_type"
    ) %>%
    tibble::as_tibble()
}

primary_subclass_cf <- run_continuous_prop(
  prop_list = primary_sub_prop_aged,
  sample_info = aged_info_cf,
  predictor = "CF51to61"
)

primary_subclass_cf %>%
  dplyr::select(
    cell_type,
    logFC,
    t,
    P.Value,
    FDR,
    spearman_rho
  ) %>%
  dplyr::arrange(FDR, P.Value) %>%
  print(n = Inf, width = Inf)

#Aging comparison
aging_info <- sample_info %>%
  dplyr::filter(
    primary_group %in%
      c(
        "young",
        "intact"
      )
  ) %>%
  dplyr::mutate(
    aging_group = dplyr::if_else(
      primary_group == "young",
      "young",
      "aged_intact"
    )
  ) %>%
  dplyr::arrange(
    mouse_id
  )

aging_cells <- cell_meta %>%
  dplyr::filter(
    mouse_id %in%
      aging_info$mouse_id
  )

table(
  aging_info$aging_group
)

primary_prop_aging <- make_prop_list(
  dat = aging_cells,
  celltype_col = "allen_subclass",
  sample_info = aging_info,
  keep_types = primary_types
)

primary_subclass_da_aging <- run_propeller_two_group(
  prop_list = primary_prop_aging,
  sample_info = aging_info,
  group_col = "aging_group",
  reference_group = "young",
  test_group = "aged_intact"
)

primary_subclass_da_aging <-
  tibble::as_tibble(
    primary_subclass_da_aging
  )

primary_subclass_da_aging %>%
  dplyr::select(
    cell_type,
    PropMean.young,
    PropMean.aged_intact,
    PropRatio,
    difference_pp,
    Tstatistic,
    P.Value,
    FDR
  ) %>%
  dplyr::arrange(FDR, P.Value) %>%
  print(n = Inf, width = Inf)

#Verify young coverage of subclasses
young_coverage <- subclass_counts %>%
  dplyr::filter(
    primary_group == "young",
    allen_subclass %in% primary_types
  ) %>%
  tidyr::complete(
    mouse_id = sample_info$mouse_id[
      sample_info$primary_group == "young"
    ],
    allen_subclass = primary_types,
    fill = list(n_cells = 0)
  ) %>%
  dplyr::group_by(allen_subclass) %>%
  dplyr::summarise(
    n_young_mice_present = sum(n_cells > 0),
    n_young_mice_ge50 = sum(n_cells >= 50),
    n_young_mice_ge100 = sum(n_cells >= 100),
    min_young_cells = min(n_cells),
    median_young_cells = median(n_cells),
    max_young_cells = max(n_cells),
    .groups = "drop"
  ) %>%
  dplyr::arrange(min_young_cells)

young_coverage %>%
  print(n = Inf, width = Inf)

#Leave-one-mouse-out sensitivity analysis
lomo_results <- purrr::map_dfr(
  aged_info$mouse_id,
  function(drop_mouse) {

    info_i <- aged_info %>%
      dplyr::filter(
        mouse_id != drop_mouse
      )

    cells_i <- aged_cells %>%
      dplyr::filter(
        mouse_id != drop_mouse
      )

    prop_i <- make_prop_list(
      dat = cells_i,
      celltype_col = "allen_subclass",
      sample_info = info_i,
      keep_types = primary_types
    )

    res_i <- run_propeller_two_group(
      prop_list = prop_i,
      sample_info = info_i,
      group_col = "primary_group",
      reference_group = "intact",
      test_group = "impaired"
    )

    res_i %>%
      dplyr::select(
        cell_type,
        difference_pp,
        P.Value,
        FDR
      ) %>%
      dplyr::mutate(
        dropped_mouse = drop_mouse
      )
  }
)

lomo_summary <- lomo_results %>%
  dplyr::left_join(
    primary_subclass_da %>%
      dplyr::select(
        cell_type,
        full_difference_pp = difference_pp
      ),
    by = "cell_type"
  ) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    full_difference_pp =
      dplyr::first(full_difference_pp),

    n_runs =
      dplyr::n(),

    direction_stability =
      mean(
        sign(difference_pp) ==
          sign(full_difference_pp)
      ),

    min_difference_pp =
      min(difference_pp),

    median_difference_pp =
      median(difference_pp),

    max_difference_pp =
      max(difference_pp),

    min_nominal_P =
      min(P.Value),

    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(direction_stability),
    cell_type
  )

lomo_summary %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

readr::write_csv(
  lomo_results,
  file.path(
    out_dir,
    "05_baiermouse_subclasses_cognition_leave_one_mouse_out_runs.csv"
  )
)

readr::write_csv(
  lomo_summary,
  file.path(
    out_dir,
    "05_baiermouse_subclasses_cognition_leave_one_mouse_out_summary.csv"
  )
)

da_robustness <- primary_subclass_da %>%
  dplyr::select(
    cell_type,
    primary_diff_pp = difference_pp,
    primary_P = P.Value,
    primary_FDR = FDR
  ) %>%

  dplyr::left_join(
    primary_subclass_da_extreme %>%
      dplyr::select(
        cell_type,
        extreme_diff_pp = difference_pp,
        extreme_P = P.Value,
        extreme_FDR = FDR
      ),
    by = "cell_type"
  ) %>%

  dplyr::left_join(
    primary_subclass_cf %>%
      dplyr::select(
        cell_type,
        CF_logFC = logFC,
        CF_P = P.Value,
        CF_FDR = FDR,
        CF_spearman_rho = spearman_rho
      ),
    by = "cell_type"
  ) %>%

  dplyr::left_join(
    lomo_summary %>%
      dplyr::select(
        cell_type,
        LOMO_direction_stability =
          direction_stability
      ),
    by = "cell_type"
  ) %>%

  dplyr::mutate(
    extreme_same_direction =
      sign(primary_diff_pp) ==
      sign(extreme_diff_pp),

    CF_direction_concordant =
      sign(primary_diff_pp) ==
      -sign(CF_spearman_rho)
  )

da_robustness %>%
  print(n = Inf, width = Inf)

# Export the supporting abundance and sensitivity results already computed above.
# These tables are used for audit and figure provenance, not new analyses.
readr::write_csv(secondary_subclass_da,
  file.path(out_dir, "05_baiermouse_secondary_subclasses_cognition_abundance.csv"))
readr::write_csv(ventricular_denominator_check,
  file.path(out_dir, "05_baiermouse_alllineages_cognition_ventricular_denominator_check.csv"))
readr::write_csv(primary_subclass_da_nonvent,
  file.path(out_dir, "05_baiermouse_subclasses_cognition_nonventricular_abundance.csv"))
readr::write_csv(denominator_sensitivity,
  file.path(out_dir, "05_baiermouse_subclasses_cognition_ventricular_denominator_sensitivity.csv"))
readr::write_csv(primary_subclass_da_extreme,
  file.path(out_dir, "05_baiermouse_subclasses_cognition_exclude_c291_abundance.csv"))
readr::write_csv(extreme_sensitivity,
  file.path(out_dir, "05_baiermouse_subclasses_cognition_exclude_c291_summary.csv"))
readr::write_csv(primary_subclass_cf,
  file.path(out_dir, "05_baiermouse_subclasses_cognition_continuous_cf_abundance.csv"))
readr::write_csv(primary_subclass_da_aging,
  file.path(out_dir, "05_baiermouse_subclasses_aging_differential_abundance.csv"))
readr::write_csv(young_coverage,
  file.path(out_dir, "05_baiermouse_subclasses_aging_young_coverage.csv"))
readr::write_csv(da_robustness,
  file.path(out_dir, "05_baiermouse_subclasses_cognition_robustness_summary.csv"))
