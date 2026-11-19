#Baier et. al
#04b: Reference-based annotation via MapMyCells

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
base_dir <- hca_output_root
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(tidyverse)
library(Seurat)
library(Matrix)

precheck <- readRDS(file.path(base_dir, "objects", "04_baiermouse_allcells_mapmycells_preparation_checkpoint.rds"))
meta <- precheck$meta
atlas <- readRDS(file.path(base_dir, "objects", "04_baiermouse_allcells_preannotation_atlas.rds"))

code_dir <- Sys.getenv("HCA_CODE_ROOT", unset = getwd())
python <- Sys.getenv("HCA_MAPMYCELLS_PYTHON", unset = "python")
status <- system2(
  python,
  args = c(shQuote(file.path(code_dir, "04_mapmycells.py")), shQuote(base_dir))
)
if (status != 0L) stop("MapMyCells mapping failed; see pipeline_results/04_baiermouse_allcells_mapmycells_run.log")

full_map_path <- file.path(base_dir, "04_annotation", "MapMyCells", "pipeline_results",
                           "04_baiermouse_allcells_mapmycells_assignments.csv")
stopifnot(ncol(atlas) == precheck$n_final_cells)

#Integrate MapMyCells result to Seurat object
full_map <- readr::read_csv(
  full_map_path,
  comment = "#",
  show_col_types = FALSE
)

dim(full_map)

stopifnot(
  nrow(full_map) == 198122,
  dplyr::n_distinct(full_map$cell_id) == 198122
)

all(
  full_map$cell_id %in%
    colnames(atlas)
)

setequal(
  full_map$cell_id,
  colnames(atlas)
)

full_class_counts <- full_map %>%
  dplyr::count(
    class_name,
    sort = TRUE
  )

full_subclass_counts <- full_map %>%
  dplyr::count(
    subclass_name,
    sort = TRUE
  )

full_cluster_counts <- full_map %>%
  dplyr::count(
    cluster_name,
    sort = TRUE
  )

full_class_counts
full_subclass_counts %>% print(n = Inf)

full_confidence <- tibble::tibble(
  level = c(
    "class",
    "subclass",
    "supertype",
    "cluster"
  ),

  median = c(
    median(full_map$class_bootstrapping_probability, na.rm = TRUE),
    median(full_map$subclass_bootstrapping_probability, na.rm = TRUE),
    median(full_map$supertype_bootstrapping_probability, na.rm = TRUE),
    median(full_map$cluster_bootstrapping_probability, na.rm = TRUE)
  ),

  p10 = c(
    quantile(full_map$class_bootstrapping_probability, 0.10, na.rm = TRUE),
    quantile(full_map$subclass_bootstrapping_probability, 0.10, na.rm = TRUE),
    quantile(full_map$supertype_bootstrapping_probability, 0.10, na.rm = TRUE),
    quantile(full_map$cluster_bootstrapping_probability, 0.10, na.rm = TRUE)
  ),

  pct_ge_0.8 = c(
    100 * mean(full_map$class_bootstrapping_probability >= 0.8, na.rm = TRUE),
    100 * mean(full_map$subclass_bootstrapping_probability >= 0.8, na.rm = TRUE),
    100 * mean(full_map$supertype_bootstrapping_probability >= 0.8, na.rm = TRUE),
    100 * mean(full_map$cluster_bootstrapping_probability >= 0.8, na.rm = TRUE)
  ),

  pct_ge_0.9 = c(
    100 * mean(full_map$class_bootstrapping_probability >= 0.9, na.rm = TRUE),
    100 * mean(full_map$subclass_bootstrapping_probability >= 0.9, na.rm = TRUE),
    100 * mean(full_map$supertype_bootstrapping_probability >= 0.9, na.rm = TRUE),
    100 * mean(full_map$cluster_bootstrapping_probability >= 0.9, na.rm = TRUE)
  )
)

full_confidence

#Save raw Allen output
write_csv(
  full_class_counts,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_allen_class_counts.csv"
  )
)

write_csv(
  full_subclass_counts,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_allen_subclass_counts.csv"
  )
)

write_csv(
  full_confidence,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_allen_mapping_confidence.csv"
  )
)

saveRDS(
  full_map,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_allen_mapping_assignments.rds"
  )
)

#Clean Allen names
stopifnot(
  nrow(full_map) == ncol(atlas),
  dplyr::n_distinct(full_map$cell_id) == nrow(full_map),
  setequal(full_map$cell_id, colnames(atlas))
)

full_map <- full_map %>%
  dplyr::mutate(
    allen_class = stringr::str_remove(
      class_name,
      "^\\d+\\s+"
    ),

    allen_subclass = stringr::str_remove(
      subclass_name,
      "^\\d+\\s+"
    ),

    allen_supertype = stringr::str_remove(
      supertype_name,
      "^\\d+\\s+"
    ),

    allen_cluster = stringr::str_remove(
      cluster_name,
      "^\\d+\\s+"
    )
  )

#Join by cell ID
allen_idx <- match(
  colnames(atlas),
  full_map$cell_id
)

stopifnot(
  !anyNA(allen_idx)
)

allen_ordered <- full_map[
  allen_idx,
]

stopifnot(
  identical(
    allen_ordered$cell_id,
    colnames(atlas)
  )
)

atlas$allen_class_label <-
  allen_ordered$class_label

atlas$allen_class <-
  allen_ordered$allen_class

atlas$allen_class_prob <-
  allen_ordered$class_bootstrapping_probability

atlas$allen_subclass_label <-
  allen_ordered$subclass_label

atlas$allen_subclass <-
  allen_ordered$allen_subclass

atlas$allen_subclass_prob <-
  allen_ordered$subclass_bootstrapping_probability

atlas$allen_supertype_label <-
  allen_ordered$supertype_label

atlas$allen_supertype <-
  allen_ordered$allen_supertype

atlas$allen_supertype_prob <-
  allen_ordered$supertype_bootstrapping_probability

atlas$allen_cluster_label <-
  allen_ordered$cluster_label

atlas$allen_cluster <-
  allen_ordered$allen_cluster

atlas$allen_cluster_prob <-
  allen_ordered$cluster_bootstrapping_probability

table(
  atlas$allen_class,
  useNA = "ifany"
)

#Reference validation
atlas_meta <- atlas[[]] %>%
  tibble::rownames_to_column(
    "cell_id"
  )

cluster_allen_full <- atlas_meta %>%
  dplyr::group_by(
    seurat_clusters,
    allen_class
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    seurat_clusters
  ) %>%
  dplyr::mutate(
    pct = 100 * n / sum(n)
  ) %>%
  dplyr::ungroup()

write_csv(
  cluster_allen_full,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_seurat_allen_cluster_concordance.csv"
  )
)

#Choroidal/ependymal identity
FeaturePlot(
  atlas,
  features = c(
    "Ttr",
    "Klotho",
    "Folr1",
    "Otx2",
    "Foxj1",
    "Pifo",
    "Aqp4",
    "Aldoc"
  ),
  ncol = 4,
  raster = TRUE
)

DotPlot(
  atlas,
  features = c(
    "Ttr",
    "Klotho",
    "Folr1",
    "Otx2",
    "Foxj1",
    "Pifo",
    "Aqp4",
    "Aldoc"
  ),
  group.by = "allen_subclass"
) +
  RotatedAxis()

ventricular_cells <- WhichCells(
  atlas,
  expression =
    allen_subclass %in% c(
      "CHOR NN",
      "Ependymal NN",
      "Astroependymal NN",
      "Astro-TE NN",
      "Astro-NT NN"
    )
)

atlas_ventricular <- subset(
  atlas,
  cells = ventricular_cells
)

DotPlot(
  atlas_ventricular,
  features = c(
    "Ttr",
    "Klotho",
    "Folr1",
    "Otx2",
    "Foxj1",
    "Pifo",
    "Aqp4",
    "Aldoc"
  ),
  group.by = "allen_subclass"
) +
  RotatedAxis()

#Define broad cell lineages
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
) ~ "Immature neuron",

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

sort(
  table(
    atlas$broad_lineage
  ),
  decreasing = TRUE
)

#Cell lineage representation by sequenced mouse
lineage_by_mouse <- atlas[[]] %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  dplyr::count(
    mouse_id,
    primary_group,
    broad_lineage
  )

lineage_by_mouse %>%
  tidyr::pivot_wider(
    names_from = mouse_id,
    values_from = n,
    values_fill = 0
  )

subclass_by_mouse <- atlas[[]] %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  dplyr::count(
    mouse_id,
    primary_group,
    allen_subclass
  )

write_csv(
  subclass_by_mouse,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_allen_subclass_counts_by_mouse.csv"
  )
)

#Save mapped checkpoint
saveRDS(
  atlas,
  file.path(
    base_dir,
    "objects",
    "04_baiermouse_allcells_allen_annotated_atlas_intermediate.rds"
  ),
  compress = FALSE
)

atlas$analysis_celltype <- dplyr::case_when(

  atlas$allen_subclass == "Microglia NN" ~
    "Microglia",

  atlas$allen_subclass == "Oligo NN" ~
    "Oligodendrocyte",

  atlas$allen_subclass == "OPC NN" ~
    "OPC",

  atlas$allen_subclass == "Astro-TE NN" ~
    "Astrocyte",

  atlas$allen_subclass == "Endo NN" ~
    "Endothelial",

  atlas$allen_subclass == "Peri NN" ~
    "Pericyte",

  atlas$allen_subclass == "SMC NN" ~
    "Smooth muscle cell",

  atlas$allen_subclass == "VLMC NN" ~
    "VLMC",

  atlas$allen_subclass == "BAM NN" ~
    "BAM",

  atlas$allen_subclass == "Lymphoid NN" ~
    "Lymphoid",

  atlas$allen_subclass == "DC NN" ~
    "Dendritic cell",

  atlas$allen_subclass == "Monocytes NN" ~
    "Monocyte",

  atlas$allen_subclass == "DG Glut" ~
    "DG glutamatergic neuron",

  atlas$allen_subclass == "HPF CR Glut" ~
    "HPF CR glutamatergic neuron",

  atlas$allen_subclass == "CA1-ProS Glut" ~
    "CA1-ProS glutamatergic neuron",

  atlas$allen_subclass == "CA3 Glut" ~
    "CA3 glutamatergic neuron",

  atlas$allen_subclass == "CHOR NN" ~
    "Choroid plexus",

  atlas$allen_subclass == "Ependymal NN" ~
    "Ependymal",

  atlas$allen_subclass == "Astroependymal NN" ~
    "Astroependymal",

  stringr::str_detect(
    atlas$allen_subclass,
    "IMN"
  ) ~
    "Immature neuron (reference)",

  TRUE ~
    "Other / descriptive"
)

aged_meta <- meta %>%
  dplyr::filter(
    primary_group %in% c(
      "intact",
      "impaired"
    )
  )

aged_mice <- aged_meta$mouse_id

all_subclasses <- sort(
  unique(atlas$allen_subclass)
)

aged_subclass_counts <- tidyr::crossing(
  mouse_id = aged_mice,
  allen_subclass = all_subclasses
) %>%

  dplyr::left_join(
    subclass_by_mouse %>%
      dplyr::select(
        mouse_id,
        allen_subclass,
        n
      ),
    by = c(
      "mouse_id",
      "allen_subclass"
    )
  ) %>%

  dplyr::mutate(
    n = tidyr::replace_na(
      n,
      0L
    )
  ) %>%

  dplyr::left_join(
    aged_meta %>%
      dplyr::select(
        mouse_id,
        primary_group
      ),
    by = "mouse_id"
  )

subclass_eligibility <- aged_subclass_counts %>%
  dplyr::group_by(
    allen_subclass
  ) %>%

  dplyr::summarise(
    total_aged_cells = sum(n),

    n_aged_mice_present =
      sum(n > 0),

    n_aged_mice_ge100 =
      sum(n >= 100),

    n_aged_mice_ge50 =
      sum(n >= 50),

    min_cells =
      min(n),

    median_cells =
      median(n),

    max_cells =
      max(n),

    .groups = "drop"
  ) %>%

  dplyr::mutate(
    statistical_tier = dplyr::case_when(

      n_aged_mice_ge100 >= 8 ~
        "Primary",

      n_aged_mice_ge50 >= 7 ~
        "Secondary",

      TRUE ~
        "Descriptive"
    )
  ) %>%

  dplyr::arrange(
    factor(
      statistical_tier,
      levels = c(
        "Primary",
        "Secondary",
        "Descriptive"
      )
    ),
    dplyr::desc(total_aged_cells)
  )

subclass_eligibility %>%
  dplyr::filter(
    statistical_tier !=
      "Descriptive"
  ) %>%
  print(n = Inf)

readr::write_csv(
  subclass_eligibility,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_allen_subclass_eligibility.csv"
  )
)

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
    "DG Glut",
    "HPF CR Glut",
    "CA1-ProS Glut",
    "CA3 Glut"
  ) ~
    "Hippocampal/parenchymal",

  TRUE ~
    "Other/reference-only"
)

#Save final annotated atlas
saveRDS(
  atlas,
  file.path(
    base_dir,
    "objects",
    "04_baiermouse_allcells_annotated_atlas.rds"
  ),
  compress = FALSE
)

annotation_final_checkpoint <- list(
  meta = meta,
  full_confidence = full_confidence,
  subclass_eligibility = subclass_eligibility,
  final_n_cells = ncol(atlas),
  allen_taxonomy = "CCN20230722",
  mapping_algorithm = "hierarchical",
  mapping_version = "cell_type_mapper 1.7.4"
)

saveRDS(
  annotation_final_checkpoint,
  file.path(
    base_dir,
    "objects",
    "04_baiermouse_allcells_annotation_checkpoint.rds"
  )
)