#Baier et. al
#04a: Reference-based annotation via MapMyCells (initial preparation)

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
base_dir <- hca_output_root
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c("04_annotation", "objects")) dir.create(file.path(base_dir,d), recursive=TRUE, showWarnings=FALSE)

#Load libraries
library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(Seurat)
library(patchwork)
library(S4Vectors)
library(zellkonverter)

s3 <- readRDS(
  file.path(
    base_dir,
    "objects",
    "03_baiermouse_allcells_final_qc_checkpoint.rds"
  )
)

meta <- s3$meta
final_qc_summary <- s3$final_qc_summary

final_sce_dir <- file.path(
  base_dir,
  "03_final_qc",
  "scDblFinder",
  "objects"
)

meta <- meta %>%
  mutate(
    final_sce_file = file.path(
      final_sce_dir,
      paste0("03_baiermouse_", mouse_id, "_final_singlets.rds")
    )
  )

stopifnot(
  all(file.exists(meta$final_sce_file))
)

final_inventory <- map_dfr(
  seq_len(nrow(meta)),
  function(i) {

    sce <- readRDS(
      meta$final_sce_file[i]
    )

    out <- tibble(
      mouse_id = meta$mouse_id[i],
      sequence_id = meta$sequence_id[i],
      n_features = nrow(sce),
      n_cells = ncol(sce),
      n_umis = sum(counts(sce))
    )

    rm(sce)
    gc()

    out
  }
)

final_inventory
sum(final_inventory$n_cells)

write_csv(
  final_inventory,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_annotation_input_inventory.csv"
  )
)

reference_sce <- readRDS(
  meta$final_sce_file[1]
)

reference_genes <- rownames(
  reference_sce
)

rm(reference_sce)
gc()

feature_check <- map_lgl(
  meta$final_sce_file,
  function(f) {

    sce <- readRDS(f)

    same <- identical(
      rownames(sce),
      reference_genes
    )

    rm(sce)
    gc()

    same
  }
)

feature_check

stopifnot(
  all(feature_check)
)

#Create Seurat object for each sequenced mouse
seurat_final_list <- vector(
  "list",
  nrow(meta)
)

names(seurat_final_list) <-
  meta$mouse_id

for (i in seq_len(nrow(meta))) {

  mouse <- meta$mouse_id[i]

  sce <- readRDS(
    meta$final_sce_file[i]
  )

  count_mat <- counts(sce)

  cell_meta <- as.data.frame(
    colData(sce)
  )

  rownames(cell_meta) <-
    colnames(count_mat)

  #Retain original 10X barcode
  cell_meta$barcode_10x <-
    colnames(count_mat)

  obj <- CreateSeuratObject(
    counts = count_mat,
    project = mouse,
    meta.data = cell_meta,
    min.cells = 0,
    min.features = 0
  )

  #Retain globally unique cell IDs
  obj <- RenameCells(
    obj,
    add.cell.id = mouse
  )

  #Retain sample metadata
  obj$mouse_id <- mouse
  obj$sequence_id <- meta$sequence_id[i]
  obj$sequence_batch <- meta$sequence_batch[i]
  obj$age_group <- meta$age_group[i]
  obj$primary_group <- meta$primary_group[i]

  seurat_final_list[[mouse]] <- obj

  rm(
    sce,
    count_mat,
    cell_meta,
    obj
  )

  gc()
}

map_int(
  seurat_final_list,
  ncol
)

sum(
  map_int(
    seurat_final_list,
    ncol
  )
)

#Merge individual Seurat objects into working atlas
atlas_raw <- merge(
  x = seurat_final_list[[1]],
  y = seurat_final_list[-1],
  add.cell.ids = NULL,
  project = "CognitiveAging_Hippocampus",
  merge.data = FALSE
)

atlas_raw
ncol(atlas_raw)

table(
  atlas_raw$mouse_id
)

table(
  atlas_raw$age_group,
  atlas_raw$primary_group
)

#Join count layers
Layers(
  atlas_raw[["RNA"]]
)

atlas_raw <- JoinLayers(
  atlas_raw
)

Layers(
  atlas_raw[["RNA"]]
)

#Calculate QC for merged atlas
atlas_raw[["percent.mt"]] <-
  PercentageFeatureSet(
    atlas_raw,
    pattern = "^mt-"
  )

atlas_raw[["percent.ribo"]] <-
  PercentageFeatureSet(
    atlas_raw,
    pattern = "^Rp[sl]"
  )

VlnPlot(
  atlas_raw,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "mouse_id",
  pt.size = 0,
  ncol = 1
)

#Save merged raw atlas object
saveRDS(
  atlas_raw,
  file.path(
    base_dir,
    "objects",
    "04_baiermouse_allcells_merged_raw_counts.rds"
  ),
  compress = FALSE
)

#Atlas normalization, PCA, UMAP visualization
atlas <- atlas_raw

atlas <- NormalizeData(
  atlas,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

atlas <- FindVariableFeatures(
  atlas,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

technical_features <- grep(
  "^mt-|^Rp[sl]|^Hb[ab]",
  VariableFeatures(atlas),
  value = TRUE
)

VariableFeatures(atlas) <-
  setdiff(
    VariableFeatures(atlas),
    technical_features
  )

length(
  VariableFeatures(atlas)
)

atlas <- ScaleData(
  atlas,
  features = VariableFeatures(atlas),
  verbose = FALSE
)

atlas <- RunPCA(
  atlas,
  features = VariableFeatures(atlas),
  npcs = 50,
  verbose = FALSE
)

ElbowPlot(
  atlas,
  ndims = 50
)

atlas <- FindNeighbors(
  atlas,
  dims = 1:30,
  verbose = FALSE
)

atlas <- FindClusters(
  atlas,
  resolution = 0.4,
  random.seed = 1234,
  verbose = FALSE
)

atlas <- RunUMAP(
  atlas,
  dims = 1:30,
  seed.use = 1234,
  verbose = FALSE
)

DimPlot(
  atlas,
  reduction = "umap",
  group.by = "mouse_id",
  raster = TRUE
)

DimPlot(
  atlas,
  reduction = "umap",
  group.by = "sequence_batch",
  raster = TRUE
)

DimPlot(
  atlas,
  reduction = "umap",
  group.by = "primary_group",
  raster = TRUE
)

DimPlot(
  atlas,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
)

DimPlot(
  atlas,
  reduction = "pca",
  group.by = "sequence_batch"
)

canonical_markers <- c(
  "Snap25",
  "Syt1",

  "Aldoc",
  "Slc1a3",
  "Aqp4",

  "P2ry12",
  "Tmem119",
  "C1qa",

  "Mog",
  "Mbp",

  "Pdgfra",
  "Cspg4",

  "Cldn5",
  "Pecam1",

  "Pdgfrb",
  "Rgs5"
)

DotPlot(
  atlas,
  features = canonical_markers
) +
  RotatedAxis()

saveRDS(
  atlas,
  file.path(
    base_dir,
    "objects",
    "04_baiermouse_allcells_preannotation_atlas.rds"
  ),
  compress = FALSE
)

#Prepare MapMyCells input
gene_id_matrix <- Read10X(
  data.dir = meta$filtered_dir[1],
  gene.column = 1,
  unique.features = TRUE
)

gene_symbol_matrix <- Read10X(
  data.dir = meta$filtered_dir[1],
  gene.column = 2,
  unique.features = TRUE
)

nrow(gene_id_matrix)
nrow(gene_symbol_matrix)
nrow(atlas_raw)

stopifnot(
  identical(
    rownames(gene_symbol_matrix),
    rownames(atlas_raw)
  )
)

gene_map <- tibble(
  gene_symbol = rownames(
    gene_symbol_matrix
  ),
  ensembl_id = rownames(
    gene_id_matrix
  )
)

head(gene_map)

sum(
  duplicated(
    gene_map$ensembl_id
  )
)

write_csv(
  gene_map,
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_gene_symbol_ensembl_map.csv"
  )
)

final_sce_dir <- file.path(
  base_dir,
  "03_final_qc",
  "scDblFinder",
  "objects"
)

meta <- meta %>%
  mutate(
    final_sce_file = file.path(
      final_sce_dir,
      paste0("03_baiermouse_", mouse_id, "_final_singlets.rds")
    )
  )

meta %>%
  select(
    mouse_id,
    final_sce_file
  )

file_check <- meta %>%
  transmute(
    mouse_id,
    final_sce_file,
    exists = file.exists(final_sce_file)
  )

file_check

stopifnot(
  all(file_check$exists)
)

#Build the MapMyCells input from the sample SCEs
map_sce_list <- vector(
  "list",
  nrow(meta)
)

names(map_sce_list) <- meta$mouse_id

for (i in seq_len(nrow(meta))) {

  mouse <- meta$mouse_id[i]

  sce <- readRDS(
    meta$final_sce_file[i]
  )

  count_mat <- counts(sce)

  stopifnot(
    identical(
      rownames(count_mat),
      reference_genes
    )
  )

  #Confirm cell IDs are globally unique
  new_cell_ids <- paste0(
    mouse,
    "_",
    colnames(count_mat)
  )

  colnames(count_mat) <- new_cell_ids

  sce_minimal <- SingleCellExperiment(
    assays = list(
      counts = count_mat
    ),
    colData = DataFrame(
      mouse_id = rep(
        mouse,
        ncol(count_mat)
      ),

      sequence_id = rep(
        as.character(meta$sequence_id[i]),
        ncol(count_mat)
      ),

      sequence_batch = rep(
        meta$sequence_batch[i],
        ncol(count_mat)
      ),

      age_group = rep(
        meta$age_group[i],
        ncol(count_mat)
      ),

      primary_group = rep(
        meta$primary_group[i],
        ncol(count_mat)
      ),

      barcode_10x = colnames(sce)
    )
  )

  colnames(sce_minimal) <- new_cell_ids

  map_sce_list[[mouse]] <- sce_minimal

  rm(
    sce,
    count_mat,
    sce_minimal
  )

  gc()
}

purrr::map_int(
  map_sce_list,
  ncol
)

sum(
  purrr::map_int(
    map_sce_list,
    ncol
  )
)

feature_identical <- purrr::map_lgl(
  map_sce_list,
  function(x) {
    identical(
      rownames(x),
      reference_genes
    )
  }
)

feature_identical
stopifnot(all(feature_identical))

map_sce <- do.call(
  cbind,
  map_sce_list
)

dim(map_sce)

table(
  colData(map_sce)$mouse_id
)

#Ensembl gene ID conversion
head(gene_map)

nrow(gene_map)

stopifnot(
  identical(
    rownames(map_sce),
    gene_map$gene_symbol
  )
)

rowData(map_sce)$gene_symbol <-
  rownames(map_sce)

rownames(map_sce) <-
  gene_map$ensembl_id

head(
  rownames(map_sce)
)

sum(
  duplicated(
    rownames(map_sce)
  )
)

map_export <- SingleCellExperiment(
  assays = list(
    counts = counts(map_sce)
  ),

  colData = DataFrame(
    mouse_id =
      colData(map_sce)$mouse_id,

    sequence_id =
      colData(map_sce)$sequence_id,

    sequence_batch =
      colData(map_sce)$sequence_batch,

    age_group =
      colData(map_sce)$age_group,

    primary_group =
      colData(map_sce)$primary_group,

    barcode_10x =
      colData(map_sce)$barcode_10x
  ),

  rowData = DataFrame(
    gene_symbol =
      rowData(map_sce)$gene_symbol
  )
)

dim(map_export)

assayNames(map_export)

class(
  counts(map_export)
)

dim(
  counts(map_export)
)

mapmycells_dir <- file.path(
  base_dir,
  "04_annotation",
  "MapMyCells"
)

dir.create(
  mapmycells_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

writeH5AD(
  map_export,
  file = file.path(
    mapmycells_dir,
    "04_baiermouse_allcells_mapmycells_input_counts.h5ad"
  ),
  X_name = "counts"
)

annotation_checkpoint <- list(
  meta = meta,
  final_inventory = final_inventory,
  gene_map = gene_map,
  n_final_cells = ncol(atlas_raw),
  n_features = nrow(atlas_raw),
  mapmycells_h5ad = file.path(
    base_dir,
    "04_annotation",
    "MapMyCells",
    "04_baiermouse_allcells_mapmycells_input_counts.h5ad"
  )
)

saveRDS(
  annotation_checkpoint,
  file.path(
    base_dir,
    "objects",
    "04_baiermouse_allcells_mapmycells_preparation_checkpoint.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    base_dir,
    "04_annotation",
    "04_baiermouse_allcells_mapmycells_preparation_sessioninfo.txt"
  )
)