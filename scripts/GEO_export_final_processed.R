#Baier et al.
#GEO processed single-cell data export

#
# Exports the four standard GEO processed-data files:
#   Baier_final_matrix.mtx.gz       sparse integer UMI counts (genes x cells)
#   Baier_final_features.tsv.gz     matrix gene features (same order as rows)
#   Baier_final_barcodes.tsv.gz     exact cell IDs (same order as columns)
#   Baier_final_cell_metadata.csv   annotations indexed by those cell IDs

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

root <- hca_output_root
out <- file.path(root, "GEO_final_processed")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

atlas_file <- file.path(root, "objects", "04_baiermouse_allcells_annotated_atlas.rds")
states_file <- file.path(root, "08_astrocyte_states", "08_frozen_cell_assignments.csv")
for (p in c(atlas_file, states_file)) {
  if (!file.exists(p)) stop("Required numbered-pipeline output missing: ", p)
}

message("Reading final annotated atlas")
atlas <- readRDS(atlas_file)
if (!inherits(atlas, "Seurat") || !"RNA" %in% names(atlas@assays))
  stop("Expected final Seurat atlas with an RNA assay")

# A single joined counts layer is necessary to avoid exporting only part of
# a Seurat v5 object. Script 04's RNA counts are the SoupX-corrected counts.
count_layers <- SeuratObject::Layers(atlas[["RNA"]], search = "^counts(\\.|$)")
if (!identical(count_layers, "counts"))
  stop("Expected a single joined RNA/counts layer, found: ",
       paste(count_layers, collapse = ", "))
counts <- SeuratObject::LayerData(atlas, assay = "RNA", layer = "counts")
if (!inherits(counts, "sparseMatrix"))
  stop("RNA/counts is not a sparse matrix")
counts <- as(counts, "dgCMatrix")

metadata_df <- atlas[[]]
if ("cell_id" %in% names(metadata_df)) {
  if (!identical(as.character(metadata_df$cell_id), rownames(metadata_df)))
    stop("Existing cell_id disagrees with atlas cell names")
  metadata_df$cell_id <- NULL
}
meta <- as.data.table(metadata_df, keep.rownames = "cell_id")
rm(metadata_df, atlas)
invisible(gc())

# Exact cell correspondence and final sample composition.
if (ncol(counts) != 198122L || nrow(meta) != 198122L)
  stop("Expected 198,122 final cells; found ", ncol(counts))
if (!identical(colnames(counts), meta$cell_id))
  stop("Count columns and metadata cell IDs are not in exactly the same order")
if (anyDuplicated(meta$cell_id) || anyDuplicated(rownames(counts)) ||
    anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
    any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts))))
  stop("Gene or cell identifiers are duplicated or missing")
if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
    any(counts@x != round(counts@x)))
  stop("Counts contain nonfinite, negative, or noninteger UMI values")

needed <- c("mouse_id", "age_group", "primary_group", "sequence_batch", "allen_subclass")
if (!all(needed %in% names(meta)))
  stop("Final atlas lacks required metadata: ",
       paste(setdiff(needed, names(meta)), collapse = ", "))
meta[, `:=`(mouse_id = as.character(mouse_id),
            age_group = as.character(age_group),
            primary_group = as.character(primary_group),
            sequence_batch = as.character(sequence_batch),
            allen_subclass = as.character(allen_subclass))]
meta[, cognition_group := fifelse(
  tolower(age_group) == "young", "Young",
  fifelse(tolower(primary_group) == "intact", "Intact",
          fifelse(tolower(primary_group) == "impaired", "Impaired", NA_character_)))]
if (anyNA(meta$mouse_id) || anyNA(meta$cognition_group) ||
    anyNA(meta$sequence_batch))
  stop("Missing mouse ID, group, or sequencing batch")
mice <- unique(meta[, .(mouse_id, cognition_group, sequence_batch)])
if (nrow(mice) != 13L || uniqueN(mice$mouse_id) != 13L ||
    !identical(as.integer(table(factor(mice$cognition_group,
                 levels = c("Young", "Intact", "Impaired")))), c(4L, 4L, 5L)))
  stop("Expected 13 distinct mice: 4 Young / 4 Intact / 5 Impaired")

# Transfer FROZEN final astrocyte states into the cell metadata only;
# non-state-assigned cells have NA, not an invented state label.
states <- fread(states_file)
if (!all(c("cell_id", "mouse_id", "state") %in% names(states)) ||
    nrow(states) != 10155L || anyDuplicated(states$cell_id) ||
    !setequal(unique(states$state), paste0("cAST", 0:5)))
  stop("Script 08 frozen astrocyte-state assignments are incomplete")
i <- match(states$cell_id, meta$cell_id)
if (anyNA(i) || !identical(as.character(meta$mouse_id[i]), as.character(states$mouse_id)))
  stop("Astrocyte-state cell IDs or mouse identities disagree with atlas")
meta[, astrocyte_state := NA_character_]
meta[i, astrocyte_state := as.character(states$state)]
rm(states, i)

# Retain a separately named original barcode if it can be verified, otherwise
# leave it blank. cell_id is always the exact column name of the count matrix.
barcode_columns <- intersect(c("original_10x_barcode", "original_barcode",
  "barcode", "cell_barcode", "barcode_10x", "barcode_raw"), names(meta))
if (length(barcode_columns)) {
  meta[, original_10x_barcode := as.character(get(barcode_columns[[1L]]))]
} else {
  prefix <- paste0(meta$mouse_id, "_")
  if (all(startsWith(meta$cell_id, prefix))) {
    meta[, original_10x_barcode := substring(cell_id, nchar(prefix) + 1L)]
  } else {
    meta[, original_10x_barcode := NA_character_]
  }
}

# A flat CSV cannot reliably encode list columns.
list_cols <- names(meta)[vapply(meta, is.list, logical(1))]
if (length(list_cols))
  stop("Metadata contains list-valued columns: ", paste(list_cols, collapse = ", "))
setcolorder(meta, c("cell_id", "original_10x_barcode", "mouse_id",
  "cognition_group", "age_group", "primary_group", "sequence_batch",
  "allen_subclass", "astrocyte_state",
  setdiff(names(meta), c("cell_id", "original_10x_barcode", "mouse_id",
    "cognition_group", "age_group", "primary_group", "sequence_batch",
    "allen_subclass", "astrocyte_state"))))

# Export four standard Matrix Market / 10x-compatible companion files.
# The atlas row names are preserved as feature identifiers. In the absence
# of verified separate stable gene IDs, feature names intentionally repeat
# the atlas identifiers; no Ensembl identifiers are inferred or invented.
# Matrix Market contains only the numeric coordinates, so the accompanying
# feature and barcode files define exactly the same row/column identities.
mtx_file <- file.path(out, "Baier_final_matrix.mtx.gz")
features_file <- file.path(out, "Baier_final_features.tsv.gz")
barcodes_file <- file.path(out, "Baier_final_barcodes.tsv.gz")
metadata_file <- file.path(out, "Baier_final_cell_metadata.csv")

message("Writing Matrix Market UMI counts (genes x cells)")
# writeMM writes an uncompressed Matrix Market file. Stream-compress it in
# bounded blocks so the full matrix is never converted to a dense array.
mtx_tmp <- tempfile(pattern = "Baier_final_matrix_", tmpdir = out,
                    fileext = ".mtx")
compress_mtx <- function(source, destination) {
  source_con <- file(source, open = "rb")
  destination_con <- gzfile(destination, open = "wb", compression = 6L)
  on.exit({
    close(source_con)
    close(destination_con)
  }, add = TRUE)
  repeat {
    chunk <- readBin(source_con, what = "raw", n = 8L * 1024L * 1024L)
    if (!length(chunk)) break
    writeBin(chunk, destination_con)
  }
}
Matrix::writeMM(counts, mtx_tmp)
compress_mtx(mtx_tmp, mtx_file)
unlink(mtx_tmp)

message("Writing feature and barcode TSV files")
data.table::fwrite(
  data.table::data.table(feature_id = rownames(counts),
                         feature_name = rownames(counts),
                         feature_type = rep("Gene Expression", nrow(counts))),
  file = features_file, sep = "\t", col.names = FALSE, quote = FALSE,
  compress = "gzip")
data.table::fwrite(
  data.table::data.table(cell_id = colnames(counts)),
  file = barcodes_file, sep = "\t", col.names = FALSE, quote = FALSE,
  compress = "gzip")

message("Writing cell-level metadata CSV")
data.table::fwrite(meta, file = metadata_file, na = "", quote = TRUE)

# Reopen the small companion files to verify that the written order agrees
# exactly with the source matrix and the exported metadata.
written_features <- data.table::fread(features_file, header = FALSE,
                                      sep = "\t", showProgress = FALSE)
written_barcodes <- data.table::fread(barcodes_file, header = FALSE,
                                      sep = "\t", showProgress = FALSE)
if (ncol(written_features) != 3L ||
    !identical(as.character(written_features[[1L]]), rownames(counts)) ||
    !identical(as.character(written_barcodes[[1L]]), colnames(counts)) ||
    !identical(as.character(written_barcodes[[1L]]), meta$cell_id))
  stop("Exported feature/barcode order does not match the count matrix/metadata")

# Validate dimensions and number of nonzero entries in the actual .mtx.gz.
mtx_con <- gzfile(mtx_file, open = "rt")
header_lines <- readLines(mtx_con, n = 12L, warn = FALSE)
size_line <- header_lines[!grepl("^%", header_lines)][1L]
size_values <- scan(text = size_line, what = integer(), quiet = TRUE)
if (length(size_values) != 3L ||
    !identical(as.integer(size_values),
               as.integer(c(nrow(counts), ncol(counts), length(counts@x)))))
  stop("Exported Matrix Market dimensions/nonzero count do not match source")
close(mtx_con)

outputs <- c(mtx_file, features_file, barcodes_file, metadata_file)
if (anyNA(file.info(outputs)$size) || any(file.info(outputs)$size <= 0))
  stop("One or more GEO export files are missing or empty")
message("GEO export complete: ", nrow(counts), " genes x ", ncol(counts),
        " cells; 13 mice; 10,155 cleaned astrocyte-state assignments")
message("GEO files written to: ", out)
for (p in outputs) message("  ", basename(p))
