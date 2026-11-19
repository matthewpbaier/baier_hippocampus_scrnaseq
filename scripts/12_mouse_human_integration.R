#Baier et al.
#12: Mouse-human integration, metabolic-program recurrence, and OXPHOS gene identity

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "12_mouse_human_integration")
ref_dir <- file.path(root, "reference")
for (p in c(out, ref_dir))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(data.table)
library(edgeR)
library(fgsea)

AXES <- c("mouse_cognition", "mouse_aging", "human_aging",
          "seaad_pathology", "seaad_cognition")
LINEAGES <- c("Astrocyte", "Microglia", "OPC", "Oligodendrocyte")
MOUSE_NAMES <- c(Astrocyte = "Astro-TE NN", Microglia = "Microglia NN",
                 OPC = "OPC NN", Oligodendrocyte = "Oligo NN")
OX <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
BLOCKS <- list(Mouse = AXES[1:2], Human_aging = AXES[3], SEAAD = AXES[4:5])
N_MATCHED <- 2000L          
N_NEAREST <- 15L                
N_CORE_PERM <- 50000L            
N_SELECTED_PERM <- 20000L        
SEED_SELECTED <- 20261026L
KEY <- function(x) toupper(trimws(as.character(x)))
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, paste0("12_", name, ".csv")))

need_file <- function(path) {
  if (!file.exists(path)) stop("Missing upstream output: ", path)
  path
}
check_cols <- function(x, columns, source) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(source, ": missing columns ", paste(missing, collapse = ", "))
}

#Load previous script outputs
mouse_dir <- file.path(root, "06_mouse_bulk")
zemke_file <- need_file(file.path(root, "09_human_normative_aging",
                                   "09_zemke_alllineages_aging_normative_aging_checkpoint.rds"))
seaad_file <- need_file(file.path(root, "10_human_SEAAD",
                                   "10_seaad_alllineages_pathology_cognition_primary_checkpoint.rds"))
mouse_cognition_file <- need_file(file.path(mouse_dir,
                              "06_baiermouse_alllineages_cognition_twogroup_gene_effects.csv"))
pb_file <- need_file(file.path(mouse_dir, "06_baiermouse_alllineages_pseudobulk_raw_counts.rds"))
mouse_meta_file <- need_file(file.path(mouse_dir, "06_baiermouse_allmice_pseudobulk_sample_metadata.csv"))
length_file <- need_file(file.path(root, "08_astrocyte_states",
                                   "08_mouse_mm39_gene_lengths_USED.csv"))
module_file <- need_file(file.path(root, "08_astrocyte_states",
                                   "08_metabolic_module_membership_USED.csv"))
rosmap_file <- need_file(file.path(root, "11_human_ROSMAP",
                                   "11_rosmap_astrocyte_pathology_regions_pathway_region_results.csv"))

human09 <- readRDS(zemke_file)
human10 <- readRDS(seaad_file)
hallmarks <- human09$human_hallmark_sets
hallmarks <- lapply(hallmarks, function(v) unique(KEY(v)))

saveRDS(hallmarks, file.path(out, "12_crosscohort_allglia_five_axes_human_hallmarks.rds"))

mouse_cog <- fread(mouse_cognition_file)
check_cols(mouse_cog, c("lineage", "gene", "logFC", "logCPM", "F", "PValue", "FDR"),
           "Mouse two-group cognition")
pb <- readRDS(pb_file)
meta <- fread(mouse_meta_file)
check_cols(meta, c("mouse_id", "primary_group", "sequence_batch"), "Mouse metadata")

#Construct two-group mouse aging input 
fit_mouse_aging <- function(lineage) {
  m <- pb[[lineage]]
  if (is.null(m)) stop("No mouse pseudobulk for ", lineage)
  inf <- copy(meta[match(colnames(m), mouse_id)])
  if (!identical(as.character(inf$mouse_id), colnames(m)))
    stop("Mouse columns and metadata disagree for ", lineage)
  idx <- which(tolower(as.character(inf$primary_group)) %in% c("young", "intact"))
  inf <- inf[idx]
  group <- factor(ifelse(tolower(inf$primary_group) == "young", "young", "aged_intact"),
                  levels = c("young", "aged_intact"))
  batch <- factor(inf$sequence_batch)
  design <- model.matrix(~ batch + group)
  if (qr(design)$rank != ncol(design)) stop("Mouse aging design is rank deficient")
  y <- edgeR::DGEList(counts = as.matrix(m[, idx, drop = FALSE]))
  keep <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  z <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = "groupaged_intact"),
                      n = Inf, sort.by = "PValue")$table
  data.table(lineage = lineage, gene = rownames(z), logFC = z$logFC,
             logCPM = z$logCPM, F = z$F, PValue = z$PValue, FDR = z$FDR)
}
mouse_age <- rbindlist(lapply(unname(MOUSE_NAMES), fit_mouse_aging))
write_table(mouse_age, "baiermouse_allglia_aging_mouse_two_group_aging_de_aged_intact_minus_young")

#Extract gene-level inputs from the two human upstream models
#Preserve raw model direction before applying favorable-state multipliers
get_seaad <- function(lineage, coefficient) {
  z <- human10$primary_de[[lineage]]$tests[[coefficient]]
  if (is.null(z)) stop("Script 10 missing SEA-AD ", lineage, " / ", coefficient)
  as.data.table(z)
}

get_zemke <- function(lineage) {
  z <- human09$zemke_age_de[[lineage]]$results
  if (is.null(z)) stop("Script 09 missing normative aging ", lineage)
  as.data.table(z)
}

native <- list()

for (ln in LINEAGES) {
  mk <- MOUSE_NAMES[[ln]]
  native[[paste(ln, "mouse_cognition", sep = "__")]] <- copy(mouse_cog[lineage == mk])
  native[[paste(ln, "mouse_aging", sep = "__")]] <- copy(mouse_age[lineage == mk])
  native[[paste(ln, "human_aging", sep = "__")]] <- get_zemke(ln)
  native[[paste(ln, "seaad_pathology", sep = "__")]] <- get_seaad(ln, "cps_z")
  native[[paste(ln, "seaad_cognition", sep = "__")]] <- get_seaad(ln, "cog_z")
}

ORIENT <- c(mouse_cognition = -1, mouse_aging = -1, human_aging = -1,
            seaad_pathology = -1, seaad_cognition = 1)
percentile <- function(effect, strength = abs(effect)) {
  ok <- is.finite(effect) & is.finite(strength)
  out <- rep(NA_real_, length(effect))
  out[ok] <- sign(effect[ok]) * rank(strength[ok], ties.method = "average") / sum(ok)
  out
}

for (ln in LINEAGES) for (axis in AXES) {
  id <- paste(ln, axis, sep = "__")
  z <- as.data.table(native[[id]])
  check_cols(z, c("gene", "logFC", "F", "PValue", "FDR"), id)
  z <- z[!is.na(gene) & nzchar(as.character(gene)) & is.finite(logFC) & is.finite(F)]
  z[, native_gene := as.character(gene)]
  if (anyDuplicated(z$native_gene)) stop("Duplicated native gene in ", id)
  z[, `:=`(lineage = ln, axis_id = axis,
           species = if (axis %in% AXES[1:2]) "mouse" else "human",
           oriented_logFC = as.numeric(logFC) * ORIENT[[axis]])]
  z[, signed_effect_percentile := percentile(oriented_logFC)]
  z[, signed_inference_percentile := percentile(oriented_logFC, sqrt(pmax(F, 0)))]
  native[[id]] <- z
}

write_table(rbindlist(native, fill = TRUE)[, .(lineage, axis_id, species,
  native_gene, logFC, F, PValue, FDR, oriented_logFC,
  signed_effect_percentile, signed_inference_percentile)], "crosscohort_allglia_five_axes_all_native_gene_statistics")

#Mouse-human ortholog mapping
previous_map_file <- file.path(ref_dir, "11_rosmap_astrocyte_pathology_regions_mouse_human_strict_1to1_orthologs.csv")
full_map_file <- file.path(out, "12_crosscohort_allglia_five_axes_full_mouse_axis_union_ortholog_map.csv")
mouse_ids <- sort(unique(unlist(lapply(native, function(z)
  if (unique(z$species) == "mouse") z$native_gene else NULL))))

#Mouse DE tables contain both gene symbols and ENSMUSG identifiers
is_ensembl <- grepl("^ENSMUSG[0-9]+(\\.[0-9]+)?$", mouse_ids,
                    ignore.case = TRUE, perl = TRUE)
mouse_symbols <- mouse_ids[!is_ensembl]
mouse_ensembl <- unique(sub("\\.[0-9]+$", "", mouse_ids[is_ensembl]))
write_table(data.table(identifier_type = c("mouse_symbol", "mouse_ensembl"),
                       n_unique = c(length(mouse_symbols), length(mouse_ensembl))),
            "crosscohort_allglia_five_axes_mouse_identifier_type_summary")
message("Ortholog mapping: ", length(mouse_symbols), " mouse symbols; ",
        length(mouse_ensembl), " mouse Ensembl IDs")

orth_symbol <- as.data.table(babelgene::orthologs(
  genes = mouse_symbols, species = "mouse", human = FALSE,
  min_support = 3, top = FALSE))
orth_ensembl <- if (length(mouse_ensembl)) as.data.table(babelgene::orthologs(
  genes = mouse_ensembl, species = "mouse", human = FALSE,
  min_support = 3, top = FALSE)) else data.table()
raw_map <- rbindlist(list(symbol = orth_symbol, ensembl = orth_ensembl),
                     idcol = "input_identifier_type", use.names = TRUE, fill = TRUE)
check_cols(raw_map, c("symbol", "human_symbol"), "Full mouse-axis ortholog mapping")
raw_map <- raw_map[!is.na(symbol) & nzchar(as.character(symbol)) &
                   !is.na(human_symbol) & nzchar(as.character(human_symbol))]
orth <- unique(raw_map[, .(mouse_gene = KEY(symbol), human_gene = KEY(human_symbol))])
orth[, n_human := uniqueN(human_gene), by = mouse_gene]
orth[, n_mouse := uniqueN(mouse_gene), by = human_gene]
orth[, strict_one_to_one := n_human == 1L & n_mouse == 1L]
fwrite(orth[, .(mouse_gene, human_gene, strict_one_to_one)], full_map_file)
orth <- orth[strict_one_to_one == TRUE, .(mouse_gene, human_gene)]
stopifnot(!anyDuplicated(orth$mouse_gene), !anyDuplicated(orth$human_gene))
write_table(orth, "crosscohort_allglia_five_axes_strict_mouse_human_orthologs")

if (file.exists(previous_map_file)) {
  prior <- fread(previous_map_file)
  check_cols(prior, c("mouse_gene", "human_gene", "strict_one_to_one"),
             "Earlier ROSMAP mapping")
  prior <- unique(prior[as.character(strict_one_to_one) %in% c("TRUE", "T", "1"),
                        .(mouse_gene = KEY(mouse_gene), human_gene = KEY(human_gene))])
  mapping_delta <- merge(prior[, .(mouse_gene, human_gene, from_11_map = TRUE)],
                         orth[, .(mouse_gene, human_gene, from_full_union = TRUE)],
                         by = c("mouse_gene", "human_gene"), all = TRUE)
  mapping_delta[is.na(from_11_map), from_11_map := FALSE]
  mapping_delta[is.na(from_full_union), from_full_union := FALSE]
  mapping_delta[, classification := fifelse(from_11_map & from_full_union,
    "shared", fifelse(from_full_union, "only_full_union", "only_script_11"))]
  write_table(mapping_delta, "crosscohort_allglia_five_axes_ortholog_map_delta_vs_script11")
}

input_alias <- unique(orth_symbol[!is.na(symbol), .(
  input_gene = KEY(symbol), mouse_gene = KEY(symbol))])
if (nrow(orth_ensembl)) {
  check_cols(orth_ensembl, "ensembl", "Mouse Ensembl ortholog mapping")
  ens_alias <- unique(orth_ensembl[!is.na(ensembl) & !is.na(symbol), .(
    input_gene = KEY(sub("\\.[0-9]+$", "", as.character(ensembl))),
    mouse_gene = KEY(symbol))])
  input_alias <- unique(rbindlist(list(input_alias, ens_alias)))
}
input_alias <- input_alias[!is.na(input_gene) & nzchar(input_gene) &
                           !is.na(mouse_gene) & nzchar(mouse_gene)]
input_alias <- input_alias[, if (uniqueN(mouse_gene) == 1L) .SD[1L] else NULL,
                           by = input_gene]
write_table(input_alias, "crosscohort_allglia_five_axes_mouse_input_to_canonical_symbol")

#Map each native axis to human-gene space
masters <- setNames(vector("list", length(LINEAGES)), LINEAGES)
for (ln in LINEAGES) {
  wide_parts <- list()
  for (ax in AXES) {
    z <- copy(native[[paste(ln, ax, sep = "__")]])
    if (unique(z$species) == "mouse") {
      z[, input_gene := KEY(sub("\\.[0-9]+$", "", native_gene))]
      z <- merge(z, input_alias, by = "input_gene", all = FALSE)
      z <- merge(z, orth, by = "mouse_gene", all = FALSE)
    } else {
      z[, human_gene := KEY(native_gene)]
    }
    col <- c("oriented_logFC", "F", "PValue", "FDR",
             "signed_effect_percentile", "signed_inference_percentile")
    v <- copy(z[, c("human_gene", col), with = FALSE])
    setnames(v, col, paste0(ax, "_", c("oriented_logFC", "F", "PValue", "FDR",
                                       "effect_pct", "inference_pct")))
    wide_parts[[ax]] <- v
  }
  coverage <- data.table(human_gene = sort(unique(unlist(
    lapply(wide_parts, function(v) v$human_gene)))))
  for (axis in AXES)
    coverage[, (axis) := human_gene %in% wide_parts[[axis]]$human_gene]
  coverage[, `:=`(lineage = ln,
                  n_axes = rowSums(as.matrix(.SD))), .SDcols = AXES]
  fwrite(coverage, file.path(out, paste0("12_crosscohort_", tolower(ln), "_five_axes_gene_universe_coverage.csv")))
  w <- Reduce(function(a, b) merge(a, b, by = "human_gene", all = FALSE), wide_parts)
  setorder(w, human_gene)
  w[, mouse_block_effect := (mouse_cognition_effect_pct + mouse_aging_effect_pct) / 2]
  w[, human_aging_block_effect := human_aging_effect_pct]
  w[, seaad_block_effect := (seaad_pathology_effect_pct + seaad_cognition_effect_pct) / 2]
  w[, consensus_effect := (mouse_block_effect + human_aging_block_effect + seaad_block_effect) / 3]
  w[, mouse_block_inference := (mouse_cognition_inference_pct + mouse_aging_inference_pct) / 2]
  w[, human_aging_block_inference := human_aging_inference_pct]
  w[, seaad_block_inference := (seaad_pathology_inference_pct + seaad_cognition_inference_pct) / 2]
  w[, consensus_inference := (mouse_block_inference + human_aging_block_inference + seaad_block_inference) / 3]
  w[, `:=`(direction_n_5axes = rowSums(as.matrix(.SD) > 0),
            all5_preserved = rowSums(as.matrix(.SD) > 0) == 5L),
    .SDcols = paste0(AXES, "_effect_pct")]
  w[, `:=`(consensus_effect_rank = frank(-consensus_effect, ties.method = "min"),
            consensus_inference_rank = frank(-consensus_inference, ties.method = "min"),
            lineage = ln)]
  setorder(w, consensus_effect_rank, human_gene)
  masters[[ln]] <- w
}
write_table(rbindlist(masters, fill = TRUE), "crosscohort_allglia_five_axes_gene_consensus_all_glia")
write_table(rbindlist(lapply(LINEAGES, function(ln)
  data.table(lineage = ln, n_common_genes = nrow(masters[[ln]])))),
  "crosscohort_allglia_five_axes_gene_coverage")

#Cohort-balanced Hallmark integration
#Each lineage/axis/fold receives its own deterministic FGSEA stream
integration_seed <- function(lineage, family, member = "") {
  j <- match(lineage, LINEAGES)
  stopifnot(!is.na(j))
  block <- switch(family, integrated = 0L, block_omission = 1000L,
                  axis = 2000L)
  k <- if (family == "block_omission") match(member, names(BLOCKS)) else
       if (family == "axis") match(member, AXES) else 0L
  stopifnot(!is.na(k))
  as.integer(190929L + 100L * j + block + k)
}

run_fgsea <- function(score, gene, seed, alphabetize = FALSE) {
  ranks <- setNames(as.numeric(score), as.character(gene))
  ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))]
  if (anyDuplicated(names(ranks))) stop("Duplicate genes in FGSEA ranks")
  if (alphabetize) ranks <- ranks[order(names(ranks))]
  # Preserve original deterministic tie-breaking and FGSEA seed.
  ranks <- sort(ranks + seq_along(ranks) * 1e-14, decreasing = TRUE)
  set.seed(seed)
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = hallmarks, stats = ranks, minSize = 10L,
    maxSize = 500L, eps = 0, nproc = 1))
  fg <- as.data.table(fg)
  fg[, `:=`(preserved_rank = frank(-NES, ties.method = "min"),
            leadingEdge_string = vapply(leadingEdge, paste,
                                        collapse = ";", FUN.VALUE = character(1)))]
  fg[, preserved_rank_percentile :=
       if (.N > 1L) 1 - (preserved_rank - 1) / (.N - 1) else 1]
  fg
}
integrated <- rbindlist(lapply(LINEAGES, function(ln) {
  m <- masters[[ln]]
  z <- run_fgsea(m$consensus_effect, m$human_gene,
                 seed = integration_seed(ln, "integrated"))
  z[, `:=`(lineage = ln, consensus_type = "effect")]
  z
}), fill = TRUE)
if (nrow(integrated) != 200L || anyDuplicated(integrated, by = c("lineage", "pathway")))
  stop("Integrated Hallmark screen must be four lineages x 50 pathways")
write_table(integrated[, .(lineage, pathway, NES, pval, padj,
                            preserved_rank, preserved_rank_percentile,
                            leadingEdge_string)], "crosscohort_allglia_five_axes_cohort_balanced_hallmarks_all200")

rank_summary <- integrated[, .(
  mean_rank = mean(preserved_rank),
  worst_rank = max(preserved_rank),
  best_rank = min(preserved_rank),
  n_top5 = sum(preserved_rank <= 5L)), by = pathway][order(mean_rank, worst_rank, pathway)]
write_table(rank_summary, "crosscohort_allglia_five_axes_integrated_cross_lineage_rank_summary")

#Cohort-omission sensitivity
block_loo <- rbindlist(lapply(names(BLOCKS), function(omit) {
  keep <- switch(omit,
    Mouse = c("human_aging_block_effect", "seaad_block_effect"),
    Human_aging = c("mouse_block_effect", "seaad_block_effect"),
    SEAAD = c("mouse_block_effect", "human_aging_block_effect"))
  rbindlist(lapply(LINEAGES, function(ln) {
    m <- masters[[ln]]
    score <- rowMeans(as.matrix(m[, ..keep]))
    z <- run_fgsea(score, m$human_gene,
                   seed = integration_seed(ln, "block_omission", omit))
    z[, `:=`(lineage = ln, omitted_block = omit)]
    z[, .(lineage, omitted_block, pathway, NES, pval, padj, preserved_rank)]
  }))
}))
write_table(block_loo, "crosscohort_allglia_five_axes_cohort_block_omission_hallmarks")

#Axis-resolved analysis
axis_results <- rbindlist(lapply(LINEAGES, function(ln) {
  m <- masters[[ln]]
  rbindlist(lapply(AXES, function(ax) {
    score <- sign(m[[paste0(ax, "_oriented_logFC")]]) *
      sqrt(pmax(m[[paste0(ax, "_F")]], 0))
    z <- run_fgsea(score, m$human_gene,
                   seed = integration_seed(ln, "axis", ax), alphabetize = TRUE)
    setnames(z, "preserved_rank", "favorable_rank")
    setnames(z, "preserved_rank_percentile", "favorable_rank_percentile")
    z[, `:=`(lineage = ln, axis_id = ax,
              statistic_type = "signed_sqrtF", universe_type = "lineage_specific")]
    z
  }), fill = TRUE)
}), fill = TRUE)
if (nrow(axis_results) != 1000L ||
    anyDuplicated(axis_results, by = c("lineage", "axis_id", "pathway")))
  stop("Axis-level Hallmark screen must be 4 x 5 x 50 = 1000 results")
axis_results[, FDR_screen1000 := p.adjust(pval, method = "BH")]
write_table(axis_results[, .(lineage, axis_id, pathway, NES, pval, padj,
  FDR_screen1000, favorable_rank, favorable_rank_percentile,
  leadingEdge_string)], "crosscohort_allglia_five_axes_hallmarks_all1000")

recurrence <- axis_results[, .(
  n_axes_favorable = sum(NES > 0),
  n_axes_favorable_FDR05 = sum(NES > 0 & padj < .05),
  all5_favorable = all(NES > 0),
  mean_NES = mean(NES),
  min_NES = min(NES),
  mean_favorable_rank_percentile = mean(favorable_rank_percentile)
), by = .(lineage, pathway)]
setorder(recurrence, lineage, -all5_favorable, -n_axes_favorable_FDR05,
         -mean_favorable_rank_percentile, -mean_NES)
recurrence[, recurrence_rank := seq_len(.N), by = lineage]
write_table(recurrence, "crosscohort_allglia_five_axes_pathway_recurrence_all200")

#Mouse-cognition versus human-only pathway ranks
rank_wide <- dcast(axis_results[, .(lineage, pathway, axis_id, favorable_rank)],
                   lineage + pathway ~ axis_id, value.var = "favorable_rank")
rank_wide[, `:=`(mouse_block = (mouse_cognition + mouse_aging) / 2,
                human_aging_block = human_aging,
                seaad_block = (seaad_pathology + seaad_cognition) / 2)]
rank_wide[, human_only_score := (human_aging_block + seaad_block) / 2]
rank_wide[, human_only_rank := frank(human_only_score, ties.method = "average"),
          by = lineage]
rank_wide[, rank_shift := human_only_rank - mouse_cognition]
write_table(rank_wide, "crosscohort_allglia_five_axes_mouse_versus_human_pathway_ranks")

#Rank correlations, block-omission stability, and NES-based mouse/human correspondence
block_cor <- rbindlist(lapply(LINEAGES, function(ln) {
  z <- integrated[lineage == ln, .(pathway, full_rank = preserved_rank)]
  rbindlist(lapply(names(BLOCKS), function(drop) {
    w <- merge(z, block_loo[lineage == ln & omitted_block == drop,
                            .(pathway, omitted_rank = preserved_rank)], by = "pathway")
    data.table(lineage = ln, omitted_block = drop, n_pathways = nrow(w),
      rho = cor(w$full_rank, w$omitted_rank, method = "spearman"))
  }))
}))
write_table(block_cor, "crosscohort_allglia_five_axes_cohort_block_omission_rank_correlations")

rank_cor <- rank_wide[, {
  z <- suppressWarnings(cor.test(mouse_cognition, human_only_rank,
                                 method = "spearman", exact = FALSE))
  .(n_pathways = .N, rho = unname(z$estimate), p_value = z$p.value)
}, by = lineage]
rank_cor[, FDR_4 := p.adjust(p_value, "BH")]
write_table(rank_cor, "crosscohort_allglia_five_axes_mouse_human_rank_correlations")

axis_nes <- axis_results[, .(lineage, pathway, axis_id, NES)]
axis_nes[, scale95 := as.numeric(quantile(abs(NES), .95)),
          by = .(lineage, axis_id)]
axis_nes[, NES_scaled := NES / scale95]
nes_wide <- dcast(axis_nes, lineage + pathway ~ axis_id,
                  value.var = "NES_scaled")
nes_wide[, human_composite := .50 * human_aging + .25 * seaad_pathology +
                             .25 * seaad_cognition]
nes_cor <- nes_wide[, {
  z <- suppressWarnings(cor.test(mouse_cognition, human_composite,
                                 method = "spearman", exact = FALSE))
  .(n_pathways = .N, rho = unname(z$estimate), p_value = z$p.value)
}, by = lineage]
nes_cor[, FDR_4 := p.adjust(p_value, "BH")]
write_table(nes_wide, "crosscohort_allglia_five_axes_mouse_human_scaled_nes_scatter_source")
write_table(nes_cor, "crosscohort_allglia_five_axes_mouse_human_scaled_nes_correlations")

separate_human_cor <- rbindlist(lapply(AXES[3:5], function(ax) {
  z <- dcast(axis_results[axis_id %in% c("mouse_cognition", ax),
                          .(lineage, pathway, axis_id, NES)],
             lineage + pathway ~ axis_id, value.var = "NES")
  z[, {
    ct <- suppressWarnings(cor.test(mouse_cognition, get(ax),
                                    method = "spearman", exact = FALSE))
    .(axis_id = ax, n_pathways = .N, rho = unname(ct$estimate),
      p_value = ct$p.value)
  }, by = lineage]
}))
write_table(separate_human_cor, "crosscohort_allglia_five_axes_mouse_cognition_vs_each_human_axis_nes_correlations")

#Five-axis gene-level OXPHOS recurrence
ox_master <- copy(masters[["Astrocyte"]][human_gene %in% hallmarks[[OX]]])
setorder(ox_master, human_gene)
O_GENES <- ox_master$human_gene

get_boolean <- function(m, ln, pathway_name) {
  gg <- sort(intersect(m$human_gene, hallmarks[[pathway_name]]))
  if (!length(gg)) stop("No shared pathway genes in ", ln, ": ", pathway_name)
  r <- m[match(gg, human_gene)]
  d <- sapply(AXES, function(ax) r[[paste0(ax, "_effect_pct")]] > 0)
  rownames(d) <- gg
  le <- matrix(FALSE, nrow = length(gg), ncol = 5L,
               dimnames = list(gg, AXES))
  hz <- axis_results[lineage == ln & pathway == pathway_name]
  if (nrow(hz) != 5L) stop("Missing one of five leading edges for ", ln, ": ", pathway_name)
  for (ax in AXES) {
    s <- hz[axis_id == ax, leadingEdge_string]
    edges <- if (length(s) == 1L && !is.na(s) && nzchar(s))
      KEY(strsplit(s, ";", fixed = TRUE)[[1]]) else character()
    le[, ax] <- gg %in% edges
  }
  list(genes = gg, direction = d, leading_edge = le, master = r,
       n5 = sum(rowSums(d) == 5L),
       n_le4 = sum(rowSums(d) == 5L & rowSums(le) >= 4L),
       n_le5 = sum(rowSums(d) == 5L & rowSums(le) == 5L))
}

ox <- get_boolean(ox_master, "Astrocyte", OX)
if (!identical(c(ox$n5, ox$n_le4, ox$n_le5), c(60L, 34L, 6L)))
  stop("Final OXPHOS recurrence differs from 185/60/34/6. ",
       "Do not interpret downstream conditional nulls until upstream identities match.")
ox_manifest <- data.table(human_gene = O_GENES,
  n_axes_favorable = rowSums(ox$direction),
  n_axes_leading_edge = rowSums(ox$leading_edge),
  all5_preserved = rowSums(ox$direction) == 5L,
  recurrent34 = rowSums(ox$direction) == 5L & rowSums(ox$leading_edge) >= 4L,
  recurrent_le5 = rowSums(ox$direction) == 5L & rowSums(ox$leading_edge) == 5L,
  consensus_effect = ox_master$consensus_effect)
write_table(ox_manifest, "crosscohort_astrocyte_recurrence_oxphos_final_185_gene_manifest")
ox_effects <- rbindlist(lapply(AXES, function(ax) data.table(
  human_gene = O_GENES, axis_id = ax,
  signed_effect_pct = ox_master[[paste0(ax, "_effect_pct")]],
  oriented_logFC = ox_master[[paste0(ax, "_oriented_logFC")]],
  LE = as.logical(ox$leading_edge[, ax]),
  recurrent34 = ox_manifest$recurrent34)))
write_table(ox_effects, "crosscohort_astrocyte_recurrence_oxphos_185_five_axis_effects_and_leading_edges")

#Donor-level score
zemke_astro <- human09$zemke_age_de[["Astrocyte"]]
expr <- edgeR::cpm(zemke_astro$dge, log = TRUE, prior.count = 2)
key <- KEY(rownames(expr))
stopifnot(!anyDuplicated(key))
use <- intersect(O_GENES, key)
expr <- expr[match(use, key), , drop = FALSE]
variable <- apply(expr, 1, stats::sd)
expr <- expr[is.finite(variable) & variable > 0, , drop = FALSE]
stopifnot(length(O_GENES) == 185L, nrow(expr) == 185L)
zs <- t(scale(t(expr)))
score <- colMeans(zs)
donor_info <- human09$zemke_sample_info
donor_score <- donor_info[match(names(score), donor_info$donor_id), , drop = FALSE]
stopifnot(identical(as.character(donor_score$donor_id), names(score)))
donor_score$oxphos_score <- unname(score)
donor_score$Gender <- factor(donor_score$Gender, levels = c("F", "M"))
score_fit <- stats::lm(oxphos_score ~ Gender + age_z, data = donor_score)
coefs <- summary(score_fit)$coefficients
ci <- confint(score_fit, "age_z")
model_summary <- data.table(
  term = "age_z", n_donors = nrow(donor_score), n_genes = nrow(expr),
  beta = unname(coefs["age_z", "Estimate"]),
  se = unname(coefs["age_z", "Std. Error"]),
  ci_low = unname(ci[1]), ci_high = unname(ci[2]),
  nominal_P = unname(coefs["age_z", "Pr(>|t|)"]))
reference_sex <- donor_score
reference_sex$Gender <- factor(rep("F", nrow(reference_sex)), levels = c("F", "M"))
donor_score$sex_adjusted_score <- as.numeric(
  predict(score_fit, newdata = reference_sex) + residuals(score_fit))
loo <- rbindlist(lapply(seq_len(nrow(donor_score)), function(i) {
  fit <- stats::lm(oxphos_score ~ Gender + age_z,
                   data = donor_score[-i, , drop = FALSE])
  ct <- summary(fit)$coefficients
  data.table(omitted_donor = donor_score$donor_id[[i]],
             beta = unname(ct["age_z", "Estimate"]),
             nominal_P = unname(ct["age_z", "Pr(>|t|)"]))
}))
write_table(donor_score, "zemke_astrocyte_aging_astrocyte_oxphos_donor_scores")
write_table(model_summary, "zemke_astrocyte_aging_astrocyte_oxphos_age_model")
write_table(loo, "zemke_astrocyte_aging_astrocyte_oxphos_leave_one_donor_out")

#OXPHOS in each of the four lineages
ox_lineages <- rbindlist(lapply(LINEAGES, function(ln) {
  z <- get_boolean(masters[[ln]], ln, OX)
  data.table(lineage = ln, n_measurable = length(z$genes),
             n_favorable_all5 = z$n5,
             n_favorable_and_LE_ge4 = z$n_le4,
             n_favorable_and_LE_all5 = z$n_le5)
}))
write_table(ox_lineages, "crosscohort_astrocyte_recurrence_oxphos_gene_recurrence_by_glial_lineage")

#Frozen mm39 length + mouse-lineage expression covariates
gl <- fread(length_file)
check_cols(gl, c("gene", "length"), "mm39 gene lengths")
gl[, `:=`(mouse_gene = KEY(gene), length = as.numeric(length))]
gl <- gl[is.finite(length) & length > 0,
         .(length = max(length)), by = mouse_gene]

build_covariates <- function(ln) {
  expr <- copy(mouse_cog[lineage == MOUSE_NAMES[[ln]]])
  expr[, input_gene := KEY(sub("\\.[0-9]+$", "", gene))]
  expr <- merge(expr, input_alias, by = "input_gene", all = FALSE)
  expr <- expr[is.finite(logCPM), .(logCPM = median(logCPM)), by = mouse_gene]
  cv <- merge(orth, gl, by = "mouse_gene", all.x = TRUE)
  cv <- merge(cv, expr, by = "mouse_gene", all.x = TRUE)
  cv[, complete := is.finite(length) & length > 0 & is.finite(logCPM)]
  cv
}
COV <- setNames(lapply(LINEAGES, build_covariates), LINEAGES)

repair_length_aliases <- function(cv) {
  missing <- cv[!is.finite(length) | length <= 0]
  db <- org.Mm.eg.db::org.Mm.eg.db
  ali <- AnnotationDbi::keys(db, keytype = "ALIAS")
  cand <- data.table(db_alias = ali, mouse_gene = KEY(ali))[
    mouse_gene %in% missing$mouse_gene]
  if (nrow(cand)) {
    ann <- as.data.table(AnnotationDbi::select(db, keys = unique(cand$db_alias),
                       keytype = "ALIAS", columns = c("SYMBOL", "ENTREZID")))
    ann <- merge(cand, ann, by.x = "db_alias", by.y = "ALIAS")
    ann[, official := KEY(SYMBOL)]
    z <- merge(missing[, .(human_gene, mouse_gene)],
               ann[, .(mouse_gene, official, ENTREZID)],
               by = "mouse_gene", allow.cartesian = TRUE)
    z <- unique(z[official == human_gene & !is.na(ENTREZID),
                  .(human_gene, mouse_gene, official, ENTREZID)])
  } else z <- data.table(human_gene = character(), mouse_gene = character(),
                        official = character(), ENTREZID = character())
  # Verified official-symbol fallback from strict human-mouse one-to-one map.
  remain <- setdiff(missing$human_gene, z$human_gene)
  if (length(remain)) {
    syms <- AnnotationDbi::keys(db, keytype = "SYMBOL")
    symtab <- data.table(symbol = syms, human_gene = KEY(syms))[
      human_gene %in% remain]
    if (nrow(symtab)) {
      a <- as.data.table(AnnotationDbi::select(db,
          keys = unique(symtab$symbol), keytype = "SYMBOL", columns = "ENTREZID"))
      a[, human_gene := KEY(SYMBOL)]
      fb <- merge(missing[, .(human_gene, mouse_gene)],
                  a[, .(human_gene, official = KEY(SYMBOL), ENTREZID)],
                  by = "human_gene", allow.cartesian = TRUE)
      z <- unique(rbindlist(list(z, fb), fill = TRUE))
    }
  }
  if (!nrow(z)) return(cv)
  z <- unique(z[!is.na(ENTREZID)])
  z <- merge(z, gl[, .(official = mouse_gene, repaired_length = length)],
             by = "official", all.x = TRUE)
  z <- z[is.finite(repaired_length)]
  if (anyDuplicated(z$human_gene)) stop("Ambiguous official gene alias for length repair")
  cv[z, on = "human_gene", length := i.repaired_length]
  cv[, complete := is.finite(length) & length > 0 & is.finite(logCPM)]
  cv
}
COV[["Astrocyte"]] <- repair_length_aliases(COV[["Astrocyte"]])
write_table(rbindlist(lapply(LINEAGES, function(ln)
  COV[[ln]][, .(lineage = ln, human_gene, mouse_gene, length, logCPM, complete)])),
  "crosscohort_allglia_five_axes_gene_length_expression_covariate_coverage")


#gene-property sensitivity gene_property_sensitivities: length/expression-controlled enrichment
l1_mouse_hallmarks <- readRDS(file.path(mouse_dir, "06_baiermouse_alllineages_hallmark_gene_sets.rds"))
l1_mouse_hallmarks <- lapply(l1_mouse_hallmarks, function(x) unique(KEY(x)))
l1_module_mouse <- fread(module_file)
l1_module_mouse[, gene := KEY(gene)]
l1_mouse_modules <- split(l1_module_mouse$gene, l1_module_mouse$module)

#Construct human modules from human MitoCarta/MSigDB membership
l1_mc_file <- hca_human_mitocarta_file
l1_mc_sheet <- grep("^A Human", readxl::excel_sheets(l1_mc_file), value = TRUE)[1]
l1_mc <- as.data.table(readxl::read_excel(l1_mc_file, sheet = l1_mc_sheet))
l1_mc_path <- grep("MitoPathways$", names(l1_mc), value = TRUE)[1]
l1_mc <- l1_mc[!is.na(get(l1_mc_path)) & get(l1_mc_path) != "0",
               .(gene = KEY(Symbol), path = as.character(get(l1_mc_path)))]
l1_mc <- l1_mc[, .(path = trimws(unlist(strsplit(path, "\\|")))), by = gene]
l1_mito_patterns <- list(
  RESP_CHAIN = "^OXPHOS", TCA = "TCA cycle", PYRUVATE = "Pyruvate",
  BETA_OX = "Fatty acid oxidation|Carnitine", KETONE = "Ketone",
  MITO_TRANSL = "central dogma|ribosome|Translation",
  MITO_IMPORT = "Protein import")
l1_mito_sets <- lapply(l1_mito_patterns, function(p)
  unique(l1_mc[grepl(p, path), gene]))
l1_reactome <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                               subcollection = "CP:REACTOME")
l1_reactome <- split(l1_reactome$gene_symbol, l1_reactome$gs_name)
l1_go <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5",
                         subcollection = "GO:BP")
l1_go <- split(l1_go$gene_symbol, l1_go$gs_name)
l1_kegg <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                           subcollection = "CP:KEGG_LEGACY")
l1_kegg <- split(l1_kegg$gene_symbol, l1_kegg$gs_name)
l1_pick <- function(collection, patterns) {
  hits <- unique(unlist(lapply(patterns, grep, x = names(collection), value = TRUE)))
  unique(KEY(unlist(collection[hits], use.names = FALSE)))
}
l1_nonmito <- list(
  PEROX_LIPID = unique(c(l1_pick(l1_reactome, "^REACTOME_PEROXISOMAL_LIPID_METABOLISM$"),
                         l1_pick(l1_go, "BETA_OXIDATION_USING_ACYL_COA_OXIDASE"))),
  PEROX_BIOGEN = l1_pick(l1_reactome, "^REACTOME_PEROXISOMAL_PROTEIN_IMPORT$"),
  LIPID_STORE = unique(c(l1_pick(l1_reactome, "^REACTOME_TRIGLYCERIDE_BIOSYNTHESIS$"),
                         l1_pick(l1_go, c("^GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT$",
                                          "^GOBP_LIPID_DROPLET_ORGANIZATION$")))),
  GLYCOLYSIS = unique(c(l1_pick(l1_reactome, "^REACTOME_GLYCOLYSIS$"),
                        l1_pick(l1_kegg, "^KEGG_GLYCOLYSIS_GLUCONEOGENESIS$"))),
  CHOLESTEROL = l1_pick(l1_reactome, "^REACTOME_CHOLESTEROL_BIOSYNTHESIS$"),
  ROS_DETOX = l1_pick(l1_reactome,
                      "^REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES$"),
  GLUTATHIONE = unique(c(l1_pick(l1_reactome, "^REACTOME_GLUTATHIONE_CONJUGATION$"),
                         l1_pick(l1_kegg, "^KEGG_GLUTATHIONE_METABOLISM$"))))
l1_mito_members <- rbindlist(lapply(names(l1_mito_sets), function(k)
  data.table(gene = l1_mito_sets[[k]], module = k,
             priority = match(k, names(l1_mito_sets)))))
setorder(l1_mito_members, gene, priority)
l1_mito_members <- unique(l1_mito_members, by = "gene")
l1_other_members <- rbindlist(lapply(names(l1_nonmito), function(k)
  data.table(gene = setdiff(l1_nonmito[[k]], l1_mito_members$gene),
             module = k)))
l1_human_members <- rbindlist(list(
  l1_mito_members[, .(gene, module)], l1_other_members))
l1_human_modules <- split(l1_human_members$gene, l1_human_members$module)

#Source gene-property sensitivity covariates: native mouse-gene span and mouse astrocyte
#cognition-model logCPM
l1_mouse_expr <- copy(mouse_cog[lineage == MOUSE_NAMES[["Astrocyte"]],
                                  .(gene = KEY(gene), logCPM)])
l1_mouse_expr <- l1_mouse_expr[, .(logCPM = median(logCPM)), by = gene]
l1_cov_mouse <- merge(gl[, .(gene = mouse_gene, length)], l1_mouse_expr,
                      by = "gene")
l1_cov_human <- COV[["Astrocyte"]][complete == TRUE,
                                      .(gene = human_gene, length, logCPM)]
l1_sets <- list(
  mouse = c(l1_mouse_modules, l1_mouse_hallmarks,
            list(CYTO_RIBOSOME = grep("^RP[SL][0-9]+[AXY]?[0-9]*$",
                                        l1_cov_mouse$gene, value = TRUE))),
  human = c(l1_human_modules, hallmarks,
            list(CYTO_RIBOSOME = grep("^RP[SL][0-9]+[AXY]?[0-9]*$",
                                        l1_cov_human$gene, value = TRUE))))
l1_sets <- lapply(l1_sets, function(s) s[lengths(s) >= 8L])
write_table(rbindlist(lapply(names(l1_sets), function(species)
  rbindlist(lapply(names(l1_sets[[species]]), function(pathway)
    data.table(species = species, pathway = pathway,
               gene = unique(l1_sets[[species]][[pathway]])))))),
  "crosscohort_astrocyte_gene_property_species_specific_program_membership")

l1_nat <- rbindlist(lapply(AXES, function(axis) {
  z <- native[[paste("Astrocyte", axis, sep = "__")]]
  data.table(axis_id = axis, species = z$species,
             native_gene = KEY(z$native_gene),
             oriented_logFC = z$oriented_logFC,
             oriented_signed_sqrtF = sign(z$oriented_logFC) *
               sqrt(pmax(z$F, 0)))
}))
l1_cov <- list(mouse = l1_cov_mouse, human = l1_cov_human)
l1_binz <- function(stat, group) {
  ave(stat, group, FUN = function(v) {
    if (length(v) > 2L && sd(v) > 0) {
      (v - mean(v)) / sd(v)
    } else {
      v * 0
    }
  })
}
l1_qb <- function(x, n) cut(rank(x, ties.method = "first"), n, labels = FALSE)
l1_gsea <- function(stat, genes, pathways, seed) {
  ranks <- sort(setNames(stat, genes)[is.finite(stat)], decreasing = TRUE)
  set.seed(seed)
  ranks <- ranks[!duplicated(names(ranks))]
  as.data.table(suppressWarnings(fgsea::fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 8L, maxSize = 1000L,
    eps = 0, BPPARAM = BiocParallel::SerialParam())))[
      , .(pathway, NES, pval, padj, size)]
}

set.seed(20261003)
l1 <- rbindlist(lapply(split(l1_nat, by = "axis_id"), function(x) {
  species <- x$species[1]
  d <- merge(x, l1_cov[[species]], by.x = "native_gene", by.y = "gene")
  d <- d[is.finite(oriented_signed_sqrtF) & is.finite(length) & length > 0 &
           is.finite(logCPM)]
  d[, `:=`(bin_length = l1_qb(log10(length), 10L),
           bin_expression = l1_qb(logCPM, 10L),
           bin_joint = paste(l1_qb(log10(length), 5L),
                             l1_qb(logCPM, 5L)))]
  rankings <- list(
    raw = d$oriented_signed_sqrtF,
    length = l1_binz(d$oriented_signed_sqrtF, d$bin_length),
    expression = l1_binz(d$oriented_signed_sqrtF, d$bin_expression),
    length_x_expression = l1_binz(d$oriented_signed_sqrtF, d$bin_joint))
  rbindlist(lapply(names(rankings), function(k) {
    l1_seed <- as.integer(20261003L +
      100L * match(x$axis_id[1], AXES) +
      match(k, c("raw", "length", "expression", "length_x_expression")))
    ans <- l1_gsea(rankings[[k]], d$native_gene, l1_sets[[species]], l1_seed)
    ans[, `:=`(axis_id = x$axis_id[1], species = species, ranking = k,
              n_covariate_genes = nrow(d),
              expression_proxy = species == "human")]
    ans
  }))
}))
write_table(l1, "crosscohort_astrocyte_gene_property_length_expression_binned_enrichment")

l1_key <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "RESP_CHAIN", "TCA",
            "BETA_OX", "HALLMARK_FATTY_ACID_METABOLISM",
            "HALLMARK_ADIPOGENESIS", "HALLMARK_PEROXISOME", "PEROX_LIPID",
            "GLYCOLYSIS", "CHOLESTEROL", "ROS_DETOX", "GLUTATHIONE",
            "MITO_IMPORT", "MITO_TRANSL", "HALLMARK_MYC_TARGETS_V1",
            "CYTO_RIBOSOME", "HALLMARK_INTERFERON_ALPHA_RESPONSE")
l1_wide <- dcast(l1[pathway %in% l1_key], pathway + axis_id ~ ranking,
                 value.var = c("NES", "padj"))
write_table(l1_wide, "crosscohort_astrocyte_gene_property_length_expression_binned_key_sets")
l1_recur <- l1[ranking == "length_x_expression",
               .(n_axes_FDR05 = sum(NES > 0 & padj < .05),
                 n_favorable = sum(NES > 0)), by = pathway]
write_table(l1_recur, "crosscohort_astrocyte_gene_property_recurrence_after_joint_binning")

#Gene-property sensitivity: 2,000 per-gene matches, 15 nearest background
# genes in normalized log10(length) and mouse-astrocyte logCPM space
l1_matched <- function(lfc, covariates, members, nperm = 2000L) {
  d <- merge(data.table(gene = names(lfc), lfc = unname(lfc)),
             covariates, by = "gene")[is.finite(lfc)]
  d[, `:=`(zl = as.numeric(scale(log10(length))),
           ze = as.numeric(scale(logCPM)))]
  selected <- d[gene %in% members]
  control <- d[!gene %in% members]
  if (nrow(selected) < 8L) return(NULL)
  nearest <- lapply(seq_len(nrow(selected)), function(i) head(order(
    (control$zl - selected$zl[i])^2 +
      (control$ze - selected$ze[i])^2), 15L))
  null <- replicate(nperm, median(control$lfc[vapply(nearest,
    function(ix) ix[sample.int(length(ix), 1L)], integer(1))]))
  observed <- median(selected$lfc)
  data.table(n = nrow(selected), obs_median = observed,
    null_median = median(null), excess = observed - median(null),
    frac_fav_set = mean(selected$lfc > 0),
    frac_fav_controls = mean(control$lfc[unique(unlist(nearest))] > 0),
    p_emp_upper = (1 + sum(null >= observed)) / (1 + nperm),
    p_emp_lower = (1 + sum(null <= observed)) / (1 + nperm))
}
l2 <- rbindlist(lapply(split(l1_nat, by = "axis_id"), function(x) {
  species <- x$species[1]
  effects <- setNames(x$oriented_logFC, x$native_gene)
  rbindlist(lapply(intersect(l1_key, names(l1_sets[[species]])), function(k) {
    ans <- l1_matched(effects, l1_cov[[species]], l1_sets[[species]][[k]])
    if (is.null(ans)) return(NULL)
    ans[, `:=`(axis_id = x$axis_id[1], species = species, pathway = k)]
    ans
  }), fill = TRUE)
}), fill = TRUE)
write_table(l2, "crosscohort_astrocyte_gene_property_15nearest_2000matched_gene_controls")

#Gene-property sensitivity: the same matched background sensitivity applied to
# the four mouse within-state cognition/aging contrasts from Script 08
l3_state <- fread(file.path(root, "08_astrocyte_states",
                            "08_within_state_cognition_aging_gene_effects.csv"))
l3 <- rbindlist(lapply(split(l3_state, by = c("state", "comparison")), function(x) {
  effects <- setNames(x$logFC_favorable, KEY(x$gene))
  rbindlist(lapply(intersect(l1_key, names(l1_sets$mouse)), function(k) {
    ans <- l1_matched(effects, l1_cov_mouse, l1_sets$mouse[[k]])
    if (is.null(ans)) return(NULL)
    ans[, `:=`(state = x$state[1], comparison = x$comparison[1],
              pathway = k)]
    ans
  }), fill = TRUE)
}), fill = TRUE)
write_table(l3, "crosscohort_astrocyte_gene_property_within_state_15nearest_matched_controls")

choose_strata <- function(cv) {
  bins <- function(v, k) pmin(k,
    ceiling(rank(v, ties.method = "first") * k / length(v)))
  for (k in c(3L, 2L)) {
    strata <- interaction(bins(log10(cv$length), k),
                          bins(cv$logCPM, k), drop = TRUE)
    if (all(table(strata) >= 3L))
      return(list(groups = split(seq_len(nrow(cv)), strata, drop = TRUE),
                  stratum = as.character(strata), k = k))
  }
  stop("Cannot form valid gene length/expression bins with >=3 genes. ",
       "Do not silently weaken the gene-matching matching specification.")
}
ox_cov <- COV[["Astrocyte"]][match(O_GENES, human_gene)]
if (!identical(ox_cov$human_gene, O_GENES) || any(!ox_cov$complete))
  stop("All 185 OXPHOS genes need verified length and logCPM. ",
       "Missing covariates: ", paste(ox_cov[complete != TRUE, human_gene],
                                     collapse = ", "))
ox_strata <- choose_strata(ox_cov)
if (ox_strata$k != 3L || length(ox_strata$groups) != 9L)
  stop("Final ED8c requires nine 3x3 length/expression strata")
write_table(data.table(human_gene = O_GENES, mouse_gene = ox_cov$mouse_gene,
                       length = ox_cov$length, logCPM = ox_cov$logCPM,
                       stratum = ox_strata$stratum),
            "crosscohort_astrocyte_recurrence_oxphos_185_verified_matching_strata")

#Conditional paired-cohort gene-identity permutations
metrics <- function(d, le) c(
  direction_all5 = sum(rowSums(d) == 5L),
  recurrent_ge4 = sum(rowSums(d) == 5L & rowSums(le) >= 4L),
  recurrent_le5 = sum(rowSums(d) == 5L & rowSums(le) == 5L))
permute_identity <- function(x, groups, n_perm, seed, mode) {
  stopifnot(mode %in% c("paired_cohort_blocks", "independent_axes"))
  d <- x$direction; le <- x$leading_edge
  observed <- metrics(d, le)
  n <- nrow(d)
  set.seed(seed)
  perm_id <- function() {
    ids <- seq_len(n)
    for (gg in groups) ids[gg] <- sample(gg, length(gg), replace = FALSE)
    ids
  }
  null <- matrix(NA_integer_, nrow = n_perm, ncol = 3L,
                 dimnames = list(NULL, names(observed)))
  for (i in seq_len(n_perm)) {
    dp <- d; lp <- le
    if (mode == "paired_cohort_blocks") {
      h <- perm_id(); s <- perm_id()
      dp[, 3] <- d[h, 3]; lp[, 3] <- le[h, 3]
      dp[, 4:5] <- d[s, 4:5, drop = FALSE]
      lp[, 4:5] <- le[s, 4:5, drop = FALSE]
    } else {
      for (j in 2:5) {
        z <- perm_id()
        dp[, j] <- d[z, j]; lp[, j] <- le[z, j]
      }
    }
    null[i, ] <- metrics(dp, lp)
  }
  data.table(metric = names(observed), observed = as.integer(observed),
    null_mean = colMeans(null), null_sd = apply(null, 2, stats::sd),
    null_q025 = apply(null, 2, stats::quantile, probs = .025),
    null_q975 = apply(null, 2, stats::quantile, probs = .975),
    empirical_p_upper = (1 + colSums(sweep(null, 2, as.integer(observed), `>=`))) /
      (n_perm + 1L), n_genes = n, null_model = mode,
    length_expression_matched = length(groups) > 1L,
    n_permutations = n_perm)
}

core_null <- rbindlist(list(
  permute_identity(ox, ox_strata$groups, N_CORE_PERM, 20261026L,
                   "independent_axes"),
  permute_identity(ox, ox_strata$groups, N_CORE_PERM, 20261027L,
                   "paired_cohort_blocks")))
write_table(core_null, "crosscohort_astrocyte_recurrence_oxphos_full185_length_expression_matched_null")

#Compare OXPHOS with competitors selected from the integrated Hallmark ranking  
selected <- rbindlist(lapply(LINEAGES, function(ln) {
  z <- integrated[lineage == ln & pathway != OX][order(preserved_rank, pathway)]
  if (nrow(z) < 2L) stop("Missing non-OXPHOS competitors for ", ln)
  data.table(lineage = ln, pathway = c(OX, z$pathway[1:2]),
    selection = c("OXPHOS reference", "Integrated competitor #1",
                  "Integrated competitor #2"))
}))
for (p in c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_FATTY_ACID_METABOLISM"))
  if (!nrow(selected[lineage == "Astrocyte" & pathway == p]))
    selected <- rbind(selected, data.table(lineage = "Astrocyte", pathway = p,
      selection = "Prespecified astrocyte alternative"))
if (!nrow(selected[lineage == "Microglia" & pathway ==
                  "HALLMARK_INTERFERON_ALPHA_RESPONSE"]))
  selected <- rbind(selected, data.table(lineage = "Microglia",
    pathway = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    selection = "Biological comparator, not selected by integrated rank"))

# Validate the selected comparison family without assuming a fixed size.
stopifnot(
  !anyDuplicated(selected, by = c("lineage", "pathway")),
  setequal(unique(selected$lineage), LINEAGES),
  all(selected[, .N, by = lineage]$N >= 3L)
)

message(
  "ED8b: ", nrow(selected),
  " distinct lineage/pathway comparisons selected."
)

print(selected[, .(lineage, pathway, selection)])

selected <- merge(selected, integrated[, .(lineage, pathway,
  preserved_rank, integrated_NES = NES, integrated_FDR = padj)],
  by = c("lineage", "pathway"), all.x = TRUE)
write_table(selected, "crosscohort_allglia_five_axes_selected_pathways_from_unbiased_integration")

selected_null <- rbindlist(lapply(seq_len(nrow(selected)), function(i) {
  ln <- selected$lineage[i]; path <- selected$pathway[i]
  x <- get_boolean(masters[[ln]], ln, path)
  cv <- COV[[ln]][match(x$genes, human_gene)]
  if (!identical(cv$human_gene, x$genes)) stop("Covariates not aligned: ", ln, "/", path)
  rbindlist(lapply(c(FALSE, TRUE), function(matched) {
    idx <- if (matched) which(cv$complete) else seq_along(x$genes)
    y <- x
    y$genes <- y$genes[idx]
    y$direction <- y$direction[idx, , drop = FALSE]
    y$leading_edge <- y$leading_edge[idx, , drop = FALSE]
    if (length(idx) < 12L) return(data.table(lineage = ln, pathway = path,
       matched = matched, reason = "Fewer than 12 genes; not tested"))
    z <- tryCatch(
      if (matched) choose_strata(cv[idx])$groups else list(all = seq_along(idx)),
      error = function(e) NULL)
    if (is.null(z)) return(data.table(lineage = ln, pathway = path,
       matched = matched, reason = "No valid matched strata"))
    rbindlist(lapply(c("paired_cohort_blocks", "independent_axes"), function(mode) {
      k <- permute_identity(y, z, N_SELECTED_PERM,
           SEED_SELECTED + 1000L * i + 10L * as.integer(matched) +
             if (mode == "paired_cohort_blocks") 1L else 2L, mode)
      k[, `:=`(lineage = ln, pathway = path, matched = matched,
               n_omitted = length(x$genes) - length(idx), reason = "OK")]
      k
    }))
  }), fill = TRUE)
}), fill = TRUE)
selected_null[!is.na(empirical_p_upper),
   BH_selected := p.adjust(empirical_p_upper, "BH"),
   by = .(matched, null_model, metric)]
write_table(selected_null, "crosscohort_allglia_five_axes_selected_pathway_conditional_nulls")

#Six-context metabolic/redox matched-gene assessment
MODULE_NAMES <- c(OX, "RESP_CHAIN", "TCA", "HALLMARK_FATTY_ACID_METABOLISM",
                  "BETA_OX", "PEROX_LIPID", "ROS_DETOX", "GLUTATHIONE",
                  "GLYCOLYSIS")
FAMILY <- c("Respiration", "Respiration", "Respiration",
            "Lipid metabolism", "Lipid metabolism", "Lipid metabolism",
            "Redox / antioxidant", "Redox / antioxidant", "Other")
module_raw <- fread(module_file)
check_cols(module_raw, c("module", "gene"), "Script 08 metabolic modules")
if ("species" %in% names(module_raw)) module_raw <- module_raw[tolower(species) == "mouse"]
module_raw[, mouse_gene := KEY(gene)]
module_human <- merge(module_raw, orth, by = "mouse_gene", all = FALSE)
module_sets <- split(module_human$human_gene, module_human$module)
SETS <- lapply(MODULE_NAMES, function(k) {
  if (startsWith(k, "HALLMARK_")) hallmarks[[k]] else module_sets[[k]]
})
names(SETS) <- MODULE_NAMES
if (any(lengths(SETS) < 5L)) stop("Original metabolic modules missing: ",
    paste(names(SETS)[lengths(SETS) < 5L], collapse = ", "))
write_table(rbindlist(lapply(names(SETS), function(k)
  data.table(set = k, gene = unique(SETS[[k]])))), "crosscohort_allglia_five_axes_metabolic_gene_membership")

matched_program <- function(fav_lfc, cov, gene_set, n_perm = N_MATCHED) {
  dt <- data.table(human_gene = names(fav_lfc), lfc = unname(fav_lfc))
  dt <- merge(dt, cov[, .(human_gene, length, logCPM)], by = "human_gene")
  dt <- dt[is.finite(lfc) & is.finite(length) & length > 0 & is.finite(logCPM)]
  if (nrow(dt) < 100L) stop("Matched background too small")
  dt[, `:=`(zl = as.numeric(scale(log10(length))),
            ze = as.numeric(scale(logCPM)))]
  inside <- dt$human_gene %in% gene_set
  a <- dt[inside]; b <- dt[!inside]
  if (nrow(a) < 5L || nrow(b) < N_NEAREST) return(NULL)
  near <- lapply(seq_len(nrow(a)), function(i) head(order(
    (b$zl - a$zl[i])^2 + (b$ze - a$ze[i])^2), N_NEAREST))
  draws <- replicate(n_perm, {
    chosen <- vapply(near, function(v) v[sample.int(length(v), 1L)], integer(1))
    median(b$lfc[chosen])
  })
  ob <- median(a$lfc)
  data.table(n = nrow(a), obs_median = ob, matched_effect = ob - median(draws),
    empirical_p = (1 + sum(draws >= ob)) / (1 + n_perm),
    frac_favorable = mean(a$lfc > 0))
}
astro_cov <- COV[["Astrocyte"]]
set.seed(20261003)
metabolic_5 <- rbindlist(lapply(AXES, function(ax) {
  m <- masters[["Astrocyte"]]
  lfc <- setNames(m[[paste0(ax, "_oriented_logFC")]], m$human_gene)
  rbindlist(lapply(MODULE_NAMES, function(k) {
    z <- matched_program(lfc, astro_cov, SETS[[k]])
    if (is.null(z)) stop("Insufficient matched genes for ", ax, " / ", k)
    z[, `:=`(axis_id = ax, set = k)]
    z
  }))
}))
ros <- fread(rosmap_file)
check_cols(ros, c("path", "region", "pathway", "ranking",
                  "excess", "p_emp_upper"), "Script 11 ROSMAP matched results")
ros <- ros[path == "nrad" & region == "allregions" & ranking == "raw" &
           pathway %in% c(OX, "BETA_OX", "GLYCOLYSIS")]
if (nrow(ros) != 3L) stop("Expected exactly three ROSMAP assessed metabolic modules")
metabolic_ros <- ros[, .(set = pathway, axis_id = "rosmap_pathology",
                         matched_effect = excess, empirical_p = p_emp_upper,
                         n, obs_median, frac_favorable)]
metabolic <- rbindlist(list(metabolic_5, metabolic_ros), fill = TRUE)
if (nrow(metabolic) != 48L || anyDuplicated(metabolic, by = c("set", "axis_id")))
  stop("Matched metabolic analysis must contain 48 unique evaluated tests")
metabolic[, FDR_48 := p.adjust(empirical_p, method = "BH")]
full_grid <- CJ(set = MODULE_NAMES, axis_id = c(AXES, "rosmap_pathology"))
metabolic_full <- merge(full_grid, metabolic, by = c("set", "axis_id"), all.x = TRUE)
metabolic_full[, `:=`(family = FAMILY[match(set, MODULE_NAMES)],
  assessed = is.finite(matched_effect),
  significant = is.finite(FDR_48) & FDR_48 < .05)]
if (nrow(metabolic_full) != 54L || sum(metabolic_full$assessed) != 48L)
  stop("Six-context coverage must be exactly 48/54 assessed comparisons")
write_table(metabolic_full, "crosscohort_allglia_five_axes_metabolic_six_context_48test_source")

# Curated MitoCarta3.0 subunit/transport classes are biological reference memberships
mitocarta_source <- file.path(out, "12_crosscohort_allglia_five_axes_broad_mito_pathways3_human_149_curated_pathways.csv")
mitopathways_html <- hca_human_mitopathways_html
doc <- rvest::read_html(mitopathways_html)
tables <- rvest::html_table(rvest::html_elements(doc, "table"), fill = TRUE)
candidates <- which(vapply(tables, function(x)
  ncol(x) >= 3L && nrow(x) >= 100L, logical(1)))
stopifnot(length(candidates) == 1L)
broad <- as.data.table(tables[[candidates]][, 1:3])
setnames(broad, c("pathway", "hierarchy", "genes"))
fwrite(broad, mitocarta_source)
broad <- fread(mitocarta_source)
broad[, pathway := trimws(as.character(pathway))]
pathway_genes <- function(name) {
  row <- broad[tolower(pathway) == tolower(name)]
  stopifnot(nrow(row) == 1L)
  unique(KEY(unlist(strsplit(as.character(row$genes), ",", fixed = TRUE))))
}
class_names <- c("Complex I", "Complex II", "Complex III", "Complex IV",
                 "Complex V", "Transport/import")
class_sets <- setNames(lapply(c("CI subunits", "CII subunits", "CIII subunits",
                                 "CIV subunits", "CV subunits"), pathway_genes),
                       class_names[1:5])
class_sets[["Transport/import"]] <- unique(c(
  pathway_genes("Small molecule transport"),
  pathway_genes("Protein import and sorting")))
classes <- copy(ox_manifest)
for (nm in class_names) classes[, (paste0("mc_", gsub("[^A-Za-z0-9]", "_", nm))) :=
                                   human_gene %in% class_sets[[nm]]]
all_mito <- unique(KEY(unlist(strsplit(paste(broad$genes, collapse = ","),
                                       ",", fixed = TRUE))))
classes[, any_MitoCarta := human_gene %in% all_mito]
write_table(classes, "crosscohort_astrocyte_recurrence_oxphos_185_with_verified_mito_carta_classes")

# Conditional matched 34-versus-185 null: original 9 length/expression strata.
flags <- grep("^mc_", names(classes), value = TRUE)
stratum <- ox_strata$stratum
selected_idx <- which(classes$recurrent34)
selected_strata <- table(stratum[selected_idx])
pools <- split(seq_len(nrow(classes)), stratum)
set.seed(20261028L + 9100L)
N_CLASS_PERM <- 10000L
null_classes <- replicate(N_CLASS_PERM, {
  ix <- unlist(lapply(names(selected_strata), function(k)
    sample(pools[[k]], as.integer(selected_strata[[k]]), replace = FALSE)))
  colSums(as.matrix(classes[ix, ..flags]))
})
class_labels <- setNames(class_names, paste0("mc_", gsub("[^A-Za-z0-9]", "_", class_names)))
class_results <- rbindlist(lapply(seq_along(flags), function(j) {
  obs <- sum(classes[[flags[j]]][selected_idx])
  elsewhere <- sum(classes[[flags[j]]]) - obs
  fisher <- stats::fisher.test(matrix(c(obs, 34L - obs, elsewhere,
                                         151L - elsewhere), nrow = 2L,
                                       byrow = TRUE), alternative = "greater")
  data.table(class = unname(class_labels[flags[j]]), n_recurrent = obs,
    n_other = elsewhere, n185_total = obs + elsewhere,
    null_mean = mean(null_classes[j, ]),
    null_q025 = unname(quantile(null_classes[j, ], .025)),
    null_q975 = unname(quantile(null_classes[j, ], .975)),
    empirical_p = (1 + sum(null_classes[j, ] >= obs)) / (1 + N_CLASS_PERM),
    fisher_one_sided_p = fisher$p.value, odds_ratio = unname(fisher$estimate))
}))
class_results[, `:=`(BH_q_matched = p.adjust(empirical_p, "BH"),
                    BH_q_Fisher = p.adjust(fisher_one_sided_p, "BH"))]
write_table(class_results, "crosscohort_astrocyte_recurrence_oxphos_mito_carta_classes_matched_34_vs185")

#checkpoint
saveRDS(list(
  lineage_states = masters,
  hallmark_effect = integrated,
  axis_native = axis_results,
  lineage_pathway_recurrence = recurrence,
  integrated_rank_summary = rank_summary,
  cohort_block_omissions = block_loo,
  ortholog_map = orth,
  OXPHOS_185 = ox_manifest,
  OXPHOS_axis_effects = ox_effects,
  Zemke_OXPHOS_donor_model = model_summary,
  Zemke_OXPHOS_donor_scores = donor_score,
  OXPHOS_paired_matched_null = core_null,
  selected_conditional_nulls = selected_null,
  metabolic_six_context = metabolic_full,
  MitoCarta_classes = classes,
  MitoCarta_matched_class_results = class_results,
  human_hallmarks = hallmarks,
  length_expression_L1 = l1,
  matched_background_L2 = l2,
  within_state_matched_L3 = l3
), file.path(out, "12_crosscohort_allglia_five_axes_mouse_human_integration_checkpoint.rds"))

