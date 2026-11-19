#Baier et al.
#13: Figure source data

code_dir <- Sys.getenv("HCA_CODE_ROOT", unset = getwd())
source(file.path(code_dir, "00_paths.R"))
library(data.table)

root <- hca_output_root
out <- file.path(root, "13_figure_source_data")
panel_dir <- file.path(out, "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

status <- list()
input_log <- list()
source_cache <- new.env(parent = emptyenv())

input_file <- function(relative) {
  path <- file.path(root, relative)
  if (!file.exists(path)) return(NULL)
  if (exists(relative, envir = source_cache, inherits = FALSE))
    return(get(relative, envir = source_cache))
  value <- if (grepl("\\.rds$", relative, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    data.table::fread(path, showProgress = FALSE)
  }
  assign(relative, value, envir = source_cache)
  input_log[[relative]] <<- data.table(
    source = relative, bytes = file.info(path)$size,
    modified_utc = format(file.info(path)$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(path)))
  value
}

mkrow <- function(panel, component, source, selector, state, rows = NA_integer_,
                  columns = NA_integer_, file = "", message = "") {
  data.table(panel, component, source, selector, status = state, rows,
             columns, output_file = file, detail = message)
}

select_source <- function(x, selector) {
  x <- as.data.table(x)
  if (selector == "all" || selector == "rds_qc") return(x)
  if (selector == "mouse_cognition_nine")
    {
    z <- x[contrast == "Cog" & !lineage %in% c("Ependymal NN", "CHOR NN")]
    selected <- unique(z[padj < .05, pathway])
    z[, display_pathway_significant_in_any_lineage := pathway %in% selected]
    return(z)
  }
  if (selector == "mouse_nine")
    return(x[!lineage %in% c("Ependymal NN", "CHOR NN")])
  if (selector == "mouse_five")
    return(x[lineage %in% c("Astro-TE NN", "Microglia NN", "OPC NN",
                          "Oligo NN", "Peri NN")])
  if (selector == "endo_cog") return(x[lineage == "Endo NN" & contrast == "Cog"])
  if (selector == "human_aging_fig5a") {
    wanted <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
      "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
      "HALLMARK_PEROXISOME", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
      "HALLMARK_GLYCOLYSIS", "HALLMARK_DNA_REPAIR",
      "HALLMARK_MYC_TARGETS_V1", "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_G2M_CHECKPOINT")
    z <- x[crossstudy_celltype == "Astrocyte" & pathway %in% wanted]
    z[, display_order := match(pathway, wanted)]
    return(z[order(display_order)])
  }
  if (selector == "human_astro") {
    field <- intersect(c("cell_type", "crossstudy_celltype", "lineage"), names(x))[1L]
    if (is.na(field)) stop("No lineage column for astrocyte selection")
    return(x[tolower(as.character(get(field))) %in% c("astrocyte", "astro-te nn")])
  }
  if (selector == "aged_AST1")
    return(x[grepl("(^|[^a-z])AST1$|cAST1$", state, ignore.case = TRUE) &
               group %in% c("Intact", "Impaired")])
  if (selector == "top10") return(head(x[order(mean_rank)], 10L))
  if (selector == "recurrent34") {
    if (!"recurrent34" %in% names(x))
      stop("The source has no recurrent34 classification; cannot select 34 genes")
    return(x[as.logical(recurrent34)])
  }
  if (startsWith(selector, "path_"))
    return(x[pathway == substring(selector, 6L)])
  stop("Unknown panel selector: ", selector)
}

special_source <- function(selector) {
  H <- "06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv"
  I <- "12_mouse_human_integration/"
  if (selector == "f1g_count") {
    h <- input_file(H)
    if (is.null(h)) return(NULL)
    h <- h[contrast == "Cog" & !lineage %in% c("Ependymal NN", "CHOR NN")]
    return(h[, .(n_intact = sum(NES > 0 & padj < .05, na.rm = TRUE),
                 n_impaired = sum(NES < 0 & padj < .05, na.rm = TRUE)), by = lineage])
  }
  if (selector == "f2b_scatter") {
    h <- input_file(H)
    if (is.null(h)) return(NULL)
    h <- h[lineage %in% c("Astro-TE NN", "Microglia NN", "OPC NN", "Oligo NN", "Peri NN") &
             contrast %in% c("Cog", "AgeP")]
    return(dcast(h, lineage + pathway ~ contrast, value.var = c("NES", "padj")))
  }
  if (selector == "f3b_top12") {
    s <- input_file(paste0(I, "12_crosscohort_allglia_five_axes_integrated_cross_lineage_rank_summary.csv"))
    h <- input_file(paste0(I, "12_crosscohort_allglia_five_axes_cohort_balanced_hallmarks_all200.csv"))
    if (is.null(s) || is.null(h)) return(NULL)
    top <- s[order(mean_rank, worst_rank)][seq_len(min(12L, .N)), pathway]
    return(merge(h[pathway %in% top], s, by = "pathway", all.x = TRUE))
  }
  if (selector == "f3f_rank") {
    x <- input_file(paste0(I, "12_crosscohort_allglia_five_axes_mouse_versus_human_pathway_ranks.csv"))
    if (is.null(x)) return(NULL)
    x <- copy(x)
    if (!all(c("mouse_cognition", "human_only_rank") %in% names(x))) return(x)
    setorder(x, lineage, mouse_cognition)
    x[, show_panel := seq_len(.N) <= 5L |
        (grepl("Microglia", lineage, ignore.case = TRUE) &
         pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), by = lineage]
    return(x[show_panel == TRUE])
  }
  if (selector == "f3g_rank") {
    h <- input_file(paste0(I, "12_crosscohort_allglia_five_axes_cohort_balanced_hallmarks_all200.csv"))
    a <- input_file(paste0(I, "12_crosscohort_allglia_five_axes_hallmarks_all1000.csv"))
    if (is.null(h) || is.null(a)) return(NULL)
    best <- h[pathway != "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]
    setorder(best, lineage, preserved_rank)
    best <- best[, head(.SD, 2L), by = lineage]
    return(merge(a, best[, .(lineage, pathway, integrated_rank = preserved_rank)],
                 by = c("lineage", "pathway")))
  }
  if (startsWith(selector, "gsea_")) {
    de <- input_file("06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_gene_effects.csv")
    sets <- input_file("06_mouse_bulk/06_baiermouse_alllineages_hallmark_gene_sets.rds")
    if (is.null(de) || is.null(sets)) return(NULL)
    path <- substring(selector, 6L)
    genes <- unique(sets[[path]])
    if (!length(genes)) stop("Absent source Hallmark: ", path)
    lineages <- c("Astro-TE NN", "Microglia NN", "OPC NN", "Oligo NN", "Peri NN")
    return(rbindlist(lapply(lineages, function(ln) {
      z <- de[lineage == ln & contrast == "Cog", .(gene, s)]
      z <- z[is.finite(s) & !is.na(gene)][order(-s, gene)]
      z <- z[!duplicated(gene)]
      hits <- z$gene %in% genes
      n_hit <- sum(hits)
      if (!n_hit || n_hit == nrow(z)) return(data.table())
      w <- abs(z$s[hits])
      denom <- sum(w)
      if (denom <= 0) return(data.table())
      v <- numeric(nrow(z))
      v[hits] <- abs(z$s[hits]) / denom
      v[!hits] <- -1 / sum(!hits)
      data.table(lineage = ln, pathway = path, rank_position = seq_len(nrow(z)),
                 gene = z$gene, signed_sqrtF = z$s, in_gene_set = hits,
                 running_ES = cumsum(v))
    }), fill = TRUE))
  }
  stop("Unknown derived source: ", selector)
}

#Resolve panel data from the declared numbered-pipeline outputs
read_panel <- function(spec) {
  selector <- spec$selector
  source <- spec$source
  if (startsWith(selector, "gsea_") ||
      selector %in% c("f1g_count", "f2b_scatter", "f3b_top12", "f3f_rank", "f3g_rank"))
    return(special_source(selector))
  if (selector %in% c("atlas_umap", "atlas_marker_dot"))
    return(NULL)  # processed once below to avoid repeated 198k-cell atlas loads
  src <- input_file(source)
  if (is.null(src)) return(NULL)
  select_source(src, selector)
}

#Figure 1C/D and ED2H use the annotated atlas from Script 04
atlas_parts <- new.env(parent = emptyenv())
load_atlas_panel_data <- function() {
  atlas_path <- file.path(root, "objects", "04_baiermouse_allcells_annotated_atlas.rds")
  if (!file.exists(atlas_path)) return(FALSE)
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Install Seurat to extract Figure 1C/D source data")
  a <- readRDS(atlas_path)
  md <- as.data.table(a[[]], keep.rownames = "cell_id")
  u <- as.data.table(Seurat::Embeddings(a[["umap"]]), keep.rownames = "cell_id")
  xy <- names(u)[-1L][seq_len(2L)]
  setnames(u, xy, c("UMAP_1", "UMAP_2"))
  fields <- intersect(c("cell_id", "mouse_id", "primary_group", "sequence_batch",
                        "broad_lineage", "analysis_celltype", "allen_subclass"), names(md))
  u <- merge(u[, .(cell_id, UMAP_1, UMAP_2)], md[, ..fields], by = "cell_id", sort = FALSE)
  assign("atlas_umap", u, envir = atlas_parts)

  #Canonical marker identities and calculation follow figure source methodology.
  markers <- list(
    Astrocyte = c("Aqp4", "Aldoc"), Microglia = c("P2ry12", "Tmem119"),
    OPC = c("Pdgfra", "Cspg4"), Oligodendrocyte = c("Mog", "Plp1"),
    `DG neuron` = c("Prox1", "Slc17a7"),
    Endothelial = c("Cldn5", "Pecam1"), Pericyte = c("Pdgfrb", "Rgs5"),
    BAM = c("Mrc1", "Pf4"), Lymphoid = c("Cd3e", "Nkg7", "Cd79a"))
  gene_names <- unique(unlist(markers, use.names = FALSE))
  if (!all(gene_names %in% rownames(a[["RNA"]])))
    stop("Figure 1D marker genes absent from annotated RNA assay")
  md[, plot_lineage := as.character(analysis_celltype)]
  md[plot_lineage == "DG glutamatergic neuron", plot_lineage := "DG neuron"]
  relevant <- md[plot_lineage %in% names(markers), cell_id]
  expr <- Seurat::FetchData(a, vars = gene_names, cells = relevant, layer = "data")
  ids <- rownames(expr)
  groups <- md$plot_lineage[match(ids, md$cell_id)]
  dots <- rbindlist(lapply(names(markers), function(ln) {
    x <- as.matrix(expr[groups == ln, gene_names, drop = FALSE])
    if (!nrow(x)) return(data.table())
    data.table(cell_type = ln, gene = gene_names,
               pct_exp = 100 * colMeans(x > 0),
               avg_exp = log1p(colMeans(expm1(x))), n_cells = nrow(x))
  }))
  dots[, avg_scaled := if (stats::sd(avg_exp) > 0) as.numeric(scale(avg_exp))
       else rep(0, .N), by = gene]
  dots[, avg_scaled := pmax(-2.5, pmin(2.5, avg_scaled))]
  assign("atlas_marker_dot", dots, envir = atlas_parts)
  input_log[["stage04_atlas_for_visualization"]] <<- data.table(
    source = normalizePath(atlas_path, winslash = "/", mustWork = TRUE),
    bytes = file.info(atlas_path)$size,
    modified_utc = format(file.info(atlas_path)$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(atlas_path)))
  rm(a, md, expr, dots, u)
  invisible(gc())
  TRUE
}

#Map each main and extended figure component to its upstream source
panel_specs <- data.table::rbindlist(list(
  list(panel="F1A", component="study_design", source="01_mouse_behavior/01_baiermouse_allmice_groups_behavior_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="F1A", component="sequencing_cohort", source="06_mouse_bulk/06_baiermouse_allmice_pseudobulk_sample_metadata.csv", selector="all", required=TRUE, note=""),
  list(panel="F1B", component="individual_behavior", source="01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv", selector="all", required=TRUE, note=""),
  list(panel="F1C", component="umap", source="", selector="atlas_umap", required=TRUE, note="Extract embedding from stage-04 atlas"),
  list(panel="F1D", component="marker_detection", source="", selector="atlas_marker_dot", required=TRUE, note="Calculate displayed marker detection/expression from normalized atlas"),
  list(panel="F1E", component="mouse_composition", source="05_differential_abundance/05_baiermouse_alllineages_allgroups_proportions_by_mouse.csv", selector="all", required=TRUE, note=""),
  list(panel="F1F", component="cognition_hallmarks", source="06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv", selector="mouse_cognition_nine", required=TRUE, note=""),
  list(panel="F1G", component="cognition_hallmarks_significant_counts", source="", selector="f1g_count", required=TRUE, note=""),
  list(panel="F1H", component="gsea_running_es", source="", selector="gsea_HALLMARK_OXIDATIVE_PHOSPHORYLATION", required=TRUE, note="Source GSEA curve computed from frozen mouse gene effects; enrichment P/FDR sourced from Script 06"),
  list(panel="F1H", component="lineage_enrichment", source="06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv", selector="path_HALLMARK_OXIDATIVE_PHOSPHORYLATION", required=TRUE, note=""),
  list(panel="F1I", component="gsea_running_es", source="", selector="gsea_HALLMARK_INTERFERON_ALPHA_RESPONSE", required=TRUE, note="Source GSEA curve computed from frozen mouse gene effects; enrichment P/FDR sourced from Script 06"),
  list(panel="F1I", component="lineage_enrichment", source="06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv", selector="path_HALLMARK_INTERFERON_ALPHA_RESPONSE", required=TRUE, note=""),
  list(panel="F1J", component="gsea_running_es", source="", selector="gsea_HALLMARK_E2F_TARGETS", required=TRUE, note="Source GSEA curve computed from frozen mouse gene effects; enrichment P/FDR sourced from Script 06"),
  list(panel="F1J", component="lineage_enrichment", source="06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv", selector="path_HALLMARK_E2F_TARGETS", required=TRUE, note=""),
  list(panel="F2A", component="design", source="06_mouse_bulk/06_baiermouse_allmice_pseudobulk_sample_metadata.csv", selector="all", required=TRUE, note=""),
  list(panel="F2B", component="aging_cognition_NES", source="", selector="f2b_scatter", required=TRUE, note=""),
  list(panel="F2C", component="pathway_patterns", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_pathway_patterns.csv", selector="mouse_nine", required=TRUE, note=""),
  list(panel="F2D", component="contrasts", source="06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv", selector="mouse_five", required=TRUE, note=""),
  list(panel="F2E", component="individual_mouse_program_scores", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_program_scores_by_mouse.csv", selector="mouse_five", required=TRUE, note=""),
  list(panel="F2E", component="group_program_summary", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_program_group_means.csv", selector="mouse_five", required=TRUE, note=""),
  list(panel="F2F", component="individual_mouse_program_scores", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_program_scores_by_mouse.csv", selector="mouse_five", required=TRUE, note=""),
  list(panel="F2F", component="group_program_summary", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_program_group_means.csv", selector="mouse_five", required=TRUE, note=""),
  list(panel="F2G", component="individual_mouse_program_scores", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_program_scores_by_mouse.csv", selector="mouse_five", required=TRUE, note=""),
  list(panel="F2G", component="group_program_summary", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_program_group_means.csv", selector="mouse_five", required=TRUE, note=""),
  list(panel="F3A", component="study_blocks", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_gene_coverage.csv", selector="all", required=TRUE, note=""),
  list(panel="F3B", component="integrated_top12", source="", selector="f3b_top12", required=TRUE, note=""),
  list(panel="F3C", component="cross_lineage_rank", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_integrated_cross_lineage_rank_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="F3D", component="OXPHOS_five_axis", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_hallmarks_all1000.csv", selector="path_HALLMARK_OXIDATIVE_PHOSPHORYLATION", required=TRUE, note=""),
  list(panel="F3E", component="cohort_block_ranks", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_mouse_versus_human_pathway_ranks.csv", selector="path_HALLMARK_OXIDATIVE_PHOSPHORYLATION", required=TRUE, note=""),
  list(panel="F3F", component="mouse_human_rank_comparison", source="", selector="f3f_rank", required=TRUE, note=""),
  list(panel="F3G", component="selected_nonOXPHOS_five_axis", source="", selector="f3g_rank", required=TRUE, note=""),
  list(panel="F3H", component="scaled_NES", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_mouse_human_scaled_nes_scatter_source.csv", selector="all", required=TRUE, note=""),
  list(panel="F3H", component="four_lineage_correlations", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_mouse_human_scaled_nes_correlations.csv", selector="all", required=TRUE, note=""),
  list(panel="F4A", component="whole_astrocyte_GO", source="06_mouse_bulk/06_baiermouse_alllineages_cognition_twogroup_gobp_enrichment.csv", selector="all", required=TRUE, note=""),
  list(panel="F4B", component="astrocyte_umap", source="08_astrocyte_states/08_UMAP_and_respiratory_UCell.csv", selector="all", required=TRUE, note=""),
  list(panel="F4B", component="state_counts", source="08_astrocyte_states/08_state_counts_fractions_by_mouse.csv", selector="all", required=TRUE, note=""),
  list(panel="F4C", component="identity_GO", source="08_astrocyte_states/08_state_identity_GO_BP_enrichment.csv", selector="all", required=TRUE, note=""),
  list(panel="F4D", component="umap_UCell", source="08_astrocyte_states/08_UMAP_and_respiratory_UCell.csv", selector="all", required=TRUE, note=""),
  list(panel="F4E", component="mouse_state_fractions", source="08_astrocyte_states/08_state_counts_fractions_by_mouse.csv", selector="all", required=TRUE, note=""),
  list(panel="F4E", component="state_means", source="08_astrocyte_states/08_state_abundance_group_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="F4E", component="state_abundance_FDR", source="08_astrocyte_states/08_state_abundance_seven_test_models.csv", selector="all", required=TRUE, note=""),
  list(panel="F4F", component="aged_mouse_AST1", source="08_astrocyte_states/08_state_counts_fractions_by_mouse.csv", selector="aged_AST1", required=TRUE, note=""),
  list(panel="F4F", component="correlation_and_permutation", source="08_astrocyte_states/08_AST1_cognitive_flexibility.csv", selector="all", required=TRUE, note=""),
  list(panel="F4G", component="paired_mouse_state_scores", source="08_astrocyte_states/08_Figure4g_three_programs_mouse_state.csv", selector="all", required=TRUE, note=""),
  list(panel="F4H", component="respiratory_complex_genes", source="08_astrocyte_states/08_respiratory_complex_gene_effects.csv", selector="all", required=TRUE, note=""),
  list(panel="F4H", component="respiratory_complex_summary", source="08_astrocyte_states/08_respiratory_complex_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="F5A", component="astrocyte_age_hallmarks", source="09_human_normative_aging/09_zemke_alllineages_aging_hallmark_younger_age_inference.csv", selector="human_aging_fig5a", required=TRUE, note=""),
  list(panel="F5B", component="donor_oxphos_scores", source="12_mouse_human_integration/12_zemke_astrocyte_aging_astrocyte_oxphos_donor_scores.csv", selector="all", required=TRUE, note=""),
  list(panel="F5B", component="sex_adjusted_model", source="12_mouse_human_integration/12_zemke_astrocyte_aging_astrocyte_oxphos_age_model.csv", selector="all", required=TRUE, note=""),
  list(panel="F5C", component="pathology_vs_cognition", source="10_human_SEAAD/10_seaad_alllineages_pathology_cognition_hallmark_cognition_vs_pathology.csv", selector="human_astro", required=TRUE, note=""),
  list(panel="F5D", component="regional_metabolic_NES", source="11_human_ROSMAP/11_rosmap_astrocyte_pathology_regions_fig5_ed7_source.csv", selector="all", required=TRUE, note=""),
  list(panel="F5E", component="matched_background_48tests", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_metabolic_six_context_48test_source.csv", selector="all", required=TRUE, note=""),
  list(panel="F5F", component="185_gene_recurrence_manifest", source="12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_final_185_gene_manifest.csv", selector="all", required=TRUE, note=""),
  list(panel="F5G", component="five_axis_gene_percentiles", source="12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_185_five_axis_effects_and_leading_edges.csv", selector="all", required=TRUE, note=""),
  list(panel="F5H", component="recurrent34_five_axis", source="12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_185_five_axis_effects_and_leading_edges.csv", selector="recurrent34", required=TRUE, note=""),
  list(panel="F5H", component="mitocarta_classification", source="12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_185_with_verified_mito_carta_classes.csv", selector="recurrent34", required=TRUE, note=""),
  list(panel="ED1A", component="behavior_observations", source="01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1B", component="behavior_observations", source="01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1B", component="spearman_statistics", source="01_mouse_behavior/01_baiermouse_allmice_behavior_spearman_correlations.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1C", component="behavior_observations", source="01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1C", component="spearman_statistics", source="01_mouse_behavior/01_baiermouse_allmice_behavior_spearman_correlations.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1D", component="behavior_observations", source="01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1D", component="spearman_statistics", source="01_mouse_behavior/01_baiermouse_allmice_behavior_spearman_correlations.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1E", component="paired_body_weights", source="01_mouse_behavior/01_baiermouse_allmice_weight_paired_observations.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1E", component="weight_time_ANOVA", source="01_mouse_behavior/01_baiermouse_allmice_weight_repeated_measures_anova.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1E", component="group_paired_tests", source="01_mouse_behavior/01_baiermouse_allmice_weight_within_group_holm_tests.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1F", component="individual_weight_loss", source="01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1F", component="welch_omnibus", source="01_mouse_behavior/01_baiermouse_allmice_weight_loss_welch_anova.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1F", component="pairwise_welch_Holm", source="01_mouse_behavior/01_baiermouse_allmice_weight_loss_welch_holm_tests.csv", selector="all", required=TRUE, note=""),
  list(panel="ED1G", component="weight_adjusted_enrichment", source="07_mouse_sensitivity/07_baiermouse_alllineages_weightadjusted_weight_sensitivity.csv", selector="all", required=TRUE, note=""),
  list(panel="ED2A", component="cellranger_counts", source="02_qc/02_baiermouse_allcells_cellranger_filtered_counts.csv", selector="all", required=TRUE, note=""),
  list(panel="ED2B", component="SoupX_contamination", source="02_qc/02_baiermouse_allcells_soupx_sample_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED2C", component="celllevel_QC", source="02_qc/02_baiermouse_allcells_presoupx_qc_metrics.rds", selector="rds_qc", required=TRUE, note=""),
  list(panel="ED2D", component="celllevel_QC", source="02_qc/02_baiermouse_allcells_presoupx_qc_metrics.rds", selector="rds_qc", required=TRUE, note=""),
  list(panel="ED2E", component="celllevel_QC", source="02_qc/02_baiermouse_allcells_presoupx_qc_metrics.rds", selector="rds_qc", required=TRUE, note=""),
  list(panel="ED2F", component="doublet_rate", source="03_final_qc/scDblFinder/03_baiermouse_allcells_doublets_sample_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED2G", component="permissive_screen", source="02_qc/02_baiermouse_allcells_presoupx_qc_retention.csv", selector="all", required=TRUE, note=""),
  list(panel="ED2G", component="sequential_QC", source="03_final_qc/03_baiermouse_allcells_final_qc_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED2H", component="umap_batch", source="", selector="atlas_umap", required=TRUE, note=""),
  list(panel="ED2I", component="allen_assignment_confidence", source="04_annotation/04_baiermouse_allcells_allen_mapping_confidence.csv", selector="all", required=TRUE, note=""),
  list(panel="ED2J", component="per_mouse_composition", source="05_differential_abundance/05_baiermouse_alllineages_allgroups_proportions_by_mouse.csv", selector="all", required=TRUE, note=""),
  list(panel="ED3A", component="mouse_hallmarks", source="06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv", selector="mouse_cognition_nine", required=TRUE, note=""),
  list(panel="ED3B", component="endothelial_all50", source="06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv", selector="endo_cog", required=TRUE, note=""),
  list(panel="ED3C", component="FGSEA_CAMERA", source="07_mouse_sensitivity/07_baiermouse_alllineages_cognition_fgsea_camera_comparison.csv", selector="all", required=TRUE, note=""),
  list(panel="ED3D", component="endothelial_mouse_program_scores", source="07_mouse_sensitivity/07_baiermouse_endothelial_cognition_endothelial_program_scores_by_mouse.csv", selector="all", required=TRUE, note=""),
  list(panel="ED3E", component="IEG_dissociation_sensitivity", source="07_mouse_sensitivity/07_baiermouse_endothelial_cognition_endothelial_immune_sensitivity.csv", selector="all", required=TRUE, note=""),
  list(panel="ED3F", component="endothelial_score_models", source="07_mouse_sensitivity/07_baiermouse_endothelial_cognition_endothelial_score_models.csv", selector="all", required=TRUE, note=""),
  list(panel="ED3G", component="pericyte_astrocyte_sensitivity", source="07_mouse_sensitivity/07_baiermouse_pericyte_cognition_pericyte_astro_transcript_sensitivity.csv", selector="all", required=TRUE, note=""),
  list(panel="ED4A", component="126_assignment_NES_empirical", source="06_mouse_bulk/06_baiermouse_alllineages_cognition_label_reassignment_pathway_p.csv", selector="all", required=TRUE, note=""),
  list(panel="ED4B", component="shared_reference_null", source="06_mouse_bulk/06_baiermouse_alllineages_cognition_label_reassignment_null_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED4C", component="astrocyte_OXPHOS_LOMO", source="07_mouse_sensitivity/07_baiermouse_astrocyte_cognition_astrocyte_oxphos_leave_one_mouse_out.csv", selector="all", required=TRUE, note=""),
  list(panel="ED4D", component="three_group_patterns", source="06_mouse_bulk/06_baiermouse_alllineages_threegroups_pathway_patterns.csv", selector="mouse_nine", required=TRUE, note=""),
  list(panel="ED5A", component="OXPHOS_rank_block_omission", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_cohort_block_omission_hallmarks.csv", selector="path_HALLMARK_OXIDATIVE_PHOSPHORYLATION", required=TRUE, note=""),
  list(panel="ED5B", component="rank_correlations", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_cohort_block_omission_rank_correlations.csv", selector="all", required=TRUE, note=""),
  list(panel="ED5C", component="pathway_mouse_human_ranks", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_mouse_versus_human_pathway_ranks.csv", selector="all", required=TRUE, note=""),
  list(panel="ED5C", component="correlations_FDR", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_mouse_human_rank_correlations.csv", selector="all", required=TRUE, note=""),
  list(panel="ED5D", component="individual_human_axis_correlations", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_mouse_cognition_vs_each_human_axis_nes_correlations.csv", selector="all", required=TRUE, note=""),
  list(panel="ED5E", component="cross_lineage_top10", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_integrated_cross_lineage_rank_summary.csv", selector="top10", required=TRUE, note=""),
  list(panel="ED6A", component="per_mouse_state_fractions", source="08_astrocyte_states/08_state_counts_fractions_by_mouse.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6A", component="state_group_means", source="08_astrocyte_states/08_state_abundance_group_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6B", component="state_abundance_model_seven", source="08_astrocyte_states/08_state_abundance_seven_test_models.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6C", component="state_coverage", source="08_astrocyte_states/08_state_coverage_by_group.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6D", component="state_leave_one_mouse_out", source="08_astrocyte_states/08_state_abundance_leave_one_mouse_out.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6E", component="AST1_QC_adjusted", source="08_astrocyte_states/08_AST1_abundance_QC_covariate_models.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6F", component="composition_reconstruction_programs", source="08_astrocyte_states/08_composition_reconstruction_program_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6G", component="within_state_program_models", source="08_astrocyte_states/08_within_state_program_group_models.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6H", component="paired_state_program_models", source="08_astrocyte_states/08_paired_state_group_models.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6I", component="matched_genes_within_state", source="08_astrocyte_states/08_matched_gene_within_state_programs.csv", selector="all", required=TRUE, note=""),
  list(panel="ED6J", component="S5_metabolic_sensitivity", source="08_astrocyte_states/08_S5_metabolic_module_sensitivity_all_models.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7A", component="Zemke_age_LODO", source="12_mouse_human_integration/12_zemke_astrocyte_aging_astrocyte_oxphos_leave_one_donor_out.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7A", component="full_donor_model", source="12_mouse_human_integration/12_zemke_astrocyte_aging_astrocyte_oxphos_age_model.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7B", component="SEAAD_metabolic_LODO", source="10_human_SEAAD/10_seaad_astrocyte_pathology_cognition_astrocyte_metabolic_loo_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7C", component="SEAAD_full_model_CI", source="10_human_SEAAD/10_seaad_astrocyte_pathology_cognition_astrocyte_metabolic_model_comparison.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7D", component="SEAAD_OXPHOS_method_comparison", source="10_human_SEAAD/10_seaad_alllineages_pathology_cognition_oxphos_method_comparison.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7E", component="ROSMAP_pooled_vs_regions", source="11_human_ROSMAP/11_rosmap_astrocyte_pathology_regions_pooled_vs_regions.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7F", component="ROSMAP_region_omission", source="11_human_ROSMAP/11_rosmap_astrocyte_pathology_regions_region_omission.csv", selector="all", required=TRUE, note=""),
  list(panel="ED7F", component="region_omission_summary", source="11_human_ROSMAP/11_rosmap_astrocyte_pathology_regions_region_omission_summary.csv", selector="all", required=TRUE, note=""),
  list(panel="ED8A", component="185_gene_recurrence_null", source="12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_full185_length_expression_matched_null.csv", selector="all", required=TRUE, note=""),
  list(panel="ED8B", component="selected_pathway_null", source="12_mouse_human_integration/12_crosscohort_allglia_five_axes_selected_pathway_conditional_nulls.csv", selector="all", required=TRUE, note=""),
  list(panel="ED8C", component="matching_strata", source="12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_185_verified_matching_strata.csv", selector="all", required=TRUE, note=""),
  list(panel="ED8D", component="four_lineage_recurrence", source="12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_gene_recurrence_by_glial_lineage.csv", selector="all", required=TRUE, note="")
), fill = TRUE)

#Extract annotated-atlas data for the panels that require it
atlas_needed <- any(panel_specs$selector %in% c("atlas_umap", "atlas_marker_dot"))
atlas_ok <- FALSE
if (atlas_needed)
  atlas_ok <- load_atlas_panel_data()

for (k in seq_len(nrow(panel_specs))) {
  sp <- panel_specs[k]
  message("Figure source: ", sp$panel, " / ", sp$component)
  x <- if (sp$selector %in% c("atlas_umap", "atlas_marker_dot")) {
    if (atlas_ok) get(sp$selector, envir = atlas_parts) else NULL
  } else read_panel(sp)
  if (is.null(x)) {
    status[[length(status) + 1L]] <- mkrow(
      sp$panel, sp$component, sp$source, sp$selector,
      "MISSING_UPSTREAM", message = "Not generated by numbered pipeline in this output root")
    next
  }
  x <- as.data.table(x)
  if (!nrow(x)) {
    status[[length(status) + 1L]] <- mkrow(
      sp$panel, sp$component, sp$source, sp$selector,
      "EMPTY_SOURCE", 0L, ncol(x), message = "Source exists but panel selection is empty")
    next
  }
  panel_file <- paste0("13_baierstudy_", tolower(sp$panel), "_", tolower(gsub("[^A-Za-z0-9]+", "_", sp$component)), "_source_data.csv")
  fwrite(x, file.path(panel_dir, panel_file))
  status[[length(status) + 1L]] <- mkrow(
    sp$panel, sp$component, sp$source, sp$selector,
    "EXPORTED", nrow(x), ncol(x), file = paste0("panels/", panel_file))
}

inventory <- rbindlist(status, fill = TRUE)
inventory <- merge(panel_specs, inventory,
                   by = c("panel", "component", "source", "selector"),
                   all.x = TRUE, sort = FALSE)
setorder(inventory, panel, component)
fwrite(inventory, file.path(out, "13_baierstudy_allpanels_figure_source_inventory.csv"))

#Display selections remain separate from the full statistical results
display_audit <- data.table(
  panel = c("F1F", "F3B", "F3F", "F3G", "F4A", "F4C", "F5A", "F5C", "F5D"),
  status = c("DETERMINISTIC", "DETERMINISTIC", "DETERMINISTIC", "DETERMINISTIC",
             "DISPLAY_SELECTION_NOT_FROZEN", "DISPLAY_SELECTION_NOT_FROZEN",
             "PREDEFINED_IN_SCRIPT", "DISPLAY_SELECTION_NOT_FROZEN",
             "DISPLAY_SELECTION_NOT_FROZEN"),
  note = c("Hallmarks with FDR<0.05 in any of nine lineages flagged",
           "Top 12 by integrated mean rank", "Top five mouse cognition ranks per lineage plus microglial OXPHOS",
           "Next two integrated non-OXPHOS pathways per lineage",
           "All GO tests exported; a 16-term figure display list must be defined explicitly",
           "All state-identity GO tests exported; representative labels require an explicit display list",
           "Eleven Hallmarks predefined in the human_aging_fig5a selector",
           "All astrocyte pathway statistics exported; selected pathways require a fixed display list",
           "All available ROSMAP program-region pairs exported; displayed row ordering requires a fixed list"))
fwrite(display_audit, file.path(out, "13_baierstudy_allpanels_figure_display_selection_audit.csv"))
fwrite(inventory[, .(n_components = .N,
                      n_exported = sum(status == "EXPORTED"),
                      n_missing = sum(status != "EXPORTED"),
                      completeness = if (all(status == "EXPORTED")) "COMPLETE" else "INCOMPLETE"),
                 by = panel][order(panel)],
       file.path(out, "13_baierstudy_allpanels_figure_completeness.csv"))
if (length(input_log))
  fwrite(rbindlist(input_log, fill = TRUE), file.path(out, "13_baierstudy_allpanels_figure_input_provenance.csv"))

#Manuscript-value audit
values <- list()
append_value <- function(label, observed, expected, tolerance, provenance) {
  values[[length(values) + 1L]] <<- data.table(
    comparison = label, observed = observed, manuscript = expected,
    delta = observed - expected, display_tolerance = tolerance,
    status = if (is.na(observed) || is.na(expected)) "NOT_EVALUATED" else
      if (abs(observed - expected) <= tolerance) "WITHIN_DISPLAY_PRECISION" else "DIFFERENT",
    upstream_file = provenance)
}

beh <- input_file("01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv")
if (!is.null(beh)) {
  append_value("behavior_n_mice", nrow(beh), 44, 0,
               "01_baiermouse_allmice_behavior_mouse_level.csv")
  append_value("behavior_n_sequenced", sum(as.logical(beh$sequenced), na.rm = TRUE),
               13, 0, "01_baiermouse_allmice_behavior_mouse_level.csv")
}
comp <- input_file("05_differential_abundance/05_baiermouse_alllineages_allgroups_proportions_by_mouse.csv")
if (!is.null(comp))
  append_value("final_singlets", sum(comp$n_cells, na.rm = TRUE), 198122, 0,
               "05_baiermouse_alllineages_allgroups_proportions_by_mouse.csv")
astro_u <- input_file("08_astrocyte_states/08_UMAP_and_respiratory_UCell.csv")
if (!is.null(astro_u))
  append_value("astrocytes_figure4", nrow(astro_u), 10155, 0,
               "08_UMAP_and_respiratory_UCell.csv")
state_da <- input_file("08_astrocyte_states/08_state_abundance_seven_test_models.csv")
if (!is.null(state_da)) {
  v <- state_da[grepl("c?AST1$", state, ignore.case = TRUE), FDR_seven_tests]
  append_value("AST1_abundance_FDR", if (length(v) == 1L) v else NA_real_,
               .015, .0005, "08_state_abundance_seven_test_models.csv")
}
ox <- input_file("12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_final_185_gene_manifest.csv")
if (!is.null(ox)) {
  append_value("eligible_OXPHOS", nrow(ox), 185, 0,
               "12_crosscohort_astrocyte_recurrence_oxphos_final_185_gene_manifest.csv")
  append_value("OXPHOS_favorable_all5", sum(as.logical(ox$all5_preserved), na.rm = TRUE),
               60, 0, "12_crosscohort_astrocyte_recurrence_oxphos_final_185_gene_manifest.csv")
  append_value("OXPHOS_recurrent_34", sum(as.logical(ox$recurrent34), na.rm = TRUE),
               34, 0, "12_crosscohort_astrocyte_recurrence_oxphos_final_185_gene_manifest.csv")
}
zd <- input_file("12_mouse_human_integration/12_zemke_astrocyte_aging_astrocyte_oxphos_donor_scores.csv")
if (!is.null(zd)) append_value("Zemke_OXPHOS_donors", nrow(zd), 40, 0,
                               "12_zemke_astrocyte_aging_astrocyte_oxphos_donor_scores.csv")
zm <- input_file("12_mouse_human_integration/12_zemke_astrocyte_aging_astrocyte_oxphos_age_model.csv")
if (!is.null(zm)) append_value("Zemke_OXPHOS_age_beta", zm$beta[1L], -.295, .0005,
                               "12_zemke_astrocyte_aging_astrocyte_oxphos_age_model.csv")
mo <- input_file("06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv")
if (!is.null(mo)) {
  v <- mo[lineage == "Astro-TE NN" & contrast == "Cog" &
            pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION", NES]
  append_value("mouse_astrocyte_OXPHOS_NES", if (length(v) == 1L) v else NA_real_,
               3.12, .005, "06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv")
}
ro <- input_file("11_human_ROSMAP/11_rosmap_astrocyte_pathology_regions_manuscript_nes_audit.csv")
if (!is.null(ro)) {
  for (j in seq_len(nrow(ro))) {
    append_value(paste0("ROSMAP_", gsub("[^A-Za-z0-9]", "_", ro$path[j]), "_NES"),
                 ro$observed_NES[j], ro$reference_NES[j], .005,
                 "11_rosmap_astrocyte_pathology_regions_manuscript_nes_audit.csv")
    append_value(paste0("ROSMAP_", gsub("[^A-Za-z0-9]", "_", ro$path[j]), "_FDR"),
                 ro$observed_FDR[j], ro$reference_FDR[j], .0005,
                 "11_rosmap_astrocyte_pathology_regions_manuscript_nes_audit.csv")
  }
}

#Supplementary manuscript checkpoints
if (!is.null(beh)) {
  aged <- beh[tolower(age_group) != "young"]
  if (nrow(aged))
    append_value("aged_not_attained", sum(!as.logical(aged$reversal_attained), na.rm = TRUE),
                 8, 0, "01_baiermouse_allmice_behavior_mouse_level.csv")
}
if (!is.null(mo)) {
  cog_astro <- mo[lineage == "Astro-TE NN" & contrast == "Cog"]
  append_value("astrocyte_intact_significant_Hallmarks",
               sum(cog_astro$padj < .05 & cog_astro$NES > 0, na.rm = TRUE),
               15, 0, "06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv")
  append_value("astrocyte_impaired_significant_Hallmarks",
               sum(cog_astro$padj < .05 & cog_astro$NES < 0, na.rm = TRUE),
               4, 0, "06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv")
}
if (!is.null(state_da)) {
  v <- state_da[grepl("cAST2_to_5$", state, ignore.case = TRUE), FDR_seven_tests]
  append_value("pooled_AST2_to_5_FDR", if (length(v) == 1L) v else NA_real_,
               .043, .0005, "08_state_abundance_seven_test_models.csv")
}
cf <- input_file("08_astrocyte_states/08_AST1_cognitive_flexibility.csv")
if (!is.null(cf)) append_value("AST1_fraction_CF_Spearman", cf$pooled_spearman_rho[1L],
                                .88, .005, "08_AST1_cognitive_flexibility.csv")
integ <- input_file("12_mouse_human_integration/12_crosscohort_allglia_five_axes_cohort_balanced_hallmarks_all200.csv")
if (!is.null(integ)) {
  oxr <- integ[pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]
  append_value("OXPHOS_number_of_glial_lineages_rank1",
               sum(oxr$preserved_rank == 1, na.rm = TRUE), 4, 0,
               "12_crosscohort_allglia_five_axes_cohort_balanced_hallmarks_all200.csv")
}

if (length(values))
  fwrite(rbindlist(values), file.path(out, "13_baierstudy_allpanels_figure_manuscript_numeric_audit.csv"))
