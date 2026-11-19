#Baier et al.
#14: Supplemental tables

#setwd("E:/baier_scrnaseq/scripts")
#source("14_supplemental_tables.R")

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))

root <- hca_output_root
out <- file.path(root, "14_supplemental_tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(openxlsx)

read_pipeline <- function(relative) {
  path <- file.path(root, relative)
  if (!file.exists(path)) stop("Required numbered-pipeline output is missing: ", relative)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  data.table::fread(path, showProgress = FALSE)
}
verify <- function(d, required, context) {
  if (!all(required %in% names(d))) {
    stop(context, " is missing columns: ",
         paste(setdiff(required, names(d)), collapse = ", "))
  }
  invisible(TRUE)
}
data_columns <- function(d, cols) {
  verify(d, cols, "Upstream table")
  d[, ..cols]
}

spec <- data.table(
  id = c(
  "S1a",
  "S1b",
  "S1c",
  "S1d",
  "S2a",
  "S2b",
  "S2c",
  "S2d",
  "S2e",
  "S3a",
  "S3b",
  "S3c",
  "S4a",
  "S4b",
  "S4c",
  "S4d",
  "S4e",
  "S4f",
  "S5a",
  "S5b",
  "S5c",
  "S5d",
  "S6a",
  "S6b",
  "S6c",
  "S7a",
  "S7b",
  "S7c",
  "S7d",
  "S7e",
  "S7f",
  "S7g",
  "S7h",
  "S7i",
  "S7j",
  "S7k",
  "S8a",
  "S8b",
  "S8c",
  "S8d",
  "S9a",
  "S9b",
  "S9c",
  "S9d",
  "S9e",
  "S9f",
  "S9g",
  "S9h",
  "S9i",
  "S9j",
  "S9k",
  "S10a",
  "S10b",
  "S10c",
  "S11a",
  "S11b",
  "S11c",
  "S11d",
  "S11e",
  "S11f"
),
  title = c(
  "Supplemental Table 1a. Sequenced mouse behavioral metadata.",
  "Supplemental Table 1b. Full mouse cohort behavioral data.",
  "Supplemental Table 1c. scRNA-seq QC.",
  "Supplemental Table 1d. Weight loss during behavioral assessment.",
  "Supplemental Table 2a. Allen class counts.",
  "Supplemental Table 2b. Allen subclass counts.",
  "Supplemental Table 2c. Cell lineage proportions by mouse.",
  "Supplemental Table 2d. Differential abundance by cell broad classification.",
  "Supplemental Table 2e. Differential abundance by cell subclass.",
  "Supplemental Table 3a. Cognition-associated differential expression across nine parenchymal cell populations.",
  "Supplemental Table 3b. Mouse pseudobulk gene-wise statistics (cognition and aging).",
  "Supplemental Table 3c. Cognitive flexibility continuous sensitivity.",
  "Supplemental Table 4a. Complete mouse Hallmark pathway screen across four contrasts.",
  "Supplemental Table 4b. GOBP cognition contrast.",
  "Supplemental Table 4c. Hallmark aging contrast.",
  "Supplemental Table 4d. CAMERA.",
  "Supplemental Table 4e. Five-axis glial Hallmark enrichment (mouse and human cohorts).",
  "Supplemental Table 4f. Hallmark, GOBP, Reactome pathways by astrocyte subclusters.",
  "Supplemental Table 5a. Three group Hallmark patterns.",
  "Supplemental Table 5b. Mouse-label reassignment correlation-null summary (126 assignments).",
  "Supplemental Table 5c. Mouse label permutation.",
  "Supplemental Table 5d. OXPHOS leave-one-mouse-out.",
  "Supplemental Table 6a. Oxidative/metabolic module assignments and gene membership.",
  "Supplemental Table 6b. Module Hallmark overlap.",
  "Supplemental Table 6c. Astrocyte metabolic modules across five species-native axes.",
  "Supplemental Table 7a. Astrocyte state abundance by mouse.",
  "Supplemental Table 7b. Astrocyte state abundance models.",
  "Supplemental Table 7c. Astrocyte mouse-level quality-control covariates.",
  "Supplemental Table 7d. Astrocyte state Hallmark, GOBP,  Reactome pathways.",
  "Supplemental Table 7e. Astrocyte state DEGs.",
  "Supplemental Table 7f. Astrocyte state differential abundance.",
  "Supplemental Table 7g. Individual mouse by astrocyte-state by transcriptional-program scores.",
  "Supplemental Table 7h. Group state means.",
  "Supplemental Table 7i. Paired differences.",
  "Supplemental Table 7j. Six batch-adjusted AST0/AST1 by program cognition effects.",
  "Supplemental Table 7k. Three paired AST1 minus AST0 cognition effects.",
  "Supplemental Table 8a. Whole-microglia Hallmark pathways across four mouse contrasts.",
  "Supplemental Table 8b. OPC and oligodendrocyte Hallmark pathways across four mouse contrasts.",
  "Supplemental Table 8c. Astrocyte-pericyte sensitivity.",
  "Supplemental Table 8d. Secondary glial lineage pathways across five cohort axes.",
  "Supplemental Table 9a. Cohort summary.",
  "Supplemental Table 9b. Cohort covariates summary.",
  "Supplemental Table 9c. Donor coverage by glial lineage.",
  "Supplemental Table 9d. Normative-aging and SEA-AD cognition/pathology human gene effects.",
  "Supplemental Table 9e. Zemke Hallmark pathways.",
  "Supplemental Table 9f. SEA-AD cognition Hallmark pathways.",
  "Supplemental Table 9g. SEA-AD pathology Hallmark pathways.",
  "Supplemental Table 9h. SEA-AD Hallmark pathways CAMERA.",
  "Supplemental Table 9i. SEA-AD astrocyte metabolic-program donor-omission sensitivity.",
  "Supplemental Table 9j. ROSMAP pathways by regions.",
  "Supplemental Table 9k. ROSMAP leave-out-regions.",
  "Supplemental Table 10a. Four-glial-lineage by five-context tested Hallmark results.",
  "Supplemental Table 10b. Cohort-balanced integrated pathway ranking by lineage.",
  "Supplemental Table 10c. OXPHOS leave-one-block-out.",
  "Supplemental Table 11a. OXPHOS gene manifest.",
  "Supplemental Table 11b. OXPHOS gene effects by axis.",
  "Supplemental Table 11c. Recurrent 34-gene leading-edge statistics.",
  "Supplemental Table 11d. Selected 13-lineage/pathway comparator family null.",
  "Supplemental Table 11e. Six-context metabolic programs (48 assessed tests).",
  "Supplemental Table 11f. Full-185-gene OXPHOS recurrence permutation null."
),
  relative = c(
  "06_mouse_bulk/06_baiermouse_allmice_pseudobulk_sample_metadata.csv;01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv",
  "01_mouse_behavior/01_baiermouse_allmice_behavior_mouse_level.csv",
  "03_final_qc/03_baiermouse_allcells_final_qc_summary.csv",
  "01_mouse_behavior/01_baiermouse_allmice_weight_group_means_sem.csv;01_mouse_behavior/01_baiermouse_allmice_weight_repeated_measures_anova.csv;01_mouse_behavior/01_baiermouse_allmice_weight_within_group_holm_tests.csv;01_mouse_behavior/01_baiermouse_allmice_weight_loss_welch_anova.csv;01_mouse_behavior/01_baiermouse_allmice_weight_loss_welch_holm_tests.csv;01_mouse_behavior/01_baiermouse_allmice_weight_loss_kruskal_sensitivity.csv;01_mouse_behavior/01_baiermouse_allmice_weight_loss_wilcoxon_sensitivity.csv",
  "04_annotation/04_baiermouse_allcells_allen_class_counts.csv",
  "04_annotation/04_baiermouse_allcells_allen_subclass_counts.csv",
  "05_differential_abundance/05_baiermouse_alllineages_allgroups_proportions_by_mouse.csv",
  "05_differential_abundance/05_baiermouse_alllineages_cognition_differential_abundance.csv",
  "05_differential_abundance/05_baiermouse_subclasses_cognition_differential_abundance.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_gene_effects.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_gene_effects.csv",
  "06_mouse_bulk/06_baiermouse_astrocyte_cognitive_flexibility_continuous_gene_effects.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_cognition_twogroup_gobp_enrichment.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv",
  "07_mouse_sensitivity/07_baiermouse_alllineages_cognition_camera_hallmarks_aged.csv",
  "12_mouse_human_integration/12_crosscohort_allglia_five_axes_hallmarks_all1000.csv",
  "08_astrocyte_states/08_state_identity_Hallmark_enrichment.csv;08_astrocyte_states/08_state_identity_GO_BP_enrichment.csv;08_astrocyte_states/08_state_identity_Reactome_enrichment.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_threegroups_pathway_patterns.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_cognition_label_reassignment_null_summary.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_cognition_label_reassignment_pathway_p.csv",
  "07_mouse_sensitivity/07_baiermouse_astrocyte_cognition_astrocyte_oxphos_leave_one_mouse_out.csv",
  "08_astrocyte_states/08_metabolic_module_membership_USED.csv;12_mouse_human_integration/12_crosscohort_astrocyte_gene_property_species_specific_program_membership.csv",
  "08_astrocyte_states/08_metabolic_module_membership_USED.csv;06_mouse_bulk/06_baiermouse_alllineages_hallmark_gene_sets.rds",
  "12_mouse_human_integration/12_crosscohort_astrocyte_gene_property_length_expression_binned_enrichment.csv",
  "08_astrocyte_states/08_state_counts_fractions_by_mouse.csv",
  "08_astrocyte_states/08_state_abundance_seven_test_models.csv",
  "08_astrocyte_states/08_astrocyte_QC_by_mouse.csv",
  "08_astrocyte_states/08_state_identity_Hallmark_enrichment.csv;08_astrocyte_states/08_state_identity_GO_BP_enrichment.csv;08_astrocyte_states/08_state_identity_Reactome_enrichment.csv",
  "08_astrocyte_states/08_state_identity_gene_effects.csv",
  "08_astrocyte_states/08_state_abundance_group_summary.csv;08_astrocyte_states/08_state_abundance_seven_test_models.csv",
  "08_astrocyte_states/08_Figure4g_three_programs_mouse_state.csv",
  "08_astrocyte_states/08_Figure4g_three_programs_mouse_state.csv",
  "08_astrocyte_states/08_mouse_paired_AST1_minus_AST0.csv",
  "08_astrocyte_states/08_within_state_program_group_models.csv;08_astrocyte_states/08_Figure4g_three_programs_mouse_state.csv",
  "08_astrocyte_states/08_paired_state_group_models.csv;08_astrocyte_states/08_mouse_paired_AST1_minus_AST0.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv",
  "06_mouse_bulk/06_baiermouse_alllineages_fourcontrasts_hallmark_enrichment.csv",
  "07_mouse_sensitivity/07_baiermouse_pericyte_cognition_pericyte_astro_transcript_all50.csv",
  "12_mouse_human_integration/12_crosscohort_allglia_five_axes_hallmarks_all1000.csv",
  "09_human_normative_aging/09_zemke_alllineages_aging_donor_design.csv;10_human_SEAAD/10_seaad_alllineages_pathology_cognition_donor_metadata.csv",
  "09_human_normative_aging/09_zemke_alllineages_aging_donor_design.csv;10_human_SEAAD/10_seaad_alllineages_pathology_cognition_donor_metadata.csv",
  "09_human_normative_aging/09_zemke_alllineages_aging_lineage_coverage.csv;10_human_SEAAD/SEAAD_inputs/10_seaad_alllineages_cognition_donor_coverage.csv",
  "09_human_normative_aging/09_zemke_astrocyte_aging_differential_expression.csv;09_human_normative_aging/09_zemke_microglia_aging_differential_expression.csv;09_human_normative_aging/09_zemke_opc_aging_differential_expression.csv;09_human_normative_aging/09_zemke_oligodendrocyte_aging_differential_expression.csv;10_human_SEAAD/10_seaad_alllineages_pathology_cognition_de_cognition_all_lineages.csv;10_human_SEAAD/10_seaad_alllineages_pathology_cognition_de_pathology_all_lineages.csv",
  "09_human_normative_aging/09_zemke_alllineages_aging_hallmark_younger_age_inference.csv",
  "10_human_SEAAD/10_seaad_alllineages_pathology_cognition_hallmark_cognition_inference.csv",
  "10_human_SEAAD/10_seaad_alllineages_pathology_cognition_hallmark_pathology_inference.csv",
  "10_human_SEAAD/10_seaad_alllineages_pathology_cognition_camera_all_hallmarks.csv",
  "10_human_SEAAD/10_seaad_astrocyte_pathology_cognition_astrocyte_metabolic_loo.csv",
  "11_human_ROSMAP/11_rosmap_astrocyte_pathology_regions_pooled_vs_regions.csv",
  "11_human_ROSMAP/11_rosmap_astrocyte_pathology_regions_region_omission_summary.csv",
  "12_mouse_human_integration/12_crosscohort_allglia_five_axes_hallmarks_all1000.csv",
  "12_mouse_human_integration/12_crosscohort_allglia_five_axes_cohort_balanced_hallmarks_all200.csv",
  "12_mouse_human_integration/12_crosscohort_allglia_five_axes_cohort_block_omission_hallmarks.csv",
  "12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_final_185_gene_manifest.csv;12_mouse_human_integration/12_crosscohort_allglia_five_axes_gene_consensus_all_glia.csv",
  "12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_185_five_axis_effects_and_leading_edges.csv;12_mouse_human_integration/12_crosscohort_allglia_five_axes_gene_consensus_all_glia.csv",
  "12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_185_five_axis_effects_and_leading_edges.csv;12_mouse_human_integration/12_crosscohort_allglia_five_axes_gene_consensus_all_glia.csv",
  "12_mouse_human_integration/12_crosscohort_allglia_five_axes_selected_pathway_conditional_nulls.csv",
  "12_mouse_human_integration/12_crosscohort_allglia_five_axes_metabolic_six_context_48test_source.csv",
  "12_mouse_human_integration/12_crosscohort_astrocyte_recurrence_oxphos_full185_length_expression_matched_null.csv"
),
  mode = c(
  "sequenced_metadata",
  "direct",
  "direct",
  "weight_loss_tests",
  "direct",
  "direct",
  "direct",
  "direct",
  "direct",
  "within_lineage_de_summary",
  "mouse_de",
  "direct",
  "wide_hallmarks",
  "direct",
  "aging",
  "direct",
  "direct",
  "state_pathways",
  "direct",
  "direct",
  "direct",
  "direct",
  "module_membership",
  "module_overlap",
  "gene_property_axes",
  "direct",
  "direct",
  "direct",
  "state_pathways",
  "direct",
  "state_summary",
  "direct",
  "state_means",
  "paired_delta",
  "state_models",
  "paired_models",
  "microglia_bulk",
  "oligolineage_bulk",
  "direct",
  "secondary_axes",
  "human_cohort",
  "human_covariates",
  "human_coverage",
  "human_de",
  "direct",
  "direct",
  "direct",
  "direct",
  "direct",
  "rosmap_regions",
  "direct",
  "direct",
  "integrated_rank",
  "oxphos_blocks",
  "oxphos_manifest",
  "oxphos_effects",
  "oxphos_recurrent",
  "direct",
  "direct",
  "direct"
),
  column_order = c(
  "mouse_id|sequence_id|sequence_batch|age_group|primary_group|gender|strain|dob|phenotyper_run_date|age|pre_weight|post_weight|pct_weight_loss_calc|cage_card|initial_entries_to_seventy|initial_entries_to_eighty|initial_entries_to_ninety|reversal_entries_to_seventy|reversal_entries_to_eighty|reversal_entries_to_ninety|initial_hours_to_seventy|initial_hours_to_eighty|initial_hours_to_ninety|reversal_hours_to_seventy|reversal_hours_to_eighty|reversal_hours_to_ninety|cf51to61|sequenced|r80_reached",
  "mouse_id|sequence_id|sequence_batch|age_group|subgroup_re80|gender|strain|DOB|phenotyper_run_date|Age|PreWeight|PostWeight|PercentWeightLost|CageCard|InitialEntriesToSeventy|InitialEntriesToEighty|InitialEntriesToNinety|ReversalEntriesToSeventy|ReversalEntriesToEighty|ReversalEntriesToNinety|InitialHoursToSeventy|InitialHoursToEighty|InitialHoursToNinety|ReversalHoursToSeventy|ReversalHoursToEighty|ReversalHoursToNinety|CF51to61",
  "mouse_id|sequence_id|sequence_batch|age_group|primary_group|n_before_miqc|n_after_miqc|n_removed_miqc|pct_retained_miqc|median_genes_before|median_mt_before|median_genes_after|median_mt_after|n_before_doublet|n_singlets|n_doublets|pct_doublets|pct_singlets|pct_final_vs_preSoup",
  "analysis|comparison|group1|group2|n_group1|n_group2|mean_pre_g|mean_post_g|mean_loss_pct|sem_loss_pct|effect|effect_unit|CI_low|CI_high|statistic|statistic_value|df_num|df_den|p_value|p_holm|source_file",
  "class_name|n",
  "subclass_name|n",
  "mouse_id|broad_lineage|n_cells|total_cells|proportion|sequence_batch|age_group|primary_group",
  "cell_type|PropMean.intact|PropMean.impaired|PropRatio|Tstatistic|P.Value|FDR|mean_reference_prop|mean_test_prop|median_reference_prop|median_test_prop|difference_pp|log2_ratio|reference_group|test_group",
  "cell_type|PropMean.intact|PropMean.impaired|PropRatio|Tstatistic|P.Value|FDR|mean_reference_prop|mean_test_prop|median_reference_prop|median_test_prop|difference_pp|log2_ratio|reference_group|test_group",
  "cell_type|contrast|n_tested|n_FDR05|n_FDR10|min_FDR",
  "gene|logFC|logCPM|F|PValue|FDR|original_source_file|source_group",
  "gene|logFC|logCPM|F|PValue|FDR",
  "cell_type|pathway|primary_NES|primary_P|primary_FDR|primary_screenwide_FDR|extreme_NES|extreme_P|extreme_FDR|CF_NES|CF_P|CF_FDR|aging_NES|aging_P|aging_FDR|CF_NES_aligned|extreme_concordant|CF_concordant|aging_relation",
  "cell_type|pathway|pval|padj|log2err|ES|NES|size|leadingEdge",
  "cell_type|pathway|pval|padj|log2err|ES|NES|size|leadingEdge",
  "pathway|NGenes|Direction|PValue|FDR|analysis_celltype|direction_sign|signed_log10P|NGenes_mapped|FDR_within|rank_P|FDR_screen",
  "cell_type|pathway|NES_native|FDR|axis|axis_label|preserve_multiplier|NES_preserved",
  "pathway_id|pval|padj|log2err|ES|NES|size|state|collection",
  "lineage|pathway|NES_AgeP|NES_AgePool|NES_Cog|NES_YvImp|padj_AgeP|padj_AgePool|padj_Cog|padj_YvImp|pattern",
  "lineage|pathway|obs_NES|p_two_sided|p_one_sided|p_one_sided_excl_1swap|n_exceed_two_sided|n_exceed_that_are_1swap|n_exceed_that_are_mirror|rho_NES_vs_overlap|mean_NES_ov0|mean_NES_ov1|mean_NES_ov2|mean_NES_ov3|p_alignment",
  "pathway|obs_NES|perm_p|lineage",
  "dropped|contrast|OX_NES|OX_padj|OX_rank|OX_median_logFC|OX_frac_pos|OX_frac_P05|other_frac_P05|core_frac_pos|core_median_logFC|group",
  "module|gene|priority|n_modules|exclusive_module|species",
  "module|hallmark|n_module|n_shared|frac_of_hallmark|jaccard",
  "pathway|pval|padj|log2err|ES|NES|size|leadingEdge|lineage|axis_id|species|set_type",
  "mouse_id|state|n|group|sequence_batch|N_astro|fraction",
  "label|estimate|p|batch_adjusted|n|FDR_seven_tests|intact_mean_fraction|impaired_mean_fraction|difference_fraction_imp_minus_intact",
  "mouse_id|group|state|metric|median_value|mean_value|n_cells",
  "pathway_id|pval|padj|log2err|ES|NES|size|state|collection",
  "gene|logFC|logCPM|F|PValue|target|signed_stat",
  "compartment|group|mean|lo|hi|n|FDR_seven_tests|segment_label|plot_x",
  "mouse_id|group|state|program|score|n_cells|source|pair_jitter|state_offset|x",
  "program|group|state|mouse_mean|mouse_SE|x",
  "mouse_id|group|program|AST0|AST1|delta",
  "program|state|family|estimate|se|ci_low|ci_high|p_value|df_residual|n_mice|FDR_BH6",
  "program|state|family|estimate|se|ci_low|ci_high|p_value|df_residual|n_mice|FDR_BH3",
  "final_state|pathway|aging_NES|aging_padj|aging_screenFDR|cognition_NES|cognition_padj|cognition_screenFDR|extreme_NES|extreme_padj|extreme_screenFDR|aging_cognition_same_direction|aging_extreme_same_direction|cognition_extreme_same_direction|aging_sig|cognition_sig|extreme_sig|significant_both_aging_cognition|significant_both_opposite|cognition_robust",
  "lineage|pathway|primary_NES_preserved|primary_FDR|primary_screenwide_FDR|extreme_NES_preserved|extreme_FDR|CF_NES_preserved|CF_FDR|aging_NES_preserved|aging_FDR|primary_sig|extreme_sig|CF_sig|aging_sig|extreme_same_direction|CF_same_direction|aging_same_direction|robust_cognition|cognition_direction|aging_relation",
  "pathway|NES|pval|padj|size|leadingEdge|model",
  "lineage|pathway|axis_id|NES|padj|FDR05",
  "dataset|n_donors|sex_distribution|release",
  "dataset|metric|source_variable|n_nonmissing|median|Q1|Q3",
  "dataset|lineage|n_donors_with_profile",
  "gene|logFC|logCPM|F|PValue|FDR|original_source_file|source_group|cell_type|model|coefficient|interpretation",
  "crossstudy_celltype|pathway|pval|padj|log2err|ES|NES|size|leadingEdge|screenwide_FDR",
  "cell_type|coefficient|direction_note|pathway|pval|padj|log2err|ES|NES|size|leadingEdge|screenwide_FDR",
  "cell_type|coefficient|direction_note|pathway|pval|padj|log2err|ES|NES|size|leadingEdge|screenwide_FDR",
  "pathway|NGenes|Correlation|Direction|PValue|FDR|model|coefficient|correlation_assumption|n_donors|n_expressed_genes|favorable_direction|FDR_recomputed",
  "model|pathway|dropped_donor|n_donors|beta_cognition|se_cognition|p_cognition|HC3_SE|HC3_p|adj_R2|max_CooksD|max_leverage|design_condition|significant_HC3_005",
  "path|pathway|ranking|combined_NES|combined_FDR|n_regions|n_favorable_NES|n_unfavorable_NES|median_NES|min_NES|max_NES|hippocampus_NES|n_region_FDR_005",
  "path|pathway|ranking|min_omission_mean_NES|max_omission_mean_NES|all_omission_means_positive",
  "pathway|pval|padj|log2err|ES|NES|size|leadingEdge|lineage|axis_id|axis_label|statistic_type|universe_type|favorable_rank|favorable_rank_percentile|direction|leadingEdge_string|FDR_screen1000|favorable|favorable_FDR05",
  "lineage|pathway|preserved_rank|integrated_NES|integrated_FDR",
  "cell_type|pathway|analysis|consensus_score|favorable_rank|lineage|analysis_label",
  "human_gene|human_aging_effect_pct|mouse_aging_effect_pct|mouse_cognition_effect_pct|seaad_cognition_effect_pct|seaad_pathology_effect_pct|human_aging_inference_pct|mouse_aging_inference_pct|mouse_cognition_inference_pct|seaad_cognition_inference_pct|seaad_pathology_inference_pct|human_aging_oriented_logFC|mouse_aging_oriented_logFC|mouse_cognition_oriented_logFC|seaad_cognition_oriented_logFC|seaad_pathology_oriented_logFC|human_aging_F|mouse_aging_F|mouse_cognition_F|seaad_cognition_F|seaad_pathology_F|human_aging_PValue|mouse_aging_PValue|mouse_cognition_PValue|seaad_cognition_PValue|seaad_pathology_PValue|human_aging_FDR|mouse_aging_FDR|mouse_cognition_FDR|seaad_cognition_FDR|seaad_pathology_FDR|mouse_gene|mouse_block_effect|human_aging_block_effect|seaad_block_effect|consensus_effect|mouse_block_inference|human_aging_block_inference|seaad_block_inference|consensus_inference|direction_n_5axes|all5_preserved|direction_n_3blocks|all3blocks_preserved|consensus_effect_rank|consensus_inference_rank|lineage|n_axes_leading_edge|recurrent34|class",
  "human_gene|human_aging_inference_pct|mouse_aging_inference_pct|mouse_cognition_inference_pct|seaad_cognition_inference_pct|seaad_pathology_inference_pct|human_aging_oriented_logFC|mouse_aging_oriented_logFC|mouse_cognition_oriented_logFC|seaad_cognition_oriented_logFC|seaad_pathology_oriented_logFC|human_aging_F|mouse_aging_F|mouse_cognition_F|seaad_cognition_F|seaad_pathology_F|human_aging_PValue|mouse_aging_PValue|mouse_cognition_PValue|seaad_cognition_PValue|seaad_pathology_PValue|human_aging_FDR|mouse_aging_FDR|mouse_cognition_FDR|seaad_cognition_FDR|seaad_pathology_FDR|mouse_gene|mouse_block_effect|human_aging_block_effect|seaad_block_effect|consensus_effect|mouse_block_inference|human_aging_block_inference|seaad_block_inference|consensus_inference|direction_n_5axes|all5_preserved|direction_n_3blocks|all3blocks_preserved|consensus_effect_rank|consensus_inference_rank|lineage|n_axes_leading_edge|recurrent34|class|axis_column|signed_effect_pct|axis_id",
  "human_gene|axis_id|human_aging_inference_pct|mouse_aging_inference_pct|mouse_cognition_inference_pct|seaad_cognition_inference_pct|seaad_pathology_inference_pct|human_aging_oriented_logFC|mouse_aging_oriented_logFC|mouse_cognition_oriented_logFC|seaad_cognition_oriented_logFC|seaad_pathology_oriented_logFC|human_aging_F|mouse_aging_F|mouse_cognition_F|seaad_cognition_F|seaad_pathology_F|human_aging_PValue|mouse_aging_PValue|mouse_cognition_PValue|seaad_cognition_PValue|seaad_pathology_PValue|human_aging_FDR|mouse_aging_FDR|mouse_cognition_FDR|seaad_cognition_FDR|seaad_pathology_FDR|mouse_gene|mouse_block_effect|human_aging_block_effect|seaad_block_effect|consensus_effect|mouse_block_inference|human_aging_block_inference|seaad_block_inference|consensus_inference|direction_n_5axes|all5_preserved|direction_n_3blocks|all3blocks_preserved|consensus_effect_rank|consensus_inference_rank|lineage|n_axes_leading_edge|recurrent34|class|axis_column|signed_effect_pct|direction_preserved|leading_edge|gene|axis|LE",
  "lineage|pathway|matched|mode|n_genes|n_genes_omitted|all5_pathway_NES_favorable|k_bins|n_strata|metric|observed|null_mean|null_q975|p_upper|reason|BH_selected",
  "set|axis_id|matched_effect|empirical_p|module|axis|family|assessed|significant|q_BH_48|BH48_q05|nominal_p05|cx|cy",
  ""
)
)

book_title <- c(
  "Supplementary Table S1 — Mouse samples, behavior, QC and weight changes",
  "Supplementary Table S2 — Cell annotation and mouse-level abundance",
  "Supplementary Table S3 — Mouse pseudobulk DEGs",
  "Supplementary Table S4 — Mouse Hallmark, GOBP, Camera results.",
  "Supplementary Table S5 — Cognition vs. aging and 126-assignment empirical null.",
  "Supplementary Table S6 — Astrocyte metabolic module gene definitions and membership",
  "Supplementary Table S7 — Six-state astrocyte abundance and within-state models.",
  "Supplementary Table S8 — Whole microglia, oligodendrocyte and pericyte results",
  "Supplementary Table S9 — Public human dataset cohort, lineage coverage, DEG and pathway results.",
  "Supplementary Table S10 — Five-axis cohort-integrated pathway analysis.",
  "Supplementary Table S11 — OXPHOS cross-context recurrence and matched comparisons"
)

make_table <- function(item) {
  files <- strsplit(item$relative, ";", fixed = TRUE)[[1L]]
  f <- function(i = 1L) read_pipeline(files[[i]])
  mode <- item$mode

  switch(mode,
    direct = as.data.table(f()),
    aging = {
      d <- as.data.table(f())
      verify(d, "contrast", "Mouse aging Hallmark output")
      d[contrast == "AgeP"]
    },
    sequenced_metadata = {
      meta <- f()
      beh <- f(2L)
      verify(meta, c("mouse_id", "sequence_batch", "age_group", "primary_group"), "Mouse metadata")
      verify(beh, c("mouse_id", "sequenced"), "Mouse behavior")
      if (nrow(meta) != 13L || anyDuplicated(meta$mouse_id))
        stop("Expected 13 distinct sequenced mice")
      # The group labels in this output originate in the 13-mouse pseudobulk ledger.
      keep <- setdiff(names(beh), c("mouse_id", "age_group"))
      z <- merge(meta, beh[, c("mouse_id", keep), with = FALSE], by = "mouse_id",
                 all.x = TRUE, sort = FALSE)
      z <- z[match(meta$mouse_id, mouse_id)]
      if (!all(as.logical(z$sequenced))) stop("Unsequenced mouse in sequenced metadata")
      if (!all(c(sum(tolower(z$primary_group) == "young"),
                 sum(tolower(z$primary_group) == "intact"),
                 sum(tolower(z$primary_group) == "impaired")) == c(4L, 4L, 5L)))
        stop("Sequenced-group counts must be 4 Young / 4 Intact / 5 Impaired")
      z
    },
    weight_loss_tests = {
      # Group summaries and all reported primary and sensitivity weight tests
      # are read from Script 01, with no models re-fitted here.
      group_means <- f(1L)
      rm_anova <- f(2L)
      within <- f(3L)
      welch_anova <- f(4L)
      welch_pairs <- f(5L)
      kruskal <- f(6L)
      wilcoxon <- f(7L)
      verify(group_means,
             c("group", "n", "mean_pre", "mean_post", "mean_loss_pct", "sem_loss_pct"),
             "Mouse weight group summaries")
      verify(rm_anova, c("effect", "numDF", "denDF", "F-value", "p-value"),
             "Mouse repeated-measures ANOVA")
      verify(within, c("group", "contrast", "estimate", "df", "CI_low", "CI_high",
                       "t", "p", "Holm_p_three_groups"),
             "Mouse within-group weight contrasts")
      verify(welch_anova, c("test", "F", "df_num", "df_den", "p"),
             "Mouse weight-loss Welch ANOVA")
      verify(welch_pairs,
             c("group1", "group2", "n1", "n2", "difference_group2_minus_group1",
               "CI_low", "CI_high", "p", "Holm_p_three_pairs"),
             "Mouse weight-loss Welch pairwise contrasts")
      verify(kruskal, c("test", "statistic", "df", "p"),
             "Mouse weight-loss Kruskal-Wallis sensitivity")
      verify(wilcoxon, c("group1", "group2", "W", "p", "Holm_p_three_pairs"),
             "Mouse weight-loss Wilcoxon sensitivity")

      rows <- list(
        group_means[, .(
          analysis = "Group descriptive statistics", comparison = group,
          group1 = group, n_group1 = n,
          mean_pre_g = mean_pre, mean_post_g = mean_post,
          mean_loss_pct, sem_loss_pct,
          source_file = basename(files[[1L]]))],
        rm_anova[, .(
          analysis = "Repeated-measures ANOVA", comparison = effect,
          statistic = "F", statistic_value = `F-value`,
          df_num = numDF, df_den = denDF, p_value = `p-value`,
          source_file = basename(files[[2L]]))],
        within[, .(
          analysis = "Pre-post within-group contrast", comparison = contrast,
          group1 = group, effect = estimate, effect_unit = "grams",
          CI_low, CI_high, statistic = "t", statistic_value = t,
          df_den = df, p_value = p, p_holm = Holm_p_three_groups,
          source_file = basename(files[[3L]]))],
        welch_anova[, .(
          analysis = "Weight-loss Welch ANOVA", comparison = test,
          statistic = "F", statistic_value = F,
          df_num, df_den, p_value = p,
          source_file = basename(files[[4L]]))],
        welch_pairs[, .(
          analysis = "Weight-loss Welch pairwise comparison",
          comparison = paste(group1, "vs", group2), group1, group2,
          n_group1 = n1, n_group2 = n2,
          effect = difference_group2_minus_group1,
          effect_unit = "percentage points", CI_low, CI_high,
          p_value = p, p_holm = Holm_p_three_pairs,
          source_file = basename(files[[5L]]))],
        kruskal[, .(
          analysis = "Weight-loss Kruskal-Wallis sensitivity", comparison = test,
          statistic = "chi-squared", statistic_value = statistic,
          df_num = df, p_value = p,
          source_file = basename(files[[6L]]))],
        wilcoxon[, .(
          analysis = "Weight-loss Wilcoxon sensitivity",
          comparison = paste(group1, "vs", group2), group1, group2,
          statistic = "W", statistic_value = W,
          p_value = p, p_holm = Holm_p_three_pairs,
          source_file = basename(files[[7L]]))]
      )
      result <- rbindlist(rows, use.names = TRUE, fill = TRUE)
      columns <- c("analysis", "comparison", "group1", "group2", "n_group1", "n_group2",
                   "mean_pre_g", "mean_post_g", "mean_loss_pct", "sem_loss_pct",
                   "effect", "effect_unit", "CI_low", "CI_high", "statistic",
                   "statistic_value", "df_num", "df_den", "p_value", "p_holm", "source_file")
      if (nrow(result) != 18L || !all(columns %in% names(result)))
        stop("Unexpected mouse weight-loss supplemental table structure")
      setcolorder(result, columns)
      result
    },
    within_lineage_de_summary = {
      d <- f()
      verify(d, c("lineage", "contrast", "gene", "FDR"),
             "Mouse cognition gene effects")
      parenchymal <- c("Astro-TE NN", "Microglia NN", "OPC NN", "Oligo NN",
                       "DG Glut", "Endo NN", "Peri NN", "BAM NN", "Lymphoid NN")
      d <- d[contrast == "Cog" & lineage %in% parenchymal]
      if (uniqueN(d$lineage) != length(parenchymal) ||
          anyDuplicated(d, by = c("lineage", "gene")))
        stop("Cognition DE source must have unique gene tests in nine parenchymal populations")
      if (anyNA(d$FDR) || any(!is.finite(d$FDR)))
        stop("Cognition DE source has invalid adjusted P values")
      result <- d[, .(n_tested = .N,
                      n_FDR05 = sum(FDR < 0.05),
                      n_FDR10 = sum(FDR < 0.10),
                      min_FDR = min(FDR)),
                  by = .(cell_type = lineage, contrast)]
      result[, lineage_order := match(cell_type, parenchymal)]
      setorder(result, lineage_order)
      result[, lineage_order := NULL]
      result
    },
    mouse_de = {
      d <- f()[contrast %in% c("Cog", "AgeP")]
      verify(d, c("gene", "lineage", "contrast", "logFC", "logCPM", "F", "PValue", "FDR"),
             "Mouse gene-wise effects")
      # Do not discard lineage and contrast: the same gene occurs in several models.
      d[, `:=`(original_source_file = basename(files[[1L]]),
               source_group = "mouse_pseudobulk")]
      d[, .(gene, logFC, logCPM, F, PValue, FDR,
            original_source_file, source_group, lineage, contrast)]
    },
    wide_hallmarks = {
      d <- f()
      verify(d, c("lineage", "contrast", "pathway", "NES", "pval", "padj"),
             "Four-contrast Hallmark screen")
      stopifnot(!anyDuplicated(d, by = c("lineage", "pathway", "contrast")))
      dcast(d, lineage + pathway ~ contrast, value.var = c("NES", "pval", "padj"))
    },
    state_pathways = {
      labels <- c("HALLMARK", "GO_BP", "REACTOME")
      rbindlist(lapply(seq_along(files), function(i) {
        d <- f(i)
        verify(d, c("state", "pathway", "NES", "pval", "padj"),
               paste("State identity collection", labels[i]))
        d[, .(state, pathway, NES, pval, padj,
              collection = labels[i])]
      }), use.names = TRUE)
    },
    module_membership = {
      mouse <- unique(f()[, .(module, gene)])
      human <- f(2L)
      verify(human, c("species", "pathway", "gene"), "Human module source")
      human <- unique(human[tolower(species) == "human" &
                               pathway %in% mouse$module,
                            .(module = pathway, gene)])
      rbindlist(list(
        mouse[, .(module, gene, species = "mouse")],
        human[, .(module, gene, species = "human")]
      ), use.names = TRUE)
    },
    module_overlap = {
      members <- unique(f()[, .(module, gene)])
      halls <- f(2L)
      wanted <- intersect(
        c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_ADIPOGENESIS"),
        names(halls))
      if (length(wanted) != 3L) stop("Frozen mouse Hallmark membership is incomplete")
      rbindlist(lapply(unique(members$module), function(mod) {
        a <- unique(toupper(as.character(members[module == mod, gene])))
        rbindlist(lapply(wanted, function(h) {
          b <- unique(toupper(as.character(halls[[h]])))
          shared <- length(intersect(a, b))
          data.table(module = mod, hallmark = h, n_module = length(a),
                     n_shared = shared, frac_of_hallmark = shared/length(b),
                     jaccard = shared/length(union(a, b)))
        }))
      }))
    },
    gene_property_axes = {
      d <- f()
      verify(d, c("axis_id", "ranking", "species", "pathway"),
             "Five-axis length/expression sensitivity")
      d[ranking == "raw"]
    },
    state_summary = {
      means <- f()
      models <- f(2L)
      verify(means, c("group","state","mean_fraction","sem_fraction","n_mice"),
             "Astrocyte state means")
      verify(models, c("state","FDR_seven_tests"), "Astrocyte seven-test models")
      if (anyDuplicated(models$state)) stop("Duplicate state in abundance models")
      merge(means, models[, .(state, FDR_seven_tests)],
            by = "state", all.x = TRUE, sort = FALSE)
    },
    state_means = {
      d <- f()
      verify(d, c("group","state","program","score","mouse_id"), "Figure 4G scores")
      if (anyDuplicated(d, by = c("mouse_id", "state", "program")))
        stop("Duplicate mouse-state-program scores")
      d[, .(mouse_mean = mean(score),
            mouse_SE = if (.N > 1L) sd(score)/sqrt(.N) else NA_real_,
            n_mice = .N), by = .(program, group, state)]
    },
    paired_delta = {
      d <- f()
      verify(d, c("mouse_id","group","program","cAST0","cAST1","AST1_minus_AST0"),
             "Paired astrocyte state scores")
      d[, .(mouse_id, group, program, AST0 = cAST0,
            AST1 = cAST1, delta = AST1_minus_AST0)]
    },
    state_models = {
      d <- f()
      scores <- f(2L)
      verify(d, c("state","program","estimate_intact_minus_impaired","SE",
                  "PValue","n_mice","FDR_six"), "Within-state model")
      z <- scores[group %in% c("Intact", "Impaired") & n_cells >= 10L,
                  .(n_observed = .N, n_batches = uniqueN(sequence_batch),
                    n_groups = uniqueN(group)), by = .(state, program)]
      d <- merge(d, z, by = c("state", "program"), all.x = TRUE, sort = FALSE)
      if (anyNA(d$n_observed) || any(d$n_observed != d$n_mice) ||
          any(d$n_groups != 2L)) stop("Within-state model degrees of freedom cannot be verified")
      #Original design was ~batch + group (full rank, two group levels).
      d[, df_residual := n_observed - n_batches - 1L]
      if (any(d$df_residual <= 0L)) stop("Nonpositive within-state residual df")
      d[, .(program, state, family = "six_within_state",
             estimate = estimate_intact_minus_impaired,
             se = SE,
             ci_low = estimate_intact_minus_impaired - qt(.975, df_residual) * SE,
             ci_high = estimate_intact_minus_impaired + qt(.975, df_residual) * SE,
             p_value = PValue, df_residual, n_mice, FDR_BH6 = FDR_six)]
    },
    paired_models = {
      d <- f()
      scores <- f(2L)
      verify(d, c("program","estimate_intact_minus_impaired","SE","PValue",
                  "n_mice","FDR_three"), "Paired-state model")
      z <- scores[group %in% c("Intact","Impaired") & is.finite(AST1_minus_AST0),
                  .(n_observed = .N, n_batches = uniqueN(sequence_batch),
                    n_groups = uniqueN(group)), by = program]
      d <- merge(d, z, by = "program", all.x = TRUE, sort = FALSE)
      if (anyNA(d$n_observed) || any(d$n_observed != d$n_mice) ||
          any(d$n_groups != 2L)) stop("Paired-state model degrees of freedom cannot be verified")
      d[, df_residual := n_observed - n_batches - 1L]
      if (any(d$df_residual <= 0L)) stop("Nonpositive paired-state residual df")
      d[, .(program, state = "AST1_minus_AST0", family = "three_paired_state",
             estimate = estimate_intact_minus_impaired,
             se = SE,
             ci_low = estimate_intact_minus_impaired - qt(.975, df_residual)*SE,
             ci_high = estimate_intact_minus_impaired + qt(.975, df_residual)*SE,
             p_value = PValue, df_residual, n_mice, FDR_BH3 = FDR_three)]
    },
    microglia_bulk = {
      d <- f()[lineage == "Microglia NN"]
      verify(d, c("pathway","contrast","NES","pval","padj"),
             "Whole-microglia Hallmark screen")
      if (uniqueN(d$pathway) != 50L || uniqueN(d$contrast) != 4L)
        stop("Whole-microglia Hallmark analysis must contain 50 pathways x 4 contrasts")
      if (anyDuplicated(d, by = c("lineage","pathway","contrast")))
        stop("Duplicate microglia pathway/contrast")
      dcast(d, lineage + pathway ~ contrast, value.var = c("NES","pval","padj"))
    },
    oligolineage_bulk = {
      d <- f()[lineage %in% c("OPC NN", "Oligo NN")]
      verify(d, c("lineage","pathway","contrast","NES","pval","padj"),
             "Whole oligodendrocyte-lineage Hallmarks")
      if (uniqueN(d$lineage) != 2L ||
          anyDuplicated(d, by = c("lineage","pathway","contrast")))
        stop("OPC/Oligo pathway specification unexpectedly changed")
      dcast(d, lineage + pathway ~ contrast, value.var = c("NES","pval","padj"))
    },
    secondary_axes = {
      d <- f()
      d[lineage %in% c("Microglia","OPC","Oligodendrocyte")]
    },
    human_cohort = {
      z <- f()
      s <- f(2L)
      verify(z, c("donor_id", "Gender"), "Zemke donor metadata")
      verify(s, c("donor_id", "Sex"), "SEA-AD donor metadata")
      sex_summary <- function(d, col) {
        v <- as.character(d[[col]])
        v <- v[!is.na(v) & nzchar(v)]
        paste(paste(names(table(v)), as.integer(table(v)), sep = ": "),
              collapse = "; ")
      }
      data.table(dataset = c("Zemke","SEA-AD HIP"),
                 n_donors = c(uniqueN(z$donor_id), uniqueN(s$donor_id)),
                 sex_distribution = c(sex_summary(z,"Gender"),
                                      sex_summary(s,"Sex")),
                 release = "Aggregate cohort data only")
    },
    human_covariates = {
      z <- f()
      s <- f(2L)
      summary_var <- function(dataset, d, variable, label) {
        verify(d, variable, dataset)
        x <- suppressWarnings(as.numeric(d[[variable]]))
        x <- x[is.finite(x)]
        data.table(dataset = dataset, metric = label, source_variable = variable,
                   n_nonmissing = length(x),
                   median = if(length(x)) median(x) else NA_real_,
                   Q1 = if(length(x)) unname(quantile(x,.25)) else NA_real_,
                   Q3 = if(length(x)) unname(quantile(x,.75)) else NA_real_)
      }
      rbindlist(list(
        summary_var("Zemke", z, "Age", "Age (years)"),
        summary_var("SEA-AD HIP", s, "age_z", "Age (standardized)"),
        summary_var("SEA-AD HIP", s, "cps_z", "Local pathology (standardized)"),
        summary_var("SEA-AD HIP", s, "cog_z", "Cognitive trajectory (standardized)")
      ))
    },
    human_coverage = {
      z <- f()
      s <- f(2L)
      verify(z, c("donor_id","crossstudy_celltype","n_nuclei"), "Zemke coverage")
      verify(s, c("donor_label","crossstudy_celltype","n_cells"), "SEA-AD coverage")
      za <- z[n_nuclei > 0, .(n_donors_with_profile = uniqueN(donor_id)),
              by = .(lineage = crossstudy_celltype)][, dataset := "Zemke"]
      sa <- s[n_cells > 0, .(n_donors_with_profile = uniqueN(donor_label)),
              by = .(lineage = crossstudy_celltype)][, dataset := "SEA-AD HIP"]
      rbindlist(list(za, sa))[, .(dataset, lineage, n_donors_with_profile)]
    },
    human_de = {
      parts <- lapply(seq_along(files), function(i) {
        d <- f(i)
        verify(d, c("gene","logFC","logCPM","F","PValue","FDR"), basename(files[[i]]))
        origin <- basename(files[[i]])
        if (i <= 4L) {
          ln <- c("Astrocyte","Microglia","OPC","Oligodendrocyte")[[i]]
          d[, `:=`(cell_type = ln, model = "normative_age",
                   coefficient = "age_z")]
        } else {
          verify(d, "cell_type", origin)
          d[, `:=`(model = if (i == 5L) "cognition" else "pathology",
                   coefficient = if (i == 5L) "cog_z" else "cps_z")]
        }
        d[, `:=`(original_source_file = origin,
                 source_group = "human_donor_pseudobulk")]
        d[, .(gene, logFC, logCPM, F, PValue, FDR,
              original_source_file, source_group, cell_type, model, coefficient)]
      })
      rbindlist(parts, use.names = TRUE, fill = TRUE)
    },
    rosmap_regions = {
      d <- f()
      verify(d, c("path","pathway","ranking","NES","padj","n_regions",
                  "n_positive","n_negative","median_NES","min_NES",
                  "max_NES","hippocampus_NES","n_FDR05"), "ROSMAP regional summaries")
      d[, .(path, pathway, ranking,
            combined_NES = NES, combined_FDR = padj,
            n_regions, n_favorable_NES = n_positive,
            n_unfavorable_NES = n_negative, median_NES, min_NES, max_NES,
            hippocampus_NES, n_region_FDR_005 = n_FDR05)]
    },
    integrated_rank = {
      d <- f()
      verify(d, c("lineage","pathway","preserved_rank","NES","padj"),
             "Integrated Hallmark rankings")
      d[, .(lineage, pathway, preserved_rank,
            integrated_NES = NES, integrated_FDR = padj)]
    },
    oxphos_blocks = {
      d <- f()
      verify(d, c("pathway","lineage"), "Cohort-block omission")
      d[pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]
    },
    oxphos_manifest = {
      man <- f()
      wide <- f(2L)[lineage == "Astrocyte" & human_gene %in% man$human_gene]
      if (nrow(man) != 185L || anyDuplicated(man$human_gene) ||
          anyDuplicated(wide$human_gene)) stop("185-gene manifest/consensus identity discrepancy")
      merge(wide, man, by = "human_gene", all.x = TRUE,
            sort = FALSE, suffixes = c("", "_manifest"))
    },
    oxphos_effects = {
      ef <- f()
      wide <- f(2L)[lineage == "Astrocyte" & human_gene %in% ef$human_gene]
      verify(ef, c("human_gene","axis_id","recurrent34"), "OXPHOS effects")
      if (nrow(ef) != 925L || anyDuplicated(ef, by = c("human_gene","axis_id")))
        stop("Expected 185 x 5 unique gene/axis effects")
      merge(ef, wide, by = "human_gene", all.x = TRUE, sort = FALSE)
    },
    oxphos_recurrent = {
      ef <- f()[as.logical(recurrent34) %in% TRUE]
      wide <- f(2L)[lineage == "Astrocyte" & human_gene %in% ef$human_gene]
      if (nrow(ef) != 170L || uniqueN(ef$human_gene) != 34L)
        stop("Expected 34 recurrent genes x five axes")
      merge(ef, wide, by = "human_gene", all.x = TRUE, sort = FALSE)
    },
    stop("Unrecognized supplemental export mode: ", mode)
  )
}

title_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 12, textDecoration = "bold",
  fontColour = "#000000")
author_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 12, fontColour = "#000000")
column_header_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 12, fontColour = "#000000",
  valign = "center", wrapText = TRUE)


excel_column <- function(n) {
  stopifnot(length(n) == 1L, is.finite(n), n >= 1L, n == as.integer(n))
  out <- ""
  while (n > 0L) {
    i <- (n - 1L) %% 26L + 1L
    out <- paste0(LETTERS[[i]], out)
    n <- (n - 1L) %/% 26L
  }
  out
}

replace_sheet_dimension <- function(path, expected_range) {
  original <- file(path, open = "rb")
  temporary <- tempfile("worksheet_", tmpdir = dirname(path), fileext = ".xml")
  revised <- file(temporary, open = "wb")
  first <- readBin(original, what = "raw", n = 4096L)
  prefix <- rawToChar(first)
  pattern <- '<dimension ref="[^"]+"\\s*/>'
  if (!grepl(pattern, prefix, perl = TRUE)) {
    close(original)
    close(revised)
    unlink(temporary)
    stop("Worksheet dimension tag missing from: ", basename(path))
  }
  prefix <- sub(pattern, paste0('<dimension ref="', expected_range, '"/>'),
                prefix, perl = TRUE)
  writeBin(charToRaw(prefix), revised)
  repeat {
    block <- readBin(original, what = "raw", n = 4L * 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, revised)
  }
  close(original)
  close(revised)
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    unlink(temporary)
    stop("Could not rewrite worksheet range metadata: ", path)
  }
  unlink(temporary)
}

repair_xlsx_package <- function(path, sheet_ranges) {
  if (!requireNamespace("zip", quietly = TRUE))
    stop("The 'zip' R package is needed for Excel compatibility; install.packages('zip')")
  directory <- tempfile("supp_xlsx_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, exdir = directory)

  for (i in seq_along(sheet_ranges)) {
    worksheet <- file.path(directory, "xl", "worksheets",
                           paste0("sheet", i, ".xml"))
    if (!file.exists(worksheet)) stop("Missing worksheet XML: ", basename(worksheet))
    replace_sheet_dimension(worksheet, unname(sheet_ranges[[i]]))
  }

  drawings_dir <- file.path(directory, "xl", "drawings")
  if (dir.exists(drawings_dir))
    stop("Unexpected drawings in a supplemental data-only workbook: ", path)
  relationship_dir <- file.path(directory, "xl", "worksheets", "_rels")
  relationship_files <- list.files(relationship_dir, pattern = "\\.rels$",
                                   full.names = TRUE)
  for (relationship_file in relationship_files) {
    xml <- paste(readLines(relationship_file, warn = FALSE), collapse = "")
    xml <- gsub('<Relationship[^>]*Type="[^"]+/(drawing|vmlDrawing)"[^>]*/>',
                "", xml, perl = TRUE)
    writeLines(xml, relationship_file, useBytes = TRUE)
  }
  content_types <- file.path(directory, "[Content_Types].xml")
  xml <- paste(readLines(content_types, warn = FALSE), collapse = "")
  xml <- gsub('<Override[^>]*PartName="/xl/drawings/[^"]+"[^>]*/>',
              "", xml, perl = TRUE)
  writeLines(xml, content_types, useBytes = TRUE)

  members <- list.files(directory, all.files = TRUE, recursive = TRUE,
                        include.dirs = FALSE, no.. = TRUE)
  repaired <- tempfile("supp_verified_", tmpdir = dirname(path),
                       fileext = ".xlsx")
  on.exit(unlink(repaired), add = TRUE)
  zip::zipr(zipfile = repaired, files = members, root = directory,
            mode = "mirror", include_directories = FALSE,
            compression_level = 6)
  archive_members <- utils::unzip(repaired, list = TRUE)$Name
  essential <- c("[Content_Types].xml", "xl/workbook.xml",
                 paste0("xl/worksheets/sheet", seq_along(sheet_ranges), ".xml"))
  if (!all(essential %in% archive_members))
    stop("Repaired XLSX is missing required package members: ", path)
  if (!file.copy(repaired, path, overwrite = TRUE))
    stop("Could not replace the exported workbook with its OOXML repair: ", path)
  invisible(path)
}

records <- list()

write_heading <- function(wb, sheet, heading) {
  openxlsx::writeData(wb, sheet, heading, startRow = 1, colNames = FALSE)
  openxlsx::writeData(wb, sheet, "BAIER et al.", startRow = 2, colNames = FALSE)
  openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1)
  openxlsx::addStyle(wb, sheet, author_style, rows = 2, cols = 1)
}

for (number in seq_len(11L)) {
  sheet_ids <- spec$id[as.integer(sub("^S([0-9]+)[a-z]$", "\\1", spec$id)) == number]
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 12,
                           fontColour = "#000000")
  openxlsx::addWorksheet(wb, "INDEX", gridLines = FALSE)
  # Range of INDEX: headings at rows 1/2, one sheet description per row 4+.
  worksheet_ranges <- c(INDEX = paste0("A1:A", 3L + length(sheet_ids)))
  write_heading(wb, "INDEX", book_title[[number]])
  openxlsx::setColWidths(wb, "INDEX", cols = 1L, widths = 90)
  openxlsx::setColWidths(wb, "INDEX", cols = 2L, widths = 18)

  for (id in sheet_ids) {
    row_ix <- match(id, spec$id)
    item <- spec[row_ix]
    if (nrow(item) != 1L) stop("Duplicate or absent sheet specification: ", id)
    message("Supplemental ", id, ": ", item$title)

    index_row <- 4L + match(id, sheet_ids) - 1L
    openxlsx::writeData(wb, "INDEX", item$title,
                        startRow = index_row, colNames = FALSE)
    openxlsx::addWorksheet(wb, id, gridLines = FALSE)
    write_heading(wb, id, item$title)

    source_files <- if (nzchar(item$relative))
      strsplit(item$relative, ";", fixed = TRUE)[[1L]] else character()
    absent <- source_files[!file.exists(file.path(root, source_files))]
    failure <- NULL

    if (length(absent)) {
      result <- NULL
      failure <- paste0("Upstream output(s) not present: ",
                        paste(absent, collapse = "; "))
      state <- "MISSING_UPSTREAM"
    } else {
      result <- tryCatch(make_table(item),
                         error = function(e) {
                           failure <<- conditionMessage(e)
                           NULL
                         })
      state <- if (is.null(result)) "EXPORT_ERROR" else "EXPORTED"
    }

    if (!is.null(result)) {
      result <- as.data.table(result)
      if (!nrow(result) || !ncol(result)) {
        state <- "EMPTY_SOURCE"
        failure <- "Upstream source had no rows or columns"
        result <- NULL
      } else if (nrow(result) > 1048572L) {
        state <- "EXCEL_ROW_LIMIT"
        failure <- paste("Excel data-row limit exceeded:", nrow(result))
        result <- NULL
      }
    }

    if (is.null(result)) {
      data_out <- data.frame(
        status = state, explanation = if (is.null(failure)) "" else failure,
        stringsAsFactors = FALSE)
      exported_rows <- NA_integer_
      exported_cols <- NA_integer_
    } else {
      ordered_columns <- if (nzchar(item$column_order))
        strsplit(item$column_order, "|", fixed = TRUE)[[1L]] else character()
      if (length(ordered_columns) && all(ordered_columns %in% names(result))) {
        data.table::setcolorder(result, c(ordered_columns,
                                          setdiff(names(result), ordered_columns)))
      }
      exported_rows <- nrow(result)
      exported_cols <- ncol(result)
      data_out <- as.data.frame(result, stringsAsFactors = FALSE)
    }

    openxlsx::writeData(wb, id, data_out, startRow = 4,
                        colNames = TRUE, rowNames = FALSE,
                        headerStyle = column_header_style,
                        withFilter = !is.null(result))
    openxlsx::setRowHeights(wb, id, rows = 4, heights = 32)
    openxlsx::setColWidths(wb, id, cols = seq_len(ncol(data_out)),
                          widths = 19)
    openxlsx::setColWidths(wb, id, cols = 1L, widths = 24)
    openxlsx::freezePane(wb, id, firstActiveRow = 5L)
    worksheet_ranges <- c(worksheet_ranges,
      setNames(paste0("A1:", excel_column(ncol(data_out)), 4L + nrow(data_out)), id))
    if (is.null(result)) message("  ", state, ": ", failure)

    records[[length(records) + 1L]] <- data.table(
      sheet = id, workbook = number, title = item$title, status = state,
      source_files = paste(source_files, collapse = "; "),
      missing_sources = paste(absent, collapse = "; "),
      data_rows = exported_rows, data_columns = exported_cols,
      error_or_missing_reason = if (is.null(failure)) "" else failure
    )
    rm(data_out, result)
    invisible(gc(verbose = FALSE))
  }

  target <- file.path(out, sprintf("BAIER_supplemental_table%d.xlsx", number))
  openxlsx::saveWorkbook(wb, target, overwrite = TRUE)
  repair_xlsx_package(target, worksheet_ranges)
  rm(wb)
  invisible(gc(verbose = FALSE))
  message("Saved ", target)
}

audit <- rbindlist(records, use.names = TRUE, fill = TRUE)
data.table::fwrite(audit, file.path(out, "14_supplemental_source_and_schema_audit.csv"))
data.table::fwrite(audit[status != "EXPORTED"],
                   file.path(out, "14_supplemental_missing_or_incomplete.csv"))

#Separate cumulative index
index_wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(index_wb, fontName = "Arial", fontSize = 12,
                         fontColour = "#000000")
openxlsx::addWorksheet(index_wb, "INDEX", gridLines = FALSE)
write_heading(index_wb, "INDEX", "Supplemental Table Index.")
openxlsx::setColWidths(index_wb, "INDEX", cols = 1L, widths = 90)
for (j in seq_len(nrow(spec))) {
  openxlsx::writeData(index_wb, "INDEX", spec$title[[j]],
                      startRow = 3L + j, colNames = FALSE)
}
index_path <- file.path(out, "BAIER_supplemental_table_index.xlsx")
openxlsx::saveWorkbook(index_wb, index_path, overwrite = TRUE)
repair_xlsx_package(index_path, c(INDEX = paste0("A1:A", 3L + nrow(spec))))
writeLines(capture.output(sessionInfo()),
           file.path(out, "14_supplemental_sessioninfo.txt"))
