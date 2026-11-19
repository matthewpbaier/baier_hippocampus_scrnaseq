#Baier et al.
#10: SEA-AD hippocampal donor-level pathology, cognition, and sensitivity analyses

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "10_human_SEAAD")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

#Derive pseudobulk from Allen HIP raw counts, original taxonomy, and SEA-AD clinical workbooks
input_dir <- file.path(out, "SEAAD_inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
primary_types <- c("Astrocyte", "Microglia", "OPC", "Oligodendrocyte", "DG neuron")
safe_types <- c(Astrocyte = "Astrocyte", Microglia = "Microglia", OPC = "OPC",
                Oligodendrocyte = "Oligodendrocyte", "DG neuron" = "DG_neuron")
input_names <- c("10_seaad_alldonors_cognition_eligible_cohort.csv",
                 paste0("10_seaad_", tolower(unname(safe_types)), "_cognition_raw_pseudobulk_counts.csv.gz"))
missing_inputs <- file.path(input_dir, input_names)
missing_inputs <- missing_inputs[!file.exists(missing_inputs)]
if (length(missing_inputs)) {
  python <- Sys.getenv("HCA_PYTHON", unset = "python")
  prep <- file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()),
                    "10_prepare_SEAAD_inputs.py")
  status <- system2(python,
                    args = c(shQuote(prep), "--input", shQuote(hca_seaad_input_root),
                             "--output", shQuote(input_dir)))
  if (status != 0L) stop("SEA-AD raw-data preparation failed (exit ", status, ")")
}
meta_file <- file.path(input_dir, "10_seaad_alldonors_cognition_eligible_cohort.csv")
pb_dir <- input_dir

#Load libraries
library(data.table)
library(edgeR)
library(limma)
library(fgsea)
library(msigdbr)
library(BiocParallel)
library(readr)

set.seed(20261002)

write_table <- function(x, name) data.table::fwrite(as.data.table(x), file.path(out, name))

#Retain donor covariates 
meta <- fread(meta_file)
meta[, donor_id := as.character(`Donor ID`)]
meta[, Sex := factor(as.character(Sex), levels = c("Female", "Male"))]
meta[, `:=`(
  age_z = as.numeric(scale(as.numeric(`Age at Death`))),
  cps_z = as.numeric(scale(as.numeric(CPS_Local))),
  cog_z = as.numeric(scale(as.numeric(slope_zcasi_irt0))),
  memory_slope_z = as.numeric(scale(as.numeric(slope_zmem0))),
  RIN_z = as.numeric(scale(as.numeric(RIN))),
  PMI_z = as.numeric(scale(as.numeric(PMI)))
)]
needed_model <- c("donor_id", "Sex", "age_z", "cps_z", "cog_z",
                  "memory_slope_z", "RIN_z", "PMI_z")
write_table(meta[, ..needed_model], "10_seaad_alllineages_pathology_cognition_donor_metadata.csv")

primary_formula <- ~ Sex + age_z + cps_z + cog_z
rin_formula <- ~ Sex + age_z + cps_z + cog_z + RIN_z
memory_formula <- ~ Sex + age_z + cps_z + memory_slope_z
check_design <- function(formula, info = meta) {
  design <- model.matrix(formula, data = info)
  rownames(design) <- info$donor_id
  design
}
design_primary <- check_design(primary_formula)
design_rin <- check_design(rin_formula)
design_memory <- check_design(memory_formula)

#Load original donor x gene pseudobulk counts for five frozen lineages
read_counts <- function(cell_type, collapse_symbols = FALSE) {
  filename <- file.path(pb_dir, paste0("10_seaad_",
                                     tolower(safe_types[[cell_type]]), "_cognition_raw_pseudobulk_counts.csv.gz"))
  if (!file.exists(filename)) stop("Missing SEA-AD pseudobulk: ", filename)
  x <- as.data.table(readr::read_csv(filename, show_col_types = FALSE))
  ids <- c("gene_identifier", "gene_symbol")
  if (!all(ids %in% names(x)) ||
      !setequal(setdiff(names(x), ids), meta$donor_id))
    stop("Unexpected SEA-AD gene/donor columns in ", filename)
  donor_cols <- meta$donor_id
  cnt <- as.matrix(x[, ..donor_cols])
  storage.mode(cnt) <- "double"
  if (any(!is.finite(cnt)) || any(cnt < 0) ||
      max(abs(cnt - round(cnt))) > 1e-6)
    stop("SEA-AD pseudobulk counts must be finite nonnegative integers: ", cell_type)
  genes <- as.character(x$gene_symbol)
  if (collapse_symbols) {
    genes <- toupper(trimws(genes))
    keep <- !is.na(genes) & nzchar(genes) & genes != "NA"
    cnt <- rowsum(cnt[keep, , drop = FALSE], group = genes[keep], reorder = TRUE)
  } else {
    # Keep original SEA-AD modeling gene-row handling for its edgeR estimates.
    rownames(cnt) <- genes
  }
  colnames(cnt) <- meta$donor_id
  cnt
}

#Donor-level edgeR QL models (age, sex, pathology, cognition)
fit_edger <- function(counts, design, coefficients) {
  if (!identical(colnames(counts), rownames(design)))
    stop("Count columns and design donor order differ")
  y <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  tests <- setNames(lapply(coefficients, function(coefficient) {
    qlf <- edgeR::glmQLFTest(fit, coef = coefficient)
    tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
    as.data.table(tab, keep.rownames = "gene")
  }), coefficients)
  list(dge = y, tests = tests, n_input = nrow(counts),
       n_tested = nrow(y))
}

counts_by_type <- setNames(lapply(primary_types, read_counts), primary_types)
primary <- setNames(lapply(primary_types, function(ct) {
  message("Primary donor-level edgeR: ", ct)
  fit_edger(counts_by_type[[ct]], design_primary, c("cog_z", "cps_z"))
}), primary_types)
write_table(rbindlist(lapply(names(primary), function(ct) {
  copy(primary[[ct]]$tests$cog_z)[, cell_type := ct]
}), fill = TRUE), "10_seaad_alllineages_pathology_cognition_de_cognition_all_lineages.csv")
write_table(rbindlist(lapply(names(primary), function(ct) {
  copy(primary[[ct]]$tests$cps_z)[, cell_type := ct]
}), fill = TRUE), "10_seaad_alllineages_pathology_cognition_de_pathology_all_lineages.csv")

# RIN-adjusted and alternative memory-slope models sensitivities
sens <- setNames(lapply(primary_types, function(ct) {
  list(RIN = fit_edger(counts_by_type[[ct]], design_rin, "cog_z")$tests$cog_z,
       memory = fit_edger(counts_by_type[[ct]], design_memory, "memory_slope_z")$tests$memory_slope_z)
}), primary_types)

#Human Hallmarks and signed-statistic/logFC-ranked FGSEA
ms <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"))
hs <- lapply(split(ms$gene_symbol, ms$gs_name), unique)
if (length(hs) != 50L) stop("Expected 50 human Hallmark gene sets")
saveRDS(hs, file.path(out, "10_seaad_alllineages_pathology_cognition_hallmark_sets.rds"))
write_table(unique(ms[, .(pathway = gs_name, gene = gene_symbol)]),
            "10_seaad_alllineages_pathology_cognition_hallmark_membership.csv")

rank_genes <- function(tab, type = c("signed_F", "logFC"), direction = 1) {
  type <- match.arg(type)
  z <- copy(as.data.table(tab))
  if (type == "signed_F") {
    if (any(z$F[is.finite(z$F)] < -1e-6)) stop("Material negative edgeR F")
    z[, statistic := direction * sign(logFC) * sqrt(pmax(F, 0))]
  } else {
    z[, statistic := direction * logFC]
  }
  z <- z[!is.na(gene) & nzchar(gene) & is.finite(statistic)]
  z <- z[order(-abs(statistic))]
  z <- z[!duplicated(gene)]
  sort(setNames(z$statistic, z$gene), decreasing = TRUE)
}

enrich <- function(tab, type = "signed_F", direction = 1) {
  ranks <- rank_genes(tab, type = type, direction = direction)
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = hs, stats = ranks, minSize = 10, maxSize = 500,
    eps = 0, BPPARAM = BiocParallel::SerialParam()))
  as.data.table(fg)[, .(pathway, NES, pval, padj)]
}

collect_gsea <- function(type, source, axis, direction = 1) {
  rbindlist(lapply(primary_types, function(ct) {
    z <- if (source == "primary") primary[[ct]]$tests[[axis]] else sens[[ct]][[source]]
    fg <- enrich(z, type, direction = direction)
    fg[, `:=`(cell_type = ct, source = source, axis = axis,
              favorable_NES = if (axis == "cps_z") -NES else NES)]
    fg
  }))
}

hall_cog <- collect_gsea("signed_F", "primary", "cog_z")
hall_path <- collect_gsea("signed_F", "primary", "cps_z")
hall_cog_fc <- collect_gsea("logFC", "primary", "cog_z")
hall_path_fc <- collect_gsea("logFC", "primary", "cps_z")
for (tab in list(hall_cog, hall_path)) {
  if (!all(tab$pathway %in% names(hs))) stop("Unknown Hallmark pathway")
}

hall_cog[, screenwide_FDR := p.adjust(pval, "BH")]
hall_path[, screenwide_FDR := p.adjust(pval, "BH")]
write_table(hall_cog, "10_seaad_alllineages_pathology_cognition_hallmark_cognition_inference.csv")
write_table(hall_path, "10_seaad_alllineages_pathology_cognition_hallmark_pathology_inference.csv")
write_table(hall_cog_fc, "10_seaad_alllineages_pathology_cognition_hallmark_cognition_log_fc_rank.csv")
write_table(hall_path_fc, "10_seaad_alllineages_pathology_cognition_hallmark_pathology_log_fc_rank.csv")

pair <- merge(hall_cog_fc[, .(cell_type, pathway, cognition_NES = NES,
                              cognition_FDR = padj)],
              hall_path_fc[, .(cell_type, pathway, pathology_NES = NES,
                               pathology_FDR = padj)],
              by = c("cell_type", "pathway"))
pair[, `:=`(pathology_favorable_NES = -pathology_NES,
            cognition_favorable_NES = cognition_NES)]
write_table(pair, "10_seaad_alllineages_pathology_cognition_hallmark_cognition_vs_pathology.csv")

hall_rin <- collect_gsea("signed_F", "RIN", "cog_z")
hall_memory <- collect_gsea("signed_F", "memory", "memory_slope_z")
write_table(hall_rin, "10_seaad_alllineages_pathology_cognition_hallmark_rin_sensitivity.csv")
write_table(hall_memory, "10_seaad_alllineages_pathology_cognition_hallmark_memory_slope_sensitivity.csv")

#Donor-level Hallmark expression scores and adjusted regressions
make_scores <- function(dge, cell_type) {
  expr <- edgeR::cpm(dge, log = TRUE, prior.count = 2,
                     normalized.lib.sizes = TRUE)
  rbindlist(lapply(names(hs), function(pathway) {
    use <- intersect(hs[[pathway]], rownames(expr))
    if (length(use) < 10L) return(NULL)
    x <- expr[use, , drop = FALSE]
    gene_sd <- apply(x, 1, stats::sd)
    x <- x[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]
    if (nrow(x) < 10L) return(NULL)
    z <- t(scale(t(x)))
    data.table(donor_id = colnames(expr), cell_type = cell_type,
               pathway = pathway, score = colMeans(z), n_genes = nrow(x))
  }), fill = TRUE)
}
score_long <- rbindlist(lapply(primary_types, function(ct)
  make_scores(primary[[ct]]$dge, ct)), fill = TRUE)
score_long <- merge(score_long, meta[, ..needed_model], by = "donor_id",
                    all.x = TRUE, sort = FALSE)
write_table(score_long, "10_seaad_alllineages_pathology_cognition_donor_hallmark_scores.csv")

score_models <- list(
  primary = score ~ Sex + age_z + cps_z + cog_z,
  plus_RIN = score ~ Sex + age_z + cps_z + cog_z + RIN_z,
  plus_PMI = score ~ Sex + age_z + cps_z + cog_z + PMI_z,
  plus_RIN_PMI = score ~ Sex + age_z + cps_z + cog_z + RIN_z + PMI_z
)
fit_score <- function(d, model_name, drop_id = NA_character_) {
  if (!is.na(drop_id)) d <- d[donor_id != drop_id]
  fit <- stats::lm(score_models[[model_name]], data = d)
  if (qr(model.matrix(fit))$rank != ncol(model.matrix(fit)))
    stop("Rank-deficient donor-score model: ", model_name)
  tab <- summary(fit)$coefficients
  data.table(model = model_name, dropped_donor = drop_id,
             n_donors = nobs(fit),
             cognition_beta = unname(tab["cog_z", "Estimate"]),
             cognition_SE = unname(tab["cog_z", "Std. Error"]),
             cognition_p = unname(tab["cog_z", "Pr(>|t|)"]),
             pathology_beta = unname(tab["cps_z", "Estimate"]),
             pathology_SE = unname(tab["cps_z", "Std. Error"]),
             pathology_p = unname(tab["cps_z", "Pr(>|t|)"]))
}

score_fits <- rbindlist(lapply(
  split(score_long, by = c("cell_type", "pathway")),
  function(d) {
    rbindlist(lapply(c("primary", "plus_RIN"), function(k) {
      z <- fit_score(d, k)
      z[, `:=`(cell_type = unique(d$cell_type), pathway = unique(d$pathway))]
      z
    }))
  }
), fill = TRUE)
score_fits[, cognition_screenwide_FDR := p.adjust(cognition_p, "BH"),
           by = model]
score_fits[, pathology_screenwide_FDR := p.adjust(pathology_p, "BH"),
           by = model]
write_table(score_fits, "10_seaad_alllineages_pathology_cognition_donor_score_regressions.csv")

#Donor omissions for each model
targets <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
             "HALLMARK_FATTY_ACID_METABOLISM")
focus_scores <- score_long[cell_type == "Astrocyte" & pathway %in% targets]
if (nrow(focus_scores) != 62L || uniqueN(focus_scores$donor_id) != 31L)
  stop("Missing 31-donor astrocyte OXPHOS/fatty-acid program scores")
focus_full <- rbindlist(lapply(split(focus_scores, by = "pathway"), function(d) {
  rbindlist(lapply(names(score_models), function(k) {
    z <- fit_score(d, k)
    z[, pathway := unique(d$pathway)]
    z
  }))
}))
focus_full[, `:=`(cognition_ci_low = cognition_beta - 1.96 * cognition_SE,
                  cognition_ci_high = cognition_beta + 1.96 * cognition_SE)]
write_table(focus_full, "10_seaad_astrocyte_pathology_cognition_astrocyte_metabolic_model_comparison.csv")

focus_loo <- rbindlist(lapply(split(focus_scores, by = "pathway"), function(d) {
  rbindlist(lapply(names(score_models), function(k) {
    rbindlist(lapply(d$donor_id, function(id) {
      z <- fit_score(d, k, drop_id = id)
      z[, pathway := unique(d$pathway)]
      z
    }))
  }))
}))
if (nrow(focus_loo) != 2L * 4L * 31L)
  stop("Incomplete SEA-AD astrocyte LOO sensitivity")
write_table(focus_loo, "10_seaad_astrocyte_pathology_cognition_astrocyte_metabolic_loo.csv")
write_table(focus_loo[, .(n_omissions = .N,
                          min_beta = min(cognition_beta),
                          max_beta = max(cognition_beta),
                          n_positive = sum(cognition_beta > 0),
                          min_p = min(cognition_p),
                          max_p = max(cognition_p)),
                      by = .(pathway, model)],
            "10_seaad_astrocyte_pathology_cognition_astrocyte_metabolic_loo_summary.csv")

#Astrocyte competitive CAMERA and self-contained FRY/MROAST
astro_counts <- read_counts("Astrocyte", collapse_symbols = TRUE)
sets_upper <- lapply(hs, function(g) unique(toupper(trimws(g))))
run_voom <- function(design) {
  y <- edgeR::DGEList(counts = astro_counts)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  vo <- limma::voom(y, design = design, plot = FALSE)
  idx <- limma::ids2indices(sets_upper, rownames(vo$E), remove.empty = TRUE)
  idx <- idx[lengths(idx) >= 10L & lengths(idx) <= 500L]
  if (length(idx) != 50L)
    stop("Expected all 50 Hallmarks in SEA-AD astrocyte voom universe; found ", length(idx))
  list(vo = vo, idx = idx)
}

camera_list <- list()
fry_list <- list()
mroast_table <- NULL
for (model in c("primary", "RIN_adjusted")) {
  design <- if (model == "primary") design_primary else design_rin
  vv <- run_voom(design)
  for (axis in c("cog_z", "cps_z")) {
    coef_index <- which(colnames(design) == axis)
    for (cor_mode in c("estimated_per_set", "fixed_0.01")) {
      cor_value <- if (cor_mode == "estimated_per_set") NA_real_ else 0.01
      z <- as.data.table(limma::camera(
        vv$vo, index = vv$idx, design = design,
        contrast = coef_index, inter.gene.cor = cor_value,
        allow.neg.cor = FALSE, use.ranks = FALSE, sort = FALSE),
        keep.rownames = "pathway")
      z[, `:=`(model = model, axis = axis, method = "CAMERA",
                correlation_assumption = cor_mode,
                n_donors = ncol(vv$vo$E),
                favorable_direction = if (axis == "cog_z")
                  Direction == "Up" else Direction == "Down",
                FDR_within_50 = p.adjust(PValue, "BH"))]
      camera_list[[length(camera_list) + 1L]] <- z
    }
    fr <- as.data.table(limma::fry(vv$vo, index = vv$idx, design = design,
                                  contrast = coef_index, sort = "none"),
                        keep.rownames = "pathway")
    fr[, `:=`(model = model, axis = axis, method = "FRY",
              n_donors = ncol(vv$vo$E),
              favorable_direction = if (axis == "cog_z")
                Direction == "Up" else Direction == "Down",
              BH_directional = p.adjust(PValue, "BH"),
              BH_mixed = p.adjust(PValue.Mixed, "BH"))]
    fry_list[[length(fry_list) + 1L]] <- fr
  }
  if (model == "primary") {
    # The final 35 sensitivity used 9,999 rotations and seed 350031.
    set.seed(350031)
    mroast_table <- as.data.table(limma::mroast(
      vv$vo, index = vv$idx, design = design,
      contrast = which(colnames(design) == "cog_z"),
      nrot = 9999, sort = "none"), keep.rownames = "pathway")
    mroast_table[, `:=`(model = "primary", axis = "cog_z",
                       method = "MROAST", nrot = 9999L,
                       n_donors = ncol(vv$vo$E),
                       BH_directional = p.adjust(PValue, "BH"),
                       BH_mixed = p.adjust(PValue.Mixed, "BH"))]
  }
}
camera <- rbindlist(camera_list, fill = TRUE)
fry_results <- rbindlist(fry_list, fill = TRUE)
write_table(camera, "10_seaad_alllineages_pathology_cognition_camera_all_hallmarks.csv")
write_table(fry_results, "10_seaad_alllineages_pathology_cognition_fry_all_hallmarks.csv")
write_table(mroast_table, "10_seaad_alllineages_pathology_cognition_mroast_primary_cognition.csv")

ox <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
method_comparison <- rbindlist(list(
  hall_cog[cell_type == "Astrocyte" & pathway == ox,
           .(method = "FGSEA", model = "primary", FDR = padj)],
  camera[pathway == ox & model == "primary" & axis == "cog_z",
         .(method = paste0("CAMERA_", correlation_assumption), model,
           FDR = FDR_within_50)],
  fry_results[pathway == ox & model == "primary" & axis == "cog_z",
              .(method = "FRY_directional", model, FDR = BH_directional)],
  mroast_table[pathway == ox, .(method = "MROAST_directional", model,
                                FDR = BH_directional)]
), fill = TRUE)
write_table(method_comparison, "10_seaad_alllineages_pathology_cognition_oxphos_method_comparison.csv")

#Checkpoint for cross-cohort integration
saveRDS(list(
  sample_info = meta,
  primary_de = primary,
  sensitivity_de = sens,
  hallmark_cognition_inference = hall_cog,
  hallmark_pathology_inference = hall_path,
  hallmark_cognition_logFC = hall_cog_fc,
  hallmark_pathology_logFC = hall_path_fc,
  hallmark_cognition_vs_pathology = pair,
  hallmark_RIN_sensitivity = hall_rin,
  hallmark_memory_sensitivity = hall_memory,
  donor_score_models = score_fits,
  astrocyte_metabolic_scores = focus_scores,
  human_hallmark_sets = hs
), file.path(out, "10_seaad_alllineages_pathology_cognition_primary_checkpoint.rds"))
