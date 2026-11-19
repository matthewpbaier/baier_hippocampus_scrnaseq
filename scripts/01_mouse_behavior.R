#Baier et al.
#01: Behavioral and body-weight analyses

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "01_mouse_behavior")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(data.table)
library(nlme)

#Load parent 44-mouse cohort and identify the sequenced subset
raw_root <- hca_mouse_input_root
source_file <- file.path(raw_root, "metadata", "behavior.csv")
meta_file <- file.path(raw_root, "metadata", "sample_metadata.csv")
if (!file.exists(source_file)) stop("Missing full-cohort behavioral file: ", source_file)
if (!file.exists(meta_file)) stop("Missing sequencing metadata: ", meta_file)
raw <- fread(source_file)
meta13 <- fread(meta_file)

b <- raw[, .(
  mouse_id = as.character(mouse_id),
  age_source = as.character(age_group),
  status_source = as.character(subgroup_re80),
  initial_entries_80 = as.numeric(InitialEntriesToEighty),
  reversal_entries_80_coded = as.numeric(ReversalEntriesToEighty),
  cf51to61 = as.numeric(CF51to61),
  pre_weight = as.numeric(PreWeight),
  post_weight = as.numeric(PostWeight)
)]

b[, age_group := fifelse(grepl("young", tolower(age_source)), "Young", "Aged")]
b[, reversal_attained := is.finite(reversal_entries_80_coded) &
                          reversal_entries_80_coded < 6000]
b[, reversal_entries_80_observed := fifelse(reversal_attained,
                                           reversal_entries_80_coded, NA_real_)]
b[, group := fifelse(age_group == "Young", "Young",
                     fifelse(reversal_attained & reversal_entries_80_coded < 1000,
                             "Intact", "Impaired"))]
b[, group := factor(group, levels = c("Young", "Intact", "Impaired"))]
b[, sequenced := mouse_id %in% as.character(meta13$mouse_id)]
b[, weight_loss_g := pre_weight - post_weight]
b[, weight_loss_pct := 100 * weight_loss_g / pre_weight]  
#Positive = weight loss

stopifnot(nrow(b) == 44L, uniqueN(b$mouse_id) == 44L,
          nrow(meta13) == 13L, uniqueN(meta13$mouse_id) == 13L,
          sum(b$sequenced) == 13L,
          all(is.finite(b$initial_entries_80)),
          all(is.finite(b$cf51to61)),
          all(is.finite(b$reversal_entries_80_coded)),
          all(is.finite(b$pre_weight)), all(is.finite(b$post_weight)),
          all(b$pre_weight > 0), all(b$post_weight > 0),
          b[, .N, by = group][group == "Young", N] == 7L,
          b[, .N, by = group][group == "Intact", N] == 25L,
          b[, .N, by = group][group == "Impaired", N] == 12L,
          b[age_group == "Aged", sum(!reversal_attained)] == 8L,
          b[age_group == "Aged", sum(reversal_attained)] == 29L,
          b[sequenced == TRUE, .N, by = group][group == "Young", N] == 4L,
          b[sequenced == TRUE, .N, by = group][group == "Intact", N] == 4L,
          b[sequenced == TRUE, .N, by = group][group == "Impaired", N] == 5L,
          b[sequenced == TRUE & group == "Impaired", sum(!reversal_attained)] == 4L)

#Cross-check original group labels
aged_b <- b[age_group == "Aged"]
stopifnot(all(tolower(aged_b$status_source) %in% c("intact", "impaired")),
          identical(tolower(aged_b$status_source), tolower(as.character(aged_b$group))))

meta13[, mouse_id := as.character(mouse_id)]
check13 <- merge(b[sequenced == TRUE, .(mouse_id, group, pre_weight, post_weight)],
                 meta13[, .(mouse_id, primary_group, pre13 = pre_weight,
                            post13 = post_weight)], by = "mouse_id", sort = FALSE)
stopifnot(nrow(check13) == 13L,
          identical(tolower(as.character(check13$group)),
                    tolower(as.character(check13$primary_group))),
          all(abs(check13$pre_weight - check13$pre13) < 1e-8),
          all(abs(check13$post_weight - check13$post13) < 1e-8))

#Mouse-level behavioral table, classification counts, and descriptive summaries
#NA for non-attainers
fwrite(b[, .(mouse_id, age_group, group, sequenced,
              initial_entries_80, reversal_entries_80_coded,
              reversal_entries_80_observed, reversal_attained, cf51to61,
              pre_weight, post_weight, weight_loss_g, weight_loss_pct)],
       file.path(out, "01_baiermouse_allmice_behavior_mouse_level.csv"))

behavior_summary <- b[, .(
  n_mice = .N,
  n_sequenced = sum(sequenced),
  n_reversal_attained = sum(reversal_attained),
  n_reversal_not_attained = sum(!reversal_attained),
  mean_initial_entries = mean(initial_entries_80),
  median_initial_entries = median(initial_entries_80),
  median_reversal_attainers = if (any(reversal_attained))
    median(reversal_entries_80_observed, na.rm = TRUE) else NA_real_,
  mean_cf51to61 = mean(cf51to61),
  sem_cf51to61 = sd(cf51to61) / sqrt(.N)
), by = group]
fwrite(behavior_summary, file.path(out, "01_baiermouse_allmice_groups_behavior_summary.csv"))

#Behavioral correlations (CF, initial entries, reversal entries)
#Reversal correlations exclude non-attainers rather than using censored value as continuous
acq_data <- copy(b[group %in% c("Intact", "Impaired")])
acq_data[, group := droplevels(group)]
acq <- wilcox.test(initial_entries_80 ~ group, data = acq_data, exact = FALSE)
fwrite(data.table(comparison = "Aged intact vs impaired: acquisition entries",
                  test = "Two-sided Wilcoxon rank-sum", n_intact = 25L,
                  n_impaired = 12L, W = unname(acq$statistic), p = acq$p.value),
       file.path(out, "01_baiermouse_allmice_acquisition_wilcoxon_test.csv"))

cor_spec <- data.table(
  comparison = c("Acquisition vs reversal", "Acquisition vs cognitive flexibility",
                 "Reversal vs cognitive flexibility"),
  x = c("initial_entries_80", "initial_entries_80", "reversal_entries_80_observed"),
  y = c("reversal_entries_80_observed", "cf51to61", "cf51to61"),
  attainers_only = c(TRUE, FALSE, TRUE)
)
cor_results <- rbindlist(lapply(seq_len(nrow(cor_spec)), function(i) {
  s <- cor_spec[i]
  z <- if (s$attainers_only) b[age_group == "Aged" & reversal_attained] else b[age_group == "Aged"]
  z <- z[is.finite(get(s$x)) & is.finite(get(s$y))]
  test <- suppressWarnings(cor.test(z[[s$x]], z[[s$y]],
                                    method = "spearman", exact = FALSE))
  data.table(comparison = s$comparison, n = nrow(z),
             scope = if (s$attainers_only) "Aged reversal attainers" else "All aged mice",
             spearman_rho = unname(test$estimate), p = test$p.value)
}))
cor_results[, FDR_BH_three := p.adjust(p, method = "BH")]
stopifnot(identical(cor_results$n, c(29L, 37L, 29L)))
fwrite(cor_results, file.path(out, "01_baiermouse_allmice_behavior_spearman_correlations.csv"))

#Paired pre/post body weight w/ random-intercept model
long <- melt(b[, .(mouse_id, group, sequenced, pre_weight, post_weight)],
             id.vars = c("mouse_id", "group", "sequenced"),
             measure.vars = c("pre_weight", "post_weight"),
             variable.name = "time", value.name = "weight")
long[, time := factor(time, levels = c("pre_weight", "post_weight"),
                      labels = c("Pre-test", "Post-test"))]
long[, group := factor(group, levels = c("Young", "Intact", "Impaired"))]
stopifnot(nrow(long) == 88L, uniqueN(long$mouse_id) == 44L,
          !anyDuplicated(long, by = c("mouse_id", "time")))

#Sum contrasts and marginal F tests
contrasts(long$group) <- contr.sum(3)
contrasts(long$time) <- contr.sum(2)
weight_model <- nlme::lme(weight ~ group * time, random = ~ 1 | mouse_id,
                          data = long, method = "REML")
weight_anova <- as.data.table(anova(weight_model, type = "marginal"),
                              keep.rownames = "effect")
fwrite(weight_anova, file.path(out, "01_baiermouse_allmice_weight_repeated_measures_anova.csv"))

cells <- expand.grid(group = c("Young", "Intact", "Impaired"),
                     time = c("Pre-test", "Post-test"))
cells$group <- factor(cells$group, levels = levels(long$group))
cells$time <- factor(cells$time, levels = levels(long$time))
contrasts(cells$group) <- contr.sum(3)
contrasts(cells$time) <- contr.sum(2)
X <- model.matrix(~ group * time, data = cells)
coef <- nlme::fixef(weight_model)
V <- vcov(weight_model)
stopifnot(identical(colnames(X), names(coef)))

weight_contrast <- function(v) {
  estimate <- as.numeric(sum(v * coef))
  SE <- as.numeric(sqrt(drop(v %*% V %*% v)))
  df <- 41
  t <- estimate / SE
  data.table(estimate = estimate, SE = SE, df = df,
             CI_low = estimate - qt(.975, df) * SE,
             CI_high = estimate + qt(.975, df) * SE,
             t = t, p = 2 * pt(-abs(t), df))
}
within <- rbindlist(lapply(levels(long$group), function(g) {
  i <- which(cells$group == g & cells$time == "Pre-test")
  j <- which(cells$group == g & cells$time == "Post-test")
  cbind(data.table(group = g, contrast = "Post-test minus pre-test"),
        weight_contrast(X[j, ] - X[i, ]))
}))
within[, Holm_p_three_groups := p.adjust(p, method = "holm")]
fwrite(within, file.path(out, "01_baiermouse_allmice_weight_within_group_holm_tests.csv"))

between <- rbindlist(lapply(levels(long$time), function(tm) {
  rbindlist(lapply(combn(levels(long$group), 2, simplify = FALSE), function(g) {
    i <- which(cells$group == g[1] & cells$time == tm)
    j <- which(cells$group == g[2] & cells$time == tm)
    cbind(data.table(time = tm, group1 = g[1], group2 = g[2]),
          weight_contrast(X[j, ] - X[i, ]))
  }))
}))
between[, Holm_p_six := p.adjust(p, method = "holm")]
fwrite(between, file.path(out, "01_baiermouse_allmice_weight_between_groups_holm_tests.csv"))

#%weight loss Welch omnibus, pairwise Welch/Holm, and nonparametric sensitivity tests
wt <- b[, .(mouse_id, group, sequenced, pre_weight, post_weight,
            weight_loss_g, weight_loss_pct)]
wt_summary <- wt[, .(
  n = .N, mean_pre = mean(pre_weight), mean_post = mean(post_weight),
  mean_loss_pct = mean(weight_loss_pct),
  sem_loss_pct = sd(weight_loss_pct) / sqrt(.N),
  mean_loss_g = mean(weight_loss_g)
), by = group]
fwrite(wt_summary, file.path(out, "01_baiermouse_allmice_weight_group_means_sem.csv"))
fwrite(long, file.path(out, "01_baiermouse_allmice_weight_paired_observations.csv"))

welch <- oneway.test(weight_loss_pct ~ group, data = wt, var.equal = FALSE)
fwrite(data.table(test = "Welch one-way ANOVA",
                  F = unname(welch$statistic),
                  df_num = unname(welch$parameter[1]),
                  df_den = unname(welch$parameter[2]), p = welch$p.value),
       file.path(out, "01_baiermouse_allmice_weight_loss_welch_anova.csv"))

pairs <- combn(levels(wt$group), 2, simplify = FALSE)
welch_pair <- rbindlist(lapply(pairs, function(g) {
  x <- wt[group == g[1], weight_loss_pct]
  y <- wt[group == g[2], weight_loss_pct]
  t <- t.test(y, x)  # Welch (unequal variances); direction = group2 - group1
  data.table(group1 = g[1], group2 = g[2], n1 = length(x), n2 = length(y),
             difference_group2_minus_group1 = mean(y) - mean(x),
             CI_low = unname(t$conf.int[1]), CI_high = unname(t$conf.int[2]),
             p = t$p.value)
}))
welch_pair[, Holm_p_three_pairs := p.adjust(p, method = "holm")]
fwrite(welch_pair, file.path(out, "01_baiermouse_allmice_weight_loss_welch_holm_tests.csv"))

kw <- kruskal.test(weight_loss_pct ~ group, data = wt)
np <- rbindlist(lapply(pairs, function(g) {
  x <- wt[group == g[1], weight_loss_pct]
  y <- wt[group == g[2], weight_loss_pct]
  t <- wilcox.test(y, x, exact = FALSE)
  data.table(group1 = g[1], group2 = g[2], W = unname(t$statistic), p = t$p.value)
}))
np[, Holm_p_three_pairs := p.adjust(p, method = "holm")]
fwrite(data.table(test = "Kruskal-Wallis", statistic = unname(kw$statistic),
                  df = unname(kw$parameter), p = kw$p.value),
       file.path(out, "01_baiermouse_allmice_weight_loss_kruskal_sensitivity.csv"))
fwrite(np, file.path(out, "01_baiermouse_allmice_weight_loss_wilcoxon_sensitivity.csv"))

fwrite(check13, file.path(out, "01_baiermouse_sequencedmice_weight_source_concordance.csv"))
print(behavior_summary)
print(cor_results)
print(weight_anova)
print(welch_pair)
