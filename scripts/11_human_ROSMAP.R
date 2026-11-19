#Baier et al.
#11: ROSMAP multiregion astrocyte pathology-associated transcription (Mathys 2024, Nature)

source(file.path(Sys.getenv("HCA_CODE_ROOT", unset = getwd()), "00_paths.R"))
root <- hca_output_root
dir.create(root, recursive = TRUE, showWarnings = FALSE)
out <- file.path(root, "11_human_ROSMAP")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

#Load libraries
library(data.table)
library(fgsea)
library(readxl)
library(msigdbr)

set.seed(20261003)
N_MATCHED <- 2000L               
N_NEAREST <- 15L                  
REGIONS <- c("HC", "EC", "MT", "AG", "PFC", "TH")
OUTCOMES <- c("nrad", "nft") 

#Published Mathys et al. astrocyte DE summary
source_name <- "aggregated_fullset.Ast_Ast.tsv.gz"
source_file <- file.path(hca_rosmap_source_root, source_name)

#Build shared gene references directly from Script 06/08/09 outputs
ref_dir <- file.path(root, "reference")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
map_file <- file.path(ref_dir, "11_rosmap_astrocyte_pathology_regions_mouse_human_strict_1to1_orthologs.csv")
length_file <- file.path(root, "08_astrocyte_states", "08_mouse_mm39_gene_lengths_USED.csv")
mouse_de_file <- file.path(root, "06_mouse_bulk", "06_baiermouse_alllineages_cognition_twogroup_gene_effects.csv")
hallmark_file <- file.path(root, "09_human_normative_aging", "09_zemke_alllineages_aging_human_hallmark_sets.rds")

mouse_de <- fread(mouse_de_file)
mouse_symbols <- unique(as.character(mouse_de$gene))
mouse_symbols <- mouse_symbols[!is.na(mouse_symbols) & nzchar(mouse_symbols) &
                               !grepl("^ENSMUSG", mouse_symbols)]
raw_map <- as.data.table(babelgene::orthologs(
  genes = mouse_symbols, species = "mouse", human = FALSE,
  min_support = 3, top = FALSE))
stopifnot(all(c("symbol", "human_symbol") %in% names(raw_map)))
orth <- unique(raw_map[!is.na(symbol) & !is.na(human_symbol), .(
  mouse_gene = toupper(trimws(symbol)),
  human_gene = toupper(trimws(human_symbol)))])
orth[, `:=`(n_human = uniqueN(human_gene)), by = mouse_gene]
orth[, n_mouse := uniqueN(mouse_gene), by = human_gene]
orth[, strict_one_to_one := n_human == 1L & n_mouse == 1L]
fwrite(orth[, .(mouse_gene, human_gene, strict_one_to_one)], map_file)
orth <- orth[strict_one_to_one == TRUE, .(mouse_gene, human_gene,
                                          strict_one_to_one)]

hallmarks <- readRDS(hallmark_file)
#Independent human metabolic modules from original MitoCarta3.0 and MSigDB,
# using human source annotations
mc_file <- hca_human_mitocarta_file
if (!file.exists(mc_file))
  stop("Missing human MitoCarta3.0 reference: ", mc_file)
mc_sheet <- grep("^A Human", readxl::excel_sheets(mc_file), value = TRUE)[1]
mc <- as.data.table(readxl::read_excel(mc_file, sheet = mc_sheet))
mc_path <- grep("MitoPathways$", names(mc), value = TRUE)[1]
mc <- mc[!is.na(get(mc_path)) & get(mc_path) != "0",
         .(gene = toupper(as.character(Symbol)), path = as.character(get(mc_path)))]
mc <- mc[, .(path = trimws(unlist(strsplit(path, "\\|")))), by = gene]
mito_patterns <- list(
  RESP_CHAIN = "^OXPHOS", TCA = "TCA cycle", PYRUVATE = "Pyruvate",
  BETA_OX = "Fatty acid oxidation|Carnitine", KETONE = "Ketone",
  MITO_TRANSL = "central dogma|ribosome|Translation",
  MITO_IMPORT = "Protein import")
mito_sets <- lapply(mito_patterns, function(pattern)
  unique(mc[grepl(pattern, path), gene]))
reactome <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                            subcollection = "CP:REACTOME")
reactome <- split(reactome$gene_symbol, reactome$gs_name)
go <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5",
                      subcollection = "GO:BP")
go <- split(go$gene_symbol, go$gs_name)
kegg <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                         subcollection = "CP:KEGG_LEGACY")
kegg <- split(kegg$gene_symbol, kegg$gs_name)
pick <- function(collection, patterns) {
  matches <- unique(unlist(lapply(patterns, grep, x = names(collection), value = TRUE)))
  unique(toupper(unlist(collection[matches], use.names = FALSE)))
}
nonmito <- list(
  PEROX_LIPID = unique(c(pick(reactome, "^REACTOME_PEROXISOMAL_LIPID_METABOLISM"),
                         pick(go, "BETA_OXIDATION_USING_ACYL_COA_OXIDASE"))),
  PEROX_BIOGEN = pick(reactome, "^REACTOME_PEROXISOMAL_PROTEIN_IMPORT"),
  LIPID_STORE = unique(c(pick(reactome, "^REACTOME_TRIGLYCERIDE_BIOSYNTHESIS"),
                         pick(go, c("^GOBP_LONG_CHAIN_FATTY_ACID_TRANSPORT$",
                                    "^GOBP_LIPID_DROPLET_ORGANIZATION$")))),
  GLYCOLYSIS = unique(c(pick(reactome, "^REACTOME_GLYCOLYSIS$"),
                        pick(kegg, "^KEGG_GLYCOLYSIS_GLUCONEOGENESIS$"))),
  CHOLESTEROL = pick(reactome, "^REACTOME_CHOLESTEROL_BIOSYNTHESIS$"),
  ROS_DETOX = pick(reactome, "^REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES$"),
  GLUTATHIONE = unique(c(pick(reactome, "^REACTOME_GLUTATHIONE_CONJUGATION$"),
                         pick(kegg, "^KEGG_GLUTATHIONE_METABOLISM$"))))
mito_members <- rbindlist(lapply(names(mito_sets), function(k)
  data.table(module = k, gene = mito_sets[[k]],
             priority = match(k, names(mito_sets)))))
setorder(mito_members, gene, priority)
mito_members <- unique(mito_members, by = "gene")[, .(module, gene)]
nonmito_members <- rbindlist(lapply(names(nonmito), function(k)
  data.table(module = k, gene = setdiff(nonmito[[k]], mito_members$gene))))
module <- unique(rbindlist(list(mito_members, nonmito_members)))
module_sets <- split(module$gene, module$module)
antiox <- toupper(c(
  "Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx4", "Prdx1", "Prdx2",
  "Prdx3", "Prdx4", "Prdx5", "Prdx6", "Txn1", "Txn2", "Txnrd1",
  "Txnrd2", "Gsr", "Gclc", "Gclm", "Gss", "Mt1", "Mt2", "Mt3",
  "Glrx", "Glrx5", "Selenow", "Park7", "Srxn1", "Nqo1", "Hmox1"
))

ms <- fread(source_file)
source_columns <- c("gene", "path", "region", "logFC_nb", "log10p_nm")
if (!all(source_columns %in% names(ms)))
  stop("Unexpected published ROSMAP file schema: missing ",
       paste(setdiff(source_columns, names(ms)), collapse = ", "))

ms <- unique(ms[is.finite(logFC_nb)], by = c("gene", "path", "region"))
if (!nrow(ms)) stop("No valid published ROSMAP astrocyte DEG rows")
ms[, `:=`(
  fav_lfc = -logFC_nb,
  fav_stat = -sign(logFC_nb) *
    pmin(fifelse(is.finite(log10p_nm), log10p_nm, 50), 50)
)]

KEY <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "RESP_CHAIN", "TCA",
         "BETA_OX", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME",
         "ROS_DETOX", "GLUTATHIONE", "ANTIOXIDANT_PANEL", "GLYCOLYSIS",
         "CHOLESTEROL", "MITO_TRANSL", "CYTO_RIBOSOME")
S <- c(lapply(hallmarks, function(v) toupper(v)), module_sets,
       list(ANTIOXIDANT_PANEL = antiox,
            CYTO_RIBOSOME = grep("^RP[SL][0-9]+[AXY]?[0-9]*$",
                                  unique(toupper(ms$gene)), value = TRUE)))
S <- lapply(S[KEY], function(x) unique(x[!is.na(x) & !grepl("^MT-", x)]))
stopifnot(all(lengths(S[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                          "RESP_CHAIN", "BETA_OX")]) >= 5L))

# Record the 13 candidate gene sets after MT-gene exclusions. FGSEA applies
# the 5-gene minimum to each observed ranking, so fewer may be tested.
fwrite(rbindlist(lapply(names(S), function(p)
  data.table(pathway = p, gene = sort(unique(S[[p]]))))),
  file.path(out, "11_rosmap_astrocyte_pathology_regions_13_tested_gene_set_memberships.csv"))

#Gene-length/expression covariates
#mm39 mouse gene span + mouse astrocyte two-group logCPM -> human 1:1 ortholog
orth <- fread(map_file)
if (!all(c("mouse_gene", "human_gene", "strict_one_to_one") %in% names(orth)))
  stop("Ortholog table requires mouse_gene, human_gene, strict_one_to_one")
orth <- orth[as.character(strict_one_to_one) %in% c("TRUE", "T", "1")]
if (anyDuplicated(orth$mouse_gene) || anyDuplicated(orth$human_gene))
  stop("1:1 ortholog mapping is not unique")

if (!file.exists(length_file)) stop("Missing Script 08 gene lengths: ", length_file)
gl <- fread(length_file)
if (!all(c("gene", "length") %in% names(gl)))
  stop("Mouse gene-length table requires gene and length")
gl[, gene := toupper(trimws(as.character(gene)))]
if (anyDuplicated(gl$gene) || any(!is.finite(gl$length) | gl$length <= 0))
  stop("Mouse gene-length table contains duplicated genes or invalid lengths")
expr <- fread(mouse_de_file)
if (!all(c("lineage", "gene", "logCPM") %in% names(expr)))
  stop("Mouse two-group DE table requires lineage, gene and logCPM")
expr <- expr[lineage == "Astro-TE NN", .(gene = toupper(trimws(as.character(gene))), logCPM)]
if (!nrow(expr)) stop("Missing original Astro-TE NN mouse logCPM")
COV <- merge(merge(gl[, .(gene, length)], expr, by = "gene"),
             orth[, .(gene = mouse_gene, human_gene)], by = "gene")
COV <- COV[, .(gene = human_gene, length = as.numeric(length),
               logCPM = as.numeric(logCPM))]
COV <- COV[is.finite(length) & length > 0 & is.finite(logCPM)]
if (anyDuplicated(COV$gene)) stop("Duplicate human gene in matched covariates")
if (nrow(COV) < 3000L) stop("Unexpectedly few matched mouse/human gene covariates")

#Rank and enrich favorable pathology-associated effects in each region
binz <- function(s, g) {
  ave(s, g, FUN = function(v) {
    if (length(v) > 2L && stats::sd(v) > 0)
      (v - mean(v)) / stats::sd(v)
    else v * 0
  })
}
#Record the actual sorted, deduplicated vectors passed to FGSEA
ranked_inputs <- list()
fgsea_call_order <- list()
run_gsea <- function(stat, sets, path_id, region_id, ranking_id) {
  ranks <- sort(stat[is.finite(stat)], decreasing = TRUE)
  ranks <- ranks[!duplicated(names(ranks))]
  call_id <- length(ranked_inputs) + 1L
  ranked_inputs[[call_id]] <<- data.table(
    call_id = call_id, path = path_id, region = region_id,
    ranking = ranking_id, rank_position = seq_along(ranks),
    gene = names(ranks), statistic = as.numeric(ranks))
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = sets, stats = ranks, minSize = 5, maxSize = 1000,
    eps = 0
  ))
  fgsea_call_order[[call_id]] <<- data.table(
    call_id = call_id, path = path_id, region = region_id,
    ranking = ranking_id, seed_at_script_start = 20261003L,
    minSize = 5L, maxSize = 1000L, eps = 0,
    gene_count = length(ranks), candidate_pathways = length(sets),
    tested_pathways = nrow(fg),
    package_version = as.character(utils::packageVersion("fgsea")))
  as.data.table(fg)[, .(pathway, NES, padj, size)]
}

#Length/expression-matched background test
matched_test <- function(log_fc, genes) {
  d <- merge(data.table(gene = names(log_fc), lfc = unname(log_fc)),
             COV, by = "gene")
  d <- d[is.finite(lfc)]
  d[, `:=`(zl = as.numeric(scale(log10(length))),
           ze = as.numeric(scale(logCPM)))]
  inside <- d$gene %in% genes
  if (sum(inside) < 5L) return(NULL)
  A <- d[inside]
  B <- d[!inside]
  if (nrow(B) < N_NEAREST) stop("Too few ROSMAP background genes")
  nearest <- lapply(seq_len(nrow(A)), function(i) {
    head(order((B$zl - A$zl[i])^2 + (B$ze - A$ze[i])^2), N_NEAREST)
  })
  draws <- replicate(N_MATCHED, {
    chosen <- vapply(nearest, function(ix) ix[sample.int(length(ix), 1L)], 1L)
    median(B$lfc[chosen])
  })
  observed <- median(A$lfc)
  data.table(n = nrow(A), obs_median = observed,
             excess = observed - median(draws),
             frac_favorable = mean(A$lfc > 0),
             p_emp_upper = (1 + sum(draws >= observed)) / (1 + N_MATCHED))
}

split_ms <- split(ms, by = c("path", "region"), keep.by = TRUE)
F1 <- rbindlist(lapply(split_ms, function(x) {
  message("ROSMAP: ", x$path[1], " / ", x$region[1])
  stat <- setNames(x$fav_stat, x$gene)
  lfc <- setNames(x$fav_lfc, x$gene)
  d <- merge(x[, .(gene, fav_stat)], COV, by = "gene")
  if (nrow(d) < 500L) stop("Too few length-matched genes for ", x$region[1])
  d[, lb := cut(rank(log10(length), ties.method = "first"),
                10, labels = FALSE)]
  paths <- rbind(
    run_gsea(stat, S, x$path[1L], x$region[1L], "raw")[,
      ranking := "raw"],
    run_gsea(setNames(binz(d$fav_stat, d$lb), d$gene), S,
             x$path[1L], x$region[1L], "length_binned")[,
      ranking := "length_binned"]
  )
  matches <- rbindlist(lapply(names(S), function(k) {
    z <- matched_test(lfc, S[[k]])
    if (is.null(z)) NULL else z[, pathway := k]
  }), fill = TRUE)
  merge(paths, matches, by = "pathway", all.x = TRUE)[,
    `:=`(path = x$path[1], region = x$region[1])]
}), fill = TRUE)
if (!nrow(F1) || anyDuplicated(F1, by = c("path", "region", "pathway", "ranking")))
  stop("ROSMAP enrichment results empty or duplicated")
fwrite(rbindlist(ranked_inputs),
  file.path(out, "11_rosmap_astrocyte_pathology_regions_exact_fgsea_ranked_inputs.csv.gz"))
fwrite(rbindlist(fgsea_call_order),
  file.path(out, "11_rosmap_astrocyte_pathology_regions_fgsea_call_order_and_parameters.csv"))
fwrite(F1, file.path(out, "11_rosmap_astrocyte_pathology_regions_pathway_region_results.csv"))
# Script 12 also consumes this short, legacy-compatible filename.
fwrite(F1, file.path(out, "11_ROSMAP_pathway_region_results.csv"))

#Pooled and individual-region outputs
SHOW <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          "HALLMARK_FATTY_ACID_METABOLISM", "BETA_OX")
plot_source <- F1[path %in% OUTCOMES & pathway %in% SHOW &
                  ranking %in% c("raw", "length_binned")]
stopifnot(nrow(plot_source) > 0L)

#Include explicit missingness (e.g., NFT in thalamus)
grid <- CJ(path = OUTCOMES, region = c(REGIONS, "allregions"),
           pathway = SHOW, ranking = c("raw", "length_binned"),
           unique = TRUE)
coverage <- merge(grid, plot_source, by = c("path", "region", "pathway", "ranking"),
                  all.x = TRUE)
coverage[, available := is.finite(NES)]
fwrite(coverage, file.path(out, "11_rosmap_astrocyte_pathology_regions_fig5_ed7_source.csv"))

#Region consistency and leave-one-region-out summaries
regional <- F1[path %in% OUTCOMES & region %in% REGIONS &
                 pathway %in% SHOW & ranking %in% c("raw", "length_binned")]
regional_summary <- regional[, .(
  n_regions = .N,
  n_positive = sum(NES > 0, na.rm = TRUE),
  n_negative = sum(NES < 0, na.rm = TRUE),
  mean_NES = mean(NES, na.rm = TRUE),
  median_NES = median(NES, na.rm = TRUE),
  min_NES = min(NES, na.rm = TRUE),
  max_NES = max(NES, na.rm = TRUE),
  n_FDR05 = sum(padj < .05, na.rm = TRUE),
  hippocampus_NES = NES[match("HC", region)]
), by = .(path, pathway, ranking)]
fwrite(regional_summary, file.path(out, "11_rosmap_astrocyte_pathology_regions_region_summary.csv"))

omissions <- rbindlist(lapply(split(regional, by = c("path", "pathway", "ranking")),
  function(x) rbindlist(lapply(intersect(REGIONS, x$region), function(drop) {
    z <- x[region != drop]
    if (!nrow(z)) return(NULL)
    data.table(path = x$path[1], pathway = x$pathway[1],
               ranking = x$ranking[1], omitted_region = drop,
               n_remaining = nrow(z), mean_NES = mean(z$NES, na.rm = TRUE),
               median_NES = median(z$NES, na.rm = TRUE),
               n_favorable = sum(z$NES > 0, na.rm = TRUE))
  }))))
fwrite(omissions, file.path(out, "11_rosmap_astrocyte_pathology_regions_region_omission.csv"))

omission_summary <- omissions[, .(
  min_omission_mean_NES = min(mean_NES),
  max_omission_mean_NES = max(mean_NES),
  all_omission_means_positive = all(mean_NES > 0)
), by = .(path, pathway, ranking)]
fwrite(omission_summary, file.path(out, "11_rosmap_astrocyte_pathology_regions_region_omission_summary.csv"))

pooled <- F1[path %in% OUTCOMES & region == "allregions" &
               pathway %in% SHOW & ranking %in% c("raw", "length_binned")]
pooled_vs_region <- merge(pooled, regional_summary, by = c("path", "pathway", "ranking"),
                          all.x = TRUE, suffixes = c("", "_regional"))
fwrite(pooled_vs_region, file.path(out, "11_rosmap_astrocyte_pathology_regions_pooled_vs_regions.csv"))

#Integration-ready favorable pathway statistics
matched_pathology <- F1[path == "nrad" & region == "allregions" &
                          pathway %in% SHOW & ranking == "raw",
                        .(pathway, NES, padj, size, n, excess,
                          p_emp_upper, frac_favorable)]
fwrite(matched_pathology, file.path(out, "11_rosmap_astrocyte_pathology_regions_pooled_pathology_for_integration.csv"))

#Rounded benchmark values
audit_expected <- data.table(
  path = c("nrad", "nft"),
  reference_NES = c(1.56, 1.81),
  reference_FDR = c(0.036, 0.019)
)
audit_observed <- F1[
  pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION" &
    region == "allregions" & ranking == "raw" & path %in% audit_expected$path,
  .(path, observed_NES = NES, observed_FDR = padj)
]
audit <- merge(audit_expected, audit_observed, by = "path", all.x = TRUE)
audit[, `:=`(NES_difference = observed_NES - reference_NES,
             FDR_difference = observed_FDR - reference_FDR)]
fwrite(audit, file.path(out, "11_rosmap_astrocyte_pathology_regions_manuscript_nes_audit.csv"))

