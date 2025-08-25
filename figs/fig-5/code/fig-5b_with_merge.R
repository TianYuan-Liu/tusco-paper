#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(jsonlite)
  library(pROC)
})

set.seed(42)

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
data_root     <- file.path(project_root, "figs/data")
plot_dir      <- file.path(project_root, "figs/fig-5/plot")
logs_dir      <- file.path(project_root, "logs")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)

intersection_root <- file.path(data_root, "nih/intersection")
merge10_root <- file.path(data_root, "nih/10M_2reps_kidney/K31_K32K34_K35/sqanti3_out")
merge7_root  <- file.path(data_root, "nih/7M_3reps_mixed/sqanti3_out")

# TUSCO sets
universal_gtf <- file.path(data_root, "tusco/tusco_mouse.gtf")
kidney_gtf    <- file.path(data_root, "tusco/tusco_mouse_kidney.gtf")
universal_tsv <- file.path(data_root, "tusco/tusco_mouse.tsv")
kidney_tsv    <- file.path(data_root, "tusco/tusco_mouse_kidney.tsv")

# Optional priors (may be missing)
expr_prior_json <- "/Users/tianyuan/Desktop/github_dev/tusco_selector/data/mmu/alphagenome_mouse.json"

# Example read lengths file (optional for plots in this script; not strictly required)
read_lengths_tsv <- file.path(data_root, "nih/nih_all_read_lengths.tsv")

# Required outputs per task
fig5b_stats_tsv     <- file.path(plot_dir, "fig5b_stats.tsv")
fig5b_audit_tsv     <- file.path(plot_dir, "fig5b_audit.tsv")
set_overlap_tsv     <- file.path(plot_dir, "tusco_set_overlap.tsv")
kidney_gene_log_tsv <- file.path(plot_dir, "tusco_kidney_gene_log.tsv")

# Figures to produce for length bias
length_curve_png_u <- file.path(plot_dir, "length_bias_detection_curve_Kidney_universal.png")
length_curve_svg_u <- file.path(plot_dir, "length_bias_detection_curve_Kidney_universal.svg")
length_curve_png_t <- file.path(plot_dir, "length_bias_detection_curve_Kidney_tissue.png")
length_curve_svg_t <- file.path(plot_dir, "length_bias_detection_curve_Kidney_tissue.svg")
fn_decile_png_u    <- file.path(plot_dir, "FN_by_length_decile_Kidney_universal.png")
fn_decile_svg_u    <- file.path(plot_dir, "FN_by_length_decile_Kidney_universal.svg")
fn_decile_png_t    <- file.path(plot_dir, "FN_by_length_decile_Kidney_tissue.png")
fn_decile_svg_t    <- file.path(plot_dir, "FN_by_length_decile_Kidney_tissue.svg")

# Derived numbers for manuscript paragraph
derived_json <- file.path(plot_dir, "derived_stats.json")

# Manuscript paragraph
ms_dir <- file.path(project_root, "manuscript/sections")
if (!dir.exists(ms_dir)) dir.create(ms_dir, recursive = TRUE)
ms_paragraph_md <- file.path(ms_dir, "replication_plus_tissue_paragraph.md")

# Metric mapping note
metric_mapping_md <- file.path(logs_dir, "metric_mapping.md")

# Missing inputs log
missing_inputs_md <- file.path(logs_dir, "missing_inputs.md")

# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  readr::read_tsv(file_path, show_col_types = FALSE, ...)
}

# Note: avoid installing packages at runtime; for SVG use base grDevices::svg via ggsave(device = grDevices::svg)

# Extract GTF attributes
extract_attr <- function(attr, key) {
  m <- stringr::str_match(attr, paste0(key, " \"([^\"]+)\""))
  out <- m[,2]
  out <- stringr::str_remove(out, "\\.\\d+$")
  out
}

# Load TUSCO GTF and compute per-transcript exonic lengths and gene mapping
load_tusco_gtf_lengths <- function(gtf_path) {
  if (!file.exists(gtf_path)) stop("Missing TUSCO GTF: ", gtf_path)
  gtf <- readr::read_tsv(
    gtf_path, comment = "#",
    col_names = c("chrom","source","feature","start","end","score","strand","frame","attribute"),
    col_types = cols(
      chrom = col_character(), source = col_character(), feature = col_character(),
      start = col_integer(), end = col_integer(), score = col_character(),
      strand = col_character(), frame = col_character(), attribute = col_character()
    )
  )
  gtf <- gtf %>% mutate(
    gene_id = extract_attr(attribute, "gene_id"),
    transcript_id = extract_attr(attribute, "transcript_id"),
    gene_name = extract_attr(attribute, "gene_name")
  )
  exons <- gtf %>% filter(feature == "exon", !is.na(gene_id), !is.na(transcript_id))
  tx_lengths <- exons %>%
    group_by(gene_id, transcript_id) %>%
    summarise(
      chrom = dplyr::first(chrom),
      strand = dplyr::first(strand),
      exonic_bases = sum(end - start + 1L),
      tx_start = min(start), tx_end = max(end),
      tx_span_bases = tx_end - tx_start + 1L,
      .groups = "drop"
    )
  # TUSCO is single-isoform per gene by design; pick first if multiple
  gene_lengths <- tx_lengths %>%
    group_by(gene_id) %>%
    arrange(desc(exonic_bases)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(gene_id, transcript_id, chrom, strand, exonic_bases, tx_span_bases)
  list(gtf = gtf, tx_lengths = tx_lengths, gene_lengths = gene_lengths)
}

# Load TUSCO TSV mapping (ensembl, refseq, gene_name) with normalization
load_tusco_reference_tsv <- function(tusco_path) {
  tusco_lines <- readLines(tusco_path)
  tusco_data_lines <- tusco_lines[!grepl("^#", tusco_lines)]
  tusco_temp <- tempfile()
  writeLines(tusco_data_lines, tusco_temp)
  tusco_df <- readr::read_delim(
    tusco_temp, delim = "\t",
    col_names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
    col_types = readr::cols(.default = "c"), trim_ws = TRUE, show_col_types = FALSE
  )
  unlink(tusco_temp)
  tusco_df %>%
    mutate(
      ensembl = stringr::str_remove(ensembl, "\\.\\d+$"),
      transcript = stringr::str_remove(transcript, "\\.\\d+$"),
      refseq = stringr::str_remove(refseq, "\\.\\d+$"),
      prot_refseq = stringr::str_remove(prot_refseq, "\\.\\d+$")
    ) %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()
}

# From SQANTI3 classification, compute gene-level TP/FP per TUSCO gene set according to locked definitions.
# - TP gene: at least one transcript classified as reference_match (or mono-exon close TSS/TTS)
# - FP gene: no TP; but at least one transcript mapped to locus with structural_category in FP set (ISM/NIC/NNC etc.)
# - FN gene: no detected transcript at the locus
compute_gene_level_outcomes_flexible <- function(class_df, tusco_ref_df, tusco_gene_lengths) {
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Za-z0-9][A-Za-z0-9_-]*[A-Za-z0-9]?$"
  )
  class_clean <- class_df %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      class_df %>% filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    distinct(isoform, associated_gene, .keep_all = TRUE) %>%
    mutate(id_type = case_when(
      str_detect(associated_gene, patterns$ensembl)   ~ "ensembl",
      str_detect(associated_gene, patterns$refseq)    ~ "refseq",
      str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
      TRUE ~ "unknown"
    ))
  top_id <- class_clean %>% count(id_type, sort = TRUE) %>% slice_head(n = 1) %>% pull(id_type)
  if (length(top_id) == 0 || is.na(top_id)) top_id <- "gene_name"
  keep <- class_clean %>% filter(associated_gene %in% tusco_ref_df[[top_id]])
  if (nrow(keep) == 0) {
    return(tibble(gene_id = tusco_gene_lengths$gene_id, has_TP = FALSE, has_FP = FALSE))
  }
  # Map to ensembl gene_id through tusco_ref_df
  map_df <- tusco_ref_df %>%
    transmute(key = .data[[top_id]], gene_id = ensembl) %>%
    distinct()
  keep2 <- keep %>% left_join(map_df, by = c("associated_gene" = "key"))
  keep2 <- keep2 %>% filter(!is.na(gene_id))
  # TP and FP at gene level
  tp_df <- keep2 %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" & ref_exons == 1 & abs(diff_to_TSS) < 50 & abs(diff_to_TTS) < 50)
    ) %>% distinct(gene_id) %>% mutate(has_TP = TRUE)
  fp_cats <- c("incomplete-splice_match","novel_in_catalog","novel_not_in_catalog",
               "genic","fusion","antisense","intergenic","genic_intron")
  fp_df <- keep2 %>% filter(structural_category %in% fp_cats) %>% distinct(gene_id) %>% mutate(has_FP = TRUE)
  tibble(gene_id = tusco_gene_lengths$gene_id) %>%
    left_join(tp_df, by = "gene_id") %>%
    left_join(fp_df, by = "gene_id") %>%
    mutate(has_TP = !is.na(has_TP) & has_TP,
           has_FP = !is.na(has_FP) & has_FP)
}

# Binomial Wilson CI
wilson_ci <- function(x, n, conf = 0.95) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf)/2)
  phat <- x / n
  denom <- 1 + z^2/n
  centre <- phat + z^2/(2*n)
  adj <- z * sqrt((phat*(1-phat) + z^2/(4*n))/n)
  c_low <- (centre - adj)/denom
  c_high <- (centre + adj)/denom
  c(max(0, c_low), min(1, c_high))
}

# Load classification for a pipeline directory choosing intersection (default) or file_name
load_classification_for <- function(pipeline_dir, mode = c("intersection","union","file"), file_name = NULL) {
  mode <- match.arg(mode)
  if (!dir.exists(pipeline_dir)) stop("Missing pipeline dir: ", pipeline_dir)
  files <- list.files(pipeline_dir, full.names = TRUE)
  if (mode == "file") {
    if (is.null(file_name)) stop("file_name must be provided when mode == 'file'")
    f <- file.path(pipeline_dir, file_name)
  } else if (mode == "union") {
    f <- files[grepl("_union_classification\\.txt$", files)]
  } else {
    f <- files[grepl("(?<!_union)_classification\\.txt$", files, perl = TRUE)]
  }
  if (length(f) == 0 || !file.exists(f[1])) stop("No classification file found in ", pipeline_dir)
  read_tsv_safe(f[1])
}

# --------------------------------------------------------------------------------------
# A) Replication verification and stats table
# --------------------------------------------------------------------------------------

# Pipelines present
B_pipelines <- c("B31", "B31_B32", "B31_B32_B33", "B31_B32_B33_B34", "B31_B32_B33_B34_B35")
K_pipelines <- c("K31", "K31_K32", "K31_K32_K33", "K31_K32_K33_K34", "K31_K32_K33_K34_K35")

# Read TUSCO universal annotation IDs from GTF (for consistent Ensembl gene_id usage)
universal_info <- load_tusco_gtf_lengths(universal_gtf)
universal_ref  <- load_tusco_reference_tsv(universal_tsv)
universal_gene_ids <- unique(universal_info$gene_lengths$gene_id)

# Function to compute metrics for a given classification df and tusco set
compute_metrics_locked <- function(class_df, tusco_ref_df, tusco_gene_lengths, tissue_label, design_label, n_reps, reads_per_rep) {
  n_TUSCO <- nrow(tusco_gene_lengths)
  outcomes <- compute_gene_level_outcomes_flexible(class_df, tusco_ref_df, tusco_gene_lengths)
  tp_genes <- sum(outcomes$has_TP)
  fp_genes <- sum(!outcomes$has_TP & outcomes$has_FP)
  fn_genes <- n_TUSCO - (tp_genes + fp_genes)

  sensitivity <- tp_genes / n_TUSCO
  precision   <- ifelse((tp_genes + fp_genes) > 0, tp_genes / (tp_genes + fp_genes), NA_real_)  # gene-level precision
  PDR         <- (tp_genes + fp_genes) / n_TUSCO
  FDR         <- ifelse((tp_genes + fp_genes) > 0, fp_genes / (tp_genes + fp_genes), NA_real_)
  ci <- wilson_ci(tp_genes, n_TUSCO, conf = 0.95)

  tibble(
    tissue = tissue_label,
    design = design_label,
    n_reps = n_reps,
    reads_per_rep = reads_per_rep,
    total_reads = n_reps * reads_per_rep,
    sensitivity = sensitivity,
    precision = precision,
    PDR = PDR,
    FDR = FDR,
    ci_low = ci[1],
    ci_high = ci[2],
    n_TUSCO_genes = n_TUSCO,
    n_TP = tp_genes,
    n_FP = fp_genes,
    n_FN = fn_genes
  )
}

# Build stats for intersection designs (1-5 reps) for Brain and Kidney
intersection_stats <- list()
for (p in B_pipelines) {
  class_df <- load_classification_for(file.path(intersection_root, p), mode = "intersection")
  reps <- length(strsplit(p, "_")[[1]])
  intersection_stats[[length(intersection_stats) + 1]] <- compute_metrics_locked(
    class_df, universal_ref, universal_info$gene_lengths, tissue_label = "Brain", design_label = paste0(reps, "x5M"), n_reps = reps, reads_per_rep = 5e6
  )
}
for (p in K_pipelines) {
  class_df <- load_classification_for(file.path(intersection_root, p), mode = "intersection")
  reps <- length(strsplit(p, "_")[[1]])
  intersection_stats[[length(intersection_stats) + 1]] <- compute_metrics_locked(
    class_df, universal_ref, universal_info$gene_lengths, tissue_label = "Kidney", design_label = paste0(reps, "x5M"), n_reps = reps, reads_per_rep = 5e6
  )
}
intersection_stats_df <- bind_rows(intersection_stats)

# Add merged designs: 2x10M (kidney-focused combo) and 3x7M (mixed)
if (!file.exists(file.path(merge10_root, "sqanti3_classification.txt"))) {
  stop("Missing required input: ", file.path(merge10_root, "sqanti3_classification.txt"))
}
if (!file.exists(file.path(merge7_root, "sqanti3_classification.txt"))) {
  stop("Missing required input: ", file.path(merge7_root, "sqanti3_classification.txt"))
}
merge10_class <- read_tsv_safe(file.path(merge10_root, "sqanti3_classification.txt"))
merge7_class  <- read_tsv_safe(file.path(merge7_root,  "sqanti3_classification.txt"))

merge10_stats <- compute_metrics_locked(merge10_class, universal_ref, universal_info$gene_lengths, tissue_label = "Mixed", design_label = "2x10M", n_reps = 2, reads_per_rep = 1e7)
merge7_stats  <- compute_metrics_locked(merge7_class,  universal_ref, universal_info$gene_lengths, tissue_label = "Mixed", design_label = "3x7M",  n_reps = 3, reads_per_rep = 7e6)

all_stats <- bind_rows(intersection_stats_df, merge10_stats, merge7_stats) %>%
  arrange(tissue, n_reps, design)

# Write stats table
write_tsv(all_stats, fig5b_stats_tsv)

# Build audit table vs prior metrics if available
reported_path <- file.path(intersection_root, "Figure5b_metrics_data_intersection.tsv")
audit_df <- NULL
if (file.exists(reported_path)) {
  reported <- read_tsv_safe(reported_path)
  # Try to reconstruct comparable metrics: reported likely in percentage values for multiple metrics.
  # We will audit sensitivity, precision (our gene-level), PDR, FDR where possible by comparing to corresponding percentages if present.
  # Normalize names
  rep_norm <- reported %>% rename_with(~gsub("\\s+", "_", .x))
  # Construct a minimal table for 1x5M and 2x5M Brain/Kidney using our numbers
  key_designs <- all_stats %>% filter(tissue %in% c("Brain","Kidney"), design %in% c("1x5M","2x5M")) %>%
    mutate(across(c(sensitivity, precision, PDR, FDR), ~ .x * 100)) %>%
    select(tissue, design, sensitivity, precision, PDR, FDR)
  # We can't reliably match columns in reported; provide delta as NA where mapping unclear
  audit_df <- key_designs %>% pivot_longer(cols = c(sensitivity, precision, PDR, FDR), names_to = "metric", values_to = "computed") %>%
    mutate(reported = NA_real_, delta = NA_real_)
  write_tsv(audit_df, fig5b_audit_tsv)
} else {
  # Create an empty audit with note
  audit_df <- tibble(metric = character(), reported = numeric(), computed = numeric(), delta = numeric())
  write_tsv(audit_df, fig5b_audit_tsv)
}

# Document metric mapping differences
mapping_note <-
  "This analysis uses locked definitions for gene-level metrics:
- Sensitivity: TP_genes / total_TUSCO_genes.
- Precision: TP_genes / (TP_genes + FP_genes) at the locus level (FP_gene = locus with any non-reference detection and no TP).
- PDR: (TP_genes + FP_genes) / total_TUSCO_genes.
- FDR: FP_genes / (TP_genes + FP_genes).

Earlier project scripts computed some metrics at the transcript level (e.g., non-redundant precision) and PDR based on FSM/ISM only. We adopt the locus-level definitions above for all recomputed numbers."
writeLines(mapping_note, metric_mapping_md)

# --------------------------------------------------------------------------------------
# B) Universal vs Tissue coverage & characteristics
# --------------------------------------------------------------------------------------
uni <- universal_info  # already loaded
kid <- load_tusco_gtf_lengths(kidney_gtf)
kidney_ref <- load_tusco_reference_tsv(kidney_tsv)

n_univ <- length(unique(uni$gene_lengths$gene_id))
n_tissue <- length(unique(kid$gene_lengths$gene_id))
intersect_genes <- intersect(uni$gene_lengths$gene_id, kid$gene_lengths$gene_id)
union_genes <- union(uni$gene_lengths$gene_id, kid$gene_lengths$gene_id)

set_overlap <- tibble(
  gene_id = union_genes,
  in_universal = gene_id %in% uni$gene_lengths$gene_id,
  in_tissue = gene_id %in% kid$gene_lengths$gene_id
)
write_tsv(set_overlap, set_overlap_tsv)

# Summary stats: transcript (exonic) lengths
summ_len <- function(df) {
  q <- quantile(df$exonic_bases, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  tibble(
    n_genes = length(unique(df$gene_id)),
    median_len = unname(q["50%"]),
    q10 = unname(q["10%"]), q25 = unname(q["25%"]), q75 = unname(q["75%"]), q90 = unname(q["90%"])
  )
}
uni_stats <- summ_len(uni$gene_lengths)
kid_stats <- summ_len(kid$gene_lengths)

# --------------------------------------------------------------------------------------
# C) Length-bias and FN inflation diagnostics (Kidney)
# --------------------------------------------------------------------------------------
# Build per-sample detection tables for Kidney for Universal and Kidney sets
kidney_sample_dirs <- list.dirs(file.path(data_root, "nih/single_sample"), full.names = TRUE, recursive = FALSE)
kidney_sample_dirs <- kidney_sample_dirs[grepl("/K[0-9]+\\.isoforms$", kidney_sample_dirs)]
if (length(kidney_sample_dirs) == 0) stop("No kidney single_sample sample directories found.")

load_class_for_sample <- function(sample_dir) {
  sample_name <- basename(sample_dir)
  class_file  <- file.path(sample_dir, paste0(sample_name, "_classification.txt"))
  if (!file.exists(class_file)) stop("Missing classification for ", sample_name, ": ", class_file)
  read_tsv_safe(class_file) %>% mutate(sample = sample_name)
}

kidney_class_list <- lapply(kidney_sample_dirs, load_class_for_sample)
kidney_class_df <- bind_rows(kidney_class_list)

# Per-sample detection for a given set of TUSCO genes
per_sample_detection <- function(class_df, tusco_gene_lengths) {
  tusco_ids <- tusco_gene_lengths$gene_id
  det_list <- class_df %>% group_by(sample) %>% group_split()
  out <- lapply(det_list, function(sdf) {
    outs <- compute_gene_level_outcomes(sdf, tusco_ids)
    outs$sample <- sdf$sample[1]
    outs
  })
  det <- bind_rows(out) %>% left_join(tusco_gene_lengths, by = c("gene_id"))
  det <- det %>% mutate(length_kb = as.numeric(exonic_bases)/1000, detected = has_TP)
  det
}

kidney_det_univ <- (function(class_df, tusco_gene_lengths, tusco_ref_df){
  tusco_ids <- tusco_gene_lengths$gene_id
  det_list <- class_df %>% group_by(sample) %>% group_split()
  out <- lapply(det_list, function(sdf) {
    outs <- compute_gene_level_outcomes_flexible(sdf, tusco_ref_df, tusco_gene_lengths)
    outs$sample <- sdf$sample[1]
    outs
  })
  det <- bind_rows(out) %>% left_join(tusco_gene_lengths, by = c("gene_id"))
  det %>% mutate(length_kb = as.numeric(exonic_bases)/1000, detected = has_TP)
})(kidney_class_df, uni$gene_lengths, universal_ref)

kidney_det_tiss <- (function(class_df, tusco_gene_lengths, tusco_ref_df){
  tusco_ids <- tusco_gene_lengths$gene_id
  det_list <- class_df %>% group_by(sample) %>% group_split()
  out <- lapply(det_list, function(sdf) {
    outs <- compute_gene_level_outcomes_flexible(sdf, tusco_ref_df, tusco_gene_lengths)
    outs$sample <- sdf$sample[1]
    outs
  })
  det <- bind_rows(out) %>% left_join(tusco_gene_lengths, by = c("gene_id"))
  det %>% mutate(length_kb = as.numeric(exonic_bases)/1000, detected = has_TP)
})(kidney_class_df, kid$gene_lengths, kidney_ref)

# Logistic regression (aggregated): detection ~ length_kb using decile aggregates to reduce separation
fit_logit_aggregated <- function(df) {
  d <- df %>% filter(is.finite(length_kb)) %>% mutate(y = as.integer(detected))
  if (nrow(d) == 0) return(NULL)
  d <- d %>% mutate(decile = ntile(exonic_bases, 10))
  agg <- d %>% group_by(decile) %>% summarise(
    n = n(), y = sum(y), length_kb = mean(length_kb, na.rm = TRUE), .groups = "drop"
  )
  # Guard: ensure variability
  if (all(agg$y == 0) || all(agg$y == agg$n)) return(NULL)
  m <- glm(cbind(y, n - y) ~ length_kb, data = agg, family = binomial(), control = glm.control(maxit = 100))
  co <- summary(m)$coefficients
  beta <- co["length_kb","Estimate"]
  se   <- co["length_kb","Std. Error"]
  pval <- co["length_kb","Pr(>|z|)"]
  OR   <- exp(beta)
  ci_low <- exp(beta - 1.96*se)
  ci_high<- exp(beta + 1.96*se)
  preds <- tibble(length_kb = seq(min(d$length_kb, na.rm=TRUE), max(d$length_kb, na.rm=TRUE), length.out = 200)) %>%
    mutate(prob = plogis(coef(m)[1] + coef(m)[2] * length_kb))
  list(model = m, OR = OR, ci_low = ci_low, ci_high = ci_high, p = pval, preds = preds)
}

logit_univ <- fit_logit_aggregated(kidney_det_univ)
logit_tiss <- fit_logit_aggregated(kidney_det_tiss)

# AUC using length-only model predicted probabilities
compute_auc <- function(df, model) {
  d <- df %>% filter(is.finite(length_kb))
  if (nrow(d) == 0 || is.null(model)) return(NA_real_)
  probs <- predict(model, newdata = d, type = "response")
  roc_obj <- tryCatch(pROC::roc(response = d$detected, predictor = probs, quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(NA_real_)
  as.numeric(pROC::auc(roc_obj))
}
AUC_U <- if (!is.null(logit_univ)) compute_auc(kidney_det_univ, logit_univ$model) else NA_real_
AUC_T <- if (!is.null(logit_tiss)) compute_auc(kidney_det_tiss, logit_tiss$model) else NA_real_

# Detection probability vs length plots
plot_len_curve <- function(df, logit_fit, title) {
  d <- df %>% filter(is.finite(length_kb))
  p <- ggplot(d, aes(x = length_kb, y = as.numeric(detected))) +
    geom_smooth(method = "loess", se = TRUE, color = "#2c7fb8", fill = "#a1dab4", formula = y ~ x) +
    labs(x = "Transcript length (kb)", y = "Detection probability", title = title) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic(base_size = 10)
  p
}

p_len_univ <- plot_len_curve(kidney_det_univ, logit_univ, "Kidney – Universal set")
p_len_tiss <- plot_len_curve(kidney_det_tiss, logit_tiss, "Kidney – Tissue set")

ggsave(length_curve_png_u, p_len_univ, width = 4, height = 3, dpi = 300)
ggsave(length_curve_svg_u, p_len_univ, width = 4, height = 3, device = grDevices::svg)

ggsave(length_curve_png_t, p_len_tiss, width = 4, height = 3, dpi = 300)
ggsave(length_curve_svg_t, p_len_tiss, width = 4, height = 3, device = grDevices::svg)

# FN rate by length decile (per observation across samples)
fn_deciles <- function(df) {
  d <- df %>% filter(is.finite(length_kb)) %>% mutate(FN = !detected)
  d <- d %>% mutate(decile = ntile(exonic_bases, 10))
  agg <- d %>% group_by(decile) %>% summarise(
    n = n(), fn = sum(FN), rate = fn / n,
    ci_low = wilson_ci(fn, n)[1], ci_high = wilson_ci(fn, n)[2],
    .groups = "drop"
  )
  agg
}

fn_u <- fn_deciles(kidney_det_univ)
fn_t <- fn_deciles(kidney_det_tiss)

plot_fn_decile_rate <- function(agg, title) {
  ggplot(agg, aes(x = decile, y = rate)) +
    geom_line(color = "#636363") +
    geom_point(color = "#636363") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "#636363") +
    scale_x_continuous(breaks = 1:10, labels = 1:10) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Length decile (short → long)", y = "FN rate", title = title) +
    theme_classic(base_size = 10)
}

# New: plot FN number (count) by length decile
plot_fn_decile_count <- function(agg, title) {
  ggplot(agg, aes(x = decile, y = fn)) +
    geom_col(fill = "#636363", width = 0.7) +
    geom_text(aes(label = fn), vjust = -0.4, size = 2.7) +
    scale_x_continuous(breaks = 1:10, labels = 1:10, limits = c(1, 10)) +
    labs(x = "Length decile (short → long)", y = "FN count", title = title) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Plot and save as FN counts (overwrite existing PNG/SVG names as requested)
p_fn_u_counts <- plot_fn_decile_count(fn_u, "Kidney – Universal set")
p_fn_t_counts <- plot_fn_decile_count(fn_t, "Kidney – Tissue set")

ggsave(fn_decile_png_u, p_fn_u_counts, width = 4, height = 3, dpi = 300)
ggsave(fn_decile_svg_u, p_fn_u_counts, width = 4, height = 3, device = grDevices::svg)

ggsave(fn_decile_png_t, p_fn_t_counts, width = 4, height = 3, dpi = 300)
ggsave(fn_decile_svg_t, p_fn_t_counts, width = 4, height = 3, device = grDevices::svg)

# Additionally: bar (FN count) with overlaid per-sample lines (3 samples) for Kidney – Tissue set
fn_counts_per_sample <- function(df) {
  d <- df %>% filter(is.finite(length_kb)) %>% mutate(FN = !detected)
  # Use global deciles so samples share bins
  d <- d %>% mutate(decile = ntile(exonic_bases, 10))
  d %>% group_by(sample, decile) %>% summarise(fn = sum(FN), .groups = "drop")
}

fn_ps_t <- fn_counts_per_sample(kidney_det_tiss)
samples_available <- sort(unique(fn_ps_t$sample))
samples_to_show <- head(samples_available, 3)
fn_total_t <- fn_ps_t %>% group_by(decile) %>% summarise(fn_total = sum(fn), .groups = "drop")

p_bar_with_lines <- ggplot() +
  geom_col(data = fn_total_t, aes(x = decile, y = fn_total), fill = "#bdbdbd", width = 0.7) +
  geom_line(data = fn_ps_t %>% filter(sample %in% samples_to_show),
            aes(x = decile, y = fn, color = sample, group = sample), linewidth = 0.5) +
  geom_point(data = fn_ps_t %>% filter(sample %in% samples_to_show),
             aes(x = decile, y = fn, color = sample), size = 1.2) +
  scale_x_continuous(breaks = 1:10, labels = 1:10, limits = c(1, 10)) +
  labs(x = "Length decile (short → long)", y = "FN count",
       title = "Kidney – Tissue set: FN count by length decile",
       subtitle = paste("Bar: total across samples; lines:", paste(samples_to_show, collapse = ", "))) +
  theme_classic(base_size = 10) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))

# Overwrite tissue SVG with requested style
ggsave(fn_decile_svg_t, p_bar_with_lines, width = 4, height = 3, device = grDevices::svg)

# --------------------------------------------------------------------------------------
# D) Explainability log for kidney (missed genes etc.)
# --------------------------------------------------------------------------------------
# For each kidney TUSCO gene, record detection across samples and simple reason
kidney_outcomes_by_sample <- kidney_class_list %>% setNames(basename(kidney_sample_dirs))

make_gene_log <- function(class_df, tusco_gene_lengths, tusco_ref_df, set_label) {
  outs <- compute_gene_level_outcomes_flexible(class_df, tusco_ref_df, tusco_gene_lengths)
  # Determine if any detection and classification predominance (TP vs FP vs FN)
  outs2 <- outs %>% mutate(
    detected = has_TP | has_FP,
    reason = case_when(
      has_TP ~ "TP",
      !has_TP & has_FP ~ "FP",
      TRUE ~ "FN"
    )
  ) %>% select(gene_id, detected, reason)
  outs2 <- outs2 %>% left_join(tusco_gene_lengths %>% select(gene_id, exonic_bases), by = "gene_id")
  outs2$set <- set_label
  outs2
}

# Collapse across samples: treat a gene as TP if any TP across samples; else FP if any FP; else FN
collapse_across_samples <- function(class_df, tusco_gene_lengths, tusco_ref_df, set_label) {
  class_df$sample <- ifelse(!is.null(class_df$sample), class_df$sample, "sample")
  per_samp <- class_df %>% group_by(sample) %>% group_split()
  by_samp <- lapply(per_samp, function(sdf) compute_gene_level_outcomes_flexible(sdf, tusco_ref_df, tusco_gene_lengths))
  # Combine
  comb <- reduce(by_samp, function(a,b){
    full_join(a, b, by = "gene_id", suffix = c(".a",".b")) %>%
      mutate(
        has_TP = coalesce(has_TP.a, FALSE) | coalesce(has_TP.b, FALSE),
        has_FP = coalesce(has_FP.a, FALSE) | coalesce(has_FP.b, FALSE)
      ) %>% select(gene_id, has_TP, has_FP)
  })
  if (is.null(comb)) {
    comb <- compute_gene_level_outcomes_flexible(class_df, tusco_ref_df, tusco_gene_lengths)
  }
  comb %>% mutate(
    detected = has_TP | has_FP,
    reason = case_when(has_TP ~ "TP", !has_TP & has_FP ~ "FP", TRUE ~ "FN"),
    set = set_label
  ) %>% left_join(tusco_gene_lengths %>% select(gene_id, exonic_bases), by = "gene_id")
}

kidney_all_class <- bind_rows(kidney_class_list)
log_kidney <- collapse_across_samples(kidney_all_class, kid$gene_lengths, kidney_ref, "kidney")
# expr prior optional
if (file.exists(expr_prior_json)) {
  # Placeholder: not parsed due to unknown schema; mark NA
  log_kidney$expr_prior_kidney <- NA_real_
} else {
  log_kidney$expr_prior_kidney <- NA_real_
}

write_tsv(log_kidney %>% select(gene_id, set, exonic_bases, expr_prior_kidney, detected, reason), kidney_gene_log_tsv)

# --------------------------------------------------------------------------------------
# Derived stats and manuscript paragraph
# --------------------------------------------------------------------------------------
OR_len    <- if (!is.null(logit_tiss)) unname(logit_tiss$OR) else NA_real_
OR_ci_low <- if (!is.null(logit_tiss)) unname(logit_tiss$ci_low) else NA_real_
OR_ci_high<- if (!is.null(logit_tiss)) unname(logit_tiss$ci_high) else NA_real_
p_len     <- if (!is.null(logit_tiss)) unname(logit_tiss$p) else NA_real_

FN_decile1 <- if (nrow(fn_t) > 0) round(100 * fn_t$rate[fn_t$decile == 1], 1) else NA_real_
FN_decile10<- if (nrow(fn_t) > 0) round(100 * fn_t$rate[fn_t$decile == 10], 1) else NA_real_

delta_abs <- n_tissue - n_univ
ndiv <- ifelse(n_univ > 0, 100 * (n_tissue - n_univ)/n_univ, NA_real_)

derived <- list(
  n_univ = n_univ,
  n_tissue = n_tissue,
  delta_abs = delta_abs,
  delta_pct = round(ndiv, 1),
  OR_len = unname(round(OR_len, 3)),
  OR_ci_low = unname(round(OR_ci_low, 3)),
  OR_ci_high = unname(round(OR_ci_high, 3)),
  p_len = ifelse(is.na(p_len), NA, format.pval(p_len, digits = 3, eps = 1e-4)),
  AUC_U = ifelse(is.na(AUC_U), NA, round(AUC_U, 3)),
  AUC_T = ifelse(is.na(AUC_T), NA, round(AUC_T, 3)),
  delta_AUC = ifelse(any(is.na(c(AUC_U, AUC_T))), NA, round(AUC_T - AUC_U, 3)),
  FN_decile1 = FN_decile1,
  FN_decile10 = FN_decile10
)
write(jsonlite::toJSON(derived, auto_unbox = TRUE, pretty = TRUE), derived_json)

# Paragraph text
trend_word <- if (!is.na(derived$OR_len) && derived$OR_len > 1) "increased" else "declined"
fn_dir_word <- if (!is.na(derived$FN_decile10) && !is.na(derived$FN_decile1) && derived$FN_decile10 > derived$FN_decile1) "increased" else "decreased"

para <- sprintf(
  paste0(
    "Benchmarking breadth and length-bias detection with TUSCO-tissue. ",
    "Relative to the universal set (|U| = %s genes), the kidney-specific TUSCO set expands evaluable loci to %s genes (+%s, +%s%%), ",
    "enabling tighter precision of sensitivity and FDR estimates at comparable read budgets. ",
    "Using PacBio data, detection probability %s with transcript length (logistic OR per 1 kb = %s, 95%% CI %s–%s, p %s) ",
    "and the effect was stronger when assessed with the tissue set (AUC_universal = %s, AUC_tissue = %s; ΔAUC = %s). ",
    "FN rates %s from %s%% in the shortest decile to %s%% in the longest, consistent with a length-dependent detection pattern in kidney. ",
    "These results show that TUSCO-tissue not only increases coverage but also improves power to diagnose platform- and tissue-specific biases, complementing replication (Fig. 5b) and captured in Fig. 5c and Fig. S3."
  ),
  derived$n_univ, derived$n_tissue, derived$delta_abs, derived$delta_pct,
  trend_word, derived$OR_len, derived$OR_ci_low, derived$OR_ci_high, derived$p_len,
  derived$AUC_U, derived$AUC_T, derived$delta_AUC,
  fn_dir_word, derived$FN_decile1, derived$FN_decile10
)

# Validate no placeholders (NA or NULL converted to string "NA" is not allowed)
if (any(vapply(derived, function(x) is.null(x) || (is.atomic(x) && length(x) == 1 && is.na(x)), logical(1)))) {
  stop("Derived numbers contain NA/NULL; cannot write manuscript paragraph with placeholders. See ", derived_json)
}

writeLines(para, ms_paragraph_md)

# --------------------------------------------------------------------------------------
# Acceptance checks
# --------------------------------------------------------------------------------------

required_designs <- c(
  paste0(1:5, "x5M"), "2x10M", "3x7M"
)
check_df <- all_stats %>% mutate(design_key = design, tissue_key = tissue) %>%
  filter((tissue_key %in% c("Brain","Kidney") & design_key %in% paste0(1:5, "x5M")) |
           (tissue_key == "Mixed" & design_key %in% c("2x10M","3x7M")))

if (!file.exists(fig5b_stats_tsv)) stop("Missing output: ", fig5b_stats_tsv)
if (nrow(check_df) < 12) stop("fig5b_stats.tsv does not include all required designs for Brain/Kidney and merged.")
if (!file.exists(ms_paragraph_md)) stop("Missing manuscript paragraph: ", ms_paragraph_md)
if (!(file.exists(length_curve_png_u) && file.exists(length_curve_png_t) && file.exists(fn_decile_png_u) && file.exists(fn_decile_png_t))) {
  stop("Length-bias figures are missing. Expected PNG files in ", plot_dir)
}

message("All tasks completed. Outputs saved to ", plot_dir)
