#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
data_root     <- file.path(project_root, "figs/data")
plot_dir      <- file.path(project_root, "figs/fig-5/plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

isoseq_root   <- file.path(data_root, "nih/single_sample")

# TUSCO references
universal_gtf <- file.path(data_root, "tusco/tusco_mouse.gtf")
kidney_gtf    <- file.path(data_root, "tusco/tusco_mouse_kidney.gtf")

universal_tsv <- file.path(data_root, "tusco/tusco_mouse.tsv")
kidney_tsv    <- file.path(data_root, "tusco/tusco_mouse_kidney.tsv")
brain_tsv     <- file.path(data_root, "tusco/tusco_mouse_brain.tsv")

# Outputs
fn_decile_png_u  <- file.path(plot_dir, "FN_by_length_decile_Kidney_universal.png")
fn_decile_png_t  <- file.path(plot_dir, "FN_by_length_decile_Kidney_tissue.png")
fn_decile_png_bu <- file.path(plot_dir, "FN_by_length_decile_Brain_universal.png")
fn_decile_png_bt <- file.path(plot_dir, "FN_by_length_decile_Brain_tissue.png")

read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  readr::read_tsv(file_path, show_col_types = FALSE, ...)
}

extract_attr <- function(attr, key) {
  m <- stringr::str_match(attr, paste0(key, " \"([^\"]+)\""))
  out <- m[,2]
  out <- stringr::str_remove(out, "\\.\\d+$")
  out
}

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
  # TUSCO is single-isoform per gene by design; pick longest if multiple
  gene_lengths <- tx_lengths %>%
    group_by(gene_id) %>% arrange(desc(exonic_bases)) %>% slice_head(n = 1) %>% ungroup() %>%
    select(gene_id, transcript_id, chrom, strand, exonic_bases, tx_span_bases)
  list(gtf = gtf, tx_lengths = tx_lengths, gene_lengths = gene_lengths)
}

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
  map_df <- tusco_ref_df %>% transmute(key = .data[[top_id]], gene_id = ensembl) %>% distinct()
  keep2 <- keep %>% left_join(map_df, by = c("associated_gene" = "key")) %>% filter(!is.na(gene_id))
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
    mutate(
      has_TP = !is.na(has_TP) & has_TP,
      has_FP = !is.na(has_FP) & has_FP
    )
}

# Sample classification tables
all_sample_dirs <- list.dirs(isoseq_root, full.names = TRUE, recursive = FALSE)
kidney_sample_dirs <- all_sample_dirs[grepl("/K[0-9]+\\.isoforms$", all_sample_dirs)]
brain_sample_dirs  <- all_sample_dirs[grepl("/B[0-9]+\\.isoforms$", all_sample_dirs)]
if (length(kidney_sample_dirs) == 0) stop("No kidney single_sample sample directories found.")
if (length(brain_sample_dirs) == 0) stop("No brain single_sample sample directories found.")

load_class_for_sample <- function(sample_dir) {
  sample_name <- basename(sample_dir)
  class_file  <- file.path(sample_dir, paste0(sample_name, "_classification.txt"))
  if (!file.exists(class_file)) stop("Missing classification for ", sample_name, ": ", class_file)
  read_tsv_safe(class_file) %>% mutate(sample = sample_name)
}

kidney_class_list <- lapply(kidney_sample_dirs, load_class_for_sample)
kidney_class_df   <- bind_rows(kidney_class_list)
brain_class_list  <- lapply(brain_sample_dirs, load_class_for_sample)
brain_class_df    <- bind_rows(brain_class_list)

# Load TUSCO sets
uni <- load_tusco_gtf_lengths(universal_gtf)
kid <- load_tusco_gtf_lengths(kidney_gtf)
uni_ref   <- load_tusco_reference_tsv(universal_tsv)
kid_ref   <- load_tusco_reference_tsv(kidney_tsv)
brain_ref <- load_tusco_reference_tsv(brain_tsv)

# Brain gene lengths derived from universal GTF (no dedicated brain GTF provided)
brain_gene_lengths <- uni$gene_lengths %>% filter(gene_id %in% brain_ref$ensembl)

# Per-sample detection enriched with lengths
per_set_detection <- function(class_df, tusco_gene_lengths, tusco_ref_df) {
  det_list <- class_df %>% group_by(sample) %>% group_split()
  out <- lapply(det_list, function(sdf) {
    outs <- compute_gene_level_outcomes_flexible(sdf, tusco_ref_df, tusco_gene_lengths)
    outs$sample <- sdf$sample[1]
    outs
  })
  det <- bind_rows(out) %>% left_join(tusco_gene_lengths, by = c("gene_id"))
  det %>% mutate(length_kb = as.numeric(exonic_bases)/1000, detected = has_TP)
}

kidney_det_univ <- per_set_detection(kidney_class_df, uni$gene_lengths, uni_ref)
kidney_det_tiss <- per_set_detection(kidney_class_df, kid$gene_lengths, kid_ref)
brain_det_univ  <- per_set_detection(brain_class_df,  uni$gene_lengths, uni_ref)
brain_det_tiss  <- per_set_detection(brain_class_df,  brain_gene_lengths, brain_ref)

# Aggregate FNs by length deciles
fn_deciles <- function(df) {
  d <- df %>% filter(is.finite(length_kb)) %>% mutate(FN = !detected)
  d <- d %>% mutate(decile = ntile(exonic_bases, 10))
  agg <- d %>% group_by(decile) %>% summarise(
    n = n(), fn = sum(FN), rate = fn / n,
    .groups = "drop"
  )
  agg
}

plot_fn_decile_count <- function(agg, title) {
  ggplot(agg, aes(x = decile, y = fn)) +
    geom_col(fill = "#636363", width = 0.7) +
    geom_text(aes(label = fn), vjust = -0.4, size = 2.7) +
    scale_x_continuous(breaks = 1:10, labels = 1:10, limits = c(1, 10)) +
    labs(x = "Length decile (short -> long)", y = "FN count", title = title) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5))
}

fn_u  <- fn_deciles(kidney_det_univ)
fn_t  <- fn_deciles(kidney_det_tiss)
fn_bu <- fn_deciles(brain_det_univ)
fn_bt <- fn_deciles(brain_det_tiss)

p_fn_u_counts  <- plot_fn_decile_count(fn_u,  "Kidney – Universal set")
p_fn_t_counts  <- plot_fn_decile_count(fn_t,  "Kidney – Tissue set")
p_fn_bu_counts <- plot_fn_decile_count(fn_bu, "Brain – Universal set")
p_fn_bt_counts <- plot_fn_decile_count(fn_bt, "Brain – Tissue set")

# Save PNGs
suppressMessages(ggsave(fn_decile_png_u,  p_fn_u_counts,  width = 4, height = 3, dpi = 300))
suppressMessages(ggsave(fn_decile_png_t,  p_fn_t_counts,  width = 4, height = 3, dpi = 300))
suppressMessages(ggsave(fn_decile_png_bu, p_fn_bu_counts, width = 4, height = 3, dpi = 300))
suppressMessages(ggsave(fn_decile_png_bt, p_fn_bt_counts, width = 4, height = 3, dpi = 300))

message("Saved: ", fn_decile_png_u)
message("Saved: ", fn_decile_png_t)
message("Saved: ", fn_decile_png_bu)
message("Saved: ", fn_decile_png_bt)

# ---------------------------------------------
# Line boxplot: per-sample FN counts per decile
# Three boxes per decile: TUSCO mouse, TUSCO mouse kidney, TUSCO mouse brain
# ---------------------------------------------

per_sample_fn_counts <- function(class_df, tusco_gene_lengths, tusco_ref_df, set_label) {
  det <- per_set_detection(class_df, tusco_gene_lengths, tusco_ref_df)
  det <- det %>% mutate(FN = !detected, decile = ntile(exonic_bases, 10))
  out <- det %>%
    group_by(sample, decile) %>%
    summarise(fn_count = sum(FN), .groups = "drop") %>%
    mutate(set = set_label)
  out
}

# Build datasets for boxplots
box_universal <- per_sample_fn_counts(
  bind_rows(kidney_class_df, brain_class_df),
  uni$gene_lengths, uni_ref, "TUSCO mouse"
)
box_kidney <- per_sample_fn_counts(
  kidney_class_df, kid$gene_lengths, kid_ref, "TUSCO mouse kidney"
)
box_brain <- per_sample_fn_counts(
  brain_class_df, brain_gene_lengths, brain_ref, "TUSCO mouse brain"
)

box_df <- bind_rows(box_universal, box_kidney, box_brain) %>%
  mutate(
    decile = as.integer(decile),
    set = recode(set,
      `TUSCO mouse` = "Universal",
      `TUSCO mouse kidney` = "Kidney",
      `TUSCO mouse brain` = "Brain"
    ),
    set = factor(set, levels = c("Universal", "Brain", "Kidney"))
  )

fn_boxplot_png <- file.path(plot_dir, "FN_by_length_decile_boxplot.png")
fn_boxplot_svg <- file.path(plot_dir, "FN_by_length_decile_boxplot.svg")

palette_sets <- c(
  "Universal" = "#1C9E77",
  "Brain"     = "#5893a4",
  "Kidney"    = "#506f68"
)

# Summaries per decile and set across samples
summary_df <- box_df %>%
  group_by(set, decile) %>%
  summarise(
    n = n(),
    mean_fn = mean(fn_count, na.rm = TRUE),
    sd_fn = sd(fn_count, na.rm = TRUE),
    se_fn = sd_fn / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(ymin = pmax(0, mean_fn - sd_fn), ymax = mean_fn + sd_fn)

pd <- position_dodge(width = 0.6)

p_box <- ggplot(summary_df, aes(x = factor(decile), y = mean_fn, color = set, group = set)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = pd, width = 0.22, linewidth = 0.5) +
  geom_line(position = pd, linewidth = 0.8) +
  geom_point(position = pd, size = 2.2) +
  scale_color_manual(values = palette_sets, drop = FALSE) +
  scale_x_discrete(breaks = as.character(1:10), labels = 1:10) +
  labs(
    x = "Length decile (short -> long)",
    y = "Mean FN per sample (± SD)",
    title = "FN by length decile (per-sample)",
    color = ""
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

suppressMessages(ggsave(fn_boxplot_png, p_box, width = 6, height = 3.8, dpi = 300))
suppressMessages(ggsave(fn_boxplot_svg, p_box, width = 6, height = 3.8))

message("Saved: ", fn_boxplot_png)
message("Saved: ", fn_boxplot_svg)

# ---------------------------------------------
# Additional plot: total transcripts per decile
# Similar style (lines + dots), counts from TUSCO sets only
# ---------------------------------------------

count_tusco_deciles <- function(gene_lengths_df, set_label) {
  gene_lengths_df %>%
    filter(is.finite(exonic_bases)) %>%
    mutate(decile = ntile(exonic_bases, 10)) %>%
    count(decile, name = "count") %>%
    mutate(set = set_label)
}

all_tx_counts <- bind_rows(
  count_tusco_deciles(uni$gene_lengths, "Universal"),
  count_tusco_deciles(brain_gene_lengths, "Brain"),
  count_tusco_deciles(kid$gene_lengths, "Kidney")
) %>% mutate(
  decile = as.integer(decile),
  set = factor(set, levels = c("Universal", "Brain", "Kidney"))
)

all_tx_png <- file.path(plot_dir, "All_transcripts_by_length_decile.png")
all_tx_svg <- file.path(plot_dir, "All_transcripts_by_length_decile.svg")

p_all <- ggplot(all_tx_counts, aes(x = factor(decile), y = count, group = set, color = set)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = palette_sets, drop = FALSE) +
  scale_x_discrete(breaks = as.character(1:10), labels = 1:10) +
  labs(
    x = "Length decile (short -> long)",
    y = "Number of transcripts",
    title = "Total transcripts by length decile",
    color = ""
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

suppressMessages(ggsave(all_tx_png, p_all, width = 6, height = 3.8, dpi = 300))
suppressMessages(ggsave(all_tx_svg, p_all, width = 6, height = 3.8))

message("Saved: ", all_tx_png)
message("Saved: ", all_tx_svg)

# ---------------------------------------------
# Combined plot: transcripts and FN in one figure
# Faceted rows, shared x, free y
# ---------------------------------------------

fn_for_plot <- summary_df %>%
  transmute(set, decile = as.integer(decile), metric = "FN mean ± SD",
            y = mean_fn, ymin = ymin, ymax = ymax)

tx_for_plot <- all_tx_counts %>%
  transmute(set, decile = as.integer(decile), metric = "Total transcripts",
            y = count, ymin = NA_real_, ymax = NA_real_)

combined_df <- bind_rows(tx_for_plot, fn_for_plot) %>%
  mutate(
    set = factor(set, levels = c("Universal", "Brain", "Kidney")),
    metric = factor(metric, levels = c("Total transcripts", "FN mean ± SD"))
  )

combined_png <- file.path(plot_dir, "Transcripts_and_FN_by_length_decile.png")
combined_svg <- file.path(plot_dir, "Transcripts_and_FN_by_length_decile.svg")

pd2 <- position_dodge(width = 0.6)

p_combined <- ggplot(combined_df, aes(x = factor(decile), y = y, color = set, group = set)) +
  geom_errorbar(
    data = dplyr::filter(combined_df, metric == "FN mean ± SD"),
    aes(ymin = ymin, ymax = ymax), position = pd2, width = 0.22, linewidth = 0.5
  ) +
  geom_line(position = pd2, linewidth = 0.8) +
  geom_point(position = pd2, size = 2.2) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = palette_sets, drop = FALSE) +
  scale_x_discrete(breaks = as.character(1:10), labels = 1:10) +
  labs(
    x = "Length decile (short -> long)",
    y = NULL,
    title = "Transcripts and FN by length decile",
    color = ""
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 6)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0)
  )

suppressMessages(ggsave(combined_png, p_combined, width = 6, height = 6.2, dpi = 300))
suppressMessages(ggsave(combined_svg, p_combined, width = 6, height = 6.2))

message("Saved: ", combined_png)
message("Saved: ", combined_svg)
