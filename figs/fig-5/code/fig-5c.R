suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(grid)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
data_root <- file.path(project_root, "figs/data")

# TUSCO reference (mouse, multi-exon)
tusco_file <- file.path(data_root, "tusco/tusco_mouse_multi_exon.tsv")

# Your classification inputs (Iso-Seq AR)
isoseq_ar_root <- file.path(data_root, "nih/isoseq_ar")

# Output directory and filename (match fig-5 widths)
plot_dir <- file.path(project_root, "figs/fig-5/plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
plot_pdf <- file.path(plot_dir, "Figure5c_TUSCO_TP_PTP_FP_FN_isoseq_ar.pdf")
plot_png <- file.path(plot_dir, "Figure5c_TUSCO_TP_PTP_FP_FN_isoseq_ar.png")
plot_tsv <- file.path(plot_dir, "Figure5c_TUSCO_TP_PTP_FP_FN_isoseq_ar.tsv")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  tryCatch(read_tsv(file_path, show_col_types = FALSE, ...),
           error = function(e) stop("Unable to read ", file_path, ": ", e$message))
}

# Load TUSCO reference into a clean table of IDs (ensembl, refseq, gene_name)
load_tusco_reference <- function(tusco_path) {
  tusco_lines <- readLines(tusco_path)
  tusco_data_lines <- tusco_lines[!grepl("^#", tusco_lines)]
  tusco_temp_file <- tempfile()
  writeLines(tusco_data_lines, tusco_temp_file)
  tusco_df <- read_delim(
    tusco_temp_file, delim = "\t",
    col_names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
    col_types = cols(.default = "c"), trim_ws = TRUE, show_col_types = FALSE
  )
  unlink(tusco_temp_file)
  tusco_df %>%
    mutate(
      ensembl = str_remove(ensembl, "\\.\\d+$"),
      transcript = str_remove(transcript, "\\.\\d+$"),
      refseq = str_remove(refseq, "\\.\\d+$"),
      prot_refseq = str_remove(prot_refseq, "\\.\\d+$")
    ) %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()
}

# Compute TP/PTP/FP/FN against TUSCO reference genes for a classification data.frame
compute_tusco_bigcats <- function(classification_data, tusco_ref_df) {
  # Patterns to detect the ID type present in associated_gene
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Za-z0-9][A-Za-z0-9_-]*[A-Za-z0-9]?$"
  )

  # Expand fusions and normalize IDs
  class_clean <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>% filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
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

  # Choose the dominant ID type
  top_id_type <- class_clean %>% count(id_type, sort = TRUE) %>% slice_head(n = 1) %>% pull(id_type)
  if (length(top_id_type) == 0 || is.na(top_id_type)) top_id_type <- "gene_name"

  # Keep only classification rows that match some TUSCO gene using the selected ID type
  class_tusco <- class_clean %>% filter(associated_gene %in% tusco_ref_df[[top_id_type]])

  # Define categories
  TP_df <- class_tusco %>% filter(subcategory == "reference_match")
  PTP_df <- class_tusco %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match"),
           !(associated_gene %in% TP_df$associated_gene))

  # False negatives: TUSCO genes of the selected ID type not seen as TP or PTP
  found_ids <- union(TP_df$associated_gene, PTP_df$associated_gene)
  ref_ids <- unique(tusco_ref_df[[top_id_type]])
  FN_count <- sum(!(ref_ids %in% found_ids))

  # False positives: remaining structural categories among matched genes
  FP_df <- class_tusco %>%
    filter(structural_category %in% c(
      "novel_in_catalog", "novel_not_in_catalog", "genic",
      "fusion", "antisense", "intergenic", "genic_intron"
    ))

  tibble(
    TP = nrow(TP_df),
    PTP = nrow(PTP_df),
    FP = nrow(FP_df),
    FN = FN_count,
    id_type = top_id_type
  )
}

# Infer tissue from sample directory name (B* -> Brain, K* -> Kidney)
infer_tissue <- function(sample_basename) {
  if (startsWith(sample_basename, "B")) return("Brain")
  if (startsWith(sample_basename, "K")) return("Kidney")
  return("Unknown")
}

# ------------------------------------------------------------
# Discover samples and compute metrics
# ------------------------------------------------------------
all_sample_dirs <- list.dirs(isoseq_ar_root, full.names = TRUE, recursive = FALSE)
# Keep only B??.isoforms or K??.isoforms
all_sample_dirs <- all_sample_dirs[grepl("/(B|K)[0-9]+\\.isoforms$", all_sample_dirs)]

if (length(all_sample_dirs) == 0) {
  stop("No sample directories found in ", isoseq_ar_root)
}

# Load TUSCO reference once
tusco_ref <- load_tusco_reference(tusco_file)

results <- list()
for (sd in all_sample_dirs) {
  sample_name <- basename(sd)                       # e.g., B31.isoforms
  class_file  <- file.path(sd, paste0(sample_name, "_classification.txt"))
  if (!file.exists(class_file)) next

  classification <- read_tsv_safe(class_file)
  metrics <- compute_tusco_bigcats(classification, tusco_ref)

  results[[length(results) + 1]] <- tibble(
    Sample = sample_name,
    Tissue = infer_tissue(sample_name),
    TP = metrics$TP,
    PTP = metrics$PTP,
    FP = metrics$FP,
    FN = metrics$FN,
    id_type = metrics$id_type
  )
}

if (length(results) == 0) stop("No classification files could be read under ", isoseq_ar_root)
results_df <- bind_rows(results)

# ------------------------------------------------------------
# Convert to percentages per sample and summarize by tissue
# ------------------------------------------------------------
long_df <- results_df %>%
  pivot_longer(cols = c(TP, PTP, FP, FN), names_to = "Metric", values_to = "Count") %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count), Percentage = 100 * Count / ifelse(Total == 0, 1, Total)) %>%
  ungroup()

summary_df <- long_df %>%
  group_by(Tissue, Metric) %>%
  summarize(
    mean_perc = mean(Percentage, na.rm = TRUE),
    sd_perc = sd(Percentage, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

# Save underlying data
write_tsv(long_df, plot_tsv)

# ------------------------------------------------------------
# Plot (match width to fig-5b_with_merge: 7.09 in)
# ------------------------------------------------------------
summary_df$Metric <- factor(summary_df$Metric, levels = c("TP", "PTP", "FP", "FN"))

p <- ggplot(summary_df, aes(x = Metric, y = mean_perc, fill = Metric)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_perc - sd_perc, ymax = mean_perc + sd_perc),
                width = 0.2, linewidth = 0.3) +
  facet_grid(. ~ Tissue, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 100), labels = scales::percent_format(scale = 1), breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = c(TP = "#2ca02c", PTP = "#9467bd", FP = "#1f77b4", FN = "#d62728")) +
  labs(x = NULL, y = "Percentage (%)", title = "TUSCO (mouse multi-exon) TP/PTP/FP/FN on Iso-Seq AR") +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(size = 7, hjust = 0.5, face = "bold", margin = margin(b = 6)),
    axis.title.y = element_text(size = 7, margin = margin(r = 4)),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7, face = "bold", margin = margin(b = 4)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = "none"
  )

# Save with similar width as fig-5b_with_merge.R
pdf(plot_pdf, width = 7.09, height = 2.8)
print(p)
dev.off()

png(plot_png, width = 2127, height = 840, res = 300)
print(p)
dev.off()

message("Saved plot to ", plot_dir)
