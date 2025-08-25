#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# --------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------
repo_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"

csv_file <- file.path(repo_root, "figs/fig-s3/FN_correlation_plot.csv")
output_dir <- file.path(repo_root, "figs/fig-s3/plots")

human_parent_dir <- file.path(repo_root, "figs/data/lrgasp/human")
mouse_parent_dir <- file.path(repo_root, "figs/data/lrgasp/mouse")

tusco_human_file <- file.path(repo_root, "figs/data/tusco/tusco_human.tsv")
tusco_mouse_file <- file.path(repo_root, "figs/data/tusco/tusco_mouse.tsv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --------------------------------------------------------------------------------------
# Helpers copied/adapted from figs/fig-3/code/figure3d-human.R
# --------------------------------------------------------------------------------------
read_tusco_genes <- function(tsv_path) {
  tusco_df <- read_delim(
    tsv_path,
    delim = "\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = cols(.default = "c"),
    trim_ws = TRUE
  )
  if (ncol(tusco_df) == 1) {
    message("Falling back to space-delimited reading for TUSCO TSV: ", tsv_path)
    tusco_df <- read_table2(
      tsv_path,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = cols(.default = "c")
    )
  }
  if (!"ensembl" %in% colnames(tusco_df)) {
    stop("Tusco TSV must contain 'ensembl' column: ", tsv_path)
  }
  tusco_df %>%
    mutate(ensembl_id = str_remove(ensembl, "\\..+$")) %>%
    distinct(ensembl_id) %>%
    pull(ensembl_id) %>%
    { .[!. %in% c("#", "") & !is.na(.)] }
}

get_tusco_fsm_ism_found <- function(class_file, tusco_set) {
  if (!file.exists(class_file)) {
    stop("Classification file not found: ", class_file)
  }
  df <- read_tsv(class_file, col_types = cols(.default = "c"))
  colnames(df)[1:7] <- c("isoform","chrom","strand","length","exons","structural_category","associated_gene")
  df_clean <- df %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      df %>%
        filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(associated_gene = str_remove(associated_gene, "\\.\\d+$"))
  found_genes <- df_clean %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    distinct(associated_gene) %>%
    pull(associated_gene)
  intersect(tusco_set, found_genes)
}

map_sample_to_classification <- function(sample_name) {
  # Samples like: WTC11-cDNA-ONT(-LS), WTC11-cDNA-PacBio(-LS), WTC11-dRNA-ONT(-LS)
  #               ES-cDNA-ONT(-LS), ES-cDNA-PacBio(-LS), ES-dRNA-ONT(-LS)
  parts <- str_split(sample_name, "-", simplify = TRUE)
  if (ncol(parts) < 3) {
    stop("Unexpected Sample format: ", sample_name)
  }
  cell <- parts[1]
  method <- tolower(parts[2])  # cdna or drna
  platform <- tolower(parts[3]) # ont or pacbio
  is_ls <- grepl("-LS$", sample_name)

  # Build folder name
  folder_prefix <- paste0(cell, "_", method, "_")
  folder <- paste0(folder_prefix, platform)
  if (is_ls) folder <- paste0(folder, "_ls")

  parent_dir <- if (cell == "WTC11") human_parent_dir else if (cell == "ES") mouse_parent_dir else stop("Unknown cell prefix: ", cell)

  # Prefer standard *_classification.txt; if not present, fall back to *_sqanti_classification.txt
  standard_path <- file.path(parent_dir, folder, paste0(folder, "_classification.txt"))
  sqanti_path <- file.path(parent_dir, folder, paste0(folder, "_sqanti_classification.txt"))

  class_file <- if (file.exists(standard_path)) standard_path else if (file.exists(sqanti_path)) sqanti_path else NA_character_

  list(cell = cell, folder = folder, class_file = class_file, parent_dir = parent_dir)
}

# --------------------------------------------------------------------------------------
# Load CSV (read metrics), recompute FN per sample using classification + TUSCO set
# --------------------------------------------------------------------------------------
if (!file.exists(csv_file)) stop("CSV not found: ", csv_file)

metrics_df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8") %>%
  mutate(
    Sample = trimws(Sample),
    `Log(total Reads Number) * Median Read Length` = suppressWarnings(as.numeric(as.character(`Log(total Reads Number) * Median Read Length`))),
    `Total Reads Number` = suppressWarnings(as.numeric(as.character(`Total Reads Number`))),
    `Median Read Length` = suppressWarnings(as.numeric(as.character(`Median Read Length`)))
  )

# Pre-load TUSCO gene sets
if (!file.exists(tusco_human_file)) stop("Missing TUSCO human TSV: ", tusco_human_file)
if (!file.exists(tusco_mouse_file)) stop("Missing TUSCO mouse TSV: ", tusco_mouse_file)

tusco_genes_human <- read_tusco_genes(tusco_human_file)
tusco_genes_mouse <- read_tusco_genes(tusco_mouse_file)

# For each sample, compute FN count following figure3d-human.R approach
compute_fn_for_sample <- function(sample_name) {
  map_info <- map_sample_to_classification(sample_name)
  class_file <- map_info$class_file
  if (is.na(class_file)) {
    warning("Classification file not found for sample ", sample_name, "; returning NA")
    return(NA_integer_)
  }
  tusco_set <- if (map_info$cell == "WTC11") tusco_genes_human else tusco_genes_mouse
  found <- get_tusco_fsm_ism_found(class_file, tusco_set)
  length(setdiff(tusco_set, found))
}

metrics_df <- metrics_df %>%
  rowwise() %>%
  mutate(`False Negatives (recomputed)` = compute_fn_for_sample(Sample)) %>%
  ungroup()

# --------------------------------------------------------------------------------------
# Correlation and plotting (using recomputed FN)
# --------------------------------------------------------------------------------------
# Filter rows with computed FN available
plot_df <- metrics_df %>% filter(!is.na(`False Negatives (recomputed)`))

cor_test <- cor.test(
  plot_df$`Log(total Reads Number) * Median Read Length`,
  plot_df$`False Negatives (recomputed)`,
  method = "pearson"
)

r_value <- as.numeric(cor_test$estimate)
p_value <- cor_test$p.value

custom_theme <- theme_classic(base_size = 7) +
  theme(
    axis.text  = element_text(size = rel(1), color = "black"),
    axis.title = element_text(size = rel(1.1)),
    axis.line  = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.title = element_text(size=rel(1.0), face="bold"),
    legend.text = element_text(size=rel(0.9))
  )

# Derive grouping aesthetics similar to prior script
plot_df <- plot_df %>%
  mutate(
    cell_line = str_extract(Sample, "^(ES|WTC11)"),
    method_simple = case_when(
      grepl("PacBio", Sample, fixed = TRUE) ~ "PacBio",
      grepl("ONT", Sample, fixed = TRUE) ~ "ONT",
      TRUE ~ "Other"
    ),
    method = case_when(
      method_simple == "PacBio" & grepl("cDNA", Sample, fixed=TRUE) ~ "PacBio-cDNA",
      method_simple == "ONT" & grepl("cDNA", Sample, fixed=TRUE) ~ "ONT-cDNA",
      method_simple == "ONT" & grepl("dRNA", Sample, fixed=TRUE) ~ "ONT-dRNA",
      TRUE ~ method_simple
    ),
    is_ls = grepl("-LS$", Sample),
    point_size = if_else(is_ls, 4.5, 3)
  )

p <- ggplot(plot_df, aes(x = `Log(total Reads Number) * Median Read Length`, y = `False Negatives (recomputed)`)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgray", fill = "gray", linetype = "solid", alpha = 0.2, linewidth = 0.4) +
  {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                               force = 1.5, segment.color = 'grey50', segment.size = 0.2, segment.alpha = 0.6, min.segment.length = 0)
    } else {
      geom_text(aes(label = Sample), size = 2.5, vjust = -0.6, color = "black")
    }
  } +
  geom_point(aes(color = cell_line, shape = method, size = point_size), alpha = 0.9, stroke = 0.5) +
  scale_color_manual(name = "Cell Line", values = c("WTC11" = "#a6dba0", "ES" = "#1C9E77")) +
  scale_shape_manual(name = "Method", values = c("PacBio-cDNA" = 16, "ONT-cDNA" = 17, "ONT-dRNA" = 18)) +
  scale_size_identity(guide = "legend", name = "Short-read Support", breaks = c(3, 4.5), labels = c("No", "Yes")) +
  guides(color = guide_legend(title="Cell Line", override.aes = list(size=4)),
         shape = guide_legend(title="Method", override.aes = list(size=4)),
         size = guide_legend(title="Short-read Support")) +
  labs(x = "Log(Total Reads) * Median Read Length",
       y = "False Negatives (recomputed)",
       subtitle = paste0("Pearson r = ", round(r_value, 2), ", p = ", format(p_value, scientific = TRUE, digits = 2))) +
  custom_theme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(5, 5, 5, 5))

print(p)

outfile_pdf <- file.path(output_dir, "fig-s3.pdf")

ggsave(filename = outfile_pdf, plot = p, width = 5, height = 3.5, device = "pdf")

message("Plot saved to: ", outfile_pdf)
