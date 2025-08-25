###############################################################
# Script for TUSCO Evaluation (Figure 4)
# Generates a single overall 8x5 radar grid PDF: fig-4b.pdf
###############################################################

message("Loading necessary libraries...")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(fmsb)
  library(cowplot)
  library(grid)
  library(png)
})

# ------------------------------------------------------------------
# Paths (updated to NIH data and repo plots directory)
# ------------------------------------------------------------------
# TUSCO references (leave as-is; required to compute metrics)
tusco_human_file <- "/Users/tianyuan/Desktop/GitHub/TUSCO/reference_dataset/tusco_human_multi_exon.tsv"
tusco_mouse_file <- "/Users/tianyuan/Desktop/GitHub/TUSCO/reference_dataset/tusco_mouse_multi_exon.tsv"

# Base data directory (lrgasp data)
base_data_dir    <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/lrgasp/tusco_novel_evl"

# Output directory (repo plots folder)
output_dir       <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/fig-4/plots"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------
# Helpers retained (minimal set required to render the overall grid)
# ------------------------------------------------------------------

# Read TSV safely
read_tsv_safe <- function(file_path, col_names = TRUE, ...) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  tryCatch(
    {
      data <- read_tsv(file_path, col_names = col_names, show_col_types = FALSE, ...)
      return(data)
    },
    error = function(e) {
      warning("Error reading file ", file_path, ": ", e$message)
      return(NULL)
    }
  )
}

# Compute TUSCO metrics (Sn, nrPre, 1/red, 1-FDR, PDR, rPre) in percent
calculate_tusco_metrics <- function(classification_data, tusco_annotation_df, rTUSCO_total_ref_transcripts) {
  if (is.null(classification_data) || nrow(classification_data) == 0) {
    return(rep(0, 6))
  }

  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)[0-9]{11}(\\.[0-9]+)?$",
    refseq    = "^(NM_|NR_|NP_)[0-9]{6,}(\\.[0-9]+)?$",
    gene_name = "^[A-Za-z0-9][A-Za-z0-9.-]*[A-Za-z0-9]$"
  )

  classification_data_cleaned <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.[0-9]+$"),
      id_type = case_when(
        str_detect(associated_gene, patterns$ensembl) ~ "ensembl",
        str_detect(associated_gene, patterns$refseq)  ~ "refseq",
        associated_gene %in% tusco_annotation_df$gene_name ~ "gene_name",
        TRUE ~ "unknown"
      )
    ) %>%
    filter(id_type != "unknown") %>%
    distinct(isoform, associated_gene, .keep_all = TRUE)

  if (nrow(classification_data_cleaned) == 0) {
    return(rep(0, 6))
  }

  id_summary <- classification_data_cleaned %>% count(id_type, sort = TRUE)
  # Choose the majority (most frequent) non-unknown id_type; fallback to gene_name
  top_id_type <- id_summary %>%
    dplyr::filter(id_type != "unknown") %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::pull(id_type)
  if (length(top_id_type) == 0 || is.na(top_id_type)) top_id_type <- "gene_name"

  tusco_transcripts <- classification_data_cleaned %>%
    filter(associated_gene %in% tusco_annotation_df[[top_id_type]])

  TUSCO_RM <- tusco_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" &
           ref_exons == 1 &
           abs(diff_to_TSS) < 50 &
           abs(diff_to_TTS) < 50)
    )
  TP_tusco <- TUSCO_RM
  TP_TUSCO_unique_genes <- unique(TP_tusco$associated_gene)

  PTP_tusco <- tusco_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") &
           !associated_gene %in% TP_TUSCO_unique_genes)

  detected_tusco_genes <- unique(tusco_transcripts$associated_gene)
  FN_tusco <- tusco_annotation_df %>%
    filter(!(!!sym(top_id_type) %in% detected_tusco_genes))

  FP_tusco <- tusco_transcripts %>%
    filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog", "genic", "fusion", "antisense", "intergenic", "genic_intron"))

  fsm_ism_count_tusco <- tusco_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    nrow()

  non_redundant_sensitivity_tusco <- if (rTUSCO_total_ref_transcripts > 0) length(TP_TUSCO_unique_genes) / rTUSCO_total_ref_transcripts else 0
  # Make nrPre transcript-based: TP transcripts / total detected transcripts
  non_redundant_precision_tusco  <- if (nrow(tusco_transcripts) > 0) nrow(TP_tusco) / nrow(tusco_transcripts) else 0

  unique_tp_ptp_tusco_genes <- length(unique(c(TP_TUSCO_unique_genes, PTP_tusco$associated_gene)))
  redundancy_tusco <- if (unique_tp_ptp_tusco_genes > 0) fsm_ism_count_tusco / unique_tp_ptp_tusco_genes else 0
  inv_redundancy_metric <- if (redundancy_tusco > 0) (1 / redundancy_tusco) else 0

  false_discovery_rate_tusco <- if (nrow(tusco_transcripts) > 0) (nrow(tusco_transcripts) - nrow(TUSCO_RM)) / nrow(tusco_transcripts) else 0
  one_minus_fdr_metric <- 1 - false_discovery_rate_tusco

  positive_detection_rate_tusco <- if (rTUSCO_total_ref_transcripts > 0) unique_tp_ptp_tusco_genes / rTUSCO_total_ref_transcripts else 0
  redundant_precision_tusco <- if (nrow(tusco_transcripts) > 0) (nrow(TP_tusco) + nrow(PTP_tusco)) / nrow(tusco_transcripts) else 0

  metrics_values <- c(
    Sn    = non_redundant_sensitivity_tusco * 100,
    nrPre = non_redundant_precision_tusco * 100,
    `1/red` = inv_redundancy_metric * 100,
    `1-FDR` = one_minus_fdr_metric * 100,
    PDR   = positive_detection_rate_tusco * 100,
    rPre  = redundant_precision_tusco * 100
  )
  metrics_values[is.na(metrics_values)] <- 0
  return(metrics_values)
}

# Render a single radar chart grob for a given pipeline/sample/species
process_and_plot_pipeline_pair <- function(pipeline_prefix, species, tusco_ref_filepath, data_dir_base, is_ml_processing = FALSE, omit_title = TRUE) {
  main_pipeline_type <- basename(data_dir_base)

  # Choose SQANTI directory and filename patterns per pipeline
  if (main_pipeline_type == "ujc_sq3") {
    sqanti_dir_name <- "sqanti3"
    ref_class_filename <- paste0(pipeline_prefix, "_sqanti3_ref_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_sqanti3_ref_classification.txt")
  } else if (main_pipeline_type == "stringtie_sq3") {
    sqanti_dir_name <- "sqanti3_out"
    ref_class_filename <- paste0(pipeline_prefix, "_stringtie_ref_sqanti_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_stringtie_tusco_sqanti_classification.txt")
  } else if (main_pipeline_type == "all-ref_sq3") {
    sqanti_dir_name <- "sqanti3"
    ref_class_filename <- paste0(pipeline_prefix, "_sqanti3_all_ref_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_sqanti3_all_ref_classification.txt")
  } else if (main_pipeline_type == "isoseq_sq3") {
    sqanti_dir_name <- "sqanti3_out"
    ref_class_filename <- "isoforms_classification.txt"
    novel_class_filename <- "isoforms_classification.txt"
  } else {
    sqanti_dir_name <- "sqanti3_out"
    ref_class_filename <- paste0(pipeline_prefix, "_sqanti_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_sqanti_classification.txt")
  }

  ref_evl_path <- file.path(data_dir_base, "ref_evl", pipeline_prefix, sqanti_dir_name)
  novel_evl_path <- file.path(data_dir_base, "novel_evl", pipeline_prefix, sqanti_dir_name)

  ref_class_file <- if (dir.exists(ref_evl_path)) file.path(ref_evl_path, ref_class_filename) else NULL
  novel_class_file <- if (dir.exists(novel_evl_path)) file.path(novel_evl_path, novel_class_filename) else NULL

  tusco_cols <- c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq")
  tusco_df <- NULL
  if (file.exists(tusco_ref_filepath)) {
    tusco_df <- tryCatch({
      read_delim(tusco_ref_filepath, delim = "\t", col_names = tusco_cols, col_types = cols(.default = "c"), trim_ws = TRUE, comment = "#", show_col_types = FALSE, guess_max = 2000)
    }, error = function(e) NULL)
    if (!is.null(tusco_df) && ncol(tusco_df) == 1) {
      tusco_df <- tryCatch({
        read_table2(tusco_ref_filepath, col_names = tusco_cols, col_types = cols(.default = "c"), comment = "#", show_col_types = FALSE, guess_max = 2000)
      }, error = function(e) NULL)
    }
  }
  if (is.null(tusco_df) || nrow(tusco_df) == 0) {
    return(ggdraw() + draw_label("No TUSCO ref", fontface = "bold"))
  }

  tusco_annotation_for_calc <- tusco_df %>% select(any_of(c("ensembl", "refseq", "gene_name"))) %>% distinct()
  rTUSCO <- nrow(tusco_annotation_for_calc)
  if (rTUSCO == 0) {
    return(ggdraw() + draw_label("Empty TUSCO ref", fontface = "bold"))
  }

  ref_classification_data <- if (!is.null(ref_class_file) && file.exists(ref_class_file)) read_tsv_safe(ref_class_file) else NULL
  novel_classification_data <- if (!is.null(novel_class_file) && file.exists(novel_class_file)) read_tsv_safe(novel_class_file) else NULL

  metrics_ref <- if (!is.null(ref_classification_data)) calculate_tusco_metrics(ref_classification_data, tusco_annotation_for_calc, rTUSCO) else NULL
  metrics_novel <- if (!is.null(novel_classification_data)) calculate_tusco_metrics(novel_classification_data, tusco_annotation_for_calc, rTUSCO) else NULL

  if (is.null(metrics_ref) && is.null(metrics_novel)) {
    return(ggdraw() + draw_label("No Valid Data", fontface = "bold"))
  }

  if (!is.null(metrics_ref) && !is.null(metrics_novel)) {
    radar_data_df <- rbind(rep(100, 6), rep(0, 6), metrics_ref, metrics_novel)
    plot_colors <- rep(if (species == "human") "#a8d5a0" else "#1b9e77", 2)
    plot_linetypes <- c(1, 2)
  } else if (!is.null(metrics_ref)) {
    radar_data_df <- rbind(rep(100, 6), rep(0, 6), metrics_ref)
    plot_colors <- if (species == "human") "#a8d5a0" else "#1b9e77"
    plot_linetypes <- 1
  } else {
    radar_data_df <- rbind(rep(100, 6), rep(0, 6), metrics_novel)
    plot_colors <- if (species == "human") "#a8d5a0" else "#1b9e77"
    plot_linetypes <- 2
  }

  colnames(radar_data_df) <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
  radar_data_df <- as.data.frame(radar_data_df)

  # Capture radarchart as a grob
  plot_grob <- NULL
  tryCatch({
    current_dev <- grDevices::dev.cur()
    temp_plot_file <- tempfile(fileext = ".png")
    png(temp_plot_file, width = 5, height = 5, units = "in", res = 200, bg = "white")
    par(mar = c(0.5, 0.5, if (omit_title) 0.5 else 1.5, 0.5), family = "sans", font = 2, cex = 7/12)
    radarchart(radar_data_df,
               axistype = 0,
               pcol = plot_colors,
               pfcol = NA,
               plwd = 8,
               plty = plot_linetypes,
               pty = 16,
               cglcol = "grey", cglty = 1, cglwd = 1.5,
               axislabcol = "black",
               vlabels = rep("", ncol(radar_data_df)),
               vlcex = 0,
               caxislabels = NULL,
               centerzero = FALSE,
               title = if (omit_title) "" else pipeline_prefix)
    dev.off()
    img_grob <- grid::rasterGrob(png::readPNG(temp_plot_file), interpolate = TRUE)
    file.remove(temp_plot_file)
    if (current_dev > 1) grDevices::dev.set(current_dev)
    plot_grob <- ggdraw() + draw_grob(img_grob)
  }, error = function(e) {
    plot_grob <<- ggdraw() + draw_label("Plot Error", fontface = "bold")
  })
  return(plot_grob)
}

# ------------------------------------------------------------------
# Build legend used underneath the grid
# ------------------------------------------------------------------
legend_data <- data.frame(
  Type = factor(c("Ref Evl", "Novel Evl"), levels = c("Ref Evl", "Novel Evl")),
  Linetype = factor(c("solid", "dashed"), levels = c("solid", "dashed")),
  Color = factor(c("human_color", "human_color"))
)
legend_plot_color <- "#a8d5a0"

p_legend_dummy <- ggplot(legend_data, aes(x = 1, y = Type, linetype = Type, color = Color)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = Type), size = 4, stroke = 1.5) +
  scale_linetype_manual(values = c("Ref Evl" = "solid", "Novel Evl" = "dashed")) +
  scale_color_manual(values = c("human_color" = legend_plot_color)) +
  scale_shape_manual(values = c("Ref Evl" = 16, "Novel Evl" = 16)) +
  guides(linetype = guide_legend(title = NULL, override.aes = list(linewidth = 1.5)),
         color = guide_legend(title = NULL, override.aes = list(color = legend_plot_color)),
         shape = guide_legend(title = NULL, override.aes = list(shape = 16))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "sans"),
        legend.key.width = unit(1.5, "cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.spacing.x = unit(0.5, "cm"))
radar_plot_legend <- cowplot::get_legend(p_legend_dummy)

# ------------------------------------------------------------------
# Generate the overall 8x5 radar grid and save as figure4.pdf
# ------------------------------------------------------------------
message("Generating overall 8x5 radar grid figure...")

# Ordered sample prefixes (columns)
overall_sample_prefixes <- c(
  "wtc11_captrap_pacbio",
  "wtc11_cdna_pacbio",
  "wtc11_cdna_ont",
  "wtc11_captrap_ont",
  "es_captrap_pacbio",
  "es_cdna_pacbio",
  "es_cdna_ont",
  "es_captrap_ont"
)

# Pipeline rows
pipeline_row_specs <- list(
  list(row_label = "Bambu",     pipeline_dir = "bambu_sq3",     is_ml = FALSE),
  list(row_label = "Stringtie", pipeline_dir = "stringtie_sq3",  is_ml = FALSE),
  list(row_label = "Flair",     pipeline_dir = "flair_sq3",     is_ml = FALSE),
  list(row_label = "All-ref",   pipeline_dir = "all-ref_sq3",   is_ml = FALSE),
  list(row_label = "IsoSeq",    pipeline_dir = "isoseq_sq3",    is_ml = FALSE)
)

# Collect grobs row-wise
overall_radar_grobs <- list()
for (row_spec in pipeline_row_specs) {
  pipeline_path <- file.path(base_data_dir, row_spec$pipeline_dir)
  for (sample_prefix in overall_sample_prefixes) {
    species <- if (grepl("^wtc11", sample_prefix, ignore.case = TRUE)) "human" else "mouse"
    tusco_file <- if (species == "human") tusco_human_file else tusco_mouse_file
    grob_cell <- tryCatch({
      process_and_plot_pipeline_pair(
        pipeline_prefix    = sample_prefix,
        species            = species,
        tusco_ref_filepath = tusco_file,
        data_dir_base      = pipeline_path,
        is_ml_processing   = row_spec$is_ml,
        omit_title         = TRUE
      )
    }, error = function(e) ggdraw() + draw_label("Error", fontface = "bold"))
    overall_radar_grobs <- append(overall_radar_grobs, list(grob_cell))
  }
}

# Column labels
column_labels <- c("WTC11 CapTrap PacBio", "WTC11 cDNA PacBio", "WTC11 cDNA ONT", "WTC11 CapTrap ONT",
                   "ES CapTrap PacBio", "ES cDNA PacBio", "ES cDNA ONT", "ES CapTrap ONT")
column_label_grobs <- lapply(column_labels, function(txt) {
  ggdraw() + draw_label(txt, fontface = "bold", size = 7, fontfamily = "sans")
})

# Build rows with row labels
row_labels <- c("Bambu", "Stringtie", "Flair", "All-ref", "IsoSeq")
row_panels <- list()
for (i in seq_along(row_labels)) {
  row_label_grob <- ggdraw() + draw_label(row_labels[i], fontfamily = "sans", fontface = "bold", size = 7, angle = 0, hjust = 1)
  row_grobs <- overall_radar_grobs[((i-1)*8 + 1):(i*8)]
  row_panel <- plot_grid(plotlist = c(list(row_label_grob), row_grobs), ncol = 9, rel_widths = c(0.5, rep(1,8)))
  row_panels[[i]] <- row_panel
}

# Header row
blank_corner <- ggdraw()
header_row <- plot_grid(plotlist = c(list(blank_corner), column_label_grobs), ncol = 9, rel_widths = c(0.5, rep(1,8)))

# Combine header + rows, add legend
overall_grid <- plot_grid(plotlist = c(list(header_row), row_panels), ncol = 1, rel_heights = c(0.6, rep(1,5)))
overall_grid_with_legend <- plot_grid(overall_grid, radar_plot_legend, ncol = 1, rel_heights = c(1, 0.05))

# Save as fig-4b.pdf (width 180 mm; height derived from 5 rows + header + legend)
cell_width_mm <- 180 / 8
overall_height_mm <- cell_width_mm * 6 + 10
  figure_path <- file.path(output_dir, "fig-4b.pdf")
message("Saving figure to ", figure_path)
ggsave(figure_path, overall_grid_with_legend, width = 180, height = overall_height_mm, units = "mm", device = "pdf", limitsize = FALSE)

message("Done.")
