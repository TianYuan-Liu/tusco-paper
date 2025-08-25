###############################################################
# Revised Code for Panel A - Mouse Version
###############################################################

# Load necessary libraries
message("Loading necessary libraries...")
suppressPackageStartupMessages({
  library(ggplot2)
  library(rmarkdown)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(rtracklayer)
  library(GenomicRanges)
  library(Gviz)
  library(fmsb)
  library(cowplot)
  library(grid)
  library(png)
})

# Define global variables
# Paths to fixed files (adjust as needed)
sirv.file <- '/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/spike-ins/lrgasp_sirvs.gtf'
tusco.file <- '/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/tusco/tusco_mouse.tsv'

data_dir <- '/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/lrgasp/mouse'

# Pipelines in the desired order:
# Row 1: ES dRNA ONT, ES cDNA ONT, ES cDNA PacBio
# Row 2: ES dRNA ONT LS, ES cDNA ONT LS, ES cDNA PacBio LS
pipelines <- c("es_drna_ont",
               "es_cdna_ont",
               "es_cdna_pacbio",
               "es_drna_ont_ls",
               "es_cdna_ont_ls",
               "es_cdna_pacbio_ls")

# Define helper functions

# Helper function to read TSV files with error handling
read_tsv_safe <- function(file_path, col_names = TRUE, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  tryCatch(
    {
      data <- read_tsv(file_path, col_names = col_names, ...)
      message("Successfully read file: ", file_path)
      return(data)
    },
    error = function(e) {
      stop("Error reading file ", file_path, ": ", e$message)
    }
  )
}

# Function to format pipeline names
format_pipeline_name <- function(prefix) {
  parts <- strsplit(prefix, "_")[[1]]
  name_map <- c(drna = "dRNA", cdna = "cDNA", ont = "ONT", pacbio = "PacBio")
  
  formatted_parts <- c("ES")
  for (p in parts[-1]) {
    if (p %in% names(name_map)) {
      formatted_parts <- c(formatted_parts, name_map[p])
    } else {
      if (p == "ls") {
        formatted_parts <- c(formatted_parts, "LS")
      } else {
        formatted_parts <- c(formatted_parts, toupper(p))
      }
    }
  }
  paste(formatted_parts, collapse = " ")
}

# Function to process each pipeline
process_pipeline <- function(pipeline_prefix) {
  # Define variables at the beginning of the function
  class.file <- file.path(data_dir, paste0(pipeline_prefix, "/", pipeline_prefix, "_classification.txt"))
  transcript_gtf_file <- file.path(data_dir, paste0(pipeline_prefix, "/", pipeline_prefix, "_corrected.gtf"))
  
  pipeline_name <- format_pipeline_name(pipeline_prefix)
  
  classification_data <- read_tsv_safe(class.file)
  
  ###############################################################
  # SIRVs Evaluation
  ###############################################################
  
  # Define SIRV-related variables
  sirv_gtf <- rtracklayer::import(sirv.file)
  sirv_gtf_df <- as.data.frame(sirv_gtf)
  
  annotation_data_sirv <- sirv_gtf_df %>%
    dplyr::filter(type == "exon") %>%
    dplyr::distinct(transcript_id) %>%
    dplyr::rename(ref_transcript_id = transcript_id)
  
  rSIRV <- nrow(annotation_data_sirv)
  
  transcript_gtf <- rtracklayer::import(transcript_gtf_file)
  transcript_gtf_df <- as.data.frame(transcript_gtf)
  
  classification_data_sirv <- classification_data %>%
    mutate(id_type = "transcript_id")
  
  classification_data_cleaned_sirv <- classification_data_sirv %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data_sirv %>%
        filter(structural_category == "fusion") %>%
        separate_rows(associated_transcript, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    distinct(isoform, associated_transcript, .keep_all = TRUE) %>%
    arrange(isoform)
  
  sirv_chromosomes <- unique(sirv_gtf_df$seqnames)
  classification_data_cleaned_sirv <- classification_data_cleaned_sirv %>%
    filter(chrom %in% sirv_chromosomes)
  
  SIRV_transcripts <- classification_data_cleaned_sirv %>%
    filter(grepl("SIRV", chrom))
  
  SIRV_RM <- SIRV_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" &
           ref_exons == 1 &
           abs(diff_to_TSS) < 50 &
           abs(diff_to_TTS) < 50)
    )
  
  TP_sirv <- SIRV_RM
  TP_SIRV <- unique(TP_sirv$associated_transcript)
  
  PTP_sirv <- SIRV_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") & 
             !associated_transcript %in% TP_sirv$associated_transcript)
  
  # Any SIRV reference transcript that is detected (in ANY structural category) 
  # should not be counted as FN. Collect all associated_transcript IDs that appear 
  # anywhere in the cleaned classification data set and subtract them from the 
  # reference list.

  detected_sirv_transcripts <- unique(classification_data_cleaned_sirv$associated_transcript)

  FN_sirv <- annotation_data_sirv %>%
    filter(!(ref_transcript_id %in% detected_sirv_transcripts))
  
  FP_sirv <- SIRV_transcripts %>%
    filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog", "genic", "fusion", "antisense"))
  
  fsm_ism_count_sirv <- SIRV_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    nrow()
  
  non_redundant_sensitivity_sirv <- length(TP_SIRV) / rSIRV
  non_redundant_precision_sirv <- if (nrow(SIRV_transcripts) > 0) nrow(TP_sirv) / nrow(SIRV_transcripts) else NA
  redundant_precision_sirv <- if (fsm_ism_count_sirv > 0) (nrow(TP_sirv) + nrow(PTP_sirv)) / nrow(SIRV_transcripts) else NA
  positive_detection_rate_sirv <- if (rSIRV > 0) length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) / rSIRV else NA
  false_discovery_rate_sirv <- if (nrow(SIRV_transcripts) > 0) (nrow(SIRV_transcripts) - nrow(SIRV_RM)) / nrow(SIRV_transcripts) else NA
  false_detection_rate_sirv <- if (nrow(SIRV_transcripts) > 0) nrow(FP_sirv) / nrow(SIRV_transcripts) else NA
  unique_tp_ptp_sirv <- length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript)))
  redundancy_sirv <- if (unique_tp_ptp_sirv > 0) fsm_ism_count_sirv / unique_tp_ptp_sirv else NA
  
  SIRV_only <- SIRV_transcripts %>%
    mutate(subcategory = case_when(
      subcategory == "reference_match" ~ "Reference match",
      # Treat qualifying mono-exon as Reference match as well
      subcategory == "mono-exon" & ref_exons == 1 & abs(diff_to_TSS) < 50 & abs(diff_to_TTS) < 50 ~ "Reference match",
      subcategory == "alternative_3end" ~ "Alternative 3'end",
      subcategory == "alternative_5end" ~ "Alternative 5'end",
      subcategory == "alternative_3end5end" ~ "Alternative 3'5'end",
      TRUE ~ subcategory
    )) %>%
    mutate(
      big_category = case_when(
        structural_category == "full-splice_match" & subcategory == "Reference match" ~ "TP",
        (structural_category == "full-splice_match" & subcategory %in% c("Alternative 3'end","Alternative 5'end","Alternative 3'5'end")) ~ "PTP",
        structural_category == "incomplete-splice_match" ~ "PTP",
        structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic_intron",
                                   "genic","antisense","fusion","intergenic") ~ "FP",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(final_label = case_when(
      big_category == "TP" ~ "RM",
      big_category == "PTP" & subcategory == "Alternative 3'end" ~ "Alternative 3'end",
      big_category == "PTP" & subcategory == "Alternative 5'end" ~ "Alternative 5'end",
      big_category == "PTP" & subcategory == "Alternative 3'5'end" ~ "Alternative 3'5'end",
      big_category == "PTP" & structural_category == "incomplete-splice_match" ~ "ISM",
      big_category == "FP" & structural_category == "novel_in_catalog" ~ "NIC",
      big_category == "FP" & structural_category == "novel_not_in_catalog" ~ "NNC",
      big_category == "FP" & structural_category == "genic_intron" ~ "Genic Intron",
      big_category == "FP" & structural_category == "genic" ~ "Genic Genomic",
      big_category == "FP" & structural_category == "antisense" ~ "Antisense",
      big_category == "FP" & structural_category == "fusion" ~ "Fusion",
      big_category == "FP" & structural_category == "intergenic" ~ "Intergenic",
      TRUE ~ NA_character_
    ))
  
  missing_df_sirv <- data.frame(
    final_label = "Missing",
    big_category = "FN",
    count = nrow(FN_sirv)
  )
  
  plot_data_sirv <- SIRV_only %>%
    filter(!is.na(final_label)) %>%
    group_by(big_category, final_label) %>%
    summarise(count = n(), .groups = "drop") %>%
    bind_rows(missing_df_sirv) %>%
    mutate(percentage = count / sum(count) * 100)
  
  plot_data_sirv$big_category <- factor(plot_data_sirv$big_category, levels = c("TP", "PTP", "FP", "FN"))
  plot_data_sirv$final_label <- factor(plot_data_sirv$final_label,
                                       levels = c("RM", "Alternative 3'end", "Alternative 5'end", "Alternative 3'5'end", "ISM",
                                                  "NIC", "NNC", "Genic Intron", "Genic Genomic", "Antisense", "Fusion", "Intergenic",
                                                  "Missing"))
  
  cat.palette <- c("FSM" = "#6BAED6", "ISM" = "#FC8D59", "NIC" = "#78C679",
                   "NNC" = "#EE6A50", "Genic Genomic" = "#969696", "Antisense" = "#66C2A4",
                   "Fusion" = "goldenrod1", "Intergenic" = "darksalmon", "Genic Intron" = "#41B6C4")
  
  subcat.palette <- c("Alternative 3'end" = '#02314d',
                      "Alternative 3'5'end" = '#0e5a87',
                      "Alternative 5'end" = '#7ccdfc',
                      "Reference match" = '#c4e1f2')
  
  mytheme <- theme_classic(base_family = "Helvetica") +
    theme(axis.line.x = element_line(color = "black", size = 0.4),
          axis.line.y = element_line(color = "black", size = 0.4),
          axis.title.x = element_text(size = 13),
          axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
          axis.title.y = element_text(size = 13),
          axis.text.y  = element_text(size = 12),
          plot.title = element_text(lineheight = .4, size = 15, hjust = 0.5, family = "Helvetica"))
  
  ###############################################################
  # TUSCO Evaluation
  ###############################################################

  # Define TUSCO-related variables
  # Robust read of header‑less Tusco annotation (falls back to whitespace if tabs not found)
  tusco_df <- read_delim(
    tusco.file,
    delim = "\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = cols(.default = "c"),
    trim_ws = TRUE
  )
  if (ncol(tusco_df) == 1) {
    tusco_df <- read_table2(
      tusco.file,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = cols(.default = "c")
    )
  }

  if (!all(c("ensembl","refseq","gene_name") %in% colnames(tusco_df))) {
    stop("Tusco TSV must contain columns: ensembl, refseq, gene_name")
  }

  annotation_data_tusco <- tusco_df %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()

  rTUSCO <- nrow(annotation_data_tusco)

  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Z0-9]+$"
  )

  classification_data_cleaned_tusco <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.\\d+$")
    ) %>%
    mutate(
      id_type = case_when(
        str_detect(associated_gene, patterns$ensembl)   ~ "ensembl",
        str_detect(associated_gene, patterns$refseq)    ~ "refseq",
        str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
        TRUE ~ "unknown"
      )
    ) %>%
    distinct(isoform, associated_gene, .keep_all = TRUE) %>%
    arrange(isoform)

  id_summary_tusco <- classification_data_cleaned_tusco %>%
    count(id_type, sort = TRUE)

  if (nrow(id_summary_tusco) > 0) {
    top_id_tusco      <- id_summary_tusco %>% slice_max(n, n = 1)
    top_id_type_tusco <- top_id_tusco$id_type
  } else {
    stop("No id_type classifications were made for TUSCO.")
  }

  if (top_id_type_tusco == "ensembl") {
    transcript_gtf_df$gene_id <- sub("\\..*", "", transcript_gtf_df$gene_id)
  }

  TUSCO_transcripts <- classification_data_cleaned_tusco %>%
    filter(associated_gene %in% annotation_data_tusco[[top_id_type_tusco]])

  TUSCO_RM <- TUSCO_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" &
           ref_exons == 1 &
           abs(diff_to_TSS) < 50 &
           abs(diff_to_TTS) < 50)
    )

  TP_tusco  <- TUSCO_RM
  TP_TUSCO  <- unique(TP_tusco$associated_gene)

  PTP_tusco <- TUSCO_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") &
           !associated_gene %in% TP_tusco$associated_gene)

  FN_tusco <- annotation_data_tusco %>%
    filter(!(!!sym(top_id_type_tusco) %in% TUSCO_transcripts$associated_gene))

  FP_tusco <- TUSCO_transcripts %>%
    filter(structural_category %in% c(
      "novel_in_catalog", "novel_not_in_catalog", "genic",
      "fusion", "antisense", "intergenic", "genic_intron"
    ))

  fsm_ism_count_tusco <- TUSCO_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    nrow()

  non_redundant_sensitivity_tusco <- length(TP_TUSCO) / rTUSCO
  non_redundant_precision_tusco  <- if (nrow(TUSCO_transcripts) > 0) nrow(TP_tusco) / nrow(TUSCO_transcripts) else NA
  redundant_precision_tusco      <- if (fsm_ism_count_tusco > 0) (nrow(TP_tusco) + nrow(PTP_tusco)) / nrow(TUSCO_transcripts) else NA
  positive_detection_rate_tusco  <- if (rTUSCO > 0) length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene))) / rTUSCO else NA
  false_discovery_rate_tusco     <- if (nrow(TUSCO_transcripts) > 0) (nrow(TUSCO_transcripts) - nrow(TUSCO_RM)) / nrow(TUSCO_transcripts) else NA
  false_detection_rate_tusco     <- if (nrow(TUSCO_transcripts) > 0) nrow(FP_tusco) / nrow(TUSCO_transcripts) else NA
  unique_tp_ptp_tusco            <- length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene)))
  redundancy_tusco               <- if (unique_tp_ptp_tusco > 0) fsm_ism_count_tusco / unique_tp_ptp_tusco else NA

  TUSCO_only <- TUSCO_transcripts %>%
    mutate(subcategory = case_when(
      subcategory == "reference_match"     ~ "Reference match",
      subcategory == "alternative_3end"    ~ "Alternative 3'end",
      subcategory == "alternative_5end"    ~ "Alternative 5'end",
      subcategory == "alternative_3end5end"~ "Alternative 3'5'end",
      TRUE ~ subcategory
    )) %>%
    mutate(
      big_category = case_when(
        structural_category == "full-splice_match" & subcategory == "Reference match"                               ~ "TP",
        (structural_category == "full-splice_match" & subcategory %in% c("Alternative 3'end","Alternative 5'end","Alternative 3'5'end")) ~ "PTP",
        structural_category == "incomplete-splice_match"                                                            ~ "PTP",
        structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic_intron",
                                   "genic","antisense","fusion","intergenic")                                       ~ "FP",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(final_label = case_when(
      big_category == "TP"  ~ "RM",
      big_category == "PTP" & subcategory == "Alternative 3'end"   ~ "Alternative 3'end",
      big_category == "PTP" & subcategory == "Alternative 5'end"   ~ "Alternative 5'end",
      big_category == "PTP" & subcategory == "Alternative 3'5'end" ~ "Alternative 3'5'end",
      big_category == "PTP" & structural_category == "incomplete-splice_match"                                      ~ "ISM",
      big_category == "FP"  & structural_category == "novel_in_catalog"                                              ~ "NIC",
      big_category == "FP"  & structural_category == "novel_not_in_catalog"                                          ~ "NNC",
      big_category == "FP"  & structural_category == "genic_intron"                                                  ~ "Genic Intron",
      big_category == "FP"  & structural_category == "genic"                                                         ~ "Genic Genomic",
      big_category == "FP"  & structural_category == "antisense"                                                     ~ "Antisense",
      big_category == "FP"  & structural_category == "fusion"                                                        ~ "Fusion",
      big_category == "FP"  & structural_category == "intergenic"                                                    ~ "Intergenic",
      TRUE ~ NA_character_
    ))

  missing_df_tusco <- data.frame(
    final_label = "Missing",
    big_category = "FN",
    count = nrow(FN_tusco)
  )

  plot_data_tusco <- TUSCO_only %>%
    filter(!is.na(final_label)) %>%
    group_by(big_category, final_label) %>%
    summarise(count = n(), .groups = "drop") %>%
    bind_rows(missing_df_tusco) %>%
    mutate(percentage = count / sum(count) * 100)

  plot_data_tusco$big_category <- factor(plot_data_tusco$big_category, levels = c("TP", "PTP", "FP", "FN"))
  plot_data_tusco$final_label  <- factor(plot_data_tusco$final_label,
                                         levels = c("RM", "Alternative 3'end", "Alternative 5'end", "Alternative 3'5'end", "ISM",
                                                    "NIC", "NNC", "Genic Intron", "Genic Genomic", "Antisense", "Fusion", "Intergenic",
                                                    "Missing"))

  metrics_labels <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")

  metrics_values_sirv <- c(
    non_redundant_sensitivity_sirv * 100,
    non_redundant_precision_sirv * 100,
    if (!is.na(redundancy_sirv) && redundancy_sirv != 0) (1 / redundancy_sirv) * 100 else 0,
    if (!is.na(false_discovery_rate_sirv)) (100 - (false_discovery_rate_sirv * 100)) else 0,
    if (!is.na(positive_detection_rate_sirv)) positive_detection_rate_sirv * 100 else 0,
    if (!is.na(redundant_precision_sirv)) redundant_precision_sirv * 100 else 0
  )

  metrics_values_tusco <- c(
    non_redundant_sensitivity_tusco * 100,
    non_redundant_precision_tusco  * 100,
    if (!is.na(redundancy_tusco) && redundancy_tusco != 0) (1 / redundancy_tusco) * 100 else 0,
    if (!is.na(false_discovery_rate_tusco)) (100 - (false_discovery_rate_tusco * 100)) else 0,
    if (!is.na(positive_detection_rate_tusco)) positive_detection_rate_tusco * 100 else 0,
    if (!is.na(redundant_precision_tusco)) redundant_precision_tusco * 100 else 0
  )

  # ---- Cosine similarity between SIRVs and TUSCO ----
  cosine_similarity <- NA
  if (!any(is.na(metrics_values_sirv)) && !any(is.na(metrics_values_tusco))) {
    denom <- sqrt(sum(metrics_values_sirv^2)) * sqrt(sum(metrics_values_tusco^2))
    if (denom != 0) {
      cosine_similarity <- sum(metrics_values_sirv * metrics_values_tusco) / denom
    }
  }
  similarity_label <- sprintf("Cosine similarity: %.3f", cosine_similarity)

  radar_df <- rbind(
    rep(100, length(metrics_labels)), # Max
    rep(0, length(metrics_labels)),   # Min
    metrics_values_sirv,
    metrics_values_tusco
  )
  colnames(radar_df) <- metrics_labels
  radar_df <- as.data.frame(radar_df)

  # Create radar grob without title and without any text or numbers
  radar_grob <- function(df) {
    tmpfile <- tempfile(fileext = ".png")
    png(tmpfile, width = 1200, height = 1200, res = 200)
    par(family = "Helvetica", mar = c(0, 0, 0, 0))
    radarchart(df,
               axistype = 0,
               cglcol = "grey", cglty = 1, cglwd = 1.5,
               plty = c(1, 1),
               axislabcol = "black",
               vlabels = rep("", length(metrics_labels)), # No variable labels
               pcol = c("#cab2d6", "#1b9e77"), plwd = 8, pty = 16,
               caxislabels = NULL,            # Remove axis labels
               vlcex = 0                        # No label size
    )
    dev.off()
    img <- png::readPNG(tmpfile)
    file.remove(tmpfile)
    grid::rasterGrob(img)
  }

  # Combine radar grob without title
  base_radar <- radar_grob(radar_df)
  # No title added

  # Return the radar plot grob
  radar_plot <- ggdraw() +
    draw_grob(base_radar, 0, 0, 1, 1) +
    draw_label(similarity_label,
               x = 0.5, y = -0.05, hjust = 0.5, vjust = 1,
               fontfamily = "Helvetica", size = 12)

  ###############################################################
  # Combined Barplot (unchanged)
  ###############################################################
  
  # Define combined barplot variables
  sirv_temp  <- plot_data_sirv  %>% dplyr::select(big_category, final_label, percentage)
  tusco_temp <- plot_data_tusco %>% dplyr::select(big_category, final_label, percentage)

  joined <- dplyr::full_join(
    sirv_temp %>% dplyr::rename(perc_sirv = percentage),
    tusco_temp %>% dplyr::rename(perc_tusco = percentage),
    by = c("big_category", "final_label")
  )

  combined_df <- joined %>%
    tidyr::replace_na(list(perc_sirv = 0, perc_tusco = 0)) %>%
    tidyr::pivot_longer(cols = c("perc_sirv", "perc_tusco"), names_to = "Type", values_to = "Percentage") %>%
    dplyr::mutate(Type = ifelse(Type == "perc_sirv", "SIRVs", "TUSCO"))

  combined_colors <- c("SIRVs" = "#cab2d6", "TUSCO" = "#1b9e77")
  
  # Keep legend for bar plot
  p_combined <- ggplot(combined_df, aes(x = final_label, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", size = 0.3, width = 0.7) +
    facet_grid(~big_category, scales = "free_x", space = "free") +
    scale_fill_manual(values = combined_colors) +
    xlab(NULL) +
    ylab("Percentage") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    mytheme +
    theme(strip.background = element_rect(fill = "grey95", colour = NA),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(pipeline_name)
  
  ###############################################################
  # Additional barplot following the style in figure7.R
  ###############################################################

  # Re-use previously defined palettes (cat.palette & subcat.palette)
  color_values <- c()
  for (label in levels(plot_data_tusco$final_label)) {
    if      (label == "RM")               color_values <- c(color_values, subcat.palette["Reference match"]) else
    if      (label == "ISM")              color_values <- c(color_values, cat.palette["ISM"]) else
    if      (label == "NIC")              color_values <- c(color_values, cat.palette["NIC"]) else
    if      (label == "NNC")              color_values <- c(color_values, cat.palette["NNC"]) else
    if      (label == "Genic Intron")     color_values <- c(color_values, cat.palette["Genic Intron"]) else
    if      (label == "Genic Genomic")    color_values <- c(color_values, cat.palette["Genic Genomic"]) else
    if      (label == "Antisense")        color_values <- c(color_values, cat.palette["Antisense"]) else
    if      (label == "Fusion")           color_values <- c(color_values, cat.palette["Fusion"]) else
    if      (label == "Intergenic")       color_values <- c(color_values, cat.palette["Intergenic"]) else
    if      (label == "Alternative 3'end")   color_values <- c(color_values, subcat.palette["Alternative 3'end"]) else
    if      (label == "Alternative 5'end")   color_values <- c(color_values, subcat.palette["Alternative 5'end"]) else
    if      (label == "Alternative 3'5'end") color_values <- c(color_values, subcat.palette["Alternative 3'5'end"]) else
    if      (label == "Missing")          color_values <- c(color_values, "gray50")
  }
  names(color_values) <- levels(plot_data_tusco$final_label)

  p_barplot <- ggplot(plot_data_tusco, aes(x = final_label, y = percentage, fill = final_label)) +
    geom_bar(stat = "identity", color = "black", size = 0.3, width = 0.7) +
    facet_grid(~big_category, scales = "free_x", space = "free") +
    scale_fill_manual(values = color_values) +
    xlab(NULL) +
    ylab("Percentage") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    mytheme +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey95", colour = NA),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(pipeline_name)

  list(
    radar    = radar_plot,
    combined = p_combined,
    barplot  = p_barplot
  )
}

# Define main processing variables

# Process all pipelines
results <- lapply(pipelines, process_pipeline)

# Extract radar plots
radar_plots <- lapply(results, function(x) x$radar)

# Radar legend with points
# Define radar legend variables
radar_legend_df <- data.frame(
  Type = c("SIRVs", "TUSCO"),
  x = c(1, 2),
  y = c(1, 2)
)

# Create a dummy radar legend plot
p_radar_legend_dummy <- ggplot(radar_legend_df, aes(x = x, y = y, color = Type, linetype = Type, shape = Type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("SIRVs" = "#cab2d6", "TUSCO" = "#1b9e77")) +
  scale_linetype_manual(values = c("SIRVs" = "solid", "TUSCO" = "solid")) +
  scale_shape_manual(values = c("SIRVs" = 16, "TUSCO" = 16)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, family = "Helvetica"),          # Increased text size
    legend.key.size = unit(1.5, "lines"),                                # Increased key size
    legend.spacing.x = unit(0.5, "cm"),                                 # Increased spacing between items
    legend.direction = "horizontal"                                     # Ensure horizontal layout
  )

# Extract the legend using cowplot::get_plot_component()
legend_component <- cowplot::get_legend(p_radar_legend_dummy)

# Convert the extracted legend component to a grob
radar_legend <- cowplot::ggdraw() +
  draw_grob(legend_component, 0, 0, 1, 1)

# Arrange the radar plots in a single row (6 columns)
panel_a <- plot_grid(plotlist = radar_plots, ncol = 6, align = "hv")

# Add the legend at the bottom of panel a
panel_a_with_legend <- plot_grid(panel_a, radar_legend, ncol = 1, rel_heights = c(1, 0.1))

 

# ------------------------------------------------------------------
# Export panels
# ------------------------------------------------------------------

# Save Panel A (radar only)
pdf(file = "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/fig-3/plots/figure3a-mouse.pdf", width = 24, height = 4) # Adjusted width and height for single row
print(panel_a_with_legend)
dev.off()