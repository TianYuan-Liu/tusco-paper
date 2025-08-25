###############################################################################
# Load libraries
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(rtracklayer)
  library(ggplot2)
  library(ggsignif)   # for significance annotations
})

###############################################################################
# 1) Define file paths and pipelines
###############################################################################
sirv_gtf_file <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/spike-ins/lrgasp_sirvs.gtf"
tusco_tsv_file <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/tusco/tusco_human.tsv" # Changed from bugsi_gtf_file

pipelines <- c(
  "WTC11_drna_ont",
  "WTC11_cdna_ont",
  "WTC11_cdna_pacbio",
  "WTC11_drna_ont_ls",
  "WTC11_cdna_ont_ls",
  "WTC11_cdna_pacbio_ls"
)

data_dir <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/lrgasp/human" # This might need to be updated if TUSCO data is in a different dir

###############################################################################
# 2) Read TSV safely
###############################################################################
read_tsv_safe <- function(file_path, col_names=TRUE, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  tryCatch({
    data <- read_tsv(file_path, col_names=col_names, ...)
    return(data)
  },
  error=function(e) {
    stop("Error reading file ", file_path, ": ", e$message)
  }
  )
}

###############################################################################
# 3) Compute SIRVs: (TP, PTP, FP, FN)
###############################################################################
compute_sirv_bigcats <- function(classification_data, transcript_gtf_file, sirv_gtf_file) {

  # Read SIRV GTF
  sirv_gtf <- rtracklayer::import(sirv_gtf_file)
  sirv_gtf_df <- as.data.frame(sirv_gtf)

  # If "transcript_id" is missing, try copying from "Parent"
  if (!"transcript_id" %in% names(sirv_gtf_df)) {
    if ("Parent" %in% names(sirv_gtf_df)) {
      sirv_gtf_df$transcript_id <- sirv_gtf_df$Parent
    } else {
      stop("No 'transcript_id' or 'Parent' attribute in SIRV GTF. Cannot continue.")
    }
  }

  # Distinct SIRV transcripts (from exons)
  annotation_data_sirv <- sirv_gtf_df %>%
    dplyr::filter(type=="exon") %>%
    dplyr::distinct(transcript_id) %>%
    dplyr::rename(ref_transcript_id=transcript_id)

  # Expand classification if "fusion"
  class_sirv <- classification_data %>%
    filter(structural_category!="fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category=="fusion") %>%
        separate_rows(associated_transcript, sep="_")
    ) %>%
    mutate(
      associated_gene=str_remove(associated_gene,"\\.\\d+$"),
      associated_transcript=str_remove(associated_transcript,"\\.\\d+$")
    ) %>%
    distinct(isoform, associated_transcript, .keep_all=TRUE)

  # Keep only SIRV chroms
  sirv_chroms <- unique(sirv_gtf_df$seqnames)
  class_sirv <- class_sirv %>% filter(chrom %in% sirv_chroms)

  SIRV_transcripts <- class_sirv %>% filter(grepl("SIRV", chrom))

  # Define categories
  TP_sirv <- SIRV_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" &
           ref_exons == 1 &
           abs(diff_to_TSS) < 50 &
           abs(diff_to_TTS) < 50)
    )
  PTP_sirv <- SIRV_transcripts %>%
    filter(structural_category %in% c("full-splice_match","incomplete-splice_match"),
           !(associated_transcript %in% TP_sirv$associated_transcript))

  ref_found <- union(TP_sirv$associated_transcript, PTP_sirv$associated_transcript)
  FN_sirv <- annotation_data_sirv %>%
    filter(!(ref_transcript_id %in% ref_found))

  FP_sirv <- SIRV_transcripts %>%
    filter(structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic",
                                      "fusion","antisense","intergenic","genic_intron"))

  summary_sirv <- data.frame(
    big_category=c("TP","PTP","FP","FN"),
    count=c(nrow(TP_sirv), nrow(PTP_sirv), nrow(FP_sirv), nrow(FN_sirv))
  )
  total_sirv <- sum(summary_sirv$count)
  summary_sirv$percentage <- 100 * summary_sirv$count / if (total_sirv==0) 1 else total_sirv

  return(summary_sirv)
}

###############################################################################
# 4) Compute TUSCO: (TP, PTP, FP, FN)
###############################################################################
compute_tusco_bigcats <- function(classification_data, tusco_tsv_file) { # Changed from compute_bugsi_bigcats, removed transcript_gtf_file

  # Read TUSCO TSV
  # Robust read of header‑less Tusco annotation (falls back to whitespace if tabs not found)
  tusco_df <- read_delim(
    tusco_tsv_file,
    delim = "\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = cols(.default = "c"),
    trim_ws = TRUE
  )
  if (ncol(tusco_df) == 1) { # Fallback if not tab-delimited
    tusco_df <- read_table2(
      tusco_tsv_file,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = cols(.default = "c")
    )
  }

  if (!all(c("ensembl","refseq","gene_name") %in% colnames(tusco_df))) {
    stop("Tusco TSV must contain columns: ensembl, refseq, gene_name")
  }

  # Distinct annotation lines for genes
  annotation_data_tusco <- tusco_df %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()

  # Patterns for ID types
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Z0-9]+$"
  )

  class_tusco <- classification_data %>%
    filter(structural_category!="fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category=="fusion") %>%
        separate_rows(associated_gene, sep="_") # Use associated_gene for TUSCO
    ) %>%
    mutate(
      associated_gene=str_remove(associated_gene,"\\.\\d+$")
    ) %>%
    mutate(
      id_type=case_when(
        str_detect(associated_gene, patterns$ensembl) ~ "ensembl",
        str_detect(associated_gene, patterns$refseq)  ~ "refseq",
        str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
        TRUE ~ "unknown"
      )
    ) %>%
    distinct(isoform, associated_gene, .keep_all=TRUE)

  # Determine top ID type used in classification data
  id_summary_tusco <- class_tusco %>% count(id_type, sort=TRUE)
  if (nrow(id_summary_tusco)==0 || id_summary_tusco$id_type[1] == "unknown") {
     # Fallback or error if no dominant known ID type
    warning("No dominant known gene ID type (ensembl, refseq, gene_name) found in classifications for TUSCO. Results might be affected.")
    # Attempt to proceed or handle as error, here we'll try to use 'gene_name' as a default if available, or skip TUSCO processing for this sample
    if("gene_name" %in% id_summary_tusco$id_type) {
        top_type <- "gene_name"
    } else if ("ensembl" %in% id_summary_tusco$id_type) {
        top_type <- "ensembl"
    } else if ("refseq" %in% id_summary_tusco$id_type) {
        top_type <- "refseq"
    } else {
        # If no valid IDs, return empty summary
        return(data.frame(big_category=c("TP","PTP","FP","FN"), count=c(0,0,0,0), percentage=c(0,0,0,0)))
    }
  } else {
    top_type <- id_summary_tusco$id_type[1]
  }
  
  cat("Using TUSCO ID type:", top_type, "\n")


  # Keep only transcripts whose associated gene is in TUSCO annotation, using the determined top_type
  TUSCO_transcripts <- class_tusco %>%
    filter(associated_gene %in% annotation_data_tusco[[top_type]])

  # Define categories
  TP_tusco <- TUSCO_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" &
           ref_exons == 1 &
           abs(diff_to_TSS) < 50 &
           abs(diff_to_TTS) < 50)
    )
  PTP_tusco <- TUSCO_transcripts %>%
    filter(structural_category %in% c("full-splice_match","incomplete-splice_match"),
           !(associated_gene %in% TP_tusco$associated_gene)) # Match by gene

  # For FN, consider genes from TUSCO annotation not found in TP or PTP
  ref_ids_tusco <- annotation_data_tusco[[top_type]]
  found_ids_tusco <- union(TP_tusco$associated_gene, PTP_tusco$associated_gene)
  
  FN_tusco_df <- annotation_data_tusco %>%
    filter(!(!!sym(top_type) %in% found_ids_tusco))
  FN_tusco_count <- nrow(FN_tusco_df)

  FP_tusco <- TUSCO_transcripts %>% # False positives are among those matched to TUSCO genes but not TP/PTP
    filter(structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic",
                                      "fusion","antisense","intergenic","genic_intron"))

  summary_tusco <- data.frame(
    big_category=c("TP","PTP","FP","FN"),
    count=c(nrow(TP_tusco), nrow(PTP_tusco), nrow(FP_tusco), FN_tusco_count) # Use FN_tusco_count
  )
  total_tusco <- sum(summary_tusco$count)
  summary_tusco$percentage <- 100 * summary_tusco$count / if (total_tusco==0) 1 else total_tusco

  return(summary_tusco)
}

###############################################################################
# 5) Process a single pipeline => returns pipeline-level results
###############################################################################
process_pipeline <- function(pipeline_prefix, data_dir, sirv_gtf, tusco_tsv) { # Changed bugsi_gtf to tusco_tsv
  class_file <- file.path(data_dir, paste0(pipeline_prefix, "/", pipeline_prefix, "_classification.txt"))
  gtf_file   <- file.path(data_dir, paste0(pipeline_prefix, "/", pipeline_prefix, "_corrected.gtf")) # Still needed for SIRV

  classification_data <- read_tsv_safe(class_file)

  # SIRVs
  sirv_res <- compute_sirv_bigcats(classification_data, gtf_file, sirv_gtf) %>%
    mutate(Type="SIRVs", pipeline=pipeline_prefix)

  # TUSCO
  tusco_res <- compute_tusco_bigcats(classification_data, tusco_tsv) %>% # Pass tusco_tsv, removed gtf_file
    mutate(Type="TUSCO", pipeline=pipeline_prefix) # Changed from BUGSI

  bind_rows(sirv_res, tusco_res) # Changed from bugsi_res
}

###############################################################################
# 6) Main
###############################################################################
cat("\n--- Gathering data for all pipelines ---\n")

all_data <- lapply(pipelines, process_pipeline,
                   data_dir   = data_dir,
                   sirv_gtf   = sirv_gtf_file,
                   tusco_tsv  = tusco_tsv_file # Changed from bugsi_gtf
) %>%
  bind_rows()

cat("\n--- Single test (one-sided paired t-test: TUSCO > SIRVs) for each category ---\n") # Changed BUGSI to TUSCO
categories <- c("TP","PTP","FP","FN")

p_values <- c()
for (cat_name in categories) {
  cat("\nCategory:", cat_name, "\n")
  df_sub <- all_data %>% filter(big_category==cat_name)
  wide <- df_sub %>%
    select(pipeline, Type, percentage) %>%
    pivot_wider(names_from=Type, values_from=percentage)

  # Check if TUSCO and SIRVs columns exist and have enough data
  if ("TUSCO" %in% colnames(wide) && "SIRVs" %in% colnames(wide) && 
      sum(!is.na(wide$TUSCO)) > 1 && sum(!is.na(wide$SIRVs)) > 1) {
    # One-sided paired t-test: TUSCO > SIRVs
    test_out <- t.test(wide$TUSCO, wide$SIRVs, paired=TRUE, alternative="greater") # Changed wide$BUGSI
    pval <- test_out$p.value
    cat("p-value = ", pval, "\n")
  } else {
    cat("Not enough data or TUSCO/SIRVs columns missing for category:", cat_name, ". Skipping t-test.\n")
    pval <- NA # Assign NA if test cannot be performed
  }
  p_values <- c(p_values, pval)
}

# Build data to show significance label
# We'll do (TUSCO vs. SIRVs) for each category
my_signifs <- data.frame(
  big_category = categories,
  p_value = p_values
)

# Define function to convert p-values to asterisk significance (original version)
# p_stars_original <- function(x) {
#   if (is.na(x)) return("NA") # Handle NA p-values
#   if (x < 0.001) return("***")
#   else if (x < 0.01) return("**")
#   else if (x < 0.05) return("*")
#   else return("ns")
# }
# my_signifs$star_label_original <- sapply(my_signifs$p_value, p_stars_original)


# Define function to convert p-values to asterisk significance (second version in original script)
p_stars <- function(x) {
  if (is.na(x)) return("NA") # Handle NA p-values
  if (x < 1e-4) {
    return("****")
  } else if (x < 1e-3) {
    return("***")
  } else if (x < 1e-2) {
    return("**")
  } else if (x < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Convert p-values to asterisk labels
my_signifs$star_label <- sapply(my_signifs$p_value, p_stars) # Using the second p_stars function

cat("\n--- Summaries: ---\n")
print(my_signifs)


# Compute mean ± SD for each category × Type
mean_data <- all_data %>%
  group_by(big_category, Type) %>%
  summarise(
    mean_perc = mean(percentage, na.rm = TRUE),
    sd_perc   = sd(percentage, na.rm = TRUE),
    .groups   = "drop"
  )

# Ensure proper ordering on the x-axis
mean_data$big_category <- factor(mean_data$big_category, levels = c("TP","PTP","FP","FN"))
mean_data$Type         <- factor(mean_data$Type, levels = c("TUSCO","SIRVs")) # Changed BUGSI to TUSCO

# Create significance annotation data
bracket_data <- my_signifs %>%
  mutate(
    cat_index  = as.numeric(factor(big_category, levels = c("TP","PTP","FP","FN"))),
    y_position = 110, # Adjust this based on your data range
    annotations = star_label
  )

# Ensure correct ordering of categories
all_data$big_category <- factor(all_data$big_category, levels = c("TP", "PTP", "FP", "FN"))
bracket_data$big_category <- factor(bracket_data$big_category, levels = c("TP", "PTP", "FP", "FN"))

# Adjust positions for significance annotations and horizontal lines
bracket_data <- bracket_data %>%
  mutate(
    y_position = 110,         # Adjust for significance text
    y_line = 105,             # Position for horizontal lines
    xmin = as.numeric(big_category) - 0.2,
    xmax = as.numeric(big_category) + 0.2
  )

# Generate the barplot with adjusted settings
p_single <- ggplot(all_data, aes(x = big_category, y = percentage, fill = Type)) +
  
  # 1. Bars for mean with reduced thickness and narrower width
  stat_summary(
    fun    = mean,
    geom   = "bar",
    width  = 0.4,             # Reduced width for narrower bars
    color  = "black",
    size   = 0.1,             # Reduced border thickness
    position = position_dodge(width = 0.6)
  ) +
  
  # 2. Error bars for SD with minimum thickness
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom     = "errorbar",
    width    = 0.2,
    size     = 0.1,             # Further reduced error bar thickness
    position = position_dodge(width = 0.6)
  ) +
  
  # 3. Jittered raw data points with further reduced size
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size     = 0.2,             # Further reduced sample point size
    alpha    = 0.7,
    shape    = 16,
    color    = "black"
  ) +
  
  # 4. Significance labels with reduced size
  geom_text(
    data       = bracket_data,
    aes(
      x          = big_category,
      y          = y_position,
      label      = annotations
    ),
    inherit.aes = FALSE,
    size        = 3,             # Reduced size for smaller font
    fontface    = "bold"
  ) +
  
  # 5. Horizontal significance lines with further reduced thickness
  geom_segment(
    data       = bracket_data,
    aes(
      x    = xmin,
      xend = xmax,
      y    = y_line,
      yend = y_line
    ),
    inherit.aes = FALSE,
    size        = 0.15,          # Further reduced line thickness
    color       = "black"
  ) +
  
  # 6. Adjust colors and theme
  scale_fill_manual(values = c("TUSCO" = "#a8d5a0", "SIRVs" = "#cab2d6")) + # Changed BUGSI color
  coord_cartesian(ylim = c(0, 100), clip = "off") +  # Set y-axis max to 100
  labs(
    x     = NULL,
    y     = "Percentage"
  ) +
  theme_classic(base_size = 7) +  # Keep base font size unchanged
  theme(
    axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 7),
    axis.title.y = element_text(size = 7, face = "bold"),
    legend.position = "none",       # Remove the legend
    legend.title = element_blank(),
    plot.margin  = margin(15, 10, 10, 10)  # Adjusted top margin
  )

# Display the plot
print(p_single)

# Define save paths
save_dir <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/fig-3/plots" # Output plots to the script's directory
pdf_path <- file.path(save_dir, "figure3b-human.pdf") # Primary output named after script

# Create the directory if it doesn't exist (optional, but good practice)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Save as high-resolution PDF with the current figure size
ggsave(
  filename = pdf_path,
  plot     = p_single,
  width    = 2.5,   # Further reduced width in inches
  height   = 1.8,   # Further reduced height in inches
  units    = "in",
  device   = cairo_pdf
)

cat("Plot saved successfully:\n")
cat(" - PDF:", pdf_path, "\n")
cat(" - PNG:", png_path, "\n")


