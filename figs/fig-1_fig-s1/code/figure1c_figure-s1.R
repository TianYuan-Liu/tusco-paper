# Function to safely load packages
safe_load <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    # Try installing if not found
    # install.packages(package_name, dependencies = TRUE)
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      warning(paste("Package", package_name, "not available and couldn't be installed. Some functionality may be limited."))
      return(FALSE)
    }
  }
  return(TRUE)
}

# Required packages
suppressPackageStartupMessages({
  # Core packages
  safe_load("rtracklayer")
  safe_load("dplyr")
  safe_load("tidyr")
  safe_load("ggplot2")
  safe_load("stringr")
  safe_load("scales")

  # Optional packages for layout
  has_cowplot <- safe_load("cowplot")
  has_gridExtra <- safe_load("gridExtra")
  has_patchwork <- safe_load("patchwork")
})

###############################################################
# Custom color palette & Define Dataset Groups
###############################################################
palette_colors <- c(
  # TUSCO group - Greens
  "Human (TUSCO)"   = "#a8d5a0",  # dark green
  "Mouse (TUSCO)"   = "#1b9e77",  # light green

  # GENCODE/MANE/RefSeq group - Oranges/Reds
  "Human (GENCODE)" = "#fdbf6f",  # light orange
  "Mouse (GENCODE)" = "#e66101",  # darker orange
  "Human (MANE)"    = "#e41a1c",  # NEW: Red
  "Human (RefSeq)"  = "#a6761d",  # NEW: Golden Brown
  "Mouse (RefSeq)"  = "#663d00",  # NEW: Dark Brown

  # SIRVs/ERCCs/Sequin group - Purples (not plotted in main figures but kept for potential use)
  "SIRVs"           = "#cab2d6",  # light purple
  "ERCCs"           = "#6a3d9a",  # darker purple
  "Sequins"         = "#413aa3"   # medium purple (User updated)
)

# Datasets for the legend (excluding RefSeq as per image)
human_datasets_legend <- c("Human (TUSCO)", "Human (GENCODE)", "Human (MANE)")
mouse_datasets_legend <- c("Mouse (TUSCO)", "Mouse (GENCODE)")
spikein_datasets_legend <- c("SIRVs", "ERCCs", "Sequins") # Keep all spike-ins

# Combined groups for plotting (still include RefSeq in plots if data exists)
human_datasets_plot <- c("Human (TUSCO)", "Human (GENCODE)", "Human (MANE)", "Human (RefSeq)")
mouse_datasets_plot <- c("Mouse (TUSCO)", "Mouse (GENCODE)", "Mouse (RefSeq)")
spikein_datasets_plot <- c("SIRVs", "ERCCs", "Sequins")

human_plus_spikeins_plot <- c(human_datasets_plot, spikein_datasets_plot)
mouse_plus_spikeins_plot <- c(mouse_datasets_plot, spikein_datasets_plot)

# Define the specific order for the legend based on the image (filled by row)
# Updated legend order for the 6-panel figure (excluding spike-ins)
# legend_order <- c(
#     "Human (BUGSI)", "Human (GENCODE)", "Human (MANE)", "Human (RefSeq)",
#     "Mouse (BUGSI)", "Mouse (GENCODE)", "Mouse (RefSeq)"
# )
# New order based on target image (2 rows)
legend_order <- c(
    "Human (TUSCO)", "Mouse (TUSCO)", "Human (GENCODE)", "Mouse (GENCODE)",
    "Human (RefSeq)", "Mouse (RefSeq)", "Human (MANE)"
)

# Combine all dataset levels relevant for plotting factor ordering (includes RefSeq)
all_plot_datasets_ordered <- unique(c(human_datasets_plot, mouse_datasets_plot, spikein_datasets_plot))

# Define the order for exon count categories
exon_levels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "10+")

###############################################################
# Define input files (switched to repo-local data under figs/data)
###############################################################
# TUSCO transcript selection lists (TSV; no header)
tusco_human_tsv   <- "figs/data/tusco/tusco_human.tsv"
tusco_mouse_tsv   <- "figs/data/tusco/tusco_mouse.tsv"

# Reference annotations (GTF/GTF.GZ)
human_gencode_gtf <- "figs/data/reference/human/gencode.v48.annotation.gtf.gz"
mouse_gencode_gtf <- "figs/data/reference/mouse/gencode.vM37.annotation.gtf.gz"
mane_gtf          <- "figs/data/reference/human/mane.v1.4.ensembl_genomic.gtf.gz"  # human only
human_refseq_gtf  <- "figs/data/reference/human/refseq/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
mouse_refseq_gtf  <- "figs/data/reference/mouse/refseq/GCF_000001635.27_GRCm39_genomic.gtf.gz"

# Spike-ins
ercc_sirv_gtf     <- "figs/data/spike-ins/lrgasp_sirvs4.gtf"
sequin_gtf        <- "figs/data/spike-ins/rnasequin_annotation_2.4.gtf"

###############################################################
# Functions for data processing (largely unchanged, minor safety checks)
###############################################################
import_and_standardize <- function(gtf_file, dataset_type, expected_species = NULL) {
  message("Importing ", dataset_type, " from ", gtf_file)

  if (!file.exists(gtf_file)) {
    warning("File not found: ", gtf_file)
    return(NULL)
  }

  tryCatch({
    gtf <- import(gtf_file)
    df <- as.data.frame(gtf)

    # Basic check for empty GTF data
    if(nrow(df) == 0) {
        warning("GTF file is empty or could not be parsed: ", gtf_file)
        return(NULL)
    }

    # Standardize gene_id and transcript_id
    if (dataset_type == "TUSCO") {
      if(!"ensembl" %in% colnames(df)) {
        warning("TUSCO GTF expected to have 'ensembl' column for gene_id: ", gtf_file)
        # Attempt to find alternative if possible, e.g., gene_id or transcript_id if present
        if("gene_id" %in% colnames(df)) df$ensembl <- df$gene_id
        else if("transcript_id" %in% colnames(df)) df$ensembl <- df$transcript_id
        else return(NULL) # Cannot proceed without an identifier
      }
      df <- df %>%
        mutate(gene_id = ensembl,
               transcript_id = ensembl) # TUSCO assumes single isoform per 'ensembl' ID
    } else {
        # Ensure gene_id exists
        if (!"gene_id" %in% colnames(df)) {
            # Try alternatives
            if ("gene" %in% colnames(df)) df$gene_id <- df$gene
            else if ("Name" %in% colnames(df) && any(grepl("gene", df$type, ignore.case = TRUE))) {
                # Heuristic: Try to extract from Name attribute if type is gene-like
                gene_entries <- df %>% filter(grepl("gene", type, ignore.case=TRUE))
                if(nrow(gene_entries) > 0 && "Name" %in% colnames(gene_entries)) {
                   # This requires more complex logic to map back to transcripts/exons, skipping for now
                   warning("Complex gene_id inference needed for ", gtf_file, ". Requires 'gene_id' column.")
                   return(NULL)
                } else {
                   warning("No 'gene_id' column or suitable alternative found in non-TUSCO GTF: ", gtf_file)
                   return(NULL)
                }
            } else {
                warning("No 'gene_id' column or suitable alternative found in non-TUSCO GTF: ", gtf_file)
                return(NULL)
            }
        }
        # Ensure transcript_id exists
        if (!"transcript_id" %in% colnames(df)) {
            if ("transcript" %in% colnames(df)) df$transcript_id <- df$transcript
            else if ("rna_id" %in% colnames(df)) df$transcript_id <- df$rna_id # Handle NCBI GTF
            else if ("transcript_name" %in% colnames(df)) df$transcript_id <- df$transcript_name
            else {
                warning("No 'transcript_id' column or suitable alternative found in non-TUSCO GTF: ", gtf_file)
                return(NULL)
            }
        }
    }

    # Minimal check for necessary columns after standardization
    if (!all(c("gene_id", "transcript_id", "type", "width", "start", "end", "seqnames") %in% colnames(df))) {
      warning("Standardized dataframe missing essential columns for ", dataset_type, " from ", gtf_file)
      return(NULL)
    }

    return(df)
  }, error = function(e) {
    warning("Error importing ", dataset_type, " from ", gtf_file, ": ", e$message)
    return(NULL)
  })
}

extract_features <- function(df, dataset_name) {
  if (is.null(df) || nrow(df) == 0) {
    warning("Empty or NULL dataframe provided for feature extraction: ", dataset_name)
    return(NULL)
  }

  # Ensure essential columns exist
  required_cols <- c("type", "transcript_id", "gene_id", "width", "start", "end")
  if (!all(required_cols %in% colnames(df))) {
      missing_cols <- setdiff(required_cols, colnames(df))
      warning("Dataframe for ", dataset_name, " is missing required columns: ", paste(missing_cols, collapse=", "))
      return(NULL)
  }

  # Filter exons
  exons <- df %>% filter(type == "exon") %>% distinct(transcript_id, start, end, .keep_all = TRUE) # Ensure unique exons per transcript

  if (nrow(exons) == 0) {
    warning("No features of type 'exon' found in dataset: ", dataset_name)
    return(NULL) # Return NULL if no exons, can't calculate features
  }

  # Calculate exon lengths
  exons <- exons %>% mutate(exon_length = width)

  # Transcript length (sum of unique exon lengths per transcript)
  transcript_lengths <- exons %>%
    group_by(transcript_id) %>%
    summarize(
      transcript_length = sum(exon_length, na.rm = TRUE),
      gene_id = first(gene_id), # Take the first gene_id associated with the transcript
      .groups = "drop"
    ) %>% filter(transcript_length > 0) # Remove transcripts with zero length

  if (nrow(transcript_lengths) == 0) {
      warning("No valid transcript lengths calculated for: ", dataset_name)
      # Don't return NULL entirely, maybe other features are valid
  }

  # Exon counts per transcript
  exon_counts <- exons %>%
    group_by(transcript_id) %>%
    summarize(
      exon_count = n(),
      gene_id = first(gene_id),
      .groups = "drop"
    )

  if (nrow(exon_counts) == 0) {
      warning("No valid exon counts calculated for: ", dataset_name)
  }

  # Introns: gaps between consecutive, ordered exons within the same transcript
  introns <- exons %>%
    filter(!is.na(transcript_id) & !is.na(start) & !is.na(end)) %>% # Ensure valid coordinates
    arrange(transcript_id, start) %>%
    group_by(transcript_id) %>%
    mutate(next_exon_start = lead(start)) %>%
    filter(!is.na(next_exon_start)) %>% # Only consider exons that have a following exon
    mutate(intron_length = (next_exon_start - 1) - end) %>% # Calculate gap
    filter(intron_length > 0) %>% # Keep only positive gaps (introns)
    ungroup() %>%
    select(transcript_id, intron_length, gene_id) # Keep relevant columns

  if (nrow(introns) == 0) {
      # This is expected for single-exon transcripts, so maybe just message, not warn
      message("No introns calculated for: ", dataset_name, " (potentially many single-exon transcripts)")
  }

  # Transcripts per gene (based on transcripts with calculated lengths)
  transcripts_per_gene <- transcript_lengths %>%
    filter(!is.na(gene_id)) %>%
    group_by(gene_id) %>%
    summarize(
      transcript_count = n(),
      .groups = "drop"
    )

  # Add dataset name to all resulting data frames
  if(nrow(transcript_lengths) > 0) transcript_lengths$dataset <- dataset_name
  if(nrow(introns) > 0) introns$dataset <- dataset_name
  if(nrow(exons) > 0) exons$dataset <- dataset_name # Keep original exons data if needed elsewhere
  if(nrow(exon_counts) > 0) exon_counts$dataset <- dataset_name
  if(nrow(transcripts_per_gene) > 0) transcripts_per_gene$dataset <- dataset_name

  # Return list, handling cases where some features might be empty
  feature_list <- list(
    transcripts = if(nrow(transcript_lengths) > 0) transcript_lengths else NULL,
    introns = if(nrow(introns) > 0) introns else NULL,
    exons = if(nrow(exons) > 0) exons else NULL, # Full exon data
    exon_counts = if(nrow(exon_counts) > 0) exon_counts else NULL,
    tpg = if(nrow(transcripts_per_gene) > 0) transcripts_per_gene else NULL
  )

  # Remove NULL elements from the list before returning
  feature_list <- feature_list[!sapply(feature_list, is.null)]
  if(length(feature_list) == 0) return(NULL) # Return NULL if nothing was extracted

  return(feature_list)
}


extract_ercc_sirv <- function(df) {
  # Simplified: Assumes df is already standardized
  if (is.null(df) || nrow(df) == 0) {
    warning("Empty dataframe provided for ERCC/SIRV extraction")
    return(list(ERCC = NULL, SIRV = NULL))
  }
  ercc_df <- df %>% filter(grepl("^ERCC", seqnames, ignore.case = TRUE))
  sirv_df <- df %>% filter(grepl("^SIRV", seqnames, ignore.case = TRUE))
  return(list(ERCC = ercc_df, SIRV = sirv_df))
}

###############################################################
# Load datasets and Extract Features
###############################################################

# Helpers to build TUSCO subsets from GENCODE using TSV transcript lists
read_tusco_transcripts <- function(tsv_path) {
  if (!file.exists(tsv_path)) {
    warning("TUSCO TSV not found: ", tsv_path)
    return(character(0))
  }
  tbl <- tryCatch({
    read.table(tsv_path, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
  }, error = function(e) {
    warning("Failed to read TUSCO TSV: ", tsv_path, "; ", e$message)
    return(NULL)
  })
  if (is.null(tbl) || nrow(tbl) == 0) return(character(0))
  # Expect column 2 to be Ensembl transcript ID (ENST*/ENSMUST*)
  tx_ids <- unique(tbl[[2]])
  tx_ids <- tx_ids[!is.na(tx_ids) & nzchar(tx_ids)]
  return(tx_ids)
}

subset_to_transcripts <- function(df, transcript_ids) {
  if (is.null(df) || length(transcript_ids) == 0) return(NULL)
  if (!"transcript_id" %in% colnames(df)) return(NULL)
  # Normalize IDs by removing version suffix (e.g., ENSTxxxx.y -> ENSTxxxx)
  normalize_id <- function(x) {
    x <- as.character(x)
    sub("\\..*$", "", x)
  }
  tx_core <- unique(normalize_id(transcript_ids))
  df_core <- df %>% mutate(transcript_id_core = normalize_id(transcript_id))
  df_core %>% filter(!is.na(transcript_id_core) & transcript_id_core %in% tx_core) %>% select(-transcript_id_core)
}

datasets_list <- list()

# Import base reference annotations
human_gencode_df <- import_and_standardize(human_gencode_gtf, "GENCODE", "human")
mouse_gencode_df <- import_and_standardize(mouse_gencode_gtf, "GENCODE", "mouse")
human_mane_df    <- import_and_standardize(mane_gtf,          "MANE",    "human")
human_refseq_df  <- import_and_standardize(human_refseq_gtf,  "RefSeq",  "human")
mouse_refseq_df  <- import_and_standardize(mouse_refseq_gtf,  "RefSeq",  "mouse")

# Build TUSCO subsets
human_tusco_ids <- read_tusco_transcripts(tusco_human_tsv)
mouse_tusco_ids <- read_tusco_transcripts(tusco_mouse_tsv)
human_tusco_df  <- subset_to_transcripts(human_gencode_df, human_tusco_ids)
mouse_tusco_df  <- subset_to_transcripts(mouse_gencode_df, mouse_tusco_ids)

# Assemble datasets list
datasets_list[["Human (TUSCO)"]]  <- extract_features(human_tusco_df,   "Human (TUSCO)")
datasets_list[["Mouse (TUSCO)"]]  <- extract_features(mouse_tusco_df,   "Mouse (TUSCO)")
datasets_list[["Human (GENCODE)"]] <- extract_features(human_gencode_df, "Human (GENCODE)")
datasets_list[["Mouse (GENCODE)"]] <- extract_features(mouse_gencode_df, "Mouse (GENCODE)")
datasets_list[["Human (MANE)"]]    <- extract_features(human_mane_df,    "Human (MANE)")
datasets_list[["Human (RefSeq)"]]  <- extract_features(human_refseq_df,  "Human (RefSeq)")
datasets_list[["Mouse (RefSeq)"]]  <- extract_features(mouse_refseq_df,  "Mouse (RefSeq)")

# Spike-ins (Load but don't assign to human/mouse groups initially)
ercc_sirv_df <- import_and_standardize(ercc_sirv_gtf, "ERCC_SIRV")
if (!is.null(ercc_sirv_df)) {
  es_split <- extract_ercc_sirv(ercc_sirv_df)
  datasets_list[["ERCCs"]] <- extract_features(es_split$ERCC, "ERCCs")
  datasets_list[["SIRVs"]] <- extract_features(es_split$SIRV, "SIRVs")
}
datasets_list[["Sequins"]] <- extract_features(import_and_standardize(sequin_gtf, "Sequins"), "Sequins")

# Remove NULL entries from the list (failed loads/extractions)
datasets_list <- datasets_list[!sapply(datasets_list, is.null)]

###############################################################
# Combine Features Separately for Human and Mouse (Main Figures - NO Spike-ins)
###############################################################
combine_features_subset <- function(datasets_list, feature_name, subset_datasets_group, all_factor_levels) {
  # Filter the main list for the subset datasets
  subset_list <- datasets_list[names(datasets_list) %in% subset_datasets_group]

  # Extract the specified feature from each dataset in the subset
  feature_list <- lapply(subset_list, function(x) {
      if (!is.null(x) && feature_name %in% names(x)) {
          return(x[[feature_name]])
      } else {
          return(NULL) # Return NULL if feature doesn't exist for this dataset
      }
  })

  # Filter out NULL feature dataframes
  valid_features <- feature_list[!sapply(feature_list, is.null)]

  if (length(valid_features) == 0) {
    warning("No valid data for feature '", feature_name, "' in the subset: ", paste(subset_datasets_group, collapse=", "))
    return(NULL)
  }

  # Combine the valid feature dataframes
  combined_data <- bind_rows(valid_features)

  # Ensure the dataset column is a factor with the specified levels
  combined_data <- combined_data %>%
    mutate(dataset = factor(dataset, levels = intersect(all_factor_levels, unique(dataset))))

  return(combined_data)
}

# Combine Features for Human ONLY (Main Plots)
transcript_lengths_human <- combine_features_subset(datasets_list, "transcripts", human_datasets_plot, all_plot_datasets_ordered)
intron_lengths_human <- combine_features_subset(datasets_list, "introns", human_datasets_plot, all_plot_datasets_ordered)
exon_counts_human <- combine_features_subset(datasets_list, "exon_counts", human_datasets_plot, all_plot_datasets_ordered)

# Combine Features for Mouse ONLY (Main Plots)
transcript_lengths_mouse <- combine_features_subset(datasets_list, "transcripts", mouse_datasets_plot, all_plot_datasets_ordered)
intron_lengths_mouse <- combine_features_subset(datasets_list, "introns", mouse_datasets_plot, all_plot_datasets_ordered)
exon_counts_mouse <- combine_features_subset(datasets_list, "exon_counts", mouse_datasets_plot, all_plot_datasets_ordered)

# Define final legend order and apply it to combined datasets
final_legend_order <- c(
    "Human (TUSCO)", "Mouse (TUSCO)", "Human (GENCODE)", "Mouse (GENCODE)",
    "Human (RefSeq)", "Mouse (RefSeq)", "Human (MANE)"
)

set_final_order <- function(df, order_levels) {
  if (!is.null(df) && nrow(df) > 0 && "dataset" %in% colnames(df)) {
    present_levels <- intersect(order_levels, unique(df$dataset))
    if (length(present_levels) > 0) {
      df$dataset <- factor(df$dataset, levels = present_levels)
    }
  }
  return(df)
}

transcript_lengths_human <- set_final_order(transcript_lengths_human, final_legend_order)
intron_lengths_human     <- set_final_order(intron_lengths_human, final_legend_order)
exon_counts_human        <- set_final_order(exon_counts_human, final_legend_order)
transcript_lengths_mouse <- set_final_order(transcript_lengths_mouse, final_legend_order)
intron_lengths_mouse     <- set_final_order(intron_lengths_mouse, final_legend_order)
exon_counts_mouse        <- set_final_order(exon_counts_mouse, final_legend_order)

###############################################################
# Nature-style Theme for Plots (Unchanged)
###############################################################
nature_theme <- theme_classic(base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0), # Titles aligned left
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Default angle 0 for x-axis text
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold"),
    axis.line = element_line(linewidth = 0.25), # Use linewidth instead of size
    axis.ticks = element_line(linewidth = 0.25),
    legend.key.size = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom", # Default legend position
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt"), # Small top margin for legend box
    legend.box.background = element_rect(color = "black", linewidth = 0.5) # Add black border
  )

###############################################################
# Define Plotting Functions
###############################################################

# 1. Transcript Length Density Plot
plot_transcript_density <- function(data, title, species_colors) {
  if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Data)")) + nature_theme)
  
  # Filter out zero-length transcripts just in case
  data <- data %>% filter(transcript_length > 0)
  if (nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Positive Length Data)")) + nature_theme)

  ggplot(data, aes(x = transcript_length, color = dataset)) +
    geom_density(linewidth = 0.5, key_glyph = "path") + # Use linewidth, path glyph better for density lines
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5), # Auto breaks
      labels = scales::trans_format("log10", scales::math_format(10^.x)), # Nicer labels
      limits = c(NA, 1e5) # Ensure upper limit consistent if needed
    ) +
    annotation_logticks(sides = "b", linewidth = 0.25, outside = TRUE, # Add log ticks below axis
                        short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    coord_cartesian(clip = "off") + # Allow ticks outside panel
    scale_color_manual(values = species_colors, name = "Dataset") +
    labs(
      title = title,
      x = "Transcript Length (nt, log scale)",
      y = "Density"
    ) +
    nature_theme +
    theme(plot.margin = margin(t = 5, r = 5, b = 15, l = 5, unit = "pt")) # Adjust margin for ticks
}

# 2. Exon Counts Bar Plot
plot_exon_counts <- function(data, title, species_colors) {
  if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Data)")) + nature_theme)

  # Prepare data for plotting: group >10, calculate frequencies
  exon_data_for_plot <- data %>%
    mutate(exon_count_group = ifelse(exon_count > 10, "10+", as.character(exon_count))) %>%
    mutate(exon_count_group = factor(exon_count_group, levels = exon_levels)) %>%
    # Calculate frequency PER DATASET
    group_by(dataset) %>%
    mutate(total_transcripts_in_dataset = n()) %>%
    ungroup() %>%
    group_by(dataset, exon_count_group, total_transcripts_in_dataset) %>%
    summarize(count = n(), .groups = "drop") %>%
    mutate(freq = (count / total_transcripts_in_dataset) * 100) %>%
    ungroup()

  if (nrow(exon_data_for_plot) == 0) return(ggplot() + ggtitle(paste(title, "(Processing Error)")) + nature_theme)

  ggplot(exon_data_for_plot,
         aes(x = exon_count_group, y = freq, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.2, color = "black") + # Add outline
    scale_fill_manual(values = species_colors, name = "Dataset") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) + # Y starts at 0
    labs(
      title = title,
      x = "Number of Exons per Transcript",
      y = "Frequency (%)"
    ) +
    nature_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate labels if many groups
}

# 3. Intron Length Density Plot
plot_intron_density <- function(data, title, species_colors) {
  if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Intron Data)")) + nature_theme)
  
  # Ensure positive intron lengths
  data <- data %>% filter(intron_length > 0)
  if (nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Positive Length Introns)")) + nature_theme)

  ggplot(data, aes(x = intron_length, color = dataset)) +
    geom_density(linewidth = 0.5, key_glyph = "path") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6), # More breaks for introns
      labels = scales::trans_format("log10", scales::math_format(10^.x))
      # Consider setting limits if needed, e.g., limits = c(10, 1e6)
    ) +
    annotation_logticks(sides = "b", linewidth = 0.25, outside = TRUE,
                        short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = species_colors, name = "Dataset") +
    labs(
      title = title,
      x = "Intron Length (nt, log scale)",
      y = "Density"
    ) +
    nature_theme +
    theme(plot.margin = margin(t = 5, r = 5, b = 15, l = 5, unit = "pt")) # Adjust margin for ticks
}


###############################################################
# Generate Plots for Human (No Spike-ins)
###############################################################
# Filter palette for Human datasets only
human_palette <- palette_colors[names(palette_colors) %in% human_datasets_plot]

plot_h_transcript_density <- plot_transcript_density(transcript_lengths_human, "", human_palette)
plot_h_exon_counts <- plot_exon_counts(exon_counts_human, "", human_palette)
plot_h_intron_density <- plot_intron_density(intron_lengths_human, "", human_palette)

###############################################################
# Generate Plots for Mouse (No Spike-ins)
###############################################################
# Filter palette for Mouse datasets only
mouse_palette <- palette_colors[names(palette_colors) %in% mouse_datasets_plot]

plot_m_transcript_density <- plot_transcript_density(transcript_lengths_mouse, "", mouse_palette)
plot_m_exon_counts <- plot_exon_counts(exon_counts_mouse, "", mouse_palette)
plot_m_intron_density <- plot_intron_density(intron_lengths_mouse, "", mouse_palette)

###############################################################
# Generate Specific Transcript Length Plot (Human/Mouse TUSCO/GENCODE/RefSeq + ERCCs + SIRVs)
###############################################################
message("Generating specific transcript length plot...")

# Define the target datasets for the specific plot
target_datasets_specific <- c(
  "Human (TUSCO)",
  "Mouse (TUSCO)",
  "Human (GENCODE)",
  "Mouse (GENCODE)",
  "SIRVs",
  "ERCCs"
)

# Extract and combine transcript length data for target datasets
transcript_lengths_specific_list <- lapply(target_datasets_specific, function(ds_name) {
  if (!is.null(datasets_list[[ds_name]]) && "transcripts" %in% names(datasets_list[[ds_name]])) {
    return(datasets_list[[ds_name]][["transcripts"]])
  } else {
    warning("Transcript data not found for: ", ds_name)
    return(NULL)
  }
})

# Filter out NULLs and bind rows
valid_transcript_lengths_specific <- transcript_lengths_specific_list[!sapply(transcript_lengths_specific_list, is.null)]
if (length(valid_transcript_lengths_specific) > 0) {
  transcript_lengths_specific <- bind_rows(valid_transcript_lengths_specific)

  # Ensure dataset is a factor with the correct levels (subset of the full list)
  transcript_lengths_specific <- transcript_lengths_specific %>%
      mutate(dataset = factor(dataset, levels = intersect(target_datasets_specific, unique(dataset))))

  # Filter the palette for the target datasets
  specific_palette <- palette_colors[names(palette_colors) %in% target_datasets_specific]

  # Create the plot
  plot_specific_transcript_density <- plot_transcript_density(
    transcript_lengths_specific,
    "", # No title
    specific_palette
  )

  # Define output path and save
  specific_plot_filename <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/fig-1/plots/transcript_length_distribution.pdf"
  # Ensure output directory exists
  try({ dir.create(dirname(specific_plot_filename), recursive = TRUE, showWarnings = FALSE) }, silent = TRUE)
  # Using a smaller size, e.g., half-width, adjusted height
  ggsave(
      specific_plot_filename,
      plot_specific_transcript_density,
      width = 85 / 25.4, # Approx half Nature column width in inches
      height = 70 / 35, # Adjusted height in inches
      dpi = 300,
      device = cairo_pdf
  )
  message("Saved specific transcript length plot: ", specific_plot_filename)

} else {
  warning("No valid transcript data found for the specific plot datasets. Skipping generation.")
}

###############################################################
# Arrange and Save Figures (Using Patchwork ONLY)
###############################################################

# Helper function using patchwork for arrangement and legend
arrange_and_save_combined <- function(plots_list, filename_base, fig_width_mm, fig_height_mm, ncol = 2, labels = c('a', 'b', 'c', 'd', 'e', 'f')) {

  # Ensure patchwork is available
  if (!requireNamespace("patchwork", quietly = TRUE)) {
     stop("patchwork package is required for arranging plots. Please install it.")
  }
 
  # Remove individual legends
  plots_no_legend <- lapply(plots_list, function(p) p + theme(legend.position = "none"))

  # Arrange using patchwork (assuming 6 plots, 2 cols)
  if (length(plots_no_legend) == 6 && ncol == 2) {
      p1 <- plots_no_legend[[1]] # Human Transcripts
      p2 <- plots_no_legend[[2]] # Mouse Transcripts
      p3 <- plots_no_legend[[3]] # Human Exons
      p4 <- plots_no_legend[[4]] # Mouse Exons
      p5 <- plots_no_legend[[5]] # Human Introns
      p6 <- plots_no_legend[[6]] # Mouse Introns

      # Arrange in 2 columns, 3 rows
      arranged_plots <- (p1 | p2) / (p3 | p4) / (p5 | p6)

      # Add common legend (collected based on plot data factors) and labels
      # Specify legend options for the COLLECTED legend via theme()
      # Use guides() to force nrow = 2 for the collected legend
      final_figure <- arranged_plots +
          plot_layout(guides = 'collect') +
          plot_annotation(tag_levels = 'a') & # Add labels
          theme(legend.position = 'bottom',       # Position collected legend at the bottom
                legend.direction = "horizontal",  # Arrange items horizontally
                legend.box = "horizontal",       # Box layout
                legend.title = element_blank(),   # No legend title
                legend.text = element_text(size = 6),
                legend.key.size = unit(0.4, "lines"), # Slightly smaller key size
                legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), # Increased margin around legend content
                legend.box.background = element_rect(color = "black", linewidth = 0.5), # Add black border
                plot.tag = element_text(size = 7, face = "bold") # Label styling
               ) &
          guides(color = "none", fill = guide_legend(nrow = 2)) # Disable COLOR guide, collect+format FILL guide

  } else {
     stop("Plot arrangement logic currently assumes 6 plots in 2 columns.")
  }

  # Save the final figure (applies to both patchwork and fallback)
  width_in <- fig_width_mm / 25.4
  height_in <- fig_height_mm / 25.4

  pdf_filename <- paste0(filename_base, ".pdf")
  png_filename <- paste0(filename_base, ".png")

  # Use ggsave for patchwork object directly
  ggsave(pdf_filename, final_figure, width = width_in, height = height_in, dpi = 300, device = cairo_pdf)
  message("Saved: ", pdf_filename)
  ggsave(png_filename, final_figure, width = width_in, height = height_in, dpi = 300, type = "cairo")
  message("Saved: ", png_filename)
}

# Define figure dimensions (e.g., half Nature width, adjustable height)
fig_width_mm <- 170 # Full Nature width for 2 columns
fig_height_mm <- 180 # Adjust height as needed for 3 rows + legend

# Combine all plots for the final figure (using plots without spike-ins)
all_plots_main <- list(
  plot_h_transcript_density, plot_m_transcript_density,
  plot_h_exon_counts, plot_m_exon_counts,
  plot_h_intron_density, plot_m_intron_density
)

# Arrange and Save Combined 6-Panel Figure (No Spike-ins)
arrange_and_save_combined(
  plots_list = all_plots_main,
  filename_base = "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/fig-1/plots/Combined_Comparative_Figure",
  fig_width_mm = fig_width_mm,
  fig_height_mm = fig_height_mm,
  ncol = 2 # Arrange in 2 columns
)