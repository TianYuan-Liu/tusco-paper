###########################
### Figure 4c (R, per-read)
###########################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(rtracklayer)
})

# Suppress data.table NSE notes for R CMD check/lint
utils::globalVariables(c(
  "exon_len", ".", "transcript_id", "associated_transcript_nov", "associated_transcript",
  "structural_category", "label", "subcategory", "pbid", "read_length",
  "exonic_length", "id", "species", "Species", "Group", "PlotGroup"
))

# ---- 0. Paths ----
BASE_DIR <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
DATA_BASE <- file.path(BASE_DIR, "figs", "data", "lrgasp", "tusco_novel_evl")
OUT_DIR_PLOTS <- file.path(BASE_DIR, "figs", "fig-4", "plots")

READ_STATS_FILES <- list(
  wtc11_cdna_pacbio = file.path(DATA_BASE, "wtc11_cdna_pacbio.transcriptome.read_stat.with_len.txt"),
  es_cdna_pacbio    = file.path(DATA_BASE, "es_cdna_pacbio.transcriptome.read_stat.with_len.txt")
)

CLASS_FILES_REF <- list(
  wtc11_cdna_pacbio = file.path(DATA_BASE, "isoseq_sq3", "ref_evl", "wtc11_cdna_pacbio", "sqanti3_out", "isoforms_classification.txt"),
  es_cdna_pacbio    = file.path(DATA_BASE, "isoseq_sq3", "ref_evl", "es_cdna_pacbio", "sqanti3_out", "isoforms_classification.txt")
)

# Add TUSCO novel evaluation classifications
CLASS_FILES_NOVEL <- list(
  wtc11_cdna_pacbio = file.path(DATA_BASE, "isoseq_sq3", "novel_evl", "wtc11_cdna_pacbio", "sqanti3_out", "isoforms_classification.txt"),
  es_cdna_pacbio    = file.path(DATA_BASE, "isoseq_sq3", "novel_evl", "es_cdna_pacbio", "sqanti3_out", "isoforms_classification.txt")
)

TUSCO_TSV <- list(
  human = file.path(BASE_DIR, "figs", "data", "tusco", "tusco_human.tsv"),
  mouse = file.path(BASE_DIR, "figs", "data", "tusco", "tusco_mouse.tsv")
)

TUSCO_GTF <- list(
  human = file.path(BASE_DIR, "figs", "data", "tusco", "tusco_human.gtf"),
  mouse = file.path(BASE_DIR, "figs", "data", "tusco", "tusco_mouse.gtf")
)

OUTPUT_PDF_COMBINED <- file.path(OUT_DIR_PLOTS, "fig-4c.pdf")
if (!dir.exists(OUT_DIR_PLOTS)) dir.create(OUT_DIR_PLOTS, recursive = TRUE)

# This script now only generates the combined fig-4c.pdf

strip_version <- function(x) {
  ifelse(is.na(x) | is.null(x), NA_character_, gsub("\\.[0-9]+$", "", as.character(x)))
}

# ---- 1. Load TUSCO transcript sets ----
load_tusco_transcripts <- function(tsv_path) {
  dt <- fread(tsv_path, sep = "\t", header = FALSE,
              col.names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
              colClasses = "character", data.table = TRUE, showProgress = FALSE)
  tx <- unique(na.omit(strip_version(dt$transcript)))
  return(tx)
}

cat("Loading TUSCO transcript lists...\n")
tusco_tx <- list(
  human = load_tusco_transcripts(TUSCO_TSV$human),
  mouse = load_tusco_transcripts(TUSCO_TSV$mouse)
)

# ---- 2. Parse GTF to compute exonic lengths per transcript ----
compute_exonic_lengths <- function(gtf_path, transcripts_of_interest) {
  if (!file.exists(gtf_path)) stop("GTF not found: ", gtf_path)
  gr <- rtracklayer::import(gtf_path)
  df <- as.data.frame(gr)
  if (!"type" %in% names(df)) {
    if (!"type" %in% names(mcols(gr))) stop("Cannot find feature type column in GTF: ", gtf_path)
    df$type <- mcols(gr)$type
  }
  if (!"transcript_id" %in% colnames(df)) {
    if ("transcript_id" %in% colnames(mcols(gr))) {
      df$transcript_id <- as.character(mcols(gr)$transcript_id)
    } else {
      stop("transcript_id not found in GTF attributes: ", gtf_path)
    }
  }
  exon_df <- df[df$type == "exon", c("seqnames","start","end","transcript_id")]
  exon_df$transcript_id <- strip_version(exon_df$transcript_id)
  toi <- unique(transcripts_of_interest[!is.na(transcripts_of_interest)])
  exon_df <- exon_df[exon_df$transcript_id %in% toi, , drop = FALSE]
  if (nrow(exon_df) == 0) return(integer(0))
  exon_dt <- as.data.table(exon_df)
  exon_dt[, exon_len := as.integer(end - start + 1L)]
  len_map <- exon_dt[, .(exonic_length = sum(exon_len, na.rm = TRUE)), by = transcript_id]
  setnames(len_map, "transcript_id", "tx")
  lens <- len_map$exonic_length
  names(lens) <- len_map$tx
  return(lens)
}

cat("Computing exonic lengths...\n")
exonic_lengths <- list(
  human = compute_exonic_lengths(TUSCO_GTF$human, tusco_tx$human),
  mouse = compute_exonic_lengths(TUSCO_GTF$mouse, tusco_tx$mouse)
)

# ---- 2b. Compute exon counts to remove mono-exon transcripts ----
compute_exon_counts <- function(gtf_path) {
  if (!file.exists(gtf_path)) return(data.table(transcript_id = character(), exon_count = integer()))
  gr <- rtracklayer::import(gtf_path)
  df <- as.data.frame(gr)
  if (!"type" %in% names(df)) {
    if (!"type" %in% names(mcols(gr))) stop("Cannot find feature type column in GTF: ", gtf_path)
    df$type <- mcols(gr)$type
  }
  if (!"transcript_id" %in% colnames(df)) {
    if ("transcript_id" %in% colnames(mcols(gr))) {
      df$transcript_id <- as.character(mcols(gr)$transcript_id)
    } else {
      stop("transcript_id not found in GTF attributes: ", gtf_path)
    }
  }
  exon_df <- df[df$type == "exon", c("transcript_id"), drop = FALSE]
  if (nrow(exon_df) == 0) return(data.table(transcript_id = character(), exon_count = integer()))
  exon_df$transcript_id <- strip_version(exon_df$transcript_id)
  as.data.table(exon_df)[, .(exon_count = .N), by = transcript_id]
}

exon_counts <- list(
  human = compute_exon_counts(TUSCO_GTF$human),
  mouse = compute_exon_counts(TUSCO_GTF$mouse)
)
multi_exon_tx <- list(
  human = unique(exon_counts$human[exon_count >= 2, transcript_id]),
  mouse = unique(exon_counts$mouse[exon_count >= 2, transcript_id])
)

# ---- 3. Load SQANTI classification, label TP/PTP for TUSCO transcripts ----
load_sqanti_classification <- function(paths, tusco_tx_set) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(data.table(pbid = character(), associated_transcript_nov = character(), label = character()))

  usecols <- c("isoform","structural_category","associated_transcript","subcategory")
  lst <- lapply(paths, function(p) fread(p, sep = "\t", select = usecols, header = TRUE, data.table = TRUE, showProgress = FALSE))
  df <- rbindlist(lst, use.names = TRUE, fill = TRUE)
  setnames(df, "isoform", "pbid")
  df[, associated_transcript_nov := strip_version(associated_transcript)]

  # Keep only mappings to TUSCO transcripts
  df <- df[associated_transcript_nov %in% tusco_tx_set]
  if (nrow(df) == 0) return(df[0])

  # Keep only FSM/ISM like Python logic, then label TP/PTP (TP if reference_match)
  df <- df[structural_category %in% c("full-splice_match", "incomplete-splice_match")]
  if (nrow(df) == 0) return(df[0])

  df[, label := NA_character_]
  df[subcategory == "reference_match", label := "TP"]
  df[is.na(label), label := "PTP"]

  # Collapse to transcript-level label: if any TP for a transcript, label TP else PTP
  tx_lab <- df[, .(label = ifelse(any(label == "TP", na.rm = TRUE), "TP", "PTP")), by = associated_transcript_nov]
  df <- unique(merge(df[, .(pbid, associated_transcript_nov)], tx_lab, by = "associated_transcript_nov", all.x = TRUE))
  df
}

## Removed variant loader with gene/TSS; not used for fig-4c

## Removed support summary loader; not used for fig-4c

## Removed support summary writer; not used for fig-4c

# ---- 4. Load read stats (keep ALL reads), compute f per read ----
load_read_stats_filter <- function(file_path, keep_pbids) {
  usecols <- c("id","pbid","read_length")
  dt <- fread(file_path, sep = "\t", header = TRUE, select = usecols, data.table = TRUE, showProgress = FALSE)
  if (length(keep_pbids)) dt <- dt[pbid %in% keep_pbids]
  dt <- dt[!is.na(pbid) & !is.na(read_length)]
  return(dt)
}

sample_species <- c(
  wtc11_cdna_pacbio = "human",
  es_cdna_pacbio    = "mouse"
)

build_reads_dt <- function(class_files_map) {
  all_reads <- list()
  for (sample in names(class_files_map)) {
    class_paths <- as.character(class_files_map[[sample]])
    class_paths <- class_paths[file.exists(class_paths)]
    if (length(class_paths) == 0) { cat("[WARN] Missing classification for:", sample, "\n"); next }
    species <- sample_species[[sample]]
    tusco_tx_set <- tusco_tx[[species]]
    tx_len_map <- exonic_lengths[[species]]

    class_df <- load_sqanti_classification(class_paths, tusco_tx_set)
    if (nrow(class_df) == 0) { cat("[WARN] No TUSCO TP/PTP in:", sample, "\n"); next }

    keep_pbids <- unique(class_df$pbid)

    read_stats_path <- READ_STATS_FILES[[sample]]
    if (!file.exists(read_stats_path)) { cat("[WARN] Missing read stats:", read_stats_path, "\n"); next }
    reads_dt <- load_read_stats_filter(read_stats_path, keep_pbids)
    if (nrow(reads_dt) == 0) { cat("[WARN] No reads for pbids in:", sample, "\n"); next }

    merged <- merge(reads_dt, class_df, by = "pbid", all.x = TRUE, allow.cartesian = TRUE)
    # Remove mono-exon transcripts for this species
    if (length(multi_exon_tx[[species]])) {
      merged <- merged[associated_transcript_nov %in% multi_exon_tx[[species]]]
    } else {
      merged <- merged[0]
    }
    merged[, exonic_length := as.numeric(tx_len_map[associated_transcript_nov])]
    merged <- merged[!is.na(exonic_length)]
    merged[, f := as.numeric(read_length) / as.numeric(exonic_length)]
    merged[, species := species]

    all_reads[[length(all_reads) + 1L]] <- merged[, .(id, pbid, f, label, species)]
  }
  if (length(all_reads) == 0) return(data.table())
  rbindlist(all_reads, use.names = TRUE, fill = TRUE)
}

## Removed TSS-filtered, with-gene, and outlier helpers; not used for fig-4c

# Build datasets separately for ref and novel
reads_dt_ref <- build_reads_dt(CLASS_FILES_REF)
reads_dt_novel <- build_reads_dt(CLASS_FILES_NOVEL)

## Support summary TSV generation removed

if (nrow(reads_dt_ref) == 0 && nrow(reads_dt_novel) == 0) stop("No per-read records for either ref or novel; cannot plot.")

prep_plot_dt <- function(reads_dt) {
  plot_dt <- reads_dt[!is.na(f) & !is.na(label) & !is.na(species)]
  plot_dt[, Species := ifelse(tolower(species) == "human", "Human", "Mouse")]
  plot_dt[, Group := factor(label, levels = c("TP","PTP"))]
  plot_dt[, PlotGroup := paste0(Species, "_", as.character(Group))]
  plot_dt[, PlotGroup := factor(PlotGroup, levels = c("Human_TP","Mouse_TP","Human_PTP","Mouse_PTP"))]
  plot_dt
}

## Species-only plot helper removed

## Ref/novel-specific datasets removed

custom_colors <- c(
  Human_TP  = "#6BAED6",   # TP blue
  Mouse_TP  = "#6BAED6",   # TP blue
  Human_PTP = "#FD8D3C",   # PTP orange
  Mouse_PTP = "#FD8D3C"     # PTP orange
)

## Outlier policy: keep all data points (no filtering)

p_make <- function(plot_dt, title_text) {
  ggplot(plot_dt, aes(x = Group, y = f, fill = PlotGroup)) +
    geom_violin(trim = TRUE, position = position_dodge(width = 0.9), alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(width = 0.9), fill = "white", alpha = 0.5, linewidth = 0.3) +
    facet_wrap(~Species, scales = "fixed") +
    scale_fill_manual(values = custom_colors, name = "Group",
                      labels = c("Human TP", "Mouse TP", "Human PTP", "Mouse PTP")) +
    coord_cartesian(ylim = c(0, 2.0)) +
    labs(x = "", y = "f (read_length / exonic_length)", title = title_text) +
    scale_x_discrete(labels = c(TP = "TP", PTP = "PTP")) +
    theme_classic(base_size = 7) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

## Species-only plotting function removed

## Ref-only plot removed

## Novel-only plot removed

# ---- Combined Ref+Novel (no TSS filter) ----
reads_dt_all <- rbindlist(list(reads_dt_ref, reads_dt_novel), use.names = TRUE, fill = TRUE)
plot_dt_combined <- if (nrow(reads_dt_all)) prep_plot_dt(reads_dt_all) else data.table()
if (nrow(plot_dt_combined)) {
  p_combined <- p_make(plot_dt_combined, "Ref+Novel combined (no TSS filter): per-read f for TUSCO TP/PTP")
  message("Saving plot to ", OUTPUT_PDF_COMBINED)
  ggsave(filename = OUTPUT_PDF_COMBINED, plot = p_combined, device = "pdf", width = 4, height = 1.8, units = "in", dpi = 300)
  message("Wrote: ", OUTPUT_PDF_COMBINED)
  # No TSV/stat outputs for fig-4c
} else {
  message("No data for combined plot; skipping fig-4c.")
}

## TSS-filtered and no-outliers sections removed
