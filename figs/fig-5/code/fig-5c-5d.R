suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(grid)
})
utils::globalVariables(c(
  "ensembl", "transcript", "refseq", "prot_refseq", "gene_name",
  "structural_category", "associated_gene", "associated_transcript",
  "isoform", "id_type", "subcategory", "ref_exons", "diff_to_TSS",
  "diff_to_TTS", "gene"
))

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
data_root <- file.path(project_root, "figs/data")

# Output panel path (new name)
plot_dir <- file.path(project_root, "figs/fig-5/plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
panel_pdf <- file.path(plot_dir, "fig-5c-5d.pdf")

# ------------------------------------------------------------
# Figure 5c: TP/PTP/FP/FN vs TUSCO sets (Iso-Seq AR)
# ------------------------------------------------------------

# TUSCO references
tusco_universal_tsv <- file.path(data_root, "tusco/tusco_mouse.tsv")
tusco_brain_tsv     <- file.path(data_root, "tusco/tusco_mouse_brain.tsv")
tusco_kidney_tsv    <- file.path(data_root, "tusco/tusco_mouse_kidney.tsv")

# Iso-Seq AR input root
isoseq_ar_root <- file.path(data_root, "nih/single_sample")

read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  tryCatch(read_tsv(file_path, show_col_types = FALSE, ...),
           error = function(e) stop("Unable to read ", file_path, ": ", e$message))
}

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

compute_tusco_bigcats <- function(classification_data, tusco_ref_df) {
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Za-z0-9][A-Za-z0-9_-]*[A-Za-z0-9]?$"
  )

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

  top_id_type <- class_clean %>% count(id_type, sort = TRUE) %>% slice_head(n = 1) %>% pull(id_type)
  if (length(top_id_type) == 0 || is.na(top_id_type)) top_id_type <- "gene_name"

  class_tusco <- class_clean %>% filter(associated_gene %in% tusco_ref_df[[top_id_type]])

  TP_df <- class_tusco %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" &
           ref_exons == 1 &
           abs(diff_to_TSS) < 50 &
           abs(diff_to_TTS) < 50)
    )
  PTP_df <- class_tusco %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match"),
           !(associated_gene %in% TP_df$associated_gene))

  found_ids <- union(TP_df$associated_gene, PTP_df$associated_gene)
  ref_ids <- unique(tusco_ref_df[[top_id_type]])
  FN_count <- sum(!(ref_ids %in% found_ids))

  tibble(
    TP = nrow(TP_df),
    PTP = nrow(PTP_df),
    FP = class_tusco %>% filter(structural_category %in% c(
      "novel_in_catalog", "novel_not_in_catalog", "genic",
      "fusion", "antisense", "intergenic", "genic_intron"
    )) %>% nrow(),
    FN = FN_count,
    id_type = top_id_type
  )
}

infer_tissue <- function(sample_basename) {
  if (startsWith(sample_basename, "B")) return("Brain")
  if (startsWith(sample_basename, "K")) return("Kidney")
  return("Unknown")
}

all_sample_dirs <- list.dirs(isoseq_ar_root, full.names = TRUE, recursive = FALSE)
all_sample_dirs <- all_sample_dirs[grepl("/(B|K)[0-9]+\\.isoforms$", all_sample_dirs)]
if (length(all_sample_dirs) == 0) stop("No sample directories found in ", isoseq_ar_root)

tusco_ref_univ   <- load_tusco_reference(tusco_universal_tsv)
tusco_ref_brain  <- load_tusco_reference(tusco_brain_tsv)
tusco_ref_kidney <- load_tusco_reference(tusco_kidney_tsv)

results <- list()
for (sd in all_sample_dirs) {
  sample_name <- basename(sd)
  class_file  <- file.path(sd, paste0(sample_name, "_classification.txt"))
  if (!file.exists(class_file)) next

  classification <- read_tsv_safe(class_file)
  sample_tissue <- infer_tissue(sample_name)

  m_univ <- compute_tusco_bigcats(classification, tusco_ref_univ)
  results[[length(results) + 1]] <- tibble(
    Sample = sample_name,
    Tissue = sample_tissue,
    Type   = "Universal",
    TP = m_univ$TP, PTP = m_univ$PTP, FP = m_univ$FP, FN = m_univ$FN,
    id_type = m_univ$id_type
  )

  if (sample_tissue == "Brain") {
    m_brain <- compute_tusco_bigcats(classification, tusco_ref_brain)
    results[[length(results) + 1]] <- tibble(
      Sample = sample_name,
      Tissue = sample_tissue,
      Type   = "Brain",
      TP = m_brain$TP, PTP = m_brain$PTP, FP = m_brain$FP, FN = m_brain$FN,
      id_type = m_brain$id_type
    )
  } else if (sample_tissue == "Kidney") {
    m_kidney <- compute_tusco_bigcats(classification, tusco_ref_kidney)
    results[[length(results) + 1]] <- tibble(
      Sample = sample_name,
      Tissue = sample_tissue,
      Type   = "Kidney",
      TP = m_kidney$TP, PTP = m_kidney$PTP, FP = m_kidney$FP, FN = m_kidney$FN,
      id_type = m_kidney$id_type
    )
  }
}

results_df <- bind_rows(results)

long_df <- results_df %>%
  pivot_longer(cols = c(TP, PTP, FP, FN), names_to = "Metric", values_to = "Count") %>%
  group_by(Sample, Type) %>%
  mutate(Total = sum(Count), Percentage = 100 * Count / ifelse(Total == 0, 1, Total)) %>%
  ungroup()

summary_df <- long_df %>%
  group_by(Tissue, Type, Metric) %>%
  summarize(
    mean_perc = mean(Percentage, na.rm = TRUE),
    sd_perc = sd(Percentage, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

summary_df$Metric <- factor(summary_df$Metric, levels = c("TP", "PTP", "FP", "FN"))
summary_df$Type <- factor(summary_df$Type, levels = c("Universal", "Brain", "Kidney"))
long_df$Type <- factor(long_df$Type, levels = c("Universal", "Brain", "Kidney"))

position_d <- position_dodge(width = 0.8)

p_c <- ggplot(summary_df, aes(x = Metric, y = mean_perc, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.2, position = position_d) +
  geom_errorbar(aes(ymin = mean_perc - sd_perc, ymax = mean_perc + sd_perc),
                width = 0.2, linewidth = 0.2, position = position_d) +
  geom_point(data = long_df, aes(x = Metric, y = Percentage, group = Type),
             position = position_d, size = 0.3, alpha = 0.5, inherit.aes = FALSE) +
  facet_grid(. ~ Tissue, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 100), labels = scales::percent_format(scale = 1), breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = c(
    "Universal" = "#1C9E77",
    "Brain"     = "#5893a4",
    "Kidney"    = "#506f68"
  )) +
  labs(x = NULL, y = "Percentage (%)",
       title = "TP/PTP/FP/FN vs TUSCO universal and tissue sets (Iso-Seq AR)") +
  theme_classic(base_size = 7) +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 7, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    strip.text = element_text(size = 7, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# ------------------------------------------------------------
# Figure 5d: Length distributions (TP transcripts, FN transcripts, PacBio reads)
# ------------------------------------------------------------

isoseq_root <- file.path(project_root, "figs/data/nih/single_sample")
fn_stats_tsv <- file.path(project_root, "figs/fig-5/plot/fig-5d_kidney_missing_gene_stats.tsv")
read_lengths_tsv <- file.path(project_root, "figs/data/nih/nih_all_read_lengths.tsv")
tusco_kidney_tsv <- file.path(project_root, "figs/data/tusco/tusco_mouse_kidney.tsv")
tp_gtf_path <- file.path(project_root, "figs/data/tusco/tusco_mouse_kidney.gtf")

load_tusco_reference_simple <- function(tusco_path) {
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
      refseq = stringr::str_remove(refseq, "\\.\\d+$")
    ) %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()
}

tusco_ref <- load_tusco_reference_simple(tusco_kidney_tsv)
refseq_to_ensembl <- setNames(
  tusco_ref$ensembl[!is.na(tusco_ref$refseq) & tusco_ref$refseq != ""],
  tusco_ref$refseq[!is.na(tusco_ref$refseq) & tusco_ref$refseq != ""]
)
gene_to_ensembl <- setNames(
  tusco_ref$ensembl[!is.na(tusco_ref$gene_name) & tusco_ref$gene_name != ""],
  tusco_ref$gene_name[!is.na(tusco_ref$gene_name) & tusco_ref$gene_name != ""]
)

get_tp_genes_for_sample <- function(sample_dir) {
  cl <- file.path(sample_dir, paste0(basename(sample_dir), "_classification.txt"))
  if (!file.exists(cl)) return(character())
  df <- read_tsv(cl, show_col_types = FALSE, progress = FALSE)
  df %>%
    filter(
      subcategory == "reference_match" |
      (subcategory == "mono-exon" & ref_exons == 1 & abs(diff_to_TSS) < 50 & abs(diff_to_TTS) < 50)
    ) %>%
    mutate(id_norm = stringr::str_remove(associated_gene, "\\.\\d+$")) %>%
    mutate(gene = dplyr::case_when(
      grepl("^(ENSMUSG|ENSG)\\d{11}$", id_norm) ~ id_norm,
      !is.na(refseq_to_ensembl[id_norm]) ~ unname(refseq_to_ensembl[id_norm]),
      !is.na(gene_to_ensembl[id_norm]) ~ unname(gene_to_ensembl[id_norm]),
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(gene) & gene != "") %>%
    distinct(gene) %>%
    pull(gene)
}

sample_dirs <- list.dirs(isoseq_root, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("/K[0-9]+\\.isoforms$", sample_dirs)]
tp_sets <- lapply(sample_dirs, get_tp_genes_for_sample)
always_tp_genes <- if (length(tp_sets) > 0) Reduce(intersect, tp_sets) else character()

# Load kidney TUSCO GTF and compute TP transcript exonic lengths
gtf <- readr::read_tsv(
  tp_gtf_path,
  comment = "#",
  col_names = c("chrom","source","feature","start","end","score","strand","frame","attribute"),
  col_types = cols(
    chrom = col_character(), source = col_character(), feature = col_character(),
    start = col_integer(), end = col_integer(), score = col_character(),
    strand = col_character(), frame = col_character(), attribute = col_character()
  )
)

extract_attr <- function(attr, key) {
  m <- stringr::str_match(attr, paste0(key, " \"([^\"]+)\""))
  out <- m[,2]
  out <- stringr::str_remove(out, "\\.\\d+$")
  out
}

gtf <- gtf %>% mutate(gene_id = extract_attr(attribute, "gene_id"), transcript_id = extract_attr(attribute, "transcript_id"))

tp_exons <- gtf %>% filter(feature == "exon", !is.na(gene_id), gene_id %in% always_tp_genes, !is.na(transcript_id))

tp_tx_stats <- tp_exons %>%
  group_by(gene_id, transcript_id) %>%
  summarise(
    chrom = dplyr::first(chrom),
    strand = dplyr::first(strand),
    n_exons = dplyr::n(),
    exonic_bases = sum(end - start + 1L),
    tx_start = min(start),
    tx_end = max(end),
    tx_span_bases = tx_end - tx_start + 1L,
    .groups = "drop"
  ) %>%
  arrange(gene_id, transcript_id)

TP_df <- tp_tx_stats %>% transmute(type = "TP", length = as.numeric(exonic_bases)) %>% distinct()

# FN transcript lengths table (precomputed)
FN_df <- read_tsv(fn_stats_tsv, show_col_types = FALSE) %>%
  transmute(type = "FN", length = as.numeric(exonic_bases)) %>%
  distinct()

# PacBio read lengths for K31-K35 only
suppressPackageStartupMessages(library(vroom))
reads <- vroom::vroom(read_lengths_tsv, col_types = list(.default = "c"), delim = "\t",
                      col_select = c("read_length","file_name","platform")) %>%
  filter(platform == "PacBio", file_name %in% paste0("K", 31:35)) %>%
  transmute(type = "PacBio reads (K31-K35)", length = as.numeric(read_length))

all_len <- bind_rows(TP_df, FN_df, reads) %>% filter(is.finite(length))

group_counts <- all_len %>% count(type) %>% tibble::deframe()
y_limits <- c("PacBio reads (K31-K35)", "TP", "FN")
y_labels <- setNames(c("Reads\n(K31−K35)", "TP", "FN"), y_limits)

p_d <- ggplot(all_len, aes(y = type, x = length, fill = type)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.8, color = NA, draw_quantiles = c(0.5)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "#222222", fill = "white", alpha = 0.9, linewidth = 0.3) +
  scale_y_discrete(limits = y_limits, labels = y_labels) +
  scale_x_continuous(breaks = c(200, 500, 1000, 2000, 3000, 5000),
                     labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 6000)) +
  scale_fill_manual(values = c(
    "PacBio reads (K31-K35)" = "#E75480",
    "TP" = "#6BAED6",
    "FN" = "#969696"
  )) +
  labs(x = "Length (nt)", y = NULL, fill = NULL) +
  theme_classic(base_size = 7) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7)
  )

# ------------------------------------------------------------
# Panel export (Figure 5c on the left, Figure 5d on the right)
# ------------------------------------------------------------
pdf(panel_pdf, width = 7.09, height = 2)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
print(p_c, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_d + labs(title = NULL), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

message("Wrote panel to ", panel_pdf)


