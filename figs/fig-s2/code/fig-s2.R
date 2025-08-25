###########################
### 0. Define Variables ###
###########################

# Input data
alphagenome_human_json <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/alphagenome/human_scores.json"
alphagenome_mouse_json <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/alphagenome/mouse_scores.json"

tusco_human_tsv <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/tusco/tusco_human.tsv"
tusco_mouse_tsv <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/tusco/tusco_mouse.tsv"

# Output
output_dir        <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/fig-s2/plots"
output_file_log   <- file.path(output_dir, "fig_s2_log.pdf")
output_file_linear<- file.path(output_dir, "fig_s2_linear.pdf")

# Visual style constants (match Fig. 1d)
LOG_OFFSET <- 1e-4

##################
### 1. Libraries ###
##################

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

############################################
### 2. Helpers: read TUSCO lists and JSON ###
############################################

read_tusco_gene_ids <- function(tsv_file) {
  # Robustly read a simple TSV with comment lines starting with '#'
  # Expect first column to be Ensembl Gene ID
  raw_lines <- readLines(tsv_file)
  raw_lines <- raw_lines[!startsWith(raw_lines, "#")]
  if (length(raw_lines) == 0) return(character(0))
  dt <- fread(text = raw_lines, sep = "\t", header = FALSE)
  gene_ids <- dt[[1]]
  gene_ids <- gsub("\\..*$", "", gene_ids) # drop version suffix
  unique(gene_ids)
}

compute_scores_from_json <- function(json_path) {
  # Robust per-record parsing to avoid column recycling/misalignment
  records <- fromJSON(json_path, simplifyVector = FALSE)
  n <- length(records)
  gene_id <- character(n)
  expression_score <- numeric(n)
  splicing_score <- numeric(n)

  for (i in seq_len(n)) {
    rec <- records[[i]]
    gid <- rec$gene_id %||% NA_character_
    gene_id[i] <- gsub("\\..*$", "", gid)

    # Expression: median across tissues
    expr_vals <- tryCatch(as.numeric(unlist(rec$tissue_expression, use.names = FALSE)), error = function(e) numeric(0))
    expression_score[i] <- if (length(expr_vals) > 0) median(expr_vals, na.rm = TRUE) else NA_real_

    # Splicing: median across tissues (may be empty for single-exon genes)
    if (!is.null(rec$tissue_splice_ratios) && length(rec$tissue_splice_ratios) > 0) {
      splice_vals <- tryCatch(as.numeric(unlist(rec$tissue_splice_ratios, use.names = FALSE)), error = function(e) numeric(0))
      splicing_score[i] <- if (length(splice_vals) > 0) suppressWarnings(median(splice_vals, na.rm = TRUE)) else NA_real_
    } else {
      splicing_score[i] <- NA_real_
    }
  }

  data.table(
    gene_id = gene_id,
    expression_score = expression_score,
    splicing_score = splicing_score
  )
}

#############################
### 3. Load and Prepare Data ###
#############################

cat("Reading TUSCO gene lists...\n")
tusco_human <- read_tusco_gene_ids(tusco_human_tsv)
tusco_mouse <- read_tusco_gene_ids(tusco_mouse_tsv)
cat("  Human TUSCO genes:", length(tusco_human), "\n")
cat("  Mouse TUSCO genes:", length(tusco_mouse), "\n\n")

cat("Reading AlphaGenome JSONs...\n")
human_scores <- compute_scores_from_json(alphagenome_human_json)
human_scores[, Species := "Human"]
mouse_scores <- compute_scores_from_json(alphagenome_mouse_json)
mouse_scores[, Species := "Mouse"]

cat("  Human genes in JSON:", nrow(human_scores), "\n")
cat("  Mouse genes in JSON:", nrow(mouse_scores), "\n\n")

# Mark groups: TUSCO vs Other (treat JSON set as single-isoform cohort)
human_scores[, Group := ifelse(gene_id %in% tusco_human, "TUSCO", "Other")]
mouse_scores[, Group := ifelse(gene_id %in% tusco_mouse, "TUSCO", "Other")]

# Combine and pivot to long for Expression and Splicing
combined <- rbindlist(list(human_scores, mouse_scores), use.names = TRUE, fill = TRUE)

long_dt <- melt(
  combined,
  id.vars = c("gene_id", "Species", "Group"),
  measure.vars = list(
    Expression = "expression_score",
    Splicing   = "splicing_score"
  ),
  variable.name = "ScoreType",
  value.name = "Score"
)

# The above "measure.vars" with list ensures two rows per gene; if it fails due to data.table version,
# fallback to manual rbind:
if (!"Score" %in% names(long_dt)) {
  long_dt <- rbindlist(list(
    combined[, .(gene_id, Species, Group, ScoreType = "Expression", Score = expression_score)],
    combined[, .(gene_id, Species, Group, ScoreType = "Splicing",   Score = splicing_score)]
  ), use.names = TRUE, fill = TRUE)
}

# Remove NA scores for plotting, but keep track of counts
cat("Non-NA counts by type:\n")
print(na.omit(long_dt)[, .N, by = .(Species, Group, ScoreType)])
long_dt <- long_dt[!is.na(Score)]

# Style: create PlotGroup for fill mapping
long_dt[, PlotGroup := paste0(Species, "_", Group)]

############################################
### 4. Plot (violin + boxplot; Fig 1d style) ###
############################################

custom_colors <- c(
  "Human_TUSCO" = "#A8D5A0",  # Light Green
  "Mouse_TUSCO" = "#1b9e77",  # Dark Green
  "Human_Other" = "#D1D3D4",  # Light Grey
  "Mouse_Other" = "#808080"   # Dark Grey
)

long_dt[, PlotGroup := factor(PlotGroup, levels = names(custom_colors))]

# Remove outliers: keep only 1st–99th percentile per Species × ScoreType
facet_limits <- long_dt[, .(
  qmin = quantile(Score, 0.01, na.rm = TRUE),
  qmax = quantile(Score, 0.99, na.rm = TRUE)
), by = .(Species, ScoreType)]

filtered_dt <- merge(long_dt, facet_limits, by = c("Species", "ScoreType"), all.x = TRUE)
filtered_dt <- filtered_dt[Score >= qmin & Score <= qmax]
filtered_dt[, c("qmin", "qmax") := NULL]

# Log-scale plot uses full data (no outlier removal)
plot_dt_log <- copy(long_dt)
plot_dt_log[, log_score := log10(Score + LOG_OFFSET)]

p_log <- ggplot(plot_dt_log, aes(x = Group, y = log_score, fill = PlotGroup)) +
  geom_violin(trim = TRUE, position = position_dodge(width = 0.9), alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(width = 0.9), fill = "white", alpha = 0.5, linewidth = 0.3,
               aes(group = interaction(Species, Group))) +
  facet_grid(ScoreType ~ Species, scales = "fixed") +
  scale_fill_manual(values = custom_colors, name = "Group",
                    labels = c("Human TUSCO", "Mouse TUSCO", "Human Other", "Mouse Other")) +
  scale_x_discrete(labels = c("TUSCO" = "TUSCO", "Other" = "Other")) +
  labs(x = "", y = bquote(log[10]~"(Score +"~.(LOG_OFFSET)~")"), title = "AlphaGenome Scores: TUSCO vs. Other Single-Isoform Genes") +
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

# Linear-scale version (force Splicing y-scale to 0–5)
plot_dt_linear <- copy(filtered_dt)
plot_dt_linear[ScoreType == "Splicing", Score := pmin(pmax(Score, 0), 5)]

# Helper blank data to clamp Splicing facets exactly to [0, 5]
splicing_bounds <- unique(plot_dt_linear[ScoreType == "Splicing", .(Species, ScoreType)])
splicing_bounds <- splicing_bounds[, .(bound = c(0, 5)), by = .(Species, ScoreType)]

p_linear <- ggplot(plot_dt_linear, aes(x = Group, y = Score, fill = PlotGroup)) +
  geom_violin(trim = TRUE, position = position_dodge(width = 0.9), alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(width = 0.9), fill = "white", alpha = 0.5, linewidth = 0.3,
               aes(group = interaction(Species, Group))) +
  # Ensure Splicing facets use 0–5 exactly
  geom_blank(data = splicing_bounds, inherit.aes = FALSE, aes(y = bound)) +
  facet_grid(ScoreType ~ Species, scales = "free_y") +
  scale_fill_manual(values = custom_colors, name = "Group",
                    labels = c("Human TUSCO", "Mouse TUSCO", "Human Other", "Mouse Other")) +
  scale_x_discrete(labels = c("TUSCO" = "TUSCO", "Other" = "Other")) +
  labs(x = "", y = "Score", title = "AlphaGenome Scores: TUSCO vs. Other Single-Isoform Genes") +
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

############################
### 5. Save figure to PDF ###
############################

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Write per-TUSCO gene scores to a log TSV (median expression, median splicing)
tusco_scores_log <- combined[Group == "TUSCO", .(gene_id, Species, expression_score, splicing_score)]
setcolorder(tusco_scores_log, c("gene_id", "Species", "expression_score", "splicing_score"))
fwrite(tusco_scores_log, file = file.path(output_dir, "fig_s2_TUSCO_scores.tsv"), sep = "\t")

ggsave(filename = output_file_log,    plot = p_log,    device = "pdf", width = 4.0, height = 2.6, units = "in", dpi = 300)
ggsave(filename = output_file_linear, plot = p_linear, device = "pdf", width = 4.0, height = 2.6, units = "in", dpi = 300)

cat("Saved log-scale figure:", output_file_log, "\n")
cat("Saved linear-scale figure:", output_file_linear, "\n")
cat("Saved TUSCO scores log:", file.path(output_dir, "fig_s2_TUSCO_scores.tsv"), "\n")
