############################################################
#  Figure S4: Better TUSCO bar-plots (combined single PDF)
#  - Uses LRGASP evaluation data under figs/data/lrgasp/tusco_novel_evl
#  - Uses TUSCO annotations under figs/data/tusco
#  - Outputs a single figure: figs/fig-s4/plot/fig-s4.pdf
############################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(glue)
})

# silence NSE notes for CRAN/lintr
utils::globalVariables(c(
  "%>%","structural_category","associated_gene","subcategory","final_label",
  "eval","sample","n","group","TP","PTP","FP","FN","cat"
))

## ------------------------------------------------------------------
## CONFIGURATION (paths relative to this script's location)
## ------------------------------------------------------------------
# Resolve script directory robustly (works with Rscript and source)
args <- commandArgs(trailingOnly = FALSE)
script_path <- NA_character_
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
}
if (is.na(script_path) || !nzchar(script_path)) {
  script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
}
if (is.na(script_path) || !nzchar(script_path)) {
  # fallback: assume repo cwd when running interactively
  script_path <- file.path(getwd(), "figs", "fig-s4", "code", "fig-s4.R")
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))

# This script is in figs/fig-s4/code → repo root is three levels up
root_dir      <- normalizePath(file.path(script_dir, "..", "..", ".."))
base_data_dir <- file.path(root_dir, "figs", "data", "lrgasp", "tusco_novel_evl")
out_dir       <- file.path(root_dir, "figs", "fig-s4", "plot")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_order <- c(
  "wtc11_captrap_pacbio","wtc11_cdna_pacbio",
  "wtc11_cdna_ont","wtc11_captrap_ont",
  "es_captrap_pacbio","es_cdna_pacbio",
  "es_cdna_ont","es_captrap_ont"
)

pipeline_specs <- tribble(
  ~label,      ~dir,
  "Bambu",     "bambu_sq3",
  "StringTie", "stringtie_sq3",
  "Flair",     "flair_sq3",
  "IsoSeq",    "isoseq_sq3"
)

cat_levels  <- c("RM","Alternative 3'end","Alternative 5'end","Alternative 3'5'end",
                 "ISM","NIC","NNC","Genic Intron","Genic Genomic","Antisense",
                 "Fusion","Intergenic","Missing")
sirv_pal <- c(
  "RM" = '#c4e1f2', "Alternative 3'end"='#02314d',
  "Alternative 5'end"='#7ccdfc',"Alternative 3'5'end"='#0e5a87',
  "ISM"="#FC8D59","NIC"="#78C679","NNC"="#EE6A50",
  "Genic Intron"="#41B6C4","Genic Genomic"="#969696",
  "Antisense"="#66C2A4","Fusion"="goldenrod1",
  "Intergenic"="darksalmon","Missing"="grey70"
)

tp_grp  <- c("RM")
ptp_grp <- c("Alternative 3'end","Alternative 5'end","Alternative 3'5'end","ISM")
fp_grp  <- c("NIC","NNC","Genic Intron","Genic Genomic","Antisense","Fusion","Intergenic")

## ------------------------------------------------------------------
## Figure S1-like theme (match label sizes and legend aesthetics)
## ------------------------------------------------------------------
theme_s4 <- ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0),
    axis.title = ggplot2::element_text(size = 7),
    axis.text = ggplot2::element_text(size = 6),
    axis.line = ggplot2::element_line(linewidth = 0.25),
    axis.ticks = ggplot2::element_line(linewidth = 0.25),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

## ------------------------------------------------------------------
## helper that loads + aggregates one classification file
## ------------------------------------------------------------------
load_one <- function(fp, tusco_ref){
  if(!file.exists(fp)) return(tibble::tibble())
  d <- readr::read_tsv(fp, show_col_types = FALSE, progress = FALSE)
  if(nrow(d)==0) return(tibble::tibble())
  d <- d %>%
        dplyr::filter(structural_category!="fusion") %>%            # unpack fusions
        dplyr::bind_rows(d %>% dplyr::filter(structural_category=="fusion") %>%
                    tidyr::separate_rows(associated_gene, sep="_")) %>%
        dplyr::mutate(associated_gene=stringr::str_remove(associated_gene,"\\.[0-9]+$")) %>%
        dplyr::filter(associated_gene %in% tusco_ref)               # keep only TUSCO genes
  
  cat_map <- function(subcat,struct){
    if(struct=="full-splice_match" && subcat=="reference_match")       "RM"
    else if(struct=="full-splice_match" && subcat=="alternative_3end") "Alternative 3'end"
    else if(struct=="full-splice_match" && subcat=="alternative_5end") "Alternative 5'end"
    else if(struct=="full-splice_match" && subcat=="alternative_3end5end") "Alternative 3'5'end"
    else if(struct=="incomplete-splice_match") "ISM"
    else if(struct=="novel_in_catalog")         "NIC"
    else if(struct=="novel_not_in_catalog")     "NNC"
    else if(struct=="genic_intron")             "Genic Intron"
    else if(struct=="genic")                    "Genic Genomic"
    else if(struct=="antisense")                "Antisense"
    else if(struct=="fusion")                   "Fusion"
    else if(struct=="intergenic")               "Intergenic"
    else NA_character_
  }
  d %>%
    dplyr::mutate(final_label = mapply(cat_map, subcategory, structural_category)) %>%
    dplyr::filter(!is.na(final_label))
}

## ------------------------------------------------------------------
## Main loop – build plots for each pipeline (Ref and Novel panels)
## ------------------------------------------------------------------
plots_for_grid <- list()
for(i in seq_len(nrow(pipeline_specs))){
  pinfo <- pipeline_specs[i,]
  message("Processing ", pinfo$label)
  
  # determine which sample prefixes really exist for this pipeline
  existing_samples <- sample_order[ sapply(sample_order, function(sp){
      ref_dir   <- file.path(base_data_dir, pinfo$dir, "ref_evl",  sp)
      novel_dir <- file.path(base_data_dir, pinfo$dir, "novel_evl", sp)
      dir.exists(ref_dir) || dir.exists(novel_dir)
  }) ]

  if(length(existing_samples)==0){
    warning("No sample directories found for ", pinfo$label)
    next
  }

  df_all <- purrr::map_dfr(existing_samples, function(sprefix){
    species <- ifelse(str_detect(sprefix,"^wtc11"),"human","mouse")
    tusco_file <- file.path(root_dir, "figs", "data", "tusco",
                            if(species=="human") "tusco_human_multi_exon.tsv" else "tusco_mouse_multi_exon.tsv")
    tusco_df <- readr::read_tsv(tusco_file, show_col_types = FALSE, progress = FALSE,
                         col_types = readr::cols(.default="c"),
                         comment = "#", col_names = FALSE) %>%
                setNames(c("ensembl","transcript","gene_name","entrez","refseq","protein"))
    tusco_ref <- tusco_df %>%
                 dplyr::select(ensembl, refseq, gene_name) %>%
                 unlist(use.names = FALSE) %>%
                 unique()
    
    # build expected paths
    dbase <- file.path(base_data_dir, pinfo$dir)
    sqanti_dir <- "sqanti3_out"
    ref_fp   <- file.path(dbase,"ref_evl",  sprefix, sqanti_dir, "isoforms_classification.txt")
    novel_fp <- file.path(dbase,"novel_evl",sprefix, sqanti_dir, "isoforms_classification.txt")
    if(pinfo$dir=="stringtie_sq3"){
      ref_fp   <- file.path(dbase,"ref_evl",  sprefix, sqanti_dir,
                            glue("{sprefix}_stringtie_ref_sqanti_classification.txt"))
      novel_fp <- file.path(dbase,"novel_evl",sprefix, sqanti_dir,
                            glue("{sprefix}_stringtie_tusco_sqanti_classification.txt"))
    }
    if(pinfo$dir %in% c("bambu_sq3","flair_sq3")){
      ref_fp   <- file.path(dbase,"ref_evl",  sprefix, sqanti_dir,
                            glue("{sprefix}_sqanti_classification.txt"))
      novel_fp <- file.path(dbase,"novel_evl",sprefix, sqanti_dir,
                            glue("{sprefix}_sqanti_classification.txt"))
    }
    
    # load both evaluation types
    ref_set   <- load_one(ref_fp,   tusco_ref) %>% dplyr::mutate(eval="Ref Evl")
    novel_set <- load_one(novel_fp, tusco_ref) %>% dplyr::mutate(eval="Novel Evl")

    comb <- bind_rows(ref_set, novel_set)

    # decide primary ID type based on what is present in data (Ensembl > RefSeq > gene_name)
    patterns <- list(ensembl="^(ENSG|ENSMUSG)", refseq="^(NM_|NR_|NP_)")
    top_id <- if(any(stringr::str_detect(comb$associated_gene, patterns$ensembl))) "ensembl" 
              else if(any(stringr::str_detect(comb$associated_gene, patterns$refseq))) "refseq" 
              else "gene_name"

    tusco_ids <- tusco_df[[top_id]] %>% stats::na.omit() %>% unique()

    # add Missing rows per evaluation
    if("associated_gene" %in% colnames(comb)){
      miss_rows <- comb %>%
                   dplyr::group_by(eval) %>%
                   dplyr::summarise(detected = dplyr::n_distinct(associated_gene), .groups="drop") %>%
                   dplyr::mutate(n_missing = pmax(length(tusco_ids) - detected, 0))
    } else {
      miss_rows <- tibble::tibble(eval = c("Ref Evl","Novel Evl"), detected = 0, n_missing = length(tusco_ids))
    }
    miss_rows <- miss_rows %>%
                   dplyr::mutate(final_label = "Missing", sample = sprefix) %>%
                   dplyr::select(sample, eval, final_label, n = n_missing)

    comb <- comb %>% dplyr::mutate(sample = sprefix, n = 1)

    # return full set rows; Missing rows will be counted later
    bind_rows(comb, miss_rows)
  })
  
  if(nrow(df_all)==0){
    warning("No data for ", pinfo$label)
    next
  }
  
  ## ----------------------------------------------------------------
  ##  If IsoSeq, print FP details for Ref/Novel to console
  ## ----------------------------------------------------------------
  if (identical(as.character(pinfo$dir), "isoseq_sq3")) {
    iso_fp <- df_all %>%
      dplyr::filter(eval %in% c("Ref Evl", "Novel Evl"),
                    final_label %in% fp_grp,
                    !is.na(associated_gene)) %>%
      dplyr::mutate(sample = as.character(sample),
                    associated_gene = as.character(associated_gene)) %>%
      dplyr::group_by(eval, sample, associated_gene, final_label) %>%
      dplyr::summarise(n = sum(n), .groups = "drop") %>%
      dplyr::arrange(eval, sample, dplyr::desc(n), associated_gene)
    message("IsoSeq – false positives (eval, sample, gene, label, n):")
    if (nrow(iso_fp) > 0) {
      print(iso_fp, n = nrow(iso_fp))
    } else {
      message("  None")
    }
  }
  
  ## summarise → stacked barplot frame + TP/PTP/FP/FN frame
  df_cnt <- df_all %>%
              dplyr::group_by(sample, eval, final_label) %>%
              dplyr::summarise(n = sum(n), .groups="drop") %>%
              dplyr::mutate(final_label=factor(final_label, levels=cat_levels))
  
  df_stacked <- df_cnt
  
  df_stats <- df_cnt %>%
      dplyr::mutate(group = dplyr::case_when(
        final_label %in% tp_grp  ~ "TP",
        final_label %in% ptp_grp ~ "PTP",
        final_label %in% fp_grp  ~ "FP",
        final_label == "Missing" ~ "FN"
      )) %>%
      dplyr::group_by(sample, eval, group) %>%
      dplyr::summarise(n = sum(n), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)
  
  # ensure all expected columns exist even if absent in data
  for(col in c("TP","PTP","FP","FN")){
    if(!(col %in% colnames(df_stats))) df_stats[[col]] <- 0
  }
  
  # pretty x-axis labels (three lines: species, prep, platform)
  format_sample_label <- function(s){
    parts <- strsplit(s, "_")[[1]]
    if(length(parts) == 3){
      species  <- parts[1]
      prep     <- parts[2]
      platform <- parts[3]
      sp  <- if (species == "wtc11") "WTC11" else toupper(species)  # ES
      pr  <- dplyr::case_when(prep == "cdna" ~ "cDNA",
                               prep == "captrap" ~ "CapTrap",
                               TRUE ~ prep)
      pl  <- dplyr::case_when(platform == "pacbio" ~ "PB",
                               platform == "pb" ~ "PB",
                               platform == "ont" ~ "ONT",
                               TRUE ~ toupper(platform))
      return(paste(sp, pr, pl, sep = "\n"))
    }
    # fallback mapping
    s <- stringr::str_replace_all(s, c("cdna" = "cDNA", "captrap" = "CapTrap", "pacbio" = "PB", "pb" = "PB", "ont" = "ONT"))
    gsub("_", "\n", s)
  }
  sample_labels <- setNames(vapply(sample_order, format_sample_label, character(1)), sample_order)

  build_one_plot <- function(eval_level){
    dst <- df_stacked %>% filter(eval==eval_level)
    if(nrow(dst)==0) return(NULL)
    ggplot2::ggplot(dst,
           ggplot2::aes(x=factor(sample, levels = sample_order),
               y=n, fill=final_label)) +
      ggplot2::geom_col(width=.7, colour="black", size=.2) +
      ggplot2::scale_fill_manual(values = sirv_pal, name = "SQANTI category") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.25))) +
      ggplot2::scale_x_discrete(labels = sample_labels) +
      theme_s4 +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.85),
                     legend.position = "none",
                     plot.margin = ggplot2::margin(5,5,18,5),
                     plot.caption = ggplot2::element_text(face = "bold", hjust = .5, size = 7)) +
      ggplot2::labs(x=NULL, y="Count",
           caption = glue::glue("{pinfo$label} – {ifelse(eval_level=='Ref Evl','Reference','Novel')}") ) +
      ggplot2::geom_text(data = df_stats %>% dplyr::filter(eval==eval_level),
                ggplot2::aes(x = sample, y = TP + PTP + FP + FN + 15,
                    label = glue::glue("TP={TP}\nPTP={PTP}\nFP={FP}\nFN={FN}")),
                size = 2.4, fontface = "bold", lineheight = .9,
                hjust = .5, inherit.aes = FALSE) +
      ggplot2::coord_cartesian(clip = "off")
  }

  plot_ref   <- build_one_plot("Ref Evl")
  plot_novel <- build_one_plot("Novel Evl")

  for(pp in list(plot_novel, plot_ref)){
    if(!is.null(pp)){
       plots_for_grid[[length(plots_for_grid)+1]] <- pp
    }
  }
}

# ------------------------------------------------------------------
#  Build composite figure with legend (single output PDF)
# ------------------------------------------------------------------
if(length(plots_for_grid)>0){
  combined_grid <- plot_grid(plotlist = plots_for_grid,
                             ncol = 2,
                             labels = letters[seq_along(plots_for_grid)],
                             label_size = 7,
                             label_fontface = "bold")

  legend_df <- tibble(cat = factor(names(sirv_pal), levels = names(sirv_pal)),
                      n   = 1)
  legend_plot <- ggplot2::ggplot(legend_df,
                        ggplot2::aes(x = cat, y = n, fill = cat)) +
                 ggplot2::geom_col(show.legend = TRUE) +
                 ggplot2::scale_fill_manual(values = sirv_pal, name = NULL) +
                 ggplot2::theme_void(base_size = 9) +
                 ggplot2::theme(legend.position   = "bottom",
                       legend.direction  = "horizontal",
                       legend.box        = "horizontal",
                       legend.key.size   = grid::unit(0.4, "lines"),
                       legend.text       = ggplot2::element_text(size = 6),
                       legend.title      = ggplot2::element_blank(),
                       legend.margin     = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
                       legend.box.background = ggplot2::element_rect(color = "black", linewidth = 0.5)) +
                 ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))
  legend_grob <- cowplot::get_legend(legend_plot)

  combined <- plot_grid(combined_grid, legend_grob, ncol = 1, rel_heights = c(1,0.08))

  out_pdf <- file.path(out_dir, "fig-s4.pdf")
  nrow_comb <- ceiling(length(plots_for_grid)/2)
  row_height <- 90
  ggsave(out_pdf, combined, width = 183, height = min(247, nrow_comb*row_height + 25), units = "mm")
  message("  → saved ", out_pdf)
} else {
  warning("No plots generated; check data paths: ", base_data_dir)
}


