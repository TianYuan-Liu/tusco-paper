# TUSCO Figure Utilities Library
# Common functions and patterns for figure generation scripts

#' Safe library loading with graceful error handling
#' @param packages Character vector of package names to load
load_required_packages <- function(packages) {
  suppressPackageStartupMessages({
    for (pkg in packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop("Required package '", pkg, "' is not installed. Please install it with: install.packages('", pkg, "')")
      }
    }
  })
}

#' Resolve file paths from multiple candidate locations
#' @param candidates Character vector of candidate file paths
#' @param is_dir Logical, whether to check for directories instead of files
#' @return First existing path, or first candidate if none exist
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (is.na(p) || is.null(p)) next
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
  }
  return(candidates[[1]])
}

#' Resolve data file paths from standard data locations
#' @param ... Path components to combine
#' @param bases Character vector of base directories to search
#' @return Resolved file path
resolve_data_path <- function(..., bases = c("data/raw", "data/processed", "../data", "../../data", "../../../data")) {
  sub_path <- file.path(...)
  for (base in bases) {
    candidate <- file.path(base, sub_path)
    if (file.exists(candidate)) return(normalizePath(candidate))
  }
  stop("Data file not found: ", sub_path, "\nSearched in: ", paste(bases, collapse = ", "))
}

#' Create standard output directories
#' @param base_dir Base directory for outputs (default: current figure directory)
#' @return List with plot_dir and table_dir paths
create_output_dirs <- function(base_dir = ".") {
  plot_dir <- file.path(base_dir, "plots")
  table_dir <- file.path(base_dir, "tables")
  
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  list(plot_dir = plot_dir, table_dir = table_dir)
}

#' Save figure with standardized naming and optional data export
#' @param plot_obj ggplot object or similar plot object
#' @param fig_id Figure identifier (e.g., "fig-3a", "fig-s2")
#' @param plot_dir Directory for plot outputs
#' @param table_dir Directory for table outputs
#' @param width Plot width in inches (default: Nature single column)
#' @param height Plot height in inches
#' @param data Optional data frame to save as TSV
#' @param data_suffix Optional suffix for data file name
save_figure <- function(plot_obj, fig_id, plot_dir = "plots", table_dir = "tables", 
                       width = 3.35, height = 4.0, data = NULL, data_suffix = NULL) {
  
  # Ensure output directories exist
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save plot
  plot_file <- file.path(plot_dir, paste0(fig_id, ".pdf"))
  ggsave(plot_file, plot_obj, width = width, height = height, 
         units = "in", device = "pdf", dpi = 300)
  message("Saved plot: ", plot_file)
  
  # Save data if provided
  if (!is.null(data)) {
    data_name <- if (is.null(data_suffix)) fig_id else paste0(fig_id, "-", data_suffix)
    data_file <- file.path(table_dir, paste0(data_name, ".tsv"))
    readr::write_tsv(data, data_file)
    message("Saved data: ", data_file)
  }
  
  invisible(plot_file)
}

#' Parse command line arguments with defaults
#' @param defaults Named list of default values
#' @return List of parsed arguments
parse_figure_args <- function(defaults = list(
  out_dir = "..",
  width = 3.35,  # Nature single column width in inches
  height = 4.0
)) {
  args <- commandArgs(trailingOnly = TRUE)
  
  result <- defaults
  if (length(args) > 0) result$out_dir <- args[1]
  if (length(args) > 1) result$width <- as.numeric(args[2])
  if (length(args) > 2) result$height <- as.numeric(args[3])
  
  return(result)
}

#' Standard TUSCO theme for consistent plot styling
#' @param base_size Base font size
#' @param base_family Font family
theme_tusco <- function(base_size = 7, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      legend.key.size = unit(0.5, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

#' Read TSV files safely with error handling
#' @param file_path Path to TSV file
#' @param ... Additional arguments passed to readr::read_tsv
read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  tryCatch(
    readr::read_tsv(file_path, show_col_types = FALSE, ...),
    error = function(e) {
      stop("Error reading file ", file_path, ": ", e$message)
    }
  )
}

#' Get figure directory from script path
#' Determines the figure directory based on script location
get_figure_dir <- function() {
  # Try to get script path from command line args
  argv <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", argv[grep("^--file=", argv)])
  
  if (length(script_path) == 1 && nzchar(script_path)) {
    # Script path available - figure dir is parent of code dir
    fig_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
  } else {
    # Fallback - assume we're in a code directory
    if (basename(getwd()) == "code") {
      fig_dir <- normalizePath("..", mustWork = FALSE)
    } else {
      fig_dir <- getwd()
    }
  }
  
  return(fig_dir)
}

# Export commonly used color palettes
TUSCO_COLORS <- list(
  human_tusco = "#a8d5a0",
  mouse_tusco = "#1b9e77", 
  human_gencode = "#fdbf6f",
  mouse_gencode = "#e66101",
  human_mane = "#e41a1c",
  sirvs = "#cab2d6",
  erccs = "#6a3d9a"
)

# Export standard dimensions
NATURE_DIMS <- list(
  single_column = 3.35,    # inches
  double_column = 7.09,    # inches
  full_page_width = 7.87,  # inches
  default_height = 4.0     # inches
)