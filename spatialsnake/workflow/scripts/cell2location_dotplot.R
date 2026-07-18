#!/usr/bin/env Rscript

# Reproduce the cell2location regional-abundance dotplot from the source TSV
# written by cell2locate_visualize.py.
# Example:
# Rscript cell2location_dotplot.R -i cell2location_relative_abundance_dotplot_source.tsv \
#   -o cell2location_relative_abundance_dotplot_R.pdf

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option(
    c("-i", "--input_tsv"),
    type = "character",
    help = "cell2location_relative_abundance_dotplot_source.tsv"
  ),
  make_option(
    c("-o", "--output_pdf"),
    type = "character",
    help = "Output PDF path"
  ),
  make_option(
    c("--max_cell_types"),
    type = "integer",
    default = 30L,
    help = "Top cell types by tissue mean; 0 displays all [default %default]"
  ),
  make_option(
    c("--enrichment_clip"),
    type = "double",
    default = 2.5,
    help = "Symmetric log2-enrichment colour limit [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_tsv) || is.null(opt$output_pdf)) {
  stop("--input_tsv and --output_pdf are required", call. = FALSE)
}
if (!file.exists(opt$input_tsv)) {
  stop("Input TSV does not exist: ", opt$input_tsv, call. = FALSE)
}
if (opt$max_cell_types < 0L) {
  stop("--max_cell_types must be >= 0", call. = FALSE)
}
if (!is.finite(opt$enrichment_clip) || opt$enrichment_clip <= 0) {
  stop("--enrichment_clip must be a positive finite number", call. = FALSE)
}

source_data <- read.delim(
  opt$input_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)
required_columns <- c(
  "cluster", "cell_type", "mean_abundance", "log2_enrichment",
  "global_mean_abundance"
)
missing_columns <- setdiff(required_columns, colnames(source_data))
if (length(missing_columns) > 0L) {
  stop(
    "Input TSV is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

source_data$cluster <- as.character(source_data$cluster)
source_data$cell_type <- as.character(source_data$cell_type)
source_data$mean_abundance <- as.numeric(source_data$mean_abundance)
source_data$log2_enrichment <- as.numeric(source_data$log2_enrichment)
source_data$global_mean_abundance <- as.numeric(source_data$global_mean_abundance)
abundance_legend <- "Mean q05 abundance\n(cells per spot)"
if (
  "abundance_key" %in% colnames(source_data) &&
  length(unique(source_data$abundance_key)) == 1L &&
  unique(source_data$abundance_key) == "means_cell_abundance_w_sf"
) {
  abundance_legend <- "Mean posterior abundance\n(cells per spot)"
}
if (
  any(!is.finite(source_data$mean_abundance)) ||
  any(source_data$mean_abundance < 0) ||
  any(!is.finite(source_data$global_mean_abundance))
) {
  stop("Abundance columns must contain finite, non-negative values", call. = FALSE)
}

cluster_order <- unique(source_data$cluster)
cell_summary <- unique(source_data[c("cell_type", "global_mean_abundance")])
cell_summary <- cell_summary[order(-cell_summary$global_mean_abundance), , drop = FALSE]
if (opt$max_cell_types > 0L && nrow(cell_summary) > opt$max_cell_types) {
  message(
    "Displaying the top ", opt$max_cell_types, " of ", nrow(cell_summary),
    " cell types by sample-balanced tissue mean"
  )
  cell_summary <- cell_summary[seq_len(opt$max_cell_types), , drop = FALSE]
}

plot_data <- source_data[source_data$cell_type %in% cell_summary$cell_type, , drop = FALSE]
peak_cluster <- vapply(
  split(plot_data, plot_data$cell_type),
  function(cell_data) {
    values <- cell_data$log2_enrichment
    if (all(is.na(values))) return(NA_character_)
    cell_data$cluster[which.max(values)]
  },
  character(1)
)
cell_summary$peak_cluster <- unname(peak_cluster[cell_summary$cell_type])
cell_summary$peak_order <- match(cell_summary$peak_cluster, cluster_order)
cell_summary$peak_order[is.na(cell_summary$peak_order)] <- length(cluster_order) + 1L
cell_summary <- cell_summary[
  order(cell_summary$peak_order, -cell_summary$global_mean_abundance, cell_summary$cell_type),
  ,
  drop = FALSE
]

plot_data$cluster <- factor(plot_data$cluster, levels = cluster_order)
plot_data$cell_type <- factor(
  plot_data$cell_type,
  levels = rev(cell_summary$cell_type)
)

positive_abundance <- plot_data$mean_abundance[plot_data$mean_abundance > 0]
if (length(positive_abundance) == 0L) {
  stop("All selected mean-abundance values are zero", call. = FALSE)
}
size_cap <- as.numeric(stats::quantile(positive_abundance, 0.99, names = FALSE))
size_breaks <- size_cap * c(1 / 3, 2 / 3, 1)
plot_data$display_abundance <- pmin(plot_data$mean_abundance, size_cap)
plot_data$display_enrichment <- pmax(
  pmin(plot_data$log2_enrichment, opt$enrichment_clip),
  -opt$enrichment_clip
)

title <- "Relative cell abundance across spatial regions"
if (nrow(cell_summary) < length(unique(source_data$cell_type))) {
  title <- sprintf(
    "%s (top %d of %d)",
    title,
    nrow(cell_summary),
    length(unique(source_data$cell_type))
  )
}

p <- ggplot(plot_data, aes(x = .data$cluster, y = .data$cell_type)) +
  geom_point(
    aes(size = .data$display_abundance, fill = .data$display_enrichment),
    shape = 21,
    colour = "#4D4D4D",
    stroke = 0.18
  ) +
  scale_size_continuous(
    range = c(0, 5.1),
    limits = c(0, size_cap),
    breaks = size_breaks,
    labels = c(
      format(size_breaks[1:2], digits = 2),
      paste0("≥", format(size_breaks[3], digits = 2))
    ),
    name = abundance_legend
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "#F7F7F7",
    high = "#B84A3A",
    midpoint = 0,
    limits = c(-opt$enrichment_clip, opt$enrichment_clip),
    breaks = c(-opt$enrichment_clip, 0, opt$enrichment_clip),
    name = "log2 enrichment\nvs tissue mean"
  ) +
  labs(
    title = title,
    x = "Unsupervised spatial cluster",
    y = "Reference cell type"
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = grid::unit(25, "mm")),
    size = guide_legend(order = 2, override.aes = list(fill = "#BDBDBD"))
  ) +
  theme_classic(base_family = "sans", base_size = 7) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = if (length(cluster_order) > 15L) 90 else 45, hjust = 1),
    axis.text.y = element_text(size = 5.8, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.key.height = grid::unit(4, "mm"),
    legend.position = "right",
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
  )

fig_width <- min(7.2, max(4.6, 2.5 + 0.32 * length(cluster_order)))
fig_height <- min(8.5, max(3.2, 1.45 + 0.22 * nrow(cell_summary)))
dir.create(dirname(normalizePath(opt$output_pdf, mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(
  filename = opt$output_pdf,
  plot = p,
  device = pdf_device,
  width = fig_width,
  height = fig_height,
  units = "in",
  bg = "white"
)
message("Saved: ", opt$output_pdf)
