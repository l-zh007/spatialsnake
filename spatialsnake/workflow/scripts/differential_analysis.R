suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggrepel)
})

option_list <- list(
  make_option("--algorithm", type = "character"),
  make_option("--counts", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--comparison", type = "character"),
  make_option("--reference", type = "character"),
  make_option("--celltype", type = "character"),
  make_option("--p_cutoff", type = "double", default = 0.05),
  make_option("--logfc_cutoff", type = "double", default = 0.5),
  make_option("--top_n", type = "integer", default = 20L),
  make_option("--result", type = "character"),
  make_option("--normalized", type = "character"),
  make_option("--gene_symbols", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c(
  "algorithm", "counts", "metadata", "output_dir", "comparison", "reference",
  "celltype", "result", "normalized", "gene_symbols"
)
missing_args <- required[vapply(required, function(key) {
  is.null(opt[[key]]) || !nzchar(trimws(as.character(opt[[key]])))
}, logical(1))]
if (length(missing_args) > 0L) {
  stop("Missing required argument(s): ", paste(missing_args, collapse = ", "), call. = FALSE)
}
if (!is.finite(opt$p_cutoff) || opt$p_cutoff <= 0 || opt$p_cutoff > 1) {
  stop("--p_cutoff must be in (0, 1]", call. = FALSE)
}
if (!is.finite(opt$logfc_cutoff) || opt$logfc_cutoff < 0) {
  stop("--logfc_cutoff must be >= 0", call. = FALSE)
}
if (opt$top_n < 1L) {
  stop("--top_n must be >= 1", call. = FALSE)
}

algorithm_key <- tolower(trimws(opt$algorithm))
algorithm <- if (algorithm_key %in% c("deseq2", "deseq")) {
  "DESeq2"
} else if (algorithm_key == "edger") {
  "edgeR"
} else {
  stop("--algorithm must be DESeq2 or edgeR", call. = FALSE)
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection), add = TRUE)
  read.delim(
    connection,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "")
  )
}

write_tsv <- function(data, path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(
    data,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

read_result_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    return(read.csv(
      path,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("NA", "NaN", "")
    ))
  }
  read_tsv(path)
}

write_result_table <- function(data, path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    write.csv(data, file = path, row.names = FALSE, na = "NA")
  } else {
    write_tsv(data, path)
  }
}

read_matrix <- function(path) {
  data <- read_tsv(path)
  if (ncol(data) < 2L || names(data)[1] != "replicate") {
    stop("Matrix file must start with a 'replicate' column: ", path, call. = FALSE)
  }
  ids <- as.character(data$replicate)
  if (anyDuplicated(ids)) {
    stop("Matrix contains duplicate replicate IDs: ", path, call. = FALSE)
  }
  values <- as.matrix(data[, -1, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "double")
  if (any(!is.finite(values))) {
    stop("Matrix contains non-finite values: ", path, call. = FALSE)
  }
  rownames(values) <- ids
  values
}

metadata <- read_tsv(opt$metadata)
if (!all(c("replicate", "condition") %in% names(metadata))) {
  stop("Metadata must contain replicate and condition columns", call. = FALSE)
}
metadata$replicate <- as.character(metadata$replicate)
metadata$condition <- as.character(metadata$condition)
if (anyDuplicated(metadata$replicate)) {
  stop("Metadata contains duplicate replicate IDs", call. = FALSE)
}
rownames(metadata) <- metadata$replicate

counts <- read_matrix(opt$counts)
if (!setequal(rownames(counts), metadata$replicate)) {
  stop("Counts and metadata contain different replicate IDs", call. = FALSE)
}
metadata <- metadata[rownames(counts), , drop = FALSE]
if (any(counts < 0) || any(abs(counts - round(counts)) > 1e-6)) {
  stop("Pseudobulk counts must be non-negative integers", call. = FALSE)
}
storage.mode(counts) <- "integer"

gene_map <- read_tsv(opt$gene_symbols)
if (!all(c("gene", "gene_symbol") %in% names(gene_map))) {
  stop("Gene-symbol mapping must contain gene and gene_symbol columns", call. = FALSE)
}
gene_map$gene <- as.character(gene_map$gene)
gene_map$gene_symbol <- as.character(gene_map$gene_symbol)
if (anyDuplicated(gene_map$gene)) {
  stop("Gene-symbol mapping contains duplicate gene IDs", call. = FALSE)
}
symbol_lookup <- setNames(gene_map$gene_symbol, gene_map$gene)

final_columns <- c(
  "gene", "gene_symbol", "base_mean", "log2FoldChange", "statistic",
  "pvalue", "padj", "comparison", "reference", "celltype", "algorithm"
)

if (algorithm == "edgeR") {
  suppressPackageStartupMessages(library(edgeR))
  metadata$condition <- factor(metadata$condition, levels = c(opt$reference, opt$comparison))
  if (anyNA(metadata$condition)) {
    stop("Metadata contains conditions outside the requested contrast", call. = FALSE)
  }
  design <- model.matrix(~condition, data = metadata)
  if (qr(design)$rank != ncol(design)) {
    stop("edgeR design matrix is not full rank", call. = FALSE)
  }
  if (nrow(design) <= ncol(design)) {
    stop("edgeR design has no residual degrees of freedom", call. = FALSE)
  }

  y <- DGEList(counts = t(counts))
  keep <- filterByExpr(y, group = metadata$condition)
  if (sum(keep) < 2L) {
    stop("edgeR filterByExpr retained fewer than two genes", call. = FALSE)
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  target_name <- make.names(paste0("condition", opt$comparison))
  coefficient <- which(make.names(colnames(design)) == target_name)
  if (length(coefficient) != 1L) {
    stop(
      "Unable to identify the edgeR comparison coefficient. Available: ",
      paste(colnames(design), collapse = ", "),
      call. = FALSE
    )
  }
  test <- glmQLFTest(fit, coef = coefficient)
  edge <- topTags(test, n = Inf, sort.by = "none")$table
  genes <- rownames(edge)
  base_mean <- rowMeans(cpm(y, normalized.lib.sizes = TRUE, log = FALSE))
  results <- data.frame(
    gene = genes,
    gene_symbol = unname(symbol_lookup[genes]),
    base_mean = as.numeric(base_mean[genes]),
    log2FoldChange = as.numeric(edge[genes, "logFC"]),
    statistic = sign(edge[genes, "logFC"]) * sqrt(pmax(edge[genes, "F"], 0)),
    pvalue = as.numeric(edge[genes, "PValue"]),
    padj = as.numeric(edge[genes, "FDR"]),
    comparison = opt$comparison,
    reference = opt$reference,
    celltype = opt$celltype,
    algorithm = "edgeR",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  results <- results[, final_columns, drop = FALSE]
  write_result_table(results, opt$result)
  normalized <- t(cpm(y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2))
} else {
  results <- read_result_table(opt$result)
  missing_result <- setdiff(final_columns, names(results))
  if (length(missing_result) > 0L) {
    stop("DESeq2 result is missing column(s): ", paste(missing_result, collapse = ", "), call. = FALSE)
  }
  results <- results[, final_columns, drop = FALSE]
  normalized <- read_matrix(opt$normalized)
}

if (!all(colnames(normalized) %in% results$gene)) {
  stop("Normalized-expression genes do not match the DEG result", call. = FALSE)
}
if (!setequal(rownames(normalized), metadata$replicate)) {
  stop("Normalized expression and metadata contain different replicate IDs", call. = FALSE)
}
normalized <- normalized[metadata$replicate, , drop = FALSE]

theme_publication <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 0.8, hjust = 0),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      panel.grid = element_blank()
    )
}

save_pdf <- function(plot, path, width_mm, height_mm) {
  grDevices::cairo_pdf(
    filename = path,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    family = "sans"
  )
  print(plot)
  grDevices::dev.off()
}

plot_volcano <- function(results, path) {
  plot_data <- results[
    is.finite(results$log2FoldChange) & is.finite(results$padj) & results$padj >= 0,
    , drop = FALSE
  ]
  if (nrow(plot_data) == 0L) {
    title <- paste0(opt$celltype, ": ", opt$comparison, " vs ", opt$reference)
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No finite FDR values available", size = 2.6) +
      xlim(-1, 1) +
      ylim(-1, 1) +
      labs(title = title) +
      theme_void(base_size = 7, base_family = "sans") +
      theme(plot.title = element_text(face = "bold", size = 7.8, hjust = 0))
    save_pdf(p, path, width_mm = 89, height_mm = 82)
    return(invisible(FALSE))
  }
  plot_data$minus_log10_fdr <- -log10(pmax(plot_data$padj, 1e-300))
  plot_data$status <- "Not significant"
  plot_data$status[
    plot_data$padj < opt$p_cutoff & plot_data$log2FoldChange >= opt$logfc_cutoff
  ] <- "Higher in comparison"
  plot_data$status[
    plot_data$padj < opt$p_cutoff & plot_data$log2FoldChange <= -opt$logfc_cutoff
  ] <- "Higher in reference"
  plot_data$status <- factor(
    plot_data$status,
    levels = c("Not significant", "Higher in comparison", "Higher in reference")
  )
  label_data <- plot_data[plot_data$status != "Not significant", , drop = FALSE]
  if (nrow(label_data) > 0L) {
    label_data <- label_data[order(label_data$padj, -abs(label_data$log2FoldChange)), , drop = FALSE]
    label_data <- head(label_data, 10L)
    label_data$display_gene <- ifelse(
      is.na(label_data$gene_symbol) | !nzchar(label_data$gene_symbol),
      label_data$gene,
      label_data$gene_symbol
    )
  }
  title <- paste0(opt$celltype, ": ", opt$comparison, " vs ", opt$reference)
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = minus_log10_fdr, colour = status)) +
    geom_point(size = 0.9, alpha = 0.72, stroke = 0) +
    geom_vline(xintercept = c(-opt$logfc_cutoff, opt$logfc_cutoff), colour = "#A7A7A7", linewidth = 0.3, linetype = 2) +
    geom_hline(yintercept = -log10(opt$p_cutoff), colour = "#A7A7A7", linewidth = 0.3, linetype = 2) +
    scale_colour_manual(
      values = c(
        "Not significant" = "#B8B8B8",
        "Higher in comparison" = "#2C7FB8",
        "Higher in reference" = "#D9822B"
      ),
      drop = FALSE
    ) +
    labs(
      title = title,
      x = expression(log[2] * " fold change"),
      y = expression(-log[10] * " FDR"),
      colour = NULL
    ) +
    theme_publication() +
    theme(legend.position = "top")
  if (nrow(label_data) > 0L) {
    p <- p + ggrepel::geom_text_repel(
      data = label_data,
      aes(label = display_gene),
      size = 2.1,
      min.segment.length = 0,
      segment.size = 0.25,
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = Inf,
      show.legend = FALSE,
      seed = 1
    )
  }
  save_pdf(p, path, width_mm = 89, height_mm = 82)
}

plot_heatmap <- function(results, normalized, metadata, path) {
  save_empty_heatmap <- function(message) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = message, size = 2.8, colour = "#4B5563") +
      xlim(-1, 1) +
      ylim(-1, 1) +
      labs(title = paste0(opt$celltype, ": top differential genes")) +
      theme_void(base_size = 7, base_family = "sans") +
      theme(plot.title = element_text(face = "bold", size = 7.8, hjust = 0))
    save_pdf(p, path, width_mm = 183, height_mm = 82)
    invisible(FALSE)
  }
  significant <- results[
    is.finite(results$padj) & results$padj < opt$p_cutoff &
      abs(results$log2FoldChange) >= opt$logfc_cutoff,
    , drop = FALSE
  ]
  if (nrow(significant) == 0L) {
    return(save_empty_heatmap("No genes pass the selected FDR and log2FC thresholds"))
  }
  significant <- significant[order(significant$padj, -abs(significant$log2FoldChange)), , drop = FALSE]
  significant <- head(significant, opt$top_n)
  genes <- intersect(significant$gene, colnames(normalized))
  if (length(genes) == 0L) {
    return(save_empty_heatmap("Significant genes are absent from normalized expression"))
  }
  expression <- t(normalized[, genes, drop = FALSE])
  row_means <- rowMeans(expression)
  row_sds <- apply(expression, 1, sd)
  row_sds[!is.finite(row_sds) | row_sds == 0] <- 1
  scaled <- sweep(sweep(expression, 1, row_means, "-"), 1, row_sds, "/")
  scaled[scaled > 2.5] <- 2.5
  scaled[scaled < -2.5] <- -2.5

  sample_order <- metadata$replicate[order(metadata$condition, metadata$replicate)]
  scaled <- scaled[, sample_order, drop = FALSE]
  heat <- as.data.frame(as.table(scaled), stringsAsFactors = FALSE)
  names(heat) <- c("gene", "replicate", "zscore")
  display_symbols <- ifelse(
    is.na(significant$gene_symbol) | !nzchar(significant$gene_symbol),
    significant$gene,
    significant$gene_symbol
  )
  symbol_map <- setNames(display_symbols, significant$gene)
  ordered_symbols <- make.unique(as.character(unname(symbol_map[genes])))
  display_map <- setNames(ordered_symbols, genes)
  heat$gene_symbol <- unname(display_map[as.character(heat$gene)])
  heat$gene_symbol <- factor(heat$gene_symbol, levels = rev(ordered_symbols))
  heat$replicate <- factor(heat$replicate, levels = sample_order)
  condition_labels <- setNames(metadata$condition, metadata$replicate)
  x_labels <- paste0(sample_order, "\n", condition_labels[sample_order])
  names(x_labels) <- sample_order
  p <- ggplot(heat, aes(x = replicate, y = gene_symbol, fill = zscore)) +
    geom_tile(linewidth = 0.18, colour = "white") +
    scale_fill_gradient2(
      low = "#3B6FB6",
      mid = "#F7F7F7",
      high = "#C94C4C",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      name = "Row z-score"
    ) +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = paste0(opt$celltype, ": top differential genes"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 7, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 6),
      axis.text.y = element_text(colour = "black", size = 6),
      plot.title = element_text(face = "bold", size = 7.8, hjust = 0),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6)
    )
  height_mm <- min(180, max(70, 6 * length(genes) + 30))
  save_pdf(p, path, width_mm = 183, height_mm = height_mm)
  invisible(TRUE)
}

plot_volcano(results, file.path(opt$output_dir, "volcano.pdf"))
plot_heatmap(results, normalized, metadata, file.path(opt$output_dir, "heatmap.pdf"))
