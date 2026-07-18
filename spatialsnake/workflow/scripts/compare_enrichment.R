#!/usr/bin/env Rscript

# Directional GO/KEGG over-representation analysis for one DEG contrast.
# GO and KEGG tables/figures are written separately, following the plotting
# logic used by the main enrichment workflow.

required_packages <- c(
  "optparse", "AnnotationDbi", "clusterProfiler", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them in the compare_stage environment.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(optparse)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(ggplot2)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--output_go_table", type = "character"),
  make_option("--output_kegg_table", type = "character"),
  make_option("--output_go_plot", type = "character"),
  make_option("--output_kegg_plot", type = "character"),
  make_option("--species", type = "character", default = "human"),
  make_option("--gene_id_type", type = "character", default = "auto"),
  make_option("--padj_cutoff", type = "double", default = 0.05),
  make_option("--logfc_cutoff", type = "double", default = 0.5),
  make_option("--go_ontology", type = "character", default = "BP"),
  make_option("--top_n", type = "integer", default = 10L)
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c(
  "input", "output_go_table", "output_kegg_table",
  "output_go_plot", "output_kegg_plot"
)
missing <- required[vapply(required, function(name) {
  is.null(opt[[name]]) || !nzchar(trimws(as.character(opt[[name]])))
}, logical(1))]
if (length(missing) > 0L) {
  stop("Missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop("Differential-expression input does not exist: ", opt$input, call. = FALSE)
}
if (!is.finite(opt$padj_cutoff) || opt$padj_cutoff <= 0 || opt$padj_cutoff > 1) {
  stop("--padj_cutoff must be in (0, 1]", call. = FALSE)
}
if (!is.finite(opt$logfc_cutoff) || opt$logfc_cutoff < 0) {
  stop("--logfc_cutoff must be >= 0", call. = FALSE)
}
if (is.na(opt$top_n) || opt$top_n < 1L) {
  stop("--top_n must be >= 1", call. = FALSE)
}

normalize_species <- function(value) {
  value <- tolower(trimws(as.character(value)))
  if (value %in% c("human", "homo_sapiens", "homo sapiens", "hs", "hsa")) return("human")
  if (value %in% c("mouse", "mus_musculus", "mus musculus", "mm", "mmu")) return("mouse")
  stop("--species must be human or mouse", call. = FALSE)
}

species <- normalize_species(opt$species)
gene_id_type <- toupper(trimws(as.character(opt$gene_id_type)))
if (gene_id_type == "AUTO") gene_id_type <- "auto"
if (!gene_id_type %in% c("auto", "SYMBOL", "ENSEMBL", "ENTREZID")) {
  stop("--gene_id_type must be auto, SYMBOL, ENSEMBL, or ENTREZID", call. = FALSE)
}
go_ontology <- toupper(trimws(as.character(opt$go_ontology)))
if (!go_ontology %in% c("BP", "CC", "MF", "ALL")) {
  stop("--go_ontology must be BP, CC, MF, or ALL", call. = FALSE)
}

orgdb_package <- if (species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
if (!requireNamespace(orgdb_package, quietly = TRUE)) {
  stop("Species '", species, "' requires R package ", orgdb_package, call. = FALSE)
}
orgdb <- get(orgdb_package, envir = asNamespace(orgdb_package))
kegg_prefix <- if (species == "human") "hsa" else "mmu"

output_paths <- c(
  opt$output_go_table, opt$output_kegg_table, opt$output_kegg_plot
)
for (path in output_paths) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  unlink(path, force = TRUE)
}
dir.create(dirname(opt$output_go_plot), recursive = TRUE, showWarnings = FALSE)
unlink(opt$output_go_plot, force = TRUE)

selected_go_ontologies <- function() {
  if (go_ontology == "ALL") c("BP", "CC", "MF") else go_ontology
}

go_plot_path <- function(ontology) {
  file.path(dirname(opt$output_go_plot), paste0("GO_", ontology, "_enrichment.pdf"))
}
for (ontology in c("BP", "CC", "MF")) {
  unlink(go_plot_path(ontology), force = TRUE)
}

empty_result <- function() {
  data.frame(
    database = character(),
    ontology = character(),
    direction = character(),
    category = character(),
    subcategory = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    FoldEnrichment = numeric(),
    Count = integer(),
    pvalue = numeric(),
    padj = numeric(),
    genes = character(),
    stringsAsFactors = FALSE
  )
}

write_csv <- function(data, path) {
  write.csv(data, file = path, row.names = FALSE, na = "")
}

theme_enrichment <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(linewidth = 0.2, colour = "#E5E7EB"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "#4B5563"),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(size = base_size - 0.5, face = "bold"),
      legend.position = "right"
    )
}

save_pdf <- function(plot, path, width_mm = 183, height_mm = 120) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white",
    limitsize = FALSE
  )
  qa_directory <- Sys.getenv("SPATIALSNAKE_FIGURE_QA_DIR", unset = "")
  if (nzchar(qa_directory) && requireNamespace("ragg", quietly = TRUE)) {
    dir.create(qa_directory, recursive = TRUE, showWarnings = FALSE)
    preview_path <- file.path(
      qa_directory,
      paste0(tools::file_path_sans_ext(basename(path)), ".png")
    )
    ragg::agg_png(
      preview_path,
      width = width_mm,
      height = height_mm,
      units = "mm",
      res = 180,
      background = "white"
    )
    print(plot)
    grDevices::dev.off()
  }
}

save_empty_plot <- function(path, title, message = "No significant enrichment terms") {
  figure <- ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 3, colour = "#4B5563") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title) +
    theme_void(base_size = 8, base_family = "sans") +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))
  save_pdf(figure, path)
}

read_deg <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    return(read.csv(
      path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
      na.strings = c("NA", "NaN", "nan", "")
    ))
  }
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection), add = TRUE)
  read.delim(
    connection, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "nan", "")
  )
}

deg <- read_deg(opt$input)
required_columns <- c(
  "gene", "gene_symbol", "log2FoldChange", "statistic", "pvalue", "padj",
  "comparison", "reference", "celltype", "algorithm"
)
missing_columns <- setdiff(required_columns, names(deg))
if (length(missing_columns) > 0L) {
  stop("DEG table is missing column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
for (column in c("log2FoldChange", "statistic", "pvalue", "padj")) {
  deg[[column]] <- suppressWarnings(as.numeric(deg[[column]]))
}
for (column in c("gene", "gene_symbol", "comparison", "reference", "celltype", "algorithm")) {
  deg[[column]] <- trimws(as.character(deg[[column]]))
}
if (any(is.na(deg$gene) | !nzchar(deg$gene)) || anyDuplicated(deg$gene)) {
  stop("DEG table requires one non-empty row per gene", call. = FALSE)
}
write_csv(empty_result(), opt$output_go_table)
write_csv(empty_result(), opt$output_kegg_table)
if (nrow(deg) == 0L) {
  for (ontology in selected_go_ontologies()) {
    save_empty_plot(go_plot_path(ontology), paste0("GO ", ontology, " enrichment"))
  }
  save_empty_plot(opt$output_kegg_plot, "KEGG enrichment")
  quit(save = "no", status = 0)
}

normalize_ids <- function(ids, id_type) {
  ids <- trimws(as.character(ids))
  if (id_type == "ENSEMBL") ids <- sub("\\.[0-9]+$", "", ids)
  ids
}

mapping_rate <- function(ids, id_type) {
  ids <- unique(normalize_ids(ids[!is.na(ids) & nzchar(ids)], id_type))
  if (length(ids) == 0L) return(0)
  valid <- AnnotationDbi::keys(orgdb, keytype = id_type)
  mean(ids %in% valid)
}

detect_id_type <- function(ids) {
  types <- c("SYMBOL", "ENSEMBL", "ENTREZID")
  rates <- setNames(vapply(types, function(type) mapping_rate(ids, type), numeric(1)), types)
  best <- max(rates)
  winners <- names(rates)[abs(rates - best) < 1e-12]
  if (!is.finite(best) || best == 0 || length(winners) != 1L) {
    stop("Unable to identify gene IDs unambiguously; set --gene_id_type explicitly", call. = FALSE)
  }
  message("Detected gene identifier type: ", winners, " (", sprintf("%.1f%%", 100 * best), " mapped)")
  winners
}

has_symbols <- any(!is.na(deg$gene_symbol) & nzchar(deg$gene_symbol))
if (has_symbols) {
  source_ids <- deg$gene_symbol
  source_type <- "SYMBOL"
} else {
  source_ids <- deg$gene
  source_type <- if (gene_id_type == "auto") detect_id_type(source_ids) else gene_id_type
}
source_ids <- normalize_ids(source_ids, source_type)

map_ids <- function(ids, id_type) {
  keys <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (length(keys) == 0L) {
    return(data.frame(source_id = character(), ENTREZID = character(), SYMBOL = character()))
  }
  mapped <- suppressMessages(
    AnnotationDbi::select(
      orgdb, keys = keys, keytype = id_type,
      columns = unique(c("ENTREZID", "SYMBOL"))
    )
  )
  if (id_type == "ENTREZID") {
    mapped$source_id <- mapped$ENTREZID
  } else if (id_type == "SYMBOL") {
    mapped$source_id <- mapped$SYMBOL
  } else {
    names(mapped)[names(mapped) == id_type] <- "source_id"
  }
  mapped <- mapped[, c("source_id", "ENTREZID", "SYMBOL"), drop = FALSE]
  mapped[] <- lapply(mapped, as.character)
  mapped <- mapped[
    !is.na(mapped$source_id) & nzchar(mapped$source_id) &
      !is.na(mapped$ENTREZID) & nzchar(mapped$ENTREZID),
    , drop = FALSE
  ]
  unique(mapped)
}

gene_map <- map_ids(source_ids, source_type)
if (nrow(gene_map) == 0L) {
  stop("No gene identifiers map to the selected species and ID type", call. = FALSE)
}
deg$row_id <- seq_len(nrow(deg))
id_frame <- data.frame(row_id = deg$row_id, source_id = source_ids, stringsAsFactors = FALSE)
mapped_deg <- merge(id_frame, gene_map, by = "source_id", all = FALSE, sort = FALSE)
mapped_deg <- merge(mapped_deg, deg, by = "row_id", all = FALSE, sort = FALSE)
universe <- unique(mapped_deg$ENTREZID[mapped_deg$row_id %in% deg$row_id[is.finite(deg$pvalue)]])
universe <- universe[!is.na(universe) & nzchar(universe)]
if (length(universe) == 0L) {
  message("No tested genes mapped to ENTREZID; enrichment is empty.")
  for (ontology in selected_go_ontologies()) {
    save_empty_plot(
      go_plot_path(ontology),
      paste0("GO ", ontology, " enrichment"),
      "No tested genes mapped to ENTREZID"
    )
  }
  save_empty_plot(opt$output_kegg_plot, "KEGG enrichment", "No tested genes mapped to ENTREZID")
  quit(save = "no", status = 0)
}

direction_rows <- list(
  higher_in_comparison = deg$row_id[
    is.finite(deg$padj) & deg$padj < opt$padj_cutoff &
      is.finite(deg$log2FoldChange) & deg$log2FoldChange >= opt$logfc_cutoff
  ],
  higher_in_reference = deg$row_id[
    is.finite(deg$padj) & deg$padj < opt$padj_cutoff &
      is.finite(deg$log2FoldChange) & deg$log2FoldChange <= -opt$logfc_cutoff
  ]
)
direction_ids <- lapply(direction_rows, function(rows) {
  intersect(unique(mapped_deg$ENTREZID[mapped_deg$row_id %in% rows]), universe)
})

symbol_map <- gene_map[
  !is.na(gene_map$SYMBOL) & nzchar(gene_map$SYMBOL), c("ENTREZID", "SYMBOL"), drop = FALSE
]
symbol_map <- symbol_map[!duplicated(symbol_map$ENTREZID), , drop = FALSE]
symbol_lookup <- setNames(symbol_map$SYMBOL, symbol_map$ENTREZID)
display_genes <- function(values) {
  vapply(as.character(values), function(value) {
    ids <- unique(strsplit(value, "/", fixed = TRUE)[[1]])
    symbols <- unique(unname(symbol_lookup[ids]))
    symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    if (length(symbols) > 0L) paste(symbols, collapse = "/") else paste(ids, collapse = "/")
  }, character(1))
}

standardize_result <- function(result, database, direction, ontology = "") {
  result <- as.data.frame(result)
  if (nrow(result) == 0L) return(empty_result())
  if (database == "GO" && "ONTOLOGY" %in% names(result)) ontology <- result$ONTOLOGY
  data.frame(
    database = database,
    ontology = as.character(rep(ontology, length.out = nrow(result))),
    direction = direction,
    category = "",
    subcategory = "",
    ID = as.character(result$ID),
    Description = as.character(result$Description),
    GeneRatio = as.character(result$GeneRatio),
    BgRatio = as.character(result$BgRatio),
    FoldEnrichment = NA_real_,
    Count = as.integer(result$Count),
    pvalue = as.numeric(result$pvalue),
    padj = as.numeric(result$p.adjust),
    genes = display_genes(result$geneID),
    stringsAsFactors = FALSE
  )
}

run_go <- function(ids, direction) {
  if (length(ids) == 0L) return(empty_result())
  result <- clusterProfiler::enrichGO(
    gene = ids, universe = universe, OrgDb = orgdb, keyType = "ENTREZID",
    ont = go_ontology, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE
  )
  standardize_result(result, "GO", direction, go_ontology)
}

kegg_categories <- local({
  data_environment <- new.env(parent = emptyenv())
  suppressWarnings(
    utils::data("kegg_category", package = "clusterProfiler", envir = data_environment)
  )
  if (!exists("kegg_category", envir = data_environment, inherits = FALSE)) {
    return(data.frame(id = character(), name = character(), category = character(), subcategory = character()))
  }
  value <- as.data.frame(get("kegg_category", envir = data_environment, inherits = FALSE))
  value[] <- lapply(value, as.character)
  value
})

get_kegg_annotation <- function() {
  if (!"PATH" %in% AnnotationDbi::columns(orgdb)) return(NULL)
  mapping <- suppressMessages(
    AnnotationDbi::select(orgdb, keys = universe, keytype = "ENTREZID", columns = "PATH")
  )
  mapping$ENTREZID <- as.character(mapping$ENTREZID)
  mapping$PATH <- as.character(mapping$PATH)
  mapping <- mapping[
    !is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID) &
      !is.na(mapping$PATH) & nzchar(mapping$PATH),
    c("PATH", "ENTREZID"), drop = FALSE
  ]
  if (nrow(mapping) == 0L) return(NULL)
  mapping$PATH <- paste0(kegg_prefix, mapping$PATH)
  term2gene <- unique(data.frame(term = mapping$PATH, gene = mapping$ENTREZID))
  pathway_ids <- sub("^[A-Za-z]+", "", unique(term2gene$term))
  pathway_names <- unname(setNames(kegg_categories$name, kegg_categories$id)[pathway_ids])
  pathway_names[is.na(pathway_names) | !nzchar(pathway_names)] <- paste(
    "KEGG pathway", unique(term2gene$term)[is.na(pathway_names) | !nzchar(pathway_names)]
  )
  term2name <- data.frame(
    term = unique(term2gene$term),
    name = pathway_names,
    stringsAsFactors = FALSE
  )
  list(term2gene = term2gene, term2name = term2name)
}

kegg_annotation <- get_kegg_annotation()
run_kegg <- function(ids, direction) {
  if (length(ids) == 0L || is.null(kegg_annotation)) return(empty_result())
  result <- clusterProfiler::enricher(
    gene = ids, universe = universe,
    TERM2GENE = kegg_annotation$term2gene,
    TERM2NAME = kegg_annotation$term2name,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500
  )
  standardize_result(result, "KEGG", direction)
}

results <- list()
for (direction in names(direction_ids)) {
  results[[paste0("GO_", direction)]] <- run_go(direction_ids[[direction]], direction)
  results[[paste0("KEGG_", direction)]] <- run_kegg(direction_ids[[direction]], direction)
}
enrichment <- do.call(rbind, c(results, list(empty_result())))
enrichment <- enrichment[
  is.finite(enrichment$padj) & enrichment$padj < opt$padj_cutoff,
  , drop = FALSE
]
if (nrow(enrichment) > 0L) {
  enrichment <- enrichment[order(enrichment$database, enrichment$direction, enrichment$padj), ]
  rownames(enrichment) <- NULL
  kegg_rows <- enrichment$database == "KEGG"
  pathway_ids <- sub("^[A-Za-z]+", "", enrichment$ID[kegg_rows])
  enrichment$category[kegg_rows] <- unname(
    setNames(kegg_categories$category, kegg_categories$id)[pathway_ids]
  )
  enrichment$subcategory[kegg_rows] <- unname(
    setNames(kegg_categories$subcategory, kegg_categories$id)[pathway_ids]
  )
  enrichment$category[is.na(enrichment$category)] <- ""
  enrichment$subcategory[is.na(enrichment$subcategory)] <- ""
}

ratio_numeric <- function(values) {
  vapply(as.character(values), function(value) {
    parts <- strsplit(value, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2L) return(NA_real_)
    numerator <- suppressWarnings(as.numeric(parts[1]))
    denominator <- suppressWarnings(as.numeric(parts[2]))
    if (!is.finite(numerator) || !is.finite(denominator) || denominator == 0) return(NA_real_)
    numerator / denominator
  }, numeric(1))
}

if (nrow(enrichment) > 0L) {
  enrichment$FoldEnrichment <- ratio_numeric(enrichment$GeneRatio) /
    ratio_numeric(enrichment$BgRatio)
}

go_data <- enrichment[enrichment$database == "GO", , drop = FALSE]
kegg_data <- enrichment[enrichment$database == "KEGG", , drop = FALSE]
write_csv(go_data, opt$output_go_table)
write_csv(kegg_data, opt$output_kegg_table)

comparison <- unique(deg$comparison)[1]
reference <- unique(deg$reference)[1]
celltype <- unique(deg$celltype)[1]
direction_labels <- c(
  higher_in_comparison = paste0("Higher in ", comparison),
  higher_in_reference = paste0("Higher in ", reference)
)

select_top_terms <- function(data, group_values) {
  if (nrow(data) == 0L) return(data)
  groups <- interaction(group_values, drop = TRUE, lex.order = TRUE)
  selected <- lapply(
    split(data, groups),
    function(part) head(part[order(part$padj, -part$Count), , drop = FALSE], opt$top_n)
  )
  result <- do.call(rbind, selected)
  rownames(result) <- NULL
  result
}

wrap_terms <- function(values, width = 42L) {
  vapply(values, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

plot_go_ontology <- function(data, ontology, path) {
  data <- data[data$ontology == ontology, , drop = FALSE]
  if (nrow(data) == 0L) {
    save_empty_plot(path, paste0(celltype, ": GO ", ontology, " enrichment"))
    return(invisible(FALSE))
  }
  plot_data <- select_top_terms(data, data$direction)
  plot_data$direction_label <- unname(direction_labels[plot_data$direction])
  plot_data$term_label <- wrap_terms(plot_data$Description)
  plot_data$row_key <- paste(
    plot_data$term_label, plot_data$direction_label, seq_len(nrow(plot_data)), sep = "|||"
  )
  plot_data$row_key <- factor(plot_data$row_key, levels = rev(plot_data$row_key))
  figure <- ggplot(plot_data, aes(x = Count, y = row_key)) +
    geom_col(width = 0.78, fill = "#4C78A8") +
    facet_grid(rows = vars(direction_label), scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(labels) sub("\\|\\|\\|.*$", "", labels)) +
    labs(
      title = paste0(celltype, ": GO ", ontology, " enrichment"),
      subtitle = paste0(comparison, " vs ", reference),
      x = "Gene count", y = "GO term"
    ) +
    theme_enrichment() +
    theme(
      axis.text.y = element_text(size = 6.2),
      legend.position = "none"
    )
  height_mm <- min(240, max(100, 6 * max(table(plot_data$direction_label)) + 45))
  save_pdf(figure, path, height_mm = height_mm)
  invisible(TRUE)
}

plot_go <- function(data) {
  data$ontology[data$ontology == ""] <- go_ontology
  for (ontology in selected_go_ontologies()) {
    plot_go_ontology(data, ontology, go_plot_path(ontology))
  }
}

plot_kegg <- function(data, path) {
  data <- data[is.finite(data$FoldEnrichment) & data$FoldEnrichment > 0, , drop = FALSE]
  if (nrow(data) == 0L) {
    save_empty_plot(path, paste0(celltype, ": KEGG enrichment"))
    return(invisible(FALSE))
  }
  plot_data <- select_top_terms(data, data$direction)
  plot_data$direction_label <- unname(direction_labels[plot_data$direction])
  plot_data$term_label <- wrap_terms(plot_data$Description)
  plot_data$row_key <- paste(
    plot_data$term_label, plot_data$direction_label, seq_len(nrow(plot_data)), sep = "|||"
  )
  plot_data$row_key <- factor(plot_data$row_key, levels = rev(plot_data$row_key))
  plot_data$plot_fdr <- pmax(plot_data$padj, .Machine$double.xmin)
  figure <- ggplot(plot_data, aes(x = FoldEnrichment, y = row_key)) +
    geom_vline(xintercept = 1, colour = "#9CA3AF", linewidth = 0.35, linetype = 2) +
    geom_point(aes(size = Count, colour = plot_fdr), alpha = 0.92) +
    facet_grid(rows = vars(direction_label), scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(labels) sub("\\|\\|\\|.*$", "", labels)) +
    scale_colour_gradient(
      low = "#CC0000", high = "#FFCCCC", trans = "log10",
      name = "Adjusted P value"
    ) +
    scale_size_continuous(range = c(2.5, 6), name = "Gene count") +
    labs(
      title = paste0(celltype, ": KEGG enrichment"),
      subtitle = paste0("Top pathways ranked by adjusted P value; ", comparison, " vs ", reference),
      x = "Fold enrichment", y = "KEGG pathway"
    ) +
    theme_enrichment() +
    theme(axis.text.y = element_text(size = 6.2))
  height_mm <- min(220, max(100, 7 * max(table(plot_data$direction_label)) + 45))
  save_pdf(figure, path, height_mm = height_mm)
  invisible(TRUE)
}

plot_go(go_data)
plot_kegg(kegg_data, opt$output_kegg_plot)

message(
  "Enrichment completed: ", nrow(go_data), " GO terms and ",
  nrow(kegg_data), " KEGG pathways written."
)
