library(optparse)
library(CellChat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(grid)

option_list <- list(
  make_option("--input_rds1", type = "character"),
  make_option("--input_rds2", type = "character"),
  make_option("--sample_name1", type = "character", default = "condition1"),
  make_option("--sample_name2", type = "character", default = "condition2"),
  make_option(c("-o", "--output_dir"), type = "character"),
  make_option("--focus_cells", type = "character", default = ""),
  make_option("--cell_pairs", type = "character", default = ""),
  make_option("--source_cells", type = "character", default = ""),
  make_option("--target_cells", type = "character", default = ""),
  make_option("--pathways", type = "character", default = ""),
  make_option("--lr_pairs", type = "character", default = ""),
  make_option("--top_cell_pairs", type = "integer", default = 3),
  make_option("--top_pathways", type = "integer", default = 3),
  make_option("--top_lr", type = "integer", default = 20),
  make_option("--bubble_angle", type = "double", default = 45),
  make_option("--bubble_remove_isolate", type = "logical", default = TRUE),
  make_option("--plot_advanced", type = "logical", default = TRUE),
  make_option("--gene_colors", type = "character", default = "white,#FEC44F,#D95F0E"),
  make_option("--gene_plot_type", type = "character", default = "dot"),
  make_option("--save_merged", type = "logical", default = TRUE),
  make_option("--receiver_cells", type = "character", default = ""),
  make_option("--do_single_bubble", type = "logical", default = FALSE),
  make_option("--pair_lr_use", type = "character", default = "")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("input_rds1", "input_rds2", "output_dir")
missing_args <- required[vapply(required, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("CellChat comparison requires exactly two RDS inputs and an output directory; missing: ",
       paste(missing_args, collapse = ", "))
}
if (!file.exists(opt$input_rds1) || !file.exists(opt$input_rds2)) {
  stop("Both CellChat RDS files must exist")
}
if (opt$top_cell_pairs < 1 || opt$top_pathways < 1 || opt$top_lr < 1) {
  stop("top_cell_pairs, top_pathways, and top_lr must be positive integers")
}
if (!opt$gene_plot_type %in% c("dot", "violin", "bar")) {
  stop("gene_plot_type must be one of: dot, violin, bar")
}

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

split_values <- function(x) {
  if (is.null(x) || !nzchar(trimws(x)) || tolower(trimws(x)) %in% c("none", "null", "na")) {
    return(character(0))
  }
  unique(trimws(unlist(strsplit(x, ",", fixed = TRUE))))
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

write_table <- function(x, filename) {
  write.csv(x, file.path(opt$output_dir, filename), row.names = FALSE, quote = TRUE)
}

message_plot <- function(filename, title, label) {
  png(file.path(opt$output_dir, filename), width = 8, height = 4.5, units = "in", res = 300)
  plot.new()
  title(main = title, cex.main = 1.1)
  text(0.5, 0.5, label, cex = 1)
  dev.off()
}

plot_difference_circle <- function(net, title_name, edge_max, n_vertices) {
  if (any(net > 0, na.rm = TRUE)) {
    netVisual_circle(net, weight.scale = TRUE, edge.weight.max = edge_max,
                     vertex.weight = rep(20, n_vertices), title.name = NULL)
    mtext(title_name, side = 3, line = 0.2, cex = 0.9)
  } else {
    plot.new()
    title(main = title_name, cex.main = 1)
    text(0.5, 0.5, "No higher interactions", cex = 0.9)
  }
}

role_heatmap_has_range <- function(cellchat, pathway) {
  if (!pathway %in% names(cellchat@netP$centr)) return(FALSE)
  centrality <- cellchat@netP$centr[[pathway]]
  values <- centrality$outdeg + centrality$indeg
  maximum <- max(values, na.rm = TRUE)
  if (!is.finite(maximum) || maximum <= 0) return(FALSE)
  normalized <- values / maximum
  length(unique(normalized[is.finite(normalized) & normalized > 0])) >= 2
}

condition_direction <- function(x, name1, name2, tolerance = 1e-15) {
  ifelse(x > tolerance, paste0(name2, "_higher"),
         ifelse(x < -tolerance, paste0(name1, "_higher"), "equal"))
}

validate_cellchat <- function(x, label) {
  if (!methods::is(x, "CellChat")) {
    stop(label, " is not a CellChat object")
  }
  if (!all(c("prob", "pval", "count", "weight") %in% names(x@net))) {
    stop(label, " does not contain a completed CellChat net slot")
  }
  if (!all(c("pathways", "prob") %in% names(x@netP))) {
    stop(label, " does not contain a completed CellChat netP slot")
  }
  if (length(levels(x@idents)) < 2) {
    stop(label, " must contain at least two cell identities")
  }
  x
}

resolve_labels <- function(values, available, label) {
  values <- split_values(values)
  if (length(values) == 0) {
    return(character(0))
  }
  matched <- available[match(tolower(values), tolower(available))]
  invalid <- values[is.na(matched)]
  if (length(invalid) > 0) {
    warning(label, " not found and ignored: ", paste(invalid, collapse = ", "))
  }
  unique(matched[!is.na(matched)])
}

parse_cell_pairs <- function(value, available) {
  items <- split_values(value)
  if (length(items) == 0) {
    return(data.frame(source = character(0), target = character(0)))
  }
  rows <- list()
  for (item in items) {
    bidirectional <- grepl("<->", item, fixed = TRUE)
    separator <- if (bidirectional) "<->" else if (grepl("->", item, fixed = TRUE)) "->" else "|"
    parts <- trimws(unlist(strsplit(item, separator, fixed = TRUE)))
    if (length(parts) != 2) {
      warning("Invalid cell-pair specification ignored: ", item)
      next
    }
    matched <- available[match(tolower(parts), tolower(available))]
    if (any(is.na(matched))) {
      warning("Cell pair contains unknown identities and was ignored: ", item)
      next
    }
    rows[[length(rows) + 1]] <- data.frame(source = matched[1], target = matched[2])
    if (bidirectional) {
      rows[[length(rows) + 1]] <- data.frame(source = matched[2], target = matched[1])
    }
  }
  if (length(rows) == 0) {
    return(data.frame(source = character(0), target = character(0)))
  }
  unique(do.call(rbind, rows))
}

cross_pairs <- function(sources, targets) {
  if (length(sources) == 0 || length(targets) == 0) {
    return(data.frame(source = character(0), target = character(0)))
  }
  expand.grid(source = sources, target = targets, stringsAsFactors = FALSE)
}

matrix_table <- function(matrix1, matrix2, value_name, name1, name2) {
  idx <- expand.grid(source = rownames(matrix1), target = colnames(matrix1), stringsAsFactors = FALSE)
  idx[[paste0(value_name, "_", name1)]] <- as.vector(matrix1)
  idx[[paste0(value_name, "_", name2)]] <- as.vector(matrix2)
  idx[[paste0(value_name, "_diff")]] <- idx[[paste0(value_name, "_", name2)]] - idx[[paste0(value_name, "_", name1)]]
  idx
}

standardize_lr <- function(x, condition) {
  d <- subsetCommunication(x)
  required_cols <- c("source", "target", "ligand", "receptor", "prob", "pval",
                     "interaction_name", "interaction_name_2", "pathway_name")
  for (column in setdiff(required_cols, names(d))) {
    d[[column]] <- NA
  }
  d <- d[, required_cols]
  d$condition <- condition
  d
}

merge_lr_tables <- function(d1, d2, name1, name2) {
  keys <- c("source", "target", "ligand", "receptor", "interaction_name",
            "interaction_name_2", "pathway_name")
  x1 <- d1[, c(keys, "prob", "pval")]
  x2 <- d2[, c(keys, "prob", "pval")]
  names(x1)[names(x1) == "prob"] <- paste0("prob_", name1)
  names(x1)[names(x1) == "pval"] <- paste0("pval_", name1)
  names(x2)[names(x2) == "prob"] <- paste0("prob_", name2)
  names(x2)[names(x2) == "pval"] <- paste0("pval_", name2)
  out <- merge(x1, x2, by = keys, all = TRUE)
  prob1 <- paste0("prob_", name1)
  prob2 <- paste0("prob_", name2)
  out[[prob1]][is.na(out[[prob1]])] <- 0
  out[[prob2]][is.na(out[[prob2]])] <- 0
  out$prob_diff <- out[[prob2]] - out[[prob1]]
  out$direction <- condition_direction(out$prob_diff, name1, name2)
  out[order(-abs(out$prob_diff), -pmax(out[[prob1]], out[[prob2]])), ]
}

pathway_summary <- function(cellchat1, cellchat2, name1, name2) {
  pathways <- union(cellchat1@netP$pathways, cellchat2@netP$pathways)
  get_flow <- function(x, pathway) {
    index <- match(pathway, x@netP$pathways)
    if (is.na(index)) 0 else sum(x@netP$prob[, , index], na.rm = TRUE)
  }
  out <- data.frame(pathway = pathways, stringsAsFactors = FALSE)
  out[[paste0("flow_", name1)]] <- vapply(pathways, function(p) get_flow(cellchat1, p), numeric(1))
  out[[paste0("flow_", name2)]] <- vapply(pathways, function(p) get_flow(cellchat2, p), numeric(1))
  out$flow_diff <- out[[paste0("flow_", name2)]] - out[[paste0("flow_", name1)]]
  out$direction <- condition_direction(out$flow_diff, name1, name2)
  out[order(-abs(out$flow_diff)), ]
}

name1 <- if (nzchar(opt$sample_name1)) opt$sample_name1 else tools::file_path_sans_ext(basename(opt$input_rds1))
name2 <- if (nzchar(opt$sample_name2)) opt$sample_name2 else tools::file_path_sans_ext(basename(opt$input_rds2))
if (name1 == name2) {
  stop("The two condition names must be different")
}

cellchat1 <- validate_cellchat(readRDS(opt$input_rds1), name1)
cellchat2 <- validate_cellchat(readRDS(opt$input_rds2), name2)

datatype1 <- as.character(cellchat1@options$datatype)
datatype2 <- as.character(cellchat2@options$datatype)
if (!identical(datatype1, datatype2)) {
  stop("CellChat datatype mismatch: ", datatype1, " versus ", datatype2)
}
db1 <- sort(as.character(cellchat1@DB$interaction$interaction_name))
db2 <- sort(as.character(cellchat2@DB$interaction$interaction_name))
if (!identical(db1, db2)) {
  stop("The two CellChat objects use incompatible interaction databases/species")
}

cell_levels <- union(levels(cellchat1@idents), levels(cellchat2@idents))
if (!identical(levels(cellchat1@idents), cell_levels)) {
  cellchat1 <- liftCellChat(cellchat1, cell_levels)
}
if (!identical(levels(cellchat2@idents), cell_levels)) {
  cellchat2 <- liftCellChat(cellchat2, cell_levels)
}

cco.list <- setNames(list(cellchat1, cellchat2), c(name1, name2))
cellchat_merged <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
joint_levels <- levels(cellchat_merged@idents$joint)
if (!setequal(joint_levels, cell_levels)) {
  stop("Merged CellChat identities do not match the aligned identity union")
}

message("CellChat comparison: ", name2, " - ", name1)
message("Cell identities: ", paste(joint_levels, collapse = ", "))
message("Datatype: ", datatype1, "; CellChat package: ", as.character(packageVersion("CellChat")))

global_summary <- data.frame(
  condition = c(name1, name2),
  interaction_count = c(sum(cellchat1@net$count), sum(cellchat2@net$count)),
  interaction_weight = c(sum(cellchat1@net$weight), sum(cellchat2@net$weight))
)
write_table(global_summary, "compare_global_summary.csv")

count_table <- matrix_table(cellchat1@net$count, cellchat2@net$count, "count", name1, name2)
weight_table <- matrix_table(cellchat1@net$weight, cellchat2@net$weight, "weight", name1, name2)
cell_pair_table <- merge(count_table, weight_table, by = c("source", "target"))
cell_pair_table$direction <- condition_direction(cell_pair_table$weight_diff, name1, name2)
cell_pair_table <- cell_pair_table[order(-abs(cell_pair_table$weight_diff), -abs(cell_pair_table$count_diff)), ]
write_table(cell_pair_table, "compare_cell_pair.csv")

pathway_table <- pathway_summary(cellchat1, cellchat2, name1, name2)
write_table(pathway_table, "compare_pathway.csv")

lr1 <- standardize_lr(cellchat1, name1)
lr2 <- standardize_lr(cellchat2, name2)
lr_table <- merge_lr_tables(lr1, lr2, name1, name2)
write_table(lr_table, "compare_lr.csv")

manual_pairs <- parse_cell_pairs(opt$cell_pairs, joint_levels)
focus_cells <- resolve_labels(opt$focus_cells, joint_levels, "focus_cells")
source_cells <- resolve_labels(opt$source_cells, joint_levels, "source_cells")
target_cells <- resolve_labels(opt$target_cells, joint_levels, "target_cells")

if (nrow(manual_pairs) > 0) {
  selected_pairs <- manual_pairs
  selection_source <- "manual_cell_pairs"
  separate_pair_plots <- TRUE
} else if (length(focus_cells) > 0) {
  selected_pairs <- cross_pairs(focus_cells, focus_cells)
  selection_source <- "focus_cells"
  separate_pair_plots <- FALSE
} else if (length(source_cells) > 0 || length(target_cells) > 0) {
  if (length(source_cells) == 0) source_cells <- joint_levels
  if (length(target_cells) == 0) target_cells <- joint_levels
  selected_pairs <- cross_pairs(source_cells, target_cells)
  selection_source <- "source_target"
  separate_pair_plots <- FALSE
} else {
  candidates <- cell_pair_table[cell_pair_table$count_diff != 0 | cell_pair_table$weight_diff != 0, ]
  if (nrow(candidates) == 0) {
    strength1 <- paste0("weight_", name1)
    strength2 <- paste0("weight_", name2)
    candidates <- cell_pair_table[(cell_pair_table[[strength1]] + cell_pair_table[[strength2]]) > 0, ]
    candidates <- candidates[order(-(candidates[[strength1]] + candidates[[strength2]])), ]
    selection_source <- "shared_strength"
  } else {
    selection_source <- "automatic_difference"
  }
  selected_pairs <- candidates[seq_len(min(opt$top_cell_pairs, nrow(candidates))), c("source", "target"), drop = FALSE]
  separate_pair_plots <- TRUE
}

if (separate_pair_plots && selection_source %in% c("automatic_difference", "shared_strength") && nrow(selected_pairs) > 1) {
  selected_keys <- paste(selected_pairs$source, selected_pairs$target, sep = "\r")
  table_keys <- paste(cell_pair_table$source, cell_pair_table$target, sep = "\r")
  selected_stats <- cell_pair_table[match(selected_keys, table_keys), ]
  weight1_name <- paste0("weight_", name1)
  weight2_name <- paste0("weight_", name2)
  both_conditions <- selected_stats[[weight1_name]] > 0 & selected_stats[[weight2_name]] > 0
  selected_pairs <- selected_pairs[order(-both_conditions, -abs(selected_stats$weight_diff)), , drop = FALSE]
}

if (nrow(selected_pairs) == 0) {
  warning("No non-zero cell-pair communication is available for focused visualization")
}

pair_key <- paste(selected_pairs$source, selected_pairs$target, sep = "\r")
lr_scope <- lr_table[paste(lr_table$source, lr_table$target, sep = "\r") %in% pair_key, , drop = FALSE]
available_pathways <- unique(lr_scope$pathway_name)
manual_pathways <- split_values(opt$pathways)
valid_manual_pathways <- available_pathways[match(tolower(manual_pathways), tolower(available_pathways))]
valid_manual_pathways <- unique(valid_manual_pathways[!is.na(valid_manual_pathways)])
if (length(manual_pathways) > length(valid_manual_pathways)) {
  warning("Some requested pathways were unavailable in the selected cell scope; automatic selection is used when necessary")
}

if (length(valid_manual_pathways) > 0) {
  selected_pathways <- valid_manual_pathways
  pathway_selection <- "manual"
} else if (nrow(lr_scope) > 0) {
  pathway_scores <- aggregate(abs(lr_scope$prob_diff), list(pathway = lr_scope$pathway_name), sum)
  names(pathway_scores)[2] <- "score"
  pathway_scores <- pathway_scores[order(-pathway_scores$score), ]
  selected_pathways <- head(pathway_scores$pathway, opt$top_pathways)
  pathway_selection <- "automatic"
} else {
  selected_pathways <- character(0)
  pathway_selection <- "none"
}

lr_scope <- lr_scope[lr_scope$pathway_name %in% selected_pathways, , drop = FALSE]
manual_lr <- split_values(opt$lr_pairs)
manual_interactions <- character(0)
if (length(manual_lr) > 0 && nrow(lr_scope) > 0) {
  for (pair in manual_lr) {
    genes <- trimws(unlist(strsplit(pair, "|", fixed = TRUE)))
    if (length(genes) == 2) {
      manual_interactions <- c(manual_interactions,
        lr_scope$interaction_name[tolower(lr_scope$ligand) == tolower(genes[1]) &
                                  tolower(lr_scope$receptor) == tolower(genes[2])])
    }
  }
  manual_interactions <- unique(manual_interactions)
  if (length(manual_interactions) == 0) {
    warning("Requested ligand-receptor pairs were unavailable in the selected cell/pathway scope; using automatic LR selection")
  }
}

if (length(manual_interactions) > 0) {
  selected_interaction_names <- manual_interactions
  lr_selection <- "manual"
} else if (nrow(lr_scope) > 0) {
  lr_scores <- aggregate(abs(lr_scope$prob_diff), list(interaction_name = lr_scope$interaction_name), max)
  names(lr_scores)[2] <- "score"
  lr_scores <- lr_scores[order(-lr_scores$score), ]
  selected_interaction_names <- head(lr_scores$interaction_name, opt$top_lr)
  lr_selection <- "automatic"
} else {
  selected_interaction_names <- character(0)
  lr_selection <- "none"
}

selected <- lr_scope[lr_scope$interaction_name %in% selected_interaction_names, , drop = FALSE]
selected$selection_source <- selection_source
selected$cell_pair <- paste(selected$source, selected$target, sep = "->")
selected$pathway_rank <- match(selected$pathway_name, selected_pathways)
selected$lr_rank <- match(selected$interaction_name, selected_interaction_names)
selected$plot_file <- ""

if (nrow(selected) > 0) {
  plot_specs <- if (separate_pair_plots) split(selected_pairs, seq_len(nrow(selected_pairs))) else list(selected_pairs)
  for (i in seq_along(plot_specs)) {
    spec <- plot_specs[[i]]
    if (is.data.frame(spec) && nrow(spec) == 1) {
      spec_sources <- spec$source
      spec_targets <- spec$target
      selected_index <- selected$source == spec_sources & selected$target == spec_targets
      pair_label <- paste0(safe_name(spec_sources), "_to_", safe_name(spec_targets))
    } else {
      spec_sources <- unique(spec$source)
      spec_targets <- unique(spec$target)
      selected_index <- paste(selected$source, selected$target, sep = "\r") %in%
        paste(spec$source, spec$target, sep = "\r")
      pair_label <- "focused_cells"
    }
    plot_data <- selected[selected_index, , drop = FALSE]
    if (nrow(plot_data) == 0) next
    filename <- if (i == 1) "compare_lr_regulated.png" else paste0("compare_lr_", pair_label, ".png")
    pair_lr <- data.frame(interaction_name = unique(plot_data$interaction_name))
    p <- netVisual_bubble(
      cellchat_merged,
      sources.use = match(spec_sources, joint_levels),
      targets.use = match(spec_targets, joint_levels),
      pairLR.use = pair_lr,
      comparison = c(1, 2),
      angle.x = opt$bubble_angle,
      remove.isolate = opt$bubble_remove_isolate
    )
    ggsave(file.path(opt$output_dir, filename), p, width = 11, height = 5.5, dpi = 300)
    selected$plot_file[selected_index] <- filename
  }
} else {
  message_plot("compare_lr_regulated.png", "Focused ligand-receptor comparison",
               "No significant ligand-receptor interaction was available in the selected scope")
}
write_table(selected, "selected_interactions.csv")

parameters <- data.frame(
  parameter = c("condition1", "condition2", "difference", "datatype", "cell_identities",
                "selection_source", "selected_cell_pairs", "pathway_selection", "selected_pathways",
                "lr_selection", "selected_lr", "cellchat_version"),
  value = c(name1, name2, paste0(name2, " - ", name1), datatype1,
            paste(joint_levels, collapse = ","), selection_source,
            paste(paste(selected_pairs$source, selected_pairs$target, sep = "->"), collapse = ","),
            pathway_selection, paste(selected_pathways, collapse = ","), lr_selection,
            paste(selected_interaction_names, collapse = ","), as.character(packageVersion("CellChat")))
)
write_table(parameters, "compare_parameters.csv")

gg_count <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "count")
gg_weight <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
ggsave(file.path(opt$output_dir, "compare_overview_number_strength.png"), gg_count + gg_weight,
       width = 6.5, height = 4, dpi = 300)

nonzero_diff <- any(cell_pair_table$count_diff != 0 | abs(cell_pair_table$weight_diff) > 1e-15)
if (nonzero_diff) {
  png(file.path(opt$output_dir, "compare_diff_number_strength_net.png"), width = 10, height = 5,
      units = "in", res = 300)
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
  count_diff <- cellchat2@net$count - cellchat1@net$count
  weight_diff <- cellchat2@net$weight - cellchat1@net$weight
  count_max <- max(abs(count_diff))
  weight_max <- max(abs(weight_diff))
  count_higher2 <- count_diff; count_higher2[count_higher2 < 0] <- 0
  count_higher1 <- -count_diff; count_higher1[count_higher1 < 0] <- 0
  weight_higher2 <- weight_diff; weight_higher2[weight_higher2 < 0] <- 0
  weight_higher1 <- -weight_diff; weight_higher1[weight_higher1 < 0] <- 0
  plot_difference_circle(count_higher2, paste("Count:", name2, "higher"), count_max, length(joint_levels))
  plot_difference_circle(count_higher1, paste("Count:", name1, "higher"), count_max, length(joint_levels))
  plot_difference_circle(weight_higher2, paste("Weight:", name2, "higher"), weight_max, length(joint_levels))
  plot_difference_circle(weight_higher1, paste("Weight:", name1, "higher"), weight_max, length(joint_levels))
  dev.off()

  heatmap_count <- netVisual_heatmap(cellchat_merged, measure = "count")
  heatmap_weight <- netVisual_heatmap(cellchat_merged, measure = "weight")
  png(file.path(opt$output_dir, "compare_heatmap_count_weight.png"), width = 10, height = 6,
      units = "in", res = 300)
  draw(heatmap_count + heatmap_weight, ht_gap = unit(0.7, "cm"))
  dev.off()
} else {
  message_plot("compare_diff_number_strength_net.png", "Differential communication network",
               "No non-zero count or weight difference was detected")
  message_plot("compare_heatmap_count_weight.png", "Differential communication heatmap",
               "No non-zero count or weight difference was detected")
}

if (nrow(pathway_table) > 0) {
  rank_stacked <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
  rank_grouped <- rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
  ggsave(file.path(opt$output_dir, "compare_pathway_strength.png"), rank_stacked + rank_grouped,
         width = 10, height = 6, dpi = 300)
} else {
  message_plot("compare_pathway_strength.png", "Pathway information flow", "No pathway was available")
}

if (opt$plot_advanced && length(selected_pathways) > 0) {
  advanced_pathways <- if (length(valid_manual_pathways) > 0) selected_pathways else selected_pathways[1]
  gene_colors <- split_values(opt$gene_colors)
  cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP")
  cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP")

  for (pathway in advanced_pathways) {
    suffix <- if (length(advanced_pathways) == 1) "" else paste0("_", safe_name(pathway))
    present1 <- pathway %in% cellchat1@netP$pathways
    present2 <- pathway %in% cellchat2@netP$pathways

    png(file.path(opt$output_dir, paste0("compare_pathway_network", suffix, ".png")),
        width = 10, height = 5, units = "in", res = 300)
    par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
    if (present1) {
      netVisual_aggregate(cellchat1, signaling = pathway, layout = "circle",
                          vertex.weight = rep(20, length(levels(cellchat1@idents))))
      mtext(name1, side = 3, line = 0.2, cex = 0.9)
    } else {
      plot.new(); title(main = name1); text(0.5, 0.5, paste(pathway, "not detected"))
    }
    if (present2) {
      netVisual_aggregate(cellchat2, signaling = pathway, layout = "circle",
                          vertex.weight = rep(20, length(levels(cellchat2@idents))))
      mtext(name2, side = 3, line = 0.2, cex = 0.9)
    } else {
      plot.new(); title(main = name2); text(0.5, 0.5, paste(pathway, "not detected"))
    }
    dev.off()

    placeholder <- function(label) ggplot() + theme_void() + annotate("text", x = 0, y = 0, label = label)
    contribution1 <- if (present1) netAnalysis_contribution(cellchat1, signaling = pathway) else placeholder(paste(pathway, "not detected in", name1))
    contribution2 <- if (present2) netAnalysis_contribution(cellchat2, signaling = pathway) else placeholder(paste(pathway, "not detected in", name2))
    ggsave(file.path(opt$output_dir, paste0("compare_pathway_contribution", suffix, ".png")),
           contribution1 + contribution2, width = 10, height = 5, dpi = 300)

    gene1 <- if (present1) plotGeneExpression(cellchat1, signaling = pathway, type = opt$gene_plot_type,
                                               color.use = gene_colors) else placeholder(paste(pathway, "not detected in", name1))
    gene2 <- if (present2) plotGeneExpression(cellchat2, signaling = pathway, type = opt$gene_plot_type,
                                               color.use = gene_colors) else placeholder(paste(pathway, "not detected in", name2))
    ggsave(file.path(opt$output_dir, paste0("compare_gene_expression", suffix, ".png")),
           gene1 + gene2, width = 11, height = 6, dpi = 300)

    role_plots <- list()
    if (present1 && role_heatmap_has_range(cellchat1, pathway)) {
      role_plots[[length(role_plots) + 1]] <- wrap_elements(full = grid.grabExpr(draw(
        netAnalysis_signalingRole_heatmap(cellchat1, signaling = pathway, pattern = "all", title = name1)
      )))
    } else {
      role_plots[[length(role_plots) + 1]] <- placeholder(paste("Insufficient signaling-role variation for", pathway, "in", name1))
    }
    if (present2 && role_heatmap_has_range(cellchat2, pathway)) {
      role_plots[[length(role_plots) + 1]] <- wrap_elements(full = grid.grabExpr(draw(
        netAnalysis_signalingRole_heatmap(cellchat2, signaling = pathway, pattern = "all", title = name2)
      )))
    } else {
      role_plots[[length(role_plots) + 1]] <- placeholder(paste("Insufficient signaling-role variation for", pathway, "in", name2))
    }
    ggsave(file.path(opt$output_dir, paste0("compare_signaling_role", suffix, ".png")),
           role_plots[[1]] + role_plots[[2]], width = 11, height = 6, dpi = 300)
  }
}

if (nzchar(opt$receiver_cells) || opt$do_single_bubble || nzchar(opt$pair_lr_use)) {
  warning("receiver_cells, do_single_bubble, and pair_lr_use are deprecated in compare-stage and were ignored")
}

if (opt$save_merged) {
  saveRDS(setNames(list(cellchat1, cellchat2), c(name1, name2)),
          file.path(opt$output_dir, "cellchat_object.list.rds"))
  saveRDS(cellchat_merged, file.path(opt$output_dir, "cellchat_merged.rds"))
}

message("Selected cell pairs: ", paste(paste(selected_pairs$source, selected_pairs$target, sep = "->"), collapse = ", "))
message("Selected pathways: ", paste(selected_pathways, collapse = ", "))
message("Selected LR interactions: ", paste(selected_interaction_names, collapse = ", "))
message("CellChat comparison completed: ", opt$output_dir)
