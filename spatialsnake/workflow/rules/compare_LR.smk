rule compare_cellchat_compute:
  input:
    rds1=downstream_file[0],
    rds2=downstream_file[1]
  output:
    global_summary=os.path.join(cellchat_compare_output_dir, "compare_global_summary.csv"),
    cell_pair=os.path.join(cellchat_compare_output_dir, "compare_cell_pair.csv"),
    pathway=os.path.join(cellchat_compare_output_dir, "compare_pathway.csv"),
    lr=os.path.join(cellchat_compare_output_dir, "compare_lr.csv"),
    cache=[os.path.join(cellchat_compare_output_dir, name) for name in ("cellchat_object.list.rds", "cellchat_merged.rds")] if cellchat_compare_save_merged else []
  threads: workflow_threads
  params:
    sample_name1=cellchat_compare_sample_name1 if cellchat_compare_sample_name1 else (samples[0] if len(samples) > 0 else "sample1"),
    sample_name2=cellchat_compare_sample_name2 if cellchat_compare_sample_name2 else (samples[1] if len(samples) > 1 else "sample2"),
    save_merged=cellchat_compare_save_merged,
    output_dir=cellchat_compare_output_dir
  shell:
      """
      Rscript "{spatialsnake_path}workflow/scripts/cellchat_compare.R" \
        --stage compute \
        --input_rds1 "{input.rds1}" \
        --input_rds2 "{input.rds2}" \
        --sample_name1 "{params.sample_name1}" \
        --sample_name2 "{params.sample_name2}" \
        --output_dir "{params.output_dir}" \
        --save_merged {params.save_merged}
      """

rule compare_cellchat_plot:
  input:
    rds1=downstream_file[0],
    rds2=downstream_file[1],
    global_summary=rules.compare_cellchat_compute.output.global_summary,
    cell_pair=rules.compare_cellchat_compute.output.cell_pair,
    pathway=rules.compare_cellchat_compute.output.pathway,
    lr=rules.compare_cellchat_compute.output.lr,
    cache=rules.compare_cellchat_compute.output.cache
  output:
    selected=os.path.join(cellchat_compare_output_dir, "selected_interactions.csv"),
    parameters=os.path.join(cellchat_compare_output_dir, "compare_parameters.csv"),
    overview=os.path.join(cellchat_compare_output_dir, "compare_overview_number_strength.png"),
    diff_net=os.path.join(cellchat_compare_output_dir, "compare_diff_number_strength_net.png"),
    heatmap=os.path.join(cellchat_compare_output_dir, "compare_heatmap_count_weight.png"),
    pathway_strength=os.path.join(cellchat_compare_output_dir, "compare_pathway_strength.png"),
    focused_lr=os.path.join(cellchat_compare_output_dir, "compare_lr_regulated.png")
  threads: workflow_threads
  params:
    sample_name1=cellchat_compare_sample_name1 if cellchat_compare_sample_name1 else (samples[0] if len(samples) > 0 else "sample1"),
    sample_name2=cellchat_compare_sample_name2 if cellchat_compare_sample_name2 else (samples[1] if len(samples) > 1 else "sample2"),
    output_dir=cellchat_compare_output_dir,
    focus_cells=cellchat_compare_focus_cells,
    cell_pairs=cellchat_compare_cell_pairs,
    pathways=cellchat_compare_pathways,
    source_cells=cellchat_compare_source_cells,
    target_cells=cellchat_compare_target_cells,
    lr_pairs=cellchat_compare_lr_pairs,
    top_cell_pairs=cellchat_compare_top_cell_pairs,
    top_pathways=cellchat_compare_top_pathways,
    top_lr=cellchat_compare_top_lr,
    plot_advanced=cellchat_compare_plot_advanced,
    receiver_cells=cellchat_compare_receiver_cells,
    bubble_angle=cellchat_compare_bubble_angle,
    bubble_remove_isolate=cellchat_compare_bubble_remove_isolate,
    do_single_bubble=cellchat_compare_do_single_bubble,
    gene_colors=cellchat_compare_gene_colors,
    gene_plot_type=cellchat_compare_gene_plot_type,
    pair_lr_use=cellchat_compare_pair_lr_use,
    save_merged=cellchat_compare_save_merged
  shell:
      """
      Rscript "{spatialsnake_path}workflow/scripts/cellchat_compare.R" \
        --stage plot \
        --input_rds1 "{input.rds1}" \
        --input_rds2 "{input.rds2}" \
        --sample_name1 "{params.sample_name1}" \
        --sample_name2 "{params.sample_name2}" \
        --output_dir "{params.output_dir}" \
        --focus_cells "{params.focus_cells}" \
        --cell_pairs "{params.cell_pairs}" \
        --pathways "{params.pathways}" \
        --source_cells "{params.source_cells}" \
        --target_cells "{params.target_cells}" \
        --lr_pairs "{params.lr_pairs}" \
        --top_cell_pairs {params.top_cell_pairs} \
        --top_pathways {params.top_pathways} \
        --top_lr {params.top_lr} \
        --plot_advanced {params.plot_advanced} \
        --receiver_cells "{params.receiver_cells}" \
        --bubble_angle {params.bubble_angle} \
        --bubble_remove_isolate {params.bubble_remove_isolate} \
        --do_single_bubble {params.do_single_bubble} \
        --gene_colors "{params.gene_colors}" \
        --gene_plot_type "{params.gene_plot_type}" \
        --pair_lr_use "{params.pair_lr_use}" \
        --save_merged {params.save_merged}
      """
