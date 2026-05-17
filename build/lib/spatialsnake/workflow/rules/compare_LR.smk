rule compare_cellchat:
  input:
    rds1 = downstream_file[0]
  output:
    outdir = directory(cellchat_compare_output_dir)
  params:
    rds2 = downstream_file[1] if len(downstream_file) > 1 else "",
    sample_name1 = cellchat_compare_sample_name1 if cellchat_compare_sample_name1 else (samples[0] if len(samples) > 0 else "sample1"),
    sample_name2 = cellchat_compare_sample_name2 if cellchat_compare_sample_name2 else (samples[1] if len(samples) > 1 else ""),
    pathways = cellchat_compare_pathways,
    source_cells = cellchat_compare_source_cells,
    target_cells = cellchat_compare_target_cells,
    receiver_cells = cellchat_compare_receiver_cells,
    bubble_angle = cellchat_compare_bubble_angle,
    bubble_remove_isolate = cellchat_compare_bubble_remove_isolate,
    do_ranknet = cellchat_compare_do_ranknet,
    do_role_heatmap = cellchat_compare_do_role_heatmap,
    do_pathway_plots = cellchat_compare_do_pathway_plots,
    do_compare_overview = cellchat_compare_do_compare_overview,
    do_compare_bubble = cellchat_compare_do_compare_bubble,
    do_single_bubble = cellchat_compare_do_single_bubble,
    do_gene_expression = cellchat_compare_do_gene_expression,
    gene_colors = cellchat_compare_gene_colors,
    gene_plot_type = cellchat_compare_gene_plot_type,
    pair_lr_use = cellchat_compare_pair_lr_use,
    save_merged = cellchat_compare_save_merged
  shell:
      """
      Rscript {spatialsnake_path}workflow/scripts/cellchat_compare.R \
        --input_rds1 {input.rds1} \
        --input_rds2 "{params.rds2}" \
        --sample_name1 "{params.sample_name1}" \
        --sample_name2 "{params.sample_name2}" \
        --output_dir {output.outdir} \
        --pathways "{params.pathways}" \
        --source_cells "{params.source_cells}" \
        --target_cells "{params.target_cells}" \
        --receiver_cells "{params.receiver_cells}" \
        --bubble_angle {params.bubble_angle} \
        --bubble_remove_isolate {params.bubble_remove_isolate} \
        --do_ranknet {params.do_ranknet} \
        --do_role_heatmap {params.do_role_heatmap} \
        --do_pathway_plots {params.do_pathway_plots} \
        --do_compare_overview {params.do_compare_overview} \
        --do_compare_bubble {params.do_compare_bubble} \
        --do_single_bubble {params.do_single_bubble} \
        --do_gene_expression {params.do_gene_expression} \
        --gene_colors "{params.gene_colors}" \
        --gene_plot_type "{params.gene_plot_type}" \
        --pair_lr_use "{params.pair_lr_use}" \
        --save_merged {params.save_merged}
      """
