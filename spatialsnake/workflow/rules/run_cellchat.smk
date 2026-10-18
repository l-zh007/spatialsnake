rule cellchat_compute:
  input:
    spatial=input_spatial
  output:
    rds=os.path.join(results_folder, "{sample}", "cellchat", "cellchat.rds"),
    stats=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_stats.csv"),
    lr=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_lr.csv")
  threads: workflow_threads
  params:
    output_dir=lambda wildcards: os.path.join(results_folder, wildcards.sample, "cellchat"),
    run_type=run_type,
    celltype_col=config.get("celltype_col", "celltype"),
    assay=config.get("cellchat_assay", config.get("assay", "Spatial")),
    species=config.get("cellchat_species", config.get("species", "human")),
    min_cells=config.get("cellchat_min_cells", config.get("min_cells", 10)),
    is_single_cell=config.get("cellchat_is_single_cell", config.get("is_single_cell", False)),
    trim=config.get("cellchat_trim", config.get("trim", 0.1)),
    interaction_length=config.get("cellchat_interaction_length", config.get("interaction_length", 250)),
    scale_factors=lambda wildcards: (
      scale_factors_files[0]
      if len(scale_factors_files) > 0 and normalize_run_type_name(run_type) in {"visium", "visiumhd", "visiumsegment", "stereoseq"}
      else ""
    ),
    scale_factors_list=",".join([str(x).strip() for x in cellchat_scale_factors if str(x).strip() != ""]) if len(cellchat_scale_factors) > 0 else "",
    sample_names=",".join(cellchat_sample_names) if len(cellchat_sample_names) > 0 else "",
    spot_size=config.get("cellchat_spot_size", 65),
    db_subset=config.get("cellchat_db_subset", "all_interactions"),
    future_max_size_gb=config.get("cellchat_future_max_size_gb", config.get("future_max_size_gb", 64))
  shell:
      """
      Rscript {spatialsnake_path}workflow/scripts/Cellchat.R \
        --stage compute \
        --input_path {input.spatial} \
        --output_dir {params.output_dir} \
        --run_type {params.run_type} \
        --sample_id {wildcards.sample} \
        --celltype_col {params.celltype_col} \
        --assay {params.assay} \
        --species {params.species} \
        --min_cells {params.min_cells} \
        --nworkers {threads} \
        --is_single_cell {params.is_single_cell} \
        --scale_factors "{params.scale_factors}" \
        --scale_factors_list "{params.scale_factors_list}" \
        --sample_names "{params.sample_names}" \
        --spot_size {params.spot_size} \
        --db_subset {params.db_subset} \
        --trim {params.trim} \
        --interaction_length {params.interaction_length} \
        --future_max_size_gb {params.future_max_size_gb}
      """

rule cellchat_plot:
  input:
    rds=rules.cellchat_compute.output.rds
  output:
    network_png=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_network.png"),
    network_pdf=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_network.pdf"),
    heatmap_png=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_heatmap.png"),
    infoflow_png=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_infoflow_bar.png"),
    selected_lr_csv=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_selected_lr.csv"),
    selected_pathway_csv=os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_selected_pathway_summary.csv")
  threads: workflow_threads
  params:
    output_dir=lambda wildcards: os.path.join(results_folder, wildcards.sample, "cellchat"),
    is_single_cell=config.get("cellchat_is_single_cell", config.get("is_single_cell", False)),
    pathways=config.get("cellchat_pathways", ""),
    top_pathways=config.get("cellchat_top_pathways", 3),
    focus_cells=config.get("cellchat_focus_cells", ""),
    source_cells=config.get("cellchat_source_cells", ""),
    target_cells=config.get("cellchat_target_cells", ""),
    pair_lr_use=config.get("cellchat_pair_lr_use", ""),
    cell_pairs=config.get("cellchat_cell_pairs", ""),
    lr_pairs=config.get("cellchat_lr_pairs", ""),
    top_cell_pairs=config.get("cellchat_top_cell_pairs", 3),
    bubble_top_lr=config.get("cellchat_bubble_top_lr", 20),
    plot_advanced=config.get("cellchat_plot_advanced", True)
  shell:
      """
      Rscript {spatialsnake_path}workflow/scripts/Cellchat.R \
        --stage plot \
        --input_path {input.rds} \
        --output_dir {params.output_dir} \
        --sample_id {wildcards.sample} \
        --is_single_cell {params.is_single_cell} \
        --pathways "{params.pathways}" \
        --top_pathways {params.top_pathways} \
        --focus_cells "{params.focus_cells}" \
        --source_cells "{params.source_cells}" \
        --target_cells "{params.target_cells}" \
        --pair_lr_use "{params.pair_lr_use}" \
        --cell_pairs "{params.cell_pairs}" \
        --lr_pairs "{params.lr_pairs}" \
        --top_cell_pairs {params.top_cell_pairs} \
        --bubble_top_lr {params.bubble_top_lr} \
        --plot_advanced {params.plot_advanced}
      """
