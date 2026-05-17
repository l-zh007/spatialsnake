rule cellchat_rule:
  input:
    spatial = input_spatial
  output:
    network_png = os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_network.png"),
    network_pdf = os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_network.pdf"),
    stats_csv = os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_stats.csv"),
    lr_csv = os.path.join(results_folder, "{sample}", "cellchat", "{sample}_cellchat_lr.csv")
  params:
    output_dir = lambda wildcards: os.path.join(results_folder, wildcards.sample, "cellchat"),
    sample_id = samples,
    run_type = run_type,
    celltype_col = config.get("celltype_col", "celltype"),
    assay = config.get("cellchat_assay", "Spatial"),
    species = config.get("cellchat_species", "human"),
    min_cells = config.get("cellchat_min_cells", 10),
    nworkers = config.get("cellchat_workers", 4),
    is_single_cell = config.get("cellchat_is_single_cell", False),
    trim = config.get("cellchat_trim", 0.1),
    interaction_length = config.get("cellchat_interaction_length", 150),
    scale_factors = lambda wildcards: (
      scale_factors_files[0]
      if len(scale_factors_files) > 0 and normalize_run_type_name(run_type) in {"visium", "visiumhd", "visiumsegment", "stereoseq"}
      else ""
    ),
    scale_factors_list = ",".join([str(x).strip() for x in cellchat_scale_factors if str(x).strip() != ""]) if len(cellchat_scale_factors) > 0 else "",
    sample_names = ",".join(cellchat_sample_names) if len(cellchat_sample_names) > 0 else "",
    spot_size = config.get("cellchat_spot_size", 65)
  shell:
      """
      Rscript {spatialsnake_path}workflow/scripts/Cellchat.R \
        --input_path {input.spatial} \
        --output_dir {params.output_dir} \
        --run_type {params.run_type} \
        --sample_id {wildcards.sample} \
        --celltype_col {params.celltype_col} \
        --assay {params.assay} \
        --species {params.species} \
        --min_cells {params.min_cells} \
        --nworkers {params.nworkers} \
        --is_single_cell {params.is_single_cell} \
        --scale_factors "{params.scale_factors}" \
        --scale_factors_list "{params.scale_factors_list}" \
        --sample_names "{params.sample_names}" \
        --spot_size {params.spot_size} \
        --trim {params.trim} \
        --interaction_length {params.interaction_length}
      """
