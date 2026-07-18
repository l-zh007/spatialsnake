def get_banksy_dir():
  if channel == "compare_analysis":
    return os.path.join(results_folder, "merge_data", "banksy")
  return os.path.join(results_folder, f"{banksy_sample_id}", "banksy")

def get_banksy_output():
  return directory(os.path.join(get_banksy_dir(), f"{banksy_sample_id}_banksy.zarr"))

rule banksy_rule:
  input:
    inputs=lambda wildcards: banksy_input
  output:
    merge=get_banksy_output(),
    results=os.path.join(get_banksy_dir(), "banksy_results", "banksy_results.csv"),
    assignments=os.path.join(get_banksy_dir(), "banksy_results", "banksy_cluster_assignments.csv"),
    parameters=os.path.join(get_banksy_dir(), "banksy_results", "banksy_parameters.csv"),
    spatial_plot=os.path.join(get_banksy_dir(), "banksy_results", "banksy_spatial_cluster.png")
  threads: workflow_threads
  params:
    k_geom = k_geom,
    max_m = max_m,
    nbr_weight_decay = nbr_weight_decay,
    banksy_n_comps = banksy_n_comps,
    lambda_list = lambda_list,
    banksy_resolution = banksy_resolution,
    banksy_num_nn = banksy_num_nn,
    banksy_max_features = banksy_max_features,
    banksy_feature_col = banksy_feature_col,
    banksy_add_umap = banksy_add_umap,
    banksy_plot_full = banksy_plot_full,
    banksy_run_nonspatial = banksy_run_nonspatial,
    banksy_plot_celltype_enrichment = banksy_plot_celltype_enrichment,
    banksy_plot_max_points = banksy_plot_max_points,
    banksy_sample_col = banksy_sample_col,
    banksy_selected_lambda = banksy_selected_lambda,
    banksy_selected_resolution = banksy_selected_resolution,
    banksy_seed = banksy_seed
  shell:
      """
      python {spatialsnake_path}workflow/scripts/run_banksy.py \
        --input_dir "{input.inputs}" \
        --output_zarr_path "{output.merge}" \
        --k_geom "{params.k_geom}" \
        --max_m "{params.max_m}" \
        --nbr_weight_decay "{params.nbr_weight_decay}" \
        --banksy_n_comps "{params.banksy_n_comps}" \
        --lambda_list "{params.lambda_list}" \
        --banksy_resolution "{params.banksy_resolution}" \
        --banksy_num_nn "{params.banksy_num_nn}" \
        --banksy_max_features "{params.banksy_max_features}" \
        --banksy_feature_col "{params.banksy_feature_col}" \
        --banksy_add_umap "{params.banksy_add_umap}" \
        --banksy_plot_full "{params.banksy_plot_full}" \
        --banksy_run_nonspatial "{params.banksy_run_nonspatial}" \
        --banksy_plot_celltype_enrichment "{params.banksy_plot_celltype_enrichment}" \
        --banksy_plot_max_points "{params.banksy_plot_max_points}" \
        --banksy_sample_col "{params.banksy_sample_col}" \
        --banksy_selected_lambda "{params.banksy_selected_lambda}" \
        --banksy_selected_resolution "{params.banksy_selected_resolution}" \
        --banksy_seed "{params.banksy_seed}"
      """
