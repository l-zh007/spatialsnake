rule liana_rule:
  input:
    liana_inputs=liana_inputs
  output:
    merge=directory(os.path.join(results_folder,f"{cpdb_sample_id}","liana_output",f"{cpdb_sample_id}.zarr")),
    results=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","liana_results.csv"),
    top=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","liana_top_interactions.csv"),
    plot_data=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","liana_plot_data.csv"),
    pair_summary=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","liana_pair_summary.csv"),
    dotplot=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","dotplot.png"),
    tileplot=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","tileplot.png"),
    heatmap=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","communication_heatmap.png"),
    gene_expression=os.path.join(results_folder,f"{cpdb_sample_id}","liana_output","gene_expression.png")
  threads: workflow_threads
  params:
    run_type = run_type,
    sample_id = cpdb_sample_id,
    liana_method = liana_method,
    liana_resource_name = liana_resource_name,
    liana_expr_prop = liana_expr_prop,
    liana_min_cells = liana_min_cells,
    liana_use_raw = liana_use_raw,
    liana_top_n = config.get("liana_top_n", 6),
    liana_pvalue = liana_pvalue,
    liana_source_celltypes = liana_source_celltypes,
    liana_target_celltypes = liana_target_celltypes,
    liana_cell_pairs = liana_cell_pairs,
    liana_pairs = liana_pairs,
    celltype=celltype
  shell:
      """
      MPLCONFIGDIR=/tmp/mplconfig NUMBA_CACHE_DIR=/tmp/numba_cache \
      python {spatialsnake_path}workflow/scripts/run_liana.py \
        --input_dir "{input.liana_inputs}" \
        --sample_id "{params.sample_id}" \
        --output_zarr_path "{output.merge}" \
        --type "{params.run_type}" \
        --method "{params.liana_method}" \
        --resource_name "{params.liana_resource_name}" \
        --celltype "{params.celltype}" \
        --expr_prop "{params.liana_expr_prop}" \
        --min_cells "{params.liana_min_cells}" \
        --use_raw "{params.liana_use_raw}" \
        --top_n "{params.liana_top_n}" \
        --pvalue "{params.liana_pvalue}" \
        --source_celltypes "{params.liana_source_celltypes}" \
        --target_celltypes "{params.liana_target_celltypes}" \
        --cell_pairs "{params.liana_cell_pairs}" \
        --pairs "{params.liana_pairs}" \
        --results_csv "{output.results}" \
        --top_csv "{output.top}" \
        --plot_data_csv "{output.plot_data}" \
        --pair_summary_csv "{output.pair_summary}"
      """
      
      
      
      
      
      
# rule visualize_liana:
#   input:
#     inputs=os.path.join(results_folder,'cellphonedb_output',f"statistical_analysis_pvalues_{output_name}.txt")
#   output:
#     merge=os.path.join(results_folder,'cellphonedb_output',"concentrate_heatmap.png")
#   params:
#     sample_id = sample_id
#     run_type = run_type,
#     counts_data=counts_data,
#     threshold= threshold,
#     threads = threads,
#     pvalue = pvalue,
#     output_name = output_name
#   shell:
#       """
#       python {spatialsnake_path}workflow/scripts/cellphoneDB_visualize.py \
#         --input_dir {input.inputs} \
#         --sample_id {params.sample_id} \
#         --output_zarr_path {output.merge} \
#         --type {params.run_type} \
#         --output_name {params.output_name}
#       """
