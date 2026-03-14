rule liana_rule:
  input:
    liana_inputs=liana_inputs
  output:
    merge=os.path.join(results_folder,"liana_output",f"{cpdb_sample_id}.zarr")
  params:
    run_type = run_type,
    sample_id = cpdb_sample_id,
    liana_method = liana_method,
    liana_resource_name = liana_resource_name,
    liana_expr_prop = liana_expr_prop,
    liana_min_cells = liana_min_cells,
    liana_use_raw = liana_use_raw,
    celltype=celltype
  shell:
      """
      python {spatialsnake_path}workflow/scripts/liana.py \
        --input_dir {input.liana_inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --method {params.liana_method} \
        --resource_name {params.liana_resource_name} \
        --celltype {params.celltype} \
        --expr_prop {params.liana_expr_prop} \
        --min_cells {params.liana_min_cells} \
        --use_raw {params.liana_use_raw}
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
