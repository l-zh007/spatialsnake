rule liana_rule:
  input:
    inputs=liana_inputs
  output:
    merge=os.path.join(results_folder,'liana_output',"liana_data.h5ad")
  params:
    run_type = run_type,
    bandwidth=bandwidth,
    cutoff=cutoff,
    resource_name=resource_name,
    celltype=celltype,
    expr_prop=expr_prop,
    sp_cellchat=sp_cellchat
  shell:
      """
      python {spatialsnake_path}workflow/scripts/liana.py \
        --input_dir {input.liana_inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --bandwidth {params.bandwidth} \
        --cutoff {params.cutoff} \
        --resource_name {params.resource_name} \
        --celltype {params.celltype} \
        --expr_prop {params.expr_prop} \
        --sp_cellchat {params.sp_cellchat}
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
