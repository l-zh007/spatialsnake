rule cellPhoneDB_rule:
  input:
    inputs=cellPhoneDB_input
  output:
    merge=os.path.join(results_folder,'cellphonedb_output',"adata.h5ad")
  params:
    sample_id = lambda wildcards: "concatenated_sdata" if channel=="compare_analysis" else wildcards.sample,
    run_type = run_type,
    counts_data=counts_data,
    threshold= threshold,
    threads = threads,
    pvalue = pvalue,
    output_name = output_name,
    degs_file_path=degs_file_path,
    microenvs_file_path=microenvs_file_path,
    active_tf_path=active_tf_path
  shell:
      """
      python {spatialsnake_path}workflow/scripts/cellphoneDB.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --counts_data {params.counts_data} \
        --threshold {params.threshold} \
        --threads {params.threads} \
        --pvalue {params.pvalue} \
        --output_name {params.output_name} \
        --active_tf_path {params.active_tf_path} \
        --microenvs_file_path {params.microenvs_file_path} \
        --degs_file_path {params.degs_file_path}
      """
rule visualize_cellPhone:
  input:
    inputs=os.path.join(results_folder,'cellphonedb_output',"adata.h5ad")
  output:
    merge=os.path.join(results_folder,'cellphonedb_output',"concentrate_heatmap.png")
  params:
    sample_id = sample_id
    run_type = run_type,
    counts_data=counts_data,
    output_name = output_name,
    cell_type1=cell_type1,
    cell_type2=cell_type2,
    gene_family=gene_family,
    celltype=celltype,
    degs_file_path=degs_file_path
  shell:
      """
      python {spatialsnake_path}workflow/scripts/cellphoneDB_visualize.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --output_name {params.output_name} \
        --degs_file_path {params.degs_file_path} \
        --celltype {params.celltype} \
        --gene_family {params.gene_family} \
        --cell_type1 {params.cell_type1} \
        --cell_type2 {params.cell_type2} \
      """
