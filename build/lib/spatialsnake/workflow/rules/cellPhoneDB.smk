import os
cpdb_out_dir = os.path.join(results_folder, cpdb_sample_id, "cellPhoneDB_results") if channel=="single_analysis" else os.path.join(results_folder, "merge_data", "cellPhoneDB_results")
cpdb_heatmap_file = f"{cpdb_sample_id}_heatmap.png" if channel=="single_analysis" else "concentrate_heatmap.png"

rule cellPhoneDB_rule:
  input:
    inputs=cellPhoneDB_input
  output:
    merge=os.path.join(cpdb_out_dir,"adata.h5ad")
  params:
    sample_id = cpdb_sample_id,
    run_type = run_type,
    counts_data=counts_data,
    threshold= threshold,
    threads = threads,
    pvalue = pvalue,
    output_name = output_name,
    degs_file_path=degs_file_path,
    microenvs_file_path=microenvs_file_path,
    active_tf_path=active_tf_path,
    celltype_col=celltype_col,
    niche_col=niche_col,
    is_singlecell=is_singlecell,
    method=cpdb_method,
    de_method=cpdb_de_method,
    iterations=iterations
  shell:
      """
      active_tf_arg=""
      microenvs_arg=""
      degs_arg=""
      if [ -n "{params.active_tf_path}" ]; then
        active_tf_arg="--active_tf_path {params.active_tf_path}"
      fi
      if [ -n "{params.microenvs_file_path}" ]; then
        microenvs_arg="--microenvs_file_path {params.microenvs_file_path}"
      fi
      if [ -n "{params.degs_file_path}" ]; then
        degs_arg="--degs_file_path {params.degs_file_path}"
      fi
      python {spatialsnake_path}workflow/scripts/cellphoneDB.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --counts_data {params.counts_data} \
        --celltype_col {params.celltype_col} \
        --niche_col {params.niche_col} \
        --is_singlecell {params.is_singlecell} \
        --method {params.method} \
        --de_method {params.de_method} \
        --iterations {params.iterations} \
        --threshold {params.threshold} \
        --threads {params.threads} \
        --pvalue {params.pvalue} \
        --output_name {params.output_name} \
        $active_tf_arg \
        $microenvs_arg \
        $degs_arg
      """
rule visualize_cellPhone:
  input:
    inputs=os.path.join(cpdb_out_dir,"adata.h5ad")
  output:
    merge=os.path.join(cpdb_out_dir,cpdb_heatmap_file)
  params:
    sample_id = cpdb_sample_id,
    run_type = run_type,
    counts_data=counts_data,
    output_name = output_name,
    cell_type1=cell_type1,
    cell_type2=cell_type2,
    celltype=celltype,
    degs_file_path=degs_file_path,
    gene_family = gene_family
  shell:
      """
      degs_arg=""
      gene_arg=""
      if [ -n "{params.gene_family}" ]; then
        gene_arg="--gene_family {params.gene_family}"
      fi
      if [ -n "{params.degs_file_path}" ]; then
        degs_arg="--degs_file_path {params.degs_file_path}"
      fi
      python {spatialsnake_path}workflow/scripts/cellphoneDB_visualize.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --output_name {params.output_name} \
        $degs_arg \
        $gene_arg \
        --celltype {params.celltype} \
        --cell_type1 {params.cell_type1} \
        --cell_type2 {params.cell_type2} \
      """
