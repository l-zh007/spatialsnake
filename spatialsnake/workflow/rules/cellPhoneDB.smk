rule cellPhoneDB_rule:
  input:
    inputs=nomal_file(run_type,"annotion")
  output:
    merge=os.path.join(results_folder,"merge_data",'cellPhoneDB',"concatenated_sdata.h5ad") if channel=="compare_analysis" else os.path.join(results_folder,"{sample}",'cellPhoneDB',"{sample}.h5ad")
  params:
    sample_id = lambda wildcards: "concatenated_sdata" if channel=="compare_analysis" else wildcards.sample,
    run_type = run_type,
    counts_data=counts_data,
    threshold= threshold,
    threads = threads,
    pvalue = pvalue,
    output_name = output_name
  shell:
      """
      python workflow/scripts/cellphoneDB.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --counts_data {params.counts_data} \
        --threshold {params.threshold} \
        --threads {params.threads} \
        --pvalue {params.pvalue} \
        --output_name {params.output_name}
      """
rule visualize_cellPhone:
  input:
    inputs=os.path.join(results_folder,"merge_data",'cellPhoneDB',"concatenated_sdata.h5ad") if channel=="compare_analysis" else os.path.join(results_folder,"{sample}",'cellPhoneDB',"{sample}.h5ad")
  output:
    merge=os.path.join(results_folder,"merge_data",'cellPhoneDB',"concentrate_heatmap.png") if channel=="compare_analysis" else os.path.join(results_folder,"{sample}",'cellPhoneDB',"heatmap.png")
  params:
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" and seg_filter==False else wildcards.sample,
    run_type = run_type,
    counts_data=counts_data,
    threshold= threshold,
    threads = threads,
    pvalue = pvalue,
    output_name = output_name
  shell:
      """
      python workflow/scripts/cellphoneDB_visualize.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --output_name {params.output_name}
      """
