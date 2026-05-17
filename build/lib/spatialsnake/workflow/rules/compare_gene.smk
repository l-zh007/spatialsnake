rule cluster_rule:
  input:
    inputs=nomal_file(run_type,"annotation")
  output:
    merge=os.path.join(results_folder,"merge_data",'compare_analysis','marker_genes_pval.csv') if compare_algorithm=="DEseq2" else os.path.join(results_folder,"merge_data",'compare_analysis','edgeR_counts.csv')
  params:
    sample_id = samples,
    run_type = run_type,
    sample_cnt = len(samples),
    compare_algorithm = compare_algorithm,
    cell_focus=cell_focus,
    group=group
  shell:
      """
      python {spatialsnake_path}workflow/scripts/DEseq2.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --algorithm {params.compare_algorithm} \
        --cell_focus {params.cell_focus} \
        --sample_cnt {params.sample_cnt} \
        --group {params.group}
      """


rule enrich_rule:
  input:
    inputs = os.path.join(results_folder,"merge_data",'compare_analysis','marker_genes_pval.csv') if compare_algorithm=="DEseq2" else os.path.join(results_folder,"merge_data",'compare_analysis','edgeR_counts.csv')
  output:
    merge = os.path.join(results_folder,"merge_data",'compare_analysis','positive','kegg_data.csv')
  params:
    spacies = spacies,
    run_type = run_type,
    compare_algorithm = compare_algorithm,
    sample_list = sample_list,
    cut_off_pvalue = config.get("cut_off_pvalue"),
    cut_off_logFC = config.get("cut_off_logFC")
  shell:
      """
      Rscript {spatialsnake_path}workflow/scripts/diffent_analysis.R \
        --input_dir {input.inputs} \
        --output_path {output.merge} \
        --type {params.run_type} \
        --spacies {spacies} \
        --algorithm {params.compare_algorithm} \
        --sample_list {params.sample_list} \
        --cut_off_pvalue {params.cut_off_pvalue} \
        --cut_off_logFC {params.cut_off_logFC}
      """
