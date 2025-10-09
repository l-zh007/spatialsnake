rule advance_analysis_rule:
  input:
    inputs=nomal_file(run_type,"annotion")
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
