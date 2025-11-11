def get_output(type):
  if seg_filter==False:
    return(os.path.join(results_folder, "merge_data", "concatenated_sdata"))
  else:
    return(os.path.join(results_folder, "merge_data", 'preprocess',"filter_concatenated_sdata"))

rule merge_in:
  input:
    outputs=parameter_output(samples,'integrate')
  output:
    merge=get_output(type) if run_type == "slide_seq" else directory(get_output(type))
  params:
    main_file = main_file,
    run_type = run_type,
    group = group,
    sample =samples,
    seg_filter = seg_filter
  shell:
      """
      python {spatialsnake_path}workflow/scripts/spatial_in_multiple.py \
        --input_path {input.outputs} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --sample_id {params.sample} \
        --group {params.group} \
        --seg_filter {params.seg_filter}
      """
    
  
