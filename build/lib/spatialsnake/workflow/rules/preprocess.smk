def preprocess_input(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{sample}_{bin}um", 'integrate',"{sample}.zarr"))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      return(os.path.join(results_folder, "{sample}", 'integrate',"{sample}.zarr"))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'integrate',"{sample}.h5ad"))
  if channel=="compare_analysis":
    return(os.path.join(results_folder, "merge_data", "concatenated_sdata"))

def preprocess_output(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{sample}_{bin}um", 'preprocess',"filter_{sample}.zarr")))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      print("correct")
      return(directory(os.path.join(results_folder, "{sample}",'preprocess', "filter_{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'preprocess',"filter_{sample}.h5ad"))
  if channel=="compare_analysis":
    return(parameter_output(samples,'preprocess'))



rule preprocess_rule:
  input:
    outputs=preprocess_input(run_type)
  output:
    merge=preprocess_output(run_type) if run_type == "slide_seq" else directory(preprocess_output(run_type))
  params:
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" else wildcards.sample,
    run_type = run_type,
    variable = variable,
    NEIGHBORS = NEIGHBORS,
    batch_method = batch_method,
    mt_threshold = mt_threshold,
    n_top_genes = n_top_genes,
    n_comps=n_comps,
    sketch=sketch,
    sample_rate=sample_rate,
    filter_dict= f"--filter_dict '{json.dumps(filter_dict)}'" if filter_list else "",
    min_cells= f"--min_cells {min_cells}" if not filter_list else "",
    min_genes= f"--min_genes {min_genes}" if not filter_list else ""
  shell:
      """
      python {spatialsnake_path}workflow/scripts/preprocessing.py \
        --input_dir {input.outputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --variable {params.variable} \
        --NEIGHBORS {params.NEIGHBORS} \
        --batch_method {params.batch_method} \
        --mt_threshold {params.mt_threshold} \
        --n_top_genes {params.n_top_genes} \
        --n_comps {params.n_comps} \
        --sketch {params.sketch} \
        --sample_rate {params.sample_rate} \
        {params.filter_dict} \
        {params.min_genes} \
        {params.min_cells}
      """


