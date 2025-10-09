def preprocess_input(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{sample}_{bin}um", "{sample}.zarr"))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      return(os.path.join(results_folder, "{sample}", "{sample}.zarr"))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", "{sample}.h5ad"))
  if channel=="compare_analysis" and seg_filter==False:
    return(os.path.join(results_folder, "merge_data", "concatenated_sdata"))
  elif channel=="compare_analysis" and seg_filter==True:
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{group}_{bin}um", "{sample}.zarr"))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      return(os.path.join(results_folder, "{group}", "{sample}.zarr"))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{group}", "{sample}.h5ad"))
  
  

def preprocess_output(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{sample}_{bin}um", 'preprocess',"filter_{sample}.zarr")))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      print("correct")
      return(directory(os.path.join(results_folder, "{sample}",'preprocess', "filter_{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'preprocess',"filter_{sample}.h5ad"))
  if channel=="compare_analysis" and seg_filter==True:
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{group}_{bin}um", "filter_{sample}.zarr")))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      return(directory(os.path.join(results_folder, "{group}","filter_{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{group}","filter_{sample}.h5ad"))
  else:
    return(parameter_output(samples,'preprocess'))
  
print(preprocess_input(run_type),preprocess_output(run_type))

rule preprocess_rule:
  input:
    outputs=preprocess_input(run_type)
  output:
    merge=preprocess_output(run_type) if run_type == "slide_seq" else directory(preprocess_output(run_type))
  params:
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" and seg_filter==False else wildcards.sample,
    run_type = run_type,
    variable = variable,
    seg_filter = seg_filter,
    NEIGHBORS = NEIGHBORS,
    harmony = harmony,
    filter_dict= f"--filter_dict '{json.dumps(filter_dict)}'" if seg_filter else "",
    min_cells= f"--min_cells {min_cells}" if not seg_filter else "",
    min_genes= f"--min_genes {min_genes}" if not seg_filter else ""
  shell:
      """
      python {spatialsnake_path}workflow/scripts/preprocessing.py \
        --input_dir {input.outputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --variable {params.variable} \
        --seg_filter {params.seg_filter} \
        --NEIGHBORS {params.NEIGHBORS} \
        --harmony {params.harmony} \
        {params.filter_dict} \
        {params.min_genes} \
        {params.min_cells}
      """


