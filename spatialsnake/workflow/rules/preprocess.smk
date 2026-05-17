STEREOSEQ_TYPES = ["stereoseq", "StereoSeq", "Stereo-seq"]

def preprocess_input(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{sample}_{bin}um", 'integrate',"{sample}.zarr"))
    elif run_type in ["visium", "xenium", "visium_segment", "Merfish", "merscope", "cosmx", "stereoseq", "StereoSeq", "Stereo-seq"]:
      return(os.path.join(results_folder, "{sample}", 'integrate',"{sample}.zarr"))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'integrate',"{sample}.h5ad"))
  if channel=="compare_analysis":
    if run_type=="slide_seq":
      return(os.path.join(results_folder, "merge_data", "integrate", "concatenated_sdata.h5ad"))
    return(os.path.join(results_folder, "merge_data", "integrate", "concatenated_sdata.zarr"))

def preprocess_output(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{sample}_{bin}um", 'preprocess',"filter_{sample}.zarr")))
    elif run_type in ["visium", "xenium", "visium_segment", "Merfish", "merscope", "cosmx", "stereoseq", "StereoSeq", "Stereo-seq"]:
      print("correct")
      return(directory(os.path.join(results_folder, "{sample}",'preprocess', "filter_{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'preprocess',"filter_{sample}.h5ad"))
  if channel=="compare_analysis":
    return(parameter_output(samples,'preprocess'))


def get_stereoseq_input_spec_param(wildcards):
  if run_type not in STEREOSEQ_TYPES:
    return ""
  sample_name = wildcards.sample if hasattr(wildcards, "sample") else None
  if not sample_name:
    return ""
  requested = globals().get("sample_stereoseq_input_spec_map", {}).get(sample_name, config.get("bin_size"))
  if requested in [None, "", False, "None", "False", "false", "NULL", "null"]:
    return ""
  return f"--input_spec '{requested}'"



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
    input_spec=get_stereoseq_input_spec_param,
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
        {params.input_spec} \
        {params.filter_dict} \
        {params.min_genes} \
        {params.min_cells}
      """
