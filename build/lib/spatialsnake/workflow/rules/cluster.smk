def cluster_input(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{sample}_{bin}um", 'preprocess',"filter_{sample}.zarr"))
    elif run_type in ["visium", "xenium", "visium_segment", "Merfish", "merscope", "cosmx", "stereoseq", "StereoSeq", "Stereo-seq"]:
      print("correct")
      return(os.path.join(results_folder, "{sample}",'preprocess', "filter_{sample}.zarr"))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'preprocess',"filter_{sample}.h5ad"))
  if channel=="compare_analysis":
    return(parameter_output(samples,'preprocess'))

def get_stereoseq_input_spec_param(wildcards):
  if run_type not in ["stereoseq", "StereoSeq", "Stereo-seq"]:
    return ""
  sample_name = wildcards.sample if hasattr(wildcards, "sample") else None
  if channel == "compare_analysis":
    requested = next(iter(globals().get("sample_stereoseq_input_spec_map", {}).values()), config.get("bin_size"))
  else:
    requested = globals().get("sample_stereoseq_input_spec_map", {}).get(sample_name, config.get("bin_size"))
  if requested in [None, "", False, "None", "False", "false", "NULL", "null"]:
    return ""
  return f"--input_spec '{requested}'"

rule cluster_rule:
  input:
    inputs=cluster_input(run_type)
  output:
    merge=nomal_file(run_type,"clustering") if run_type == "slide_seq" else directory(nomal_file(run_type,"clustering"))
  params:
    sample_id = lambda wildcards: "concatenated_sdata" if channel=="compare_analysis"  else wildcards.sample,
    run_type = run_type,
    tsene = tsene,
    MIN_DIST = MIN_DIST,
    SPREAD = SPREAD,
    RES = RES,
    cluster_algorithm = cluster_algorithm,
    n_clusters = n_clusters,
    k_geom = k_geom,
    max_m = max_m,
    nbr_weight_decay = nbr_weight_decay,
    n_comps = n_comps,
    lambda_list = lambda_list,
    sketch = sketch,
    pcs = pcs,
    NEIGHBORS = NEIGHBORS,
    input_spec = get_stereoseq_input_spec_param
  shell:
      """
        python {spatialsnake_path}workflow/scripts/clustering.py \
              --input_dir {input.inputs} \
              --sample_id {params.sample_id} \
              --output_zarr_path {output.merge} \
              --type {params.run_type} \
              --tsene {params.tsene} \
              --MIN_DIST {params.MIN_DIST} \
              --SPREAD {params.SPREAD} \
              --RES {params.RES} \
              --cluster_algorithm {params.cluster_algorithm} \
              --n_clusters {params.n_clusters} \
              --sketch {params.sketch} \
              --NEIGHBORS {params.NEIGHBORS} \
              --pcs {params.pcs} \
              {params.input_spec}
      """
      
