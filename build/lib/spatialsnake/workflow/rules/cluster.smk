def cluster_input(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{sample}_{bin}um", 'preprocess',"filter_{sample}.zarr"))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      print("correct")
      return(os.path.join(results_folder, "{sample}",'preprocess', "filter_{sample}.zarr"))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'preprocess',"filter_{sample}.h5ad"))
  if channel=="compare_analysis":
    return(parameter_output(samples,'preprocess'))
# def cluster_output(run_type):
#   if channel == 'single_analysis':
#     if run_type=="visium_HD":
#       return(os.path.join(results_folder, "{sample}_{bin}um", 'clustering',"{sample}.zarr"))
#     elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
#       print("correct")
#       return(os.path.join(results_folder, "{sample}",'clustering', "{sample}.zarr"))
#     elif run_type=="slide_seq":
#       return(os.path.join(results_folder, "{sample}", 'clustering',"{sample}.h5ad"))
#   if channel=="compare_analysis":
#     return(parameter_output(samples,'clustering'))
  

rule cluster_rule:
  input:
    inputs=cluster_input(run_type)
  output:
    merge=nomal_file(run_type,"clustering") if run_type == "slide_seq" else directory(nomal_file(run_type,"clustering"))
  params:
    sample_id = lambda wildcards: "concatenated_sdata" if channel=="compare_analysis" and seg_filter==False else wildcards.sample,
    run_type = run_type,
    tsene = tsene,
    MIN_DIST = MIN_DIST,
    SPREAD = SPREAD,
    RES = RES,
    harmony = harmony,
    cluster_algorithm = cluster_algorithm
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
        --cluster_algorithm {params.cluster_algorithm}
      """
