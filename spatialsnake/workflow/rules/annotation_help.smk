def annotation_help_output(run_type,file_name):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{sample}_{bin}um", 'clustering',file_name))
    elif run_type in ["visium", "xenium", "visium_segment", "Merfish", "merscope", "cosmx", "stereoseq", "StereoSeq", "Stereo-seq"]:
      print("correct")
      return(os.path.join(results_folder, "{sample}",'clustering', file_name))
  if channel=="compare_analysis":
    return(os.path.join(results_folder,"merge_data",'clustering',file_name))

def enrich_output_map(run_type):
  return {
    "kegg_csv": annotation_help_output(run_type, 'kegg_data.csv'),
    "go_csv": annotation_help_output(run_type, 'GO_data.csv'),
    "bp_png": annotation_help_output(run_type, 'BP_GO_cluster.png'),
    "cc_png": annotation_help_output(run_type, 'CC_GO_cluster.png'),
    "mf_png": annotation_help_output(run_type, 'MF_GO_cluster.png'),
  }

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


rule annotation_help:
  input:
    inputs=normal_file(run_type,"clustering")
  output:
    merge=annotation_help_output(run_type,'marker_genes_pval.csv')
    # os.path.join(results_folder,"{sample}",'clustering','marker_genes_pval.csv') if channel == 'single_analysis' else os.path.join(results_folder,"merge_data",'clustering','marker_genes_pval.csv') 
  threads: workflow_threads
  params:
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" else wildcards.sample,
    run_type = run_type,
    image_type = image_type,
    image_slice = image_slice,
    vis_mode = vis_mode,
    point_size = config.get("point_size"),
    coord = f"--coord {coord[0]} {coord[1]} {coord[2]} {coord[3]}" if image_slice==True else "",
    sample_cnt = len(samples) if channel=="compare_analysis" else 1,
    markers_algorithm = markers_algorithm,
    shape_type=shape_type,
    input_spec=get_stereoseq_input_spec_param
  shell:
      """
      python {spatialsnake_path}workflow/scripts/cluster_visualize.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --image_type {params.image_type} \
        --image_slice {params.image_slice} \
        --vis_mode {params.vis_mode} \
        --point_size {params.point_size} \
        --markers_algorithm {params.markers_algorithm} \
        --sample_cnt {params.sample_cnt} \
        --shape_type {params.shape_type} \
        --threads {threads} \
        {params.input_spec} \
        {params.coord}
      """

rule enrich_rule:
  input:
    inputs = annotation_help_output(run_type,'marker_genes_pval.csv')
  output:
    kegg_csv = enrich_output_map(run_type)['kegg_csv'],
    go_csv = enrich_output_map(run_type)['go_csv'],
    bp_png = enrich_output_map(run_type)['bp_png'],
    cc_png = enrich_output_map(run_type)['cc_png'],
    mf_png = enrich_output_map(run_type)['mf_png']
  threads: workflow_threads
  params:
    species = species,
    run_type = run_type,
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" else wildcards.sample,
  shell:
      """
      Rscript {spatialsnake_path}workflow/scripts/enrichment.R \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_path {output.kegg_csv} \
        --type {params.run_type} \
        --species {params.species}
      """

    










