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

rule reannotation_rule:
  input:
    inputs=lambda wildcards: subset_task_lookup[(wildcards.sample, wildcards.subset)]["input_path"]
  output:
    merge=directory(os.path.join(results_folder, "{sample}", "reannotation", "{subset}", "{subset}.zarr")),
    csv=os.path.join(results_folder, "{sample}", "reannotation", "{subset}", "celltype_annotations.csv"),
    proportion_png=os.path.join(results_folder, "{sample}", "reannotation", "{subset}", "celltype_proportion.png"),
    umap_png=os.path.join(results_folder, "{sample}", "reannotation", "{subset}", "umap_recluster.png"),
    spatial_png=os.path.join(results_folder, "{sample}", "reannotation", "{subset}", "spatial_clusters.png")
  threads: workflow_threads
  params:
    sample_id=lambda wildcards: wildcards.subset,
    run_type=run_type,
    anno_data=lambda wildcards: f"--anno_data '{json.dumps(anno_data)}'",
    vis_mode=vis_mode,
    point_size=config.get("point_size"),
    image_slice=image_slice,
    coord=f"--coord {coord[0]} {coord[1]} {coord[2]} {coord[3]}" if image_slice == True else "",
    input_spec=get_stereoseq_input_spec_param
  shell:
      """
      python {spatialsnake_path}workflow/scripts/reannotation.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --output_csv {output.csv} \
        --type {params.run_type} \
        --image_slice {params.image_slice} \
        --vis_mode {params.vis_mode} \
        --point_size {params.point_size} \
        {params.input_spec} \
        {params.coord} \
        {params.anno_data}
      """
