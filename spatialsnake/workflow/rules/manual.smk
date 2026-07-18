def annotation_output(run_type, file_name, compare_file_name=None):
  if channel == 'single_analysis':
    if run_type == "visium_HD":
      return os.path.join(results_folder, "{sample}_{bin}um", "annotation", file_name)
    elif run_type in ["visium", "xenium", "visium_segment", "Merfish", "merscope", "cosmx", "stereoseq", "StereoSeq", "Stereo-seq"]:
      return os.path.join(results_folder, "{sample}", "annotation", file_name)
  if channel == "compare_analysis":
    if compare_file_name is None:
      compare_file_name = file_name.replace("{sample}", "concatenated_sdata")
    return os.path.join(results_folder, "merge_data", "annotation", compare_file_name)

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

rule annotation_manual:
  input:
    inputs=normal_file(run_type,"clustering")
  output:
    main=directory(annotation_output(run_type, "{sample}.zarr", "concatenated_sdata.zarr")),
    cell_clusters=annotation_output(run_type, "{sample}_cell_clusters.csv", "concatenated_sdata_cell_clusters.csv"),
    proportion_png=annotation_output(run_type, "celltype_proportion.png"),
    umap_png=annotation_output(run_type, "{sample}_UMAP.png", "concatenated_sdata_UMAP.png"),
    gene_enrich_png=annotation_output(run_type, "{sample}_gene_enrich.png", "concatenated_sdata_gene_enrich.png")
  threads: workflow_threads
  params:
    sample_id = lambda wildcards: "concatenated_sdata" if channel=="compare_analysis" else wildcards.sample,
    run_type = run_type,
    anno_data = lambda wildcards: f"--anno_data '{json.dumps(anno_data)}'",
    vis_mode = vis_mode,
    point_size = config.get("point_size"),
    image_slice = image_slice,
    coord = f"--coord {coord[0]} {coord[1]} {coord[2]} {coord[3]}" if image_slice==True else "",
    sample_cnt = len(samples) if channel=="compare_analysis" else 1,
    input_spec = get_stereoseq_input_spec_param
  shell:
      """
      python {spatialsnake_path}workflow/scripts/manual_annotation.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.main} \
        --type {params.run_type} \
        --image_slice {params.image_slice} \
        --vis_mode {params.vis_mode} \
        --point_size {params.point_size} \
        {params.input_spec} \
        {params.anno_data} \
        {params.coord}
      """
