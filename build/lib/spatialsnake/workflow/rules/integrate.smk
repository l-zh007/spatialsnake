STEREOSEQ_TYPES = ["stereoseq", "StereoSeq", "Stereo-seq"]
IMAGING_ZARR_TYPES = ["visium", "xenium", "visium_segment", "Merfish", "merscope", "cosmx"] + STEREOSEQ_TYPES

def resolve_nested_input_dir(resolved_dir, default_dir, terminal_dirname, *suffix_parts):
  if not resolved_dir:
    return default_dir
  normalized_dir = os.path.normpath(resolved_dir)
  if os.path.basename(normalized_dir) == terminal_dirname:
    return normalized_dir
  return os.path.join(normalized_dir, *suffix_parts)

def resolve_input_path(wildcards):
  sample_name = wildcards.sample
  sample_dir_map = globals().get("sample_input_dir_map", {})
  compare_key = f"{wildcards.group}::{sample_name}" if channel == "compare_analysis" and hasattr(wildcards, "group") else None
  resolved_dir = sample_dir_map.get(compare_key, sample_dir_map.get(sample_name))
  if run_type in ["visium", "xenium", "slide_seq"]:
    base_dir = resolved_dir if resolved_dir else os.path.join(data_fold, sample_name)
    return os.path.join(base_dir, main_file)
  elif run_type in ["Merfish", "merscope", "cosmx"] or run_type in STEREOSEQ_TYPES:
    return resolved_dir if resolved_dir else os.path.join(data_fold, sample_name)
  elif run_type == "visium_HD":
    square_dir = f"square_{wildcards.bin}um"
    base_dir = resolve_nested_input_dir(
      resolved_dir,
      os.path.join(data_fold, sample_name, "binned_outputs", "square_{bin}um"),
      square_dir,
      "binned_outputs",
      square_dir,
    )
    return os.path.join(base_dir, main_file)
  elif run_type == "visium_segment":
    base_dir = resolve_nested_input_dir(
      resolved_dir,
      os.path.join(data_fold, sample_name, "segmented_outputs"),
      "segmented_outputs",
      "segmented_outputs",
    )
    return os.path.join(base_dir, main_file)

def get_bin_size_param(wildcards):
  if run_type == "visium_HD":
    return f"--bin_size {wildcards.bin}"
  if run_type in STEREOSEQ_TYPES:
    compare_key = f"{wildcards.group}::{wildcards.sample}" if channel == "compare_analysis" and hasattr(wildcards, "group") else wildcards.sample
    requested = globals().get("sample_stereoseq_input_spec_map", {}).get(compare_key, globals().get("sample_stereoseq_input_spec_map", {}).get(wildcards.sample, config.get("bin_size")))
    if requested in [None, "", False, "None", "False", "false"]:
      return ""
    return f"--input_spec '{requested}'"
  return ""
def get_output(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{sample}_{bin}um",'integrate', "{sample}.zarr")))
    elif run_type in IMAGING_ZARR_TYPES:
      return(directory(os.path.join(results_folder, "{sample}",'integrate', "{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}",'integrate',"{sample}.h5ad"))
  elif channel=="compare_analysis":
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{group}_{bin}um", "{sample}.zarr")))
    elif run_type in IMAGING_ZARR_TYPES:
      return(directory(os.path.join(results_folder, "{group}", "{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{group}", "{sample}.h5ad"))
  


rule get_zarr:
    input:
        path = resolve_input_path
    output:
        outputs = get_output(run_type)
    params:
        main_file = main_file,
        run_type = run_type,
        bin_size = get_bin_size_param,
        cells_boundaries=cells_boundaries if run_type == "xenium" else "False",
        nucleus_boundaries=nucleus_boundaries if run_type == "xenium" else "False",
        nucleus_labels=nucleus_labels if run_type == "xenium" else "False",
        morphology_mip=morphology_mip if run_type == "xenium" else "False",
        channel=channel,
        scale_factors=scale_factors,
        image=image,
        geojson=geojson,
        coor_file=coor_file,
        merscope_z_layers_arg=(
            f"--merscope_z_layers {merscope_z_layers}"
            if run_type == "Merfish" and merscope_z_layers not in [None, "", "None", "null"]
            else ""
        ),
        merscope_region_name_arg=(
            f"--merscope_region_name {merscope_region_name}"
            if run_type == "Merfish" and merscope_region_name not in [None, "", "None", "null"]
            else ""
        ),
        merscope_transcripts=merscope_transcripts if run_type == "Merfish" and channel == "single_analysis" else "True",
        merscope_cells_boundaries=merscope_cells_boundaries if run_type == "Merfish" and channel == "single_analysis" else "True",
        merscope_cells_table=merscope_cells_table if run_type == "Merfish" and channel == "single_analysis" else "True",
        merscope_mosaic_images=merscope_mosaic_images if run_type == "Merfish" and channel == "single_analysis" else "True"
    shell:
      """
      python {spatialsnake_path}workflow/scripts/spatial_in.py \
        --input_dir {input.path} \
        --output_zarr_path {output.outputs} \
        --sample_id {wildcards.sample} \
        --count_file {params.main_file} \
        --type {params.run_type} \
        --cells_boundaries {params.cells_boundaries} \
        --nucleus_boundaries {params.nucleus_boundaries} \
        --nucleus_labels {params.nucleus_labels} \
        --morphology_mip {params.morphology_mip} \
        --scale_factors {params.scale_factors} \
        --image {params.image} \
        --geojson {params.geojson} \
        --coor_file {params.coor_file} \
        --channel {params.channel} \
        {params.merscope_z_layers_arg} \
        {params.merscope_region_name_arg} \
        --merscope_transcripts {params.merscope_transcripts} \
        --merscope_cells_boundaries {params.merscope_cells_boundaries} \
        --merscope_cells_table {params.merscope_cells_table} \
        --merscope_mosaic_images {params.merscope_mosaic_images} \
        {params.bin_size}
      """
