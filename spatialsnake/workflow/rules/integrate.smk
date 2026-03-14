def input_file(run_type):
  if run_type=="visium" or run_type=="xenium" or run_type=="slide_seq":
    return(os.path.join(data_fold,'{sample}',main_file))
  elif run_type=="Merfish":
    return(os.path.join(data_fold,'{sample}'))
  elif run_type=="visium_HD":
    return(os.path.join(data_fold,'{sample}',"binned_outputs","square_{bin}um",main_file))
  elif run_type=="visium_segment":
    return os.path.join(data_fold, "{sample}", "segmented_outputs", main_file)
def get_output(run_type):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{sample}_{bin}um",'integrate', "{sample}.zarr")))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment" or run_type=="Merfish":
      return(directory(os.path.join(results_folder, "{sample}",'integrate', "{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}",'integrate',"{sample}.h5ad"))
  elif channel=="compare_analysis":
    if run_type=="visium_HD":
      return(directory(os.path.join(results_folder, "{group}_{bin}um", "{sample}.zarr")))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment" or run_type=="Merfish":
      return(directory(os.path.join(results_folder, "{group}", "{sample}.zarr")))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{group}", "{sample}.h5ad"))
  


rule get_zarr:
    input:
        path = input_file(run_type)
    output:
        outputs = get_output(run_type)
    params:
        main_file = main_file,
        run_type = run_type,
        bin_size = "--bin_size {bin}" if run_type == "visium_HD" else "",
        cells_boundaries=cells_boundaries if run_type == "xenium" else "False",
        nucleus_boundaries=nucleus_boundaries if run_type == "xenium" else "False",
        nucleus_labels=nucleus_labels if run_type == "xenium" else "False",
        morphology_mip=morphology_mip if run_type == "xenium" else "False",
        channel=channel,
        scale_factors=scale_factors,
        image=image,
        geojson=geojson,
        coor_file=coor_file
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
        {params.bin_size}
      """
