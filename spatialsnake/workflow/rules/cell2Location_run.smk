rule cell2Location_rule:
  input:
    input_spatial=input_spatial,
    input_singlecell = input_singlecell
  output:
    output_dir_zarr=temp(directory(os.path.join(results_folder,"{sample}",'cell2Location','tem.zarr')))
  params:
    sample_id = samples,
    run_type = run_type,
    max_epochs_reference=config.get('max_epochs_reference',250),
    remove_mt=config.get('remove_mt',True),
    N_cells_per_location=config.get('N_cells_per_location',30),
    max_epochs_st = config.get('max_epochs_st',30000),
    device = device
  shell:
      """
      python {spatialsnake_path}workflow/scripts/cell2Location.py \
        --input_spatial {input.input_spatial} \
        --sample_id {params.sample_id} \
        --output_dir_zarr {output.output_dir_zarr} \
        --input_singlecell {input.input_singlecell} \
        --type {params.run_type} \
        --max_epochs_reference {params.max_epochs_reference} \
        --remove_mt {params.remove_mt} \
        --N_cells_per_location {params.N_cells_per_location} \
        --max_epochs_st {params.max_epochs_st} \
        --device {params.device}
      """

rule cell2Location_visualize_rule:
  input:
    inputs=os.path.join(results_folder,"{sample}",'cell2Location','tem.zarr')
  output:
    merge = directory(os.path.join(results_folder,"{sample}",'cell2Location','{sample}.zarr'))
  params:
    sample_id = samples,
    run_type = run_type,
    image_type = image_type,
    image_slice = False,
    sample_cnt = len(samples) if channel=="comparision_analysis" else 1,
    shape_type=shape_type
  shell:
      """
      python {spatialsnake_path}workflow/scripts/cell2locate_visualize.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --image_type {params.image_type} \
        --image_slice {params.image_slice} \
        --sample_cnt {params.sample_cnt} \
        --shape_type {params.shape_type}
      """