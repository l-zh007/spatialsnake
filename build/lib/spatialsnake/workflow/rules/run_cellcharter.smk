rule cellcharter_rule:
  input:
    inputs=cellcharter_input
  output:
    merge=os.path.join(results_folder,"cellcharter",f'{sample_id}_cellcharter.zarr') 
  params:
    sample_id = sample_id,
    run_type = run_type,
    image_type = image_type,
    channal = channal,
    shape_type=shape_type,
    significance=significance,
    max_cluster=max_cluster,
    condition_col=condition_col,
    sample_col=sample_col,
    celltype_col=celltype_col,
    cellcharter_col=cellcharter_col
  shell:
      """
      python {spatialsnake_path}workflow/scripts/cellcharter.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --image_type {params.image_type} \
        --channal {params.channal} \
        --shape_type {params.shape_type}
      """
