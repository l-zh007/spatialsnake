def get_cellcharter_output():
  if channel == "compare_analysis":
    return directory(os.path.join(results_folder, "merge_data", "cellcharter", f"{cellcharter_sample_id}_cellcharter.zarr"))
  return directory(os.path.join(results_folder, f"{cellcharter_sample_id}", "cellcharter", f"{cellcharter_sample_id}_cellcharter.zarr"))

rule cellcharter_rule:
  input:
    inputs=lambda wildcards: cellcharter_input
  output:
    merge=get_cellcharter_output()
  params:
    sample_id = cellcharter_sample_id,
    run_type = run_type,
    image_type = image_type,
    channal = channel,
    shape_type=shape_type,
    significance=significance,
    max_cluster=max_cluster,
    condition_col=condition_col,
    sample_col=sample_col,
    celltype_col=celltype_col,
    cellcharter_col=cellcharter_col
  shell:
      """
      python {spatialsnake_path}workflow/scripts/run_cellcharter.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --image_type {params.image_type} \
        --channal {params.channal} \
        --shape_type {params.shape_type} \
        --significance {params.significance} \
        --max_cluster {params.max_cluster} \
        --condition_col {params.condition_col} \
        --sample_col {params.sample_col} \
        --celltype_col {params.celltype_col} \
        --cellcharter_col {params.cellcharter_col}
      """
