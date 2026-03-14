sample_path_map = {sample: path for sample, path in zip(samples, downstream_file)}

rule reannotation_rule:
  input:
    inputs=lambda wildcards: sample_path_map[wildcards.sample]
  output:
    merge=os.path.join(results_folder, "{sample}", "reannotation", "{sample}.zarr") if run_type != "slide_seq" else os.path.join(results_folder, "{sample}", "reannotation", "{sample}.h5ad"),
    csv=os.path.join(results_folder, "{sample}", "reannotation", "celltype_annotations.csv")
  params:
    sample_id=lambda wildcards: wildcards.sample,
    run_type=run_type,
    anno_data=lambda wildcards: f"--anno_data '{json.dumps(anno_data)}'"
  shell:
      """
      python {spatialsnake_path}workflow/scripts/reannotation.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --output_csv {output.csv} \
        --type {params.run_type} \
        {params.anno_data}
      """
