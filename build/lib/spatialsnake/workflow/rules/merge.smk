def get_output(type):
    if run_type == "slide_seq":
        return(os.path.join(results_folder, "merge_data", "integrate", "concatenated_sdata.h5ad"))
    return(os.path.join(results_folder, "merge_data", "integrate", "concatenated_sdata.zarr"))

def get_merge_stereoseq_input_specs():
    if run_type not in ["stereoseq", "StereoSeq", "Stereo-seq"]:
        return ""
    requested_specs = []
    sample_spec_map = globals().get("sample_stereoseq_input_spec_map", {})
    for sample_name, group_name in zip(samples, group):
        compare_key = f"{group_name}::{sample_name}"
        requested = sample_spec_map.get(compare_key, sample_spec_map.get(sample_name, config.get("bin_size")))
        if requested in [None, "", False, "None", "False", "false", "NULL", "null"]:
            continue
        requested_specs.append(str(requested))
    if not requested_specs:
        return ""
    return "--input_spec " + " ".join(f"'{spec}'" for spec in requested_specs)

rule merge_in:
  input:
    outputs=parameter_output(samples,'integrate')
  output:
    merge=get_output(type) if run_type == "slide_seq" else directory(get_output(type))
  params:
    main_file = main_file,
    run_type = run_type,
    group = group,
    sample =samples,
    input_spec = get_merge_stereoseq_input_specs()
  shell:
      """
      python {spatialsnake_path}workflow/scripts/spatial_in_multiple.py \
        --input_path {input.outputs} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --sample_id {params.sample} \
        --group {params.group} \
        {params.input_spec}
      """
    
  
