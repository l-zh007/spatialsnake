rule annotion_mannel:
  input:
    inputs=nomal_file(run_type,"clustering")
  output:
    merge=directory(nomal_file(run_type,"annotion")) if run_type!="slide_seq" else nomal_file(run_type,"annotion")
  params:
    sample_id = lambda wildcards: "concatenated_sdata" if channel=="compare_analysis" and seg_filter==False else wildcards.sample,
    run_type = run_type,
    anno_data = lambda wildcards: f"--anno_data '{json.dumps(anno_data)}'",
    image_type = image_type,
    image_slice = image_slice,
    coord = coord,
    sample_cnt = len(samples) if channel=="comparision_analysis" else 1,
    shape_type=shape_type
  shell:
      """
      python workflow/scripts/mannul_annotion.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --image_type {params.image_type} \
        --image_slice {params.image_slice} \
        --sample_cnt {params.sample_cnt} \
        --coord {params.coord} \
        --shape_type {params.shape_type} \
        {params.anno_data}
      """
