rule run_slice_rule:
  input:
    input_path=INPUT_FIlE
  output:
    merge=INPUT_FIlE if os.path.splitext(INPUT_FIlE)[1]==".h5ad" else directory(get_output(type))
  params:
    barcode=barcode,
    max_x=max_x,
    min_x=min_x,
    max_y=max_y,
    min_y=min_y
  shell:
      """
      python {spatialsnake_path}workflow/scripts/slice_out.py \
        --input_path {input.input_path} \
        --output_zarr_path {output.merge} \
        --barcode {params.barcode} \
        --max_x {params.max_x} \
        --min_x {params.min_x} \
        --max_y {params.max_y} \
        --min_y {params.min_y}
      """
