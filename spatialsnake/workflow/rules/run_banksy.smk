rule banksy_rule:
  input:
    inputs=lambda wildcards: banksy_input
  output:
    merge=directory(os.path.join(results_folder,"banksy",f'{banksy_sample_id}_banksy.zarr'))
  params:
    k_geom = k_geom,
    max_m = max_m,
    nbr_weight_decay = nbr_weight_decay,
    n_comps = n_comps,
    lambda_list = lambda_list,
    RES = RES
  shell:
      """
      python {spatialsnake_path}workflow/scripts/banksy.py \
        --input_dir {input.inputs} \
        --output_zarr_path {output.merge} \
        --k_geom {params.k_geom} \
        --max_m {params.max_m} \
        --nbr_weight_decay {params.nbr_weight_decay} \
        --n_comps {params.n_comps} \
        --lambda_list {params.lambda_list} \
        --resolution {params.RES}
      """
