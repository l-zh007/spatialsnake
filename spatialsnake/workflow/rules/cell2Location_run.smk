rule cell2Location_rule:
  input:
    input_spatial=cell2location_input_spatial if channel=="compare_analysis" else lambda wildcards: cell2location_spatial_by_sample[wildcards.sample],
    input_singlecell=cell2location_input_singlecell if channel=="compare_analysis" else lambda wildcards: cell2location_reference_by_sample[wildcards.sample]
  output:
    raw_zarr=temp(directory(os.path.join(results_folder, "merge_data", "cell2Location", "concatenated_sdata_raw.zarr"))) if channel=="compare_analysis" else temp(directory(os.path.join(results_folder,"{sample}",'cell2Location','{sample}_raw.zarr'))),
    figure_dir=directory(os.path.join(results_folder, "merge_data", "cell2Location", "figure")) if channel=="compare_analysis" else directory(os.path.join(results_folder,"{sample}",'cell2Location','figure'))
  threads: workflow_threads
  params:
    sample_id = "concatenated_sdata" if channel=="compare_analysis" else "{sample}",
    run_type = run_type,
    max_epochs_reference=config.get('max_epochs_reference',250),
    remove_mt=config.get('remove_mt',True),
    N_cells_per_location=config.get('N_cells_per_location',30),
    max_epochs_st = config.get('max_epochs_st',30000),
    labels_key_reference = config.get('labels_key_reference',"celltype"),
    batch_key_reference = config.get('batch_key_reference',"sample"),
    cell_count_cutoff = config.get('cell_count_cutoff',15),
    cell_percentage_cutoff2 = config.get('cell_percentage_cutoff2',0.05),
    nonz_mean_cutoff = config.get('nonz_mean_cutoff',1.12),
    batch_key_st = config.get('batch_key_st',"sample"),
    detection_alpha = config.get('detection_alpha',20),
    save_models = config.get('save_models',True),
    device = device
  shell:
      """
      OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
      python {spatialsnake_path}workflow/scripts/cell2Location.py \
        --input_spatial {input.input_spatial} \
        --sample_id {params.sample_id} \
        --output_dir_zarr {output.raw_zarr} \
        --input_singlecell {input.input_singlecell} \
        --type {params.run_type} \
        --max_epochs_reference {params.max_epochs_reference} \
        --remove_mt {params.remove_mt} \
        --N_cells_per_location {params.N_cells_per_location} \
        --max_epochs_st {params.max_epochs_st} \
        --labels_key_reference "{params.labels_key_reference}" \
        --batch_key_reference "{params.batch_key_reference}" \
        --cell_count_cutoff {params.cell_count_cutoff} \
        --cell_percentage_cutoff2 {params.cell_percentage_cutoff2} \
        --nonz_mean_cutoff {params.nonz_mean_cutoff} \
        --batch_key_st "{params.batch_key_st}" \
        --detection_alpha {params.detection_alpha} \
        --save_models {params.save_models} \
        --device {params.device}
      """

rule cell2Location_visualize_rule:
  input:
    raw_zarr=rules.cell2Location_rule.output.raw_zarr,
    figure_dir=rules.cell2Location_rule.output.figure_dir
  output:
    zarr=directory(os.path.join(results_folder, "merge_data", "cell2Location", "concatenated_sdata.zarr")) if channel=="compare_analysis" else directory(os.path.join(results_folder,"{sample}",'cell2Location','{sample}.zarr')),
    coloc_dir=directory(os.path.join(results_folder, "merge_data", "cell2Location", "CoLocatedComb")) if channel=="compare_analysis" else directory(os.path.join(results_folder,"{sample}",'cell2Location','CoLocatedComb')),
    microenvironment=os.path.join(results_folder, "merge_data", "cell2Location", "cellphonedb_microenvironments.tsv") if channel=="compare_analysis" else os.path.join(results_folder,"{sample}",'cell2Location','cellphonedb_microenvironments.tsv')
  threads: workflow_threads
  params:
    sample_id = "concatenated_sdata" if channel=="compare_analysis" else "{sample}",
    run_type = run_type,
    image_type = image_type,
    image_slice = config.get("image_slice", False),
    sample_cnt = len(samples) if channel=="compare_analysis" else 1,
    shape_type=shape_type,
    celltype_col = config.get("celltype_col", "celltype"),
    batch_key_st = config.get("batch_key_st", "sample"),
    microenvironment_threshold = config.get("cell2location_microenvironment_threshold", 0.10),
    dotplot_max_cell_types = config.get("cell2location_dotplot_max_cell_types", 30),
    dotplot_enrichment_clip = config.get("cell2location_dotplot_enrichment_clip", 2.5)
  shell:
      """
      OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
      python {spatialsnake_path}workflow/scripts/cell2locate_visualize.py \
        --input_dir {input.raw_zarr} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.zarr} \
        --type {params.run_type} \
        --image_type {params.image_type} \
        --image_slice {params.image_slice} \
        --sample_cnt {params.sample_cnt} \
        --shape_type {params.shape_type} \
        --celltype_col "{params.celltype_col}" \
        --batch_key_st "{params.batch_key_st}" \
        --microenvironment_output {output.microenvironment} \
        --microenvironment_threshold {params.microenvironment_threshold} \
        --dotplot_max_cell_types {params.dotplot_max_cell_types} \
        --dotplot_enrichment_clip {params.dotplot_enrichment_clip}
      """
