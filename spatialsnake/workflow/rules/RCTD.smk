rule RCTD_run:
    input:
        spatial = input_spatial,
        single_cell = input_singlecell
    output:
        results = os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_results.csv"),
        weights = os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_weights.csv")
    params:
        output_dir = lambda wildcards: os.path.join(results_folder, wildcards.sample, "RCTD"),
        sample_id = "{sample}",
        mode = config.get("RCTD_mode", "doublet"),
        sc_cell_type_col = config.get("sc_cell_type_col", config.get("cell_type_col", "celltype")),
        spatial_cell_type_col = config.get("spatial_cell_type_col", "celltype"),
        group_by = config.get("group_by", "sample"),
        max_cores = config.get("max_cores", 8)
    shell:
        """
        Rscript {spatialsnake_path}workflow/scripts/RCTD.R \
            --spatial_input {input.spatial} \
            --sc_input {input.single_cell} \
            --output_dir {params.output_dir} \
            --sample_id {params.sample_id} \
            --mode {params.mode} \
            --max_cores {params.max_cores} \
            --sc_cell_type_col {params.sc_cell_type_col} \
            --spatial_cell_type_col {params.spatial_cell_type_col} \
            --group_by {params.group_by}
        """

rule RCTD_visualize:
    input:
        results = rules.RCTD_run.output.results,
        weights = rules.RCTD_run.output.weights,
        zarr_input = config.get("zarr_input")
    output:
        output_zarr = directory(os.path.join(results_folder, "{sample}", "RCTD", "{sample}.zarr")),
        plot = os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_spatial_plot.png")
    params:
        sample_id = "{sample}",
        mode = config.get("RCTD_mode", "doublet")
    shell:
        """
        python {spatialsnake_path}workflow/scripts/RCTD_visualize.py \
            --zarr_input {input.zarr_input} \
            --rctd_weights {input.weights} \
            --rctd_results {input.results} \
            --output_zarr {output.output_zarr} \
            --sample_id {params.sample_id} \
            --mode {params.mode} \
            --output_plot {output.plot}
        """
