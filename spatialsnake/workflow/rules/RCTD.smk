RCTD_MODE = str(config.get("RCTD_mode", "doublet")).strip().lower()
if RCTD_MODE not in {"full", "doublet"}:
    raise ValueError("RCTD_mode must be either 'full' or 'doublet'")

if "zarr_input" in config:
    raise ValueError(
        "RCTD no longer accepts the zarr_input parameter. Move that SpatialData "
        "Zarr path to the second column of sample.txt: "
        "sample_id spatial_zarr sc_reference."
    )

if len(samples) != len(set(samples)):
    duplicate_samples = sorted(
        {sample for sample in samples if samples.count(sample) > 1}
    )
    raise ValueError(
        "RCTD sample IDs must be unique in sample.txt; duplicate ID(s): "
        + ", ".join(duplicate_samples)
    )


def _rctd_spatial_zarr(wildcards):
    path = (
        input_spatial
        if channel == "compare_analysis"
        else rctd_spatial_by_sample[wildcards.sample]
    )
    path = os.path.normpath(str(path).strip())
    if not path.lower().endswith(".zarr"):
        raise ValueError(
            f"RCTD sample '{wildcards.sample}' requires an existing SpatialData "
            f".zarr directory in the second column of sample.txt; got: {path}"
        )
    if not os.path.isdir(path):
        raise ValueError(
            f"RCTD SpatialData Zarr directory does not exist for sample "
            f"'{wildcards.sample}': {path}"
        )
    return path


def _rctd_single_cell_reference(wildcards):
    path = (
        input_singlecell
        if channel == "compare_analysis"
        else rctd_reference_by_sample[wildcards.sample]
    )
    path = os.path.normpath(str(path).strip())
    if os.path.splitext(path)[1].lower() not in {".h5ad", ".rds"}:
        raise ValueError(
            f"RCTD sample '{wildcards.sample}' requires an .h5ad or .rds "
            f"single-cell reference in the third column of sample.txt; got: {path}"
        )
    if not os.path.isfile(path):
        raise ValueError(
            f"RCTD single-cell reference does not exist for sample "
            f"'{wildcards.sample}': {path}"
        )
    return path


def _rctd_temp_h5ad(wildcards):
    zarr_stem = os.path.splitext(os.path.basename(_rctd_spatial_zarr(wildcards)))[0]
    return os.path.join(
        results_folder,
        wildcards.sample,
        "RCTD",
        ".rctd_h5ad",
        f"{zarr_stem}.h5ad",
    )


RCTD_TEMP_DIR = os.path.join(results_folder, "{sample}", "RCTD", ".rctd_h5ad")
RCTD_PLOT_NAME = (
    "{sample}_RCTD_full_dotplot.pdf"
    if RCTD_MODE == "full"
    else "{sample}_RCTD_spatial_plot.pdf"
)
RCTD_SOURCE_NAME = (
    "{sample}_RCTD_full_dotplot_source.tsv"
    if RCTD_MODE == "full"
    else "{sample}_RCTD_spatial_plot_source.tsv"
)
RCTD_EXTRA_OUTPUTS = (
    []
    if RCTD_MODE == "full"
    else [
        os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_doublet_proportion_dotplot.pdf"),
        os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_doublet_proportion_dotplot_source.tsv"),
        os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_spot_class_bar.pdf"),
        os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_spot_class_bar_source.tsv"),
    ]
)


rule RCTD_prepare_spatial:
    input:
        spatial_zarr = _rctd_spatial_zarr
    output:
        converted_dir = temp(directory(RCTD_TEMP_DIR))
    threads: workflow_threads
    params:
        converted_h5ad = _rctd_temp_h5ad,
        validator = os.path.join(spatialsnake_path, "workflow", "scripts", "RCTD_validate_input.py")
    shell:
        """
        OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
        spatialsnake useful_tool --option=transform {input.spatial_zarr:q} \
            --transform_from=zarr \
            --transform_to=h5ad \
            --save_image=True \
            --output_dir={output.converted_dir:q}
        python {params.validator:q} \
            --spatial_zarr {input.spatial_zarr:q} \
            --converted_h5ad {params.converted_h5ad:q}
        """


rule RCTD_run:
    input:
        converted_dir = rules.RCTD_prepare_spatial.output.converted_dir,
        single_cell = _rctd_single_cell_reference
    output:
        results = os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_results.csv"),
        weights = os.path.join(results_folder, "{sample}", "RCTD", "{sample}_RCTD_weights.csv")
    threads: workflow_threads
    params:
        spatial_h5ad = _rctd_temp_h5ad,
        output_dir = lambda wildcards: os.path.join(results_folder, wildcards.sample, "RCTD"),
        sample_id = "{sample}",
        mode = RCTD_MODE,
        sc_cell_type_col = config.get("sc_cell_type_col", config.get("cell_type_col", "celltype")),
        spatial_cell_type_col = config.get("spatial_cell_type_col", "celltype"),
        group_by = config.get("group_by", "sample")
    shell:
        """
        Rscript {spatialsnake_path}workflow/scripts/RCTD.R \
            --spatial_input {params.spatial_h5ad:q} \
            --sc_input {input.single_cell:q} \
            --output_dir {params.output_dir:q} \
            --sample_id {params.sample_id:q} \
            --mode {params.mode:q} \
            --max_cores {threads} \
            --sc_cell_type_col {params.sc_cell_type_col:q} \
            --spatial_cell_type_col {params.spatial_cell_type_col:q} \
            --group_by {params.group_by:q}
        """


rule RCTD_visualize:
    input:
        results = rules.RCTD_run.output.results,
        weights = rules.RCTD_run.output.weights,
        spatial_zarr = _rctd_spatial_zarr
    output:
        output_zarr = directory(os.path.join(results_folder, "{sample}", "RCTD", "{sample}.zarr")),
        plot = os.path.join(results_folder, "{sample}", "RCTD", RCTD_PLOT_NAME),
        source = os.path.join(results_folder, "{sample}", "RCTD", RCTD_SOURCE_NAME),
        extras = RCTD_EXTRA_OUTPUTS
    threads: workflow_threads
    params:
        sample_id = "{sample}",
        mode = RCTD_MODE,
        run_type = run_type,
        cluster_col = config.get("spatial_cell_type_col", "celltype"),
        sample_col = config.get("group_by", "sample"),
        max_cell_types = config.get("rctd_dotplot_max_cell_types", 30),
        enrichment_clip = config.get("rctd_dotplot_enrichment_clip", 2.5)
    shell:
        """
        OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
        python {spatialsnake_path}workflow/scripts/RCTD_visualize.py \
            --spatial_zarr {input.spatial_zarr:q} \
            --rctd_weights {input.weights:q} \
            --rctd_results {input.results:q} \
            --output_zarr {output.output_zarr:q} \
            --sample_id {params.sample_id:q} \
            --mode {params.mode:q} \
            --run_type {params.run_type:q} \
            --output_plot {output.plot:q} \
            --source_data {output.source:q} \
            --cluster_col {params.cluster_col:q} \
            --sample_col {params.sample_col:q} \
            --max_cell_types {params.max_cell_types} \
            --enrichment_clip {params.enrichment_clip}
        """
