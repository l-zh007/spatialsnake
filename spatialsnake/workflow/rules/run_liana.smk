liana_dir = os.path.join(results_folder, f"{cpdb_sample_id}", "liana_output")

rule liana_compute:
  input:
    data=liana_inputs
  output:
    results=os.path.join(liana_dir, "liana_results.csv")
  threads: workflow_threads
  params:
    sample_id=cpdb_sample_id,
    output_zarr=os.path.join(liana_dir, f"{cpdb_sample_id}.zarr"),
    run_type=run_type,
    method=liana_method,
    resource=liana_resource_name,
    expr_prop=liana_expr_prop,
    min_cells=liana_min_cells,
    use_raw=liana_use_raw,
    celltype=celltype
  shell:
      """
      MPLCONFIGDIR=/tmp/mplconfig NUMBA_CACHE_DIR=/tmp/numba_cache \
      python {spatialsnake_path}workflow/scripts/run_liana.py \
        --stage compute \
        --input_dir "{input.data}" \
        --sample_id "{params.sample_id}" \
        --output_zarr_path "{params.output_zarr}" \
        --type "{params.run_type}" \
        --method "{params.method}" \
        --resource_name "{params.resource}" \
        --celltype "{params.celltype}" \
        --expr_prop "{params.expr_prop}" \
        --min_cells "{params.min_cells}" \
        --use_raw "{params.use_raw}" \
        --results_csv "{output.results}"
      """

rule liana_plot:
  input:
    data=liana_inputs,
    results=rules.liana_compute.output.results
  output:
    merge=directory(os.path.join(liana_dir, f"{cpdb_sample_id}.zarr")),
    top=os.path.join(liana_dir, "liana_top_interactions.csv"),
    plot_data=os.path.join(liana_dir, "liana_plot_data.csv"),
    pair_summary=os.path.join(liana_dir, "liana_pair_summary.csv"),
    dotplot=os.path.join(liana_dir, "dotplot.png"),
    tileplot=os.path.join(liana_dir, "tileplot.png"),
    heatmap=os.path.join(liana_dir, "communication_heatmap.png"),
    gene_expression=os.path.join(liana_dir, "gene_expression.png")
  threads: workflow_threads
  params:
    sample_id=cpdb_sample_id,
    run_type=run_type,
    method=liana_method,
    celltype=celltype,
    use_raw=liana_use_raw,
    top_n=config.get("liana_top_n", 6),
    pvalue=liana_pvalue,
    source_celltypes=liana_source_celltypes,
    target_celltypes=liana_target_celltypes,
    cell_pairs=liana_cell_pairs,
    pairs=liana_pairs
  shell:
      """
      MPLCONFIGDIR=/tmp/mplconfig NUMBA_CACHE_DIR=/tmp/numba_cache \
      python {spatialsnake_path}workflow/scripts/run_liana.py \
        --stage plot \
        --input_dir "{input.data}" \
        --sample_id "{params.sample_id}" \
        --output_zarr_path "{output.merge}" \
        --type "{params.run_type}" \
        --method "{params.method}" \
        --celltype "{params.celltype}" \
        --use_raw "{params.use_raw}" \
        --top_n "{params.top_n}" \
        --pvalue "{params.pvalue}" \
        --source_celltypes "{params.source_celltypes}" \
        --target_celltypes "{params.target_celltypes}" \
        --cell_pairs "{params.cell_pairs}" \
        --pairs "{params.pairs}" \
        --results_csv "{input.results}" \
        --top_csv "{output.top}" \
        --plot_data_csv "{output.plot_data}" \
        --pair_summary_csv "{output.pair_summary}"
      """
