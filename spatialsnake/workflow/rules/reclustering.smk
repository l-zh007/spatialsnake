rule reclustering_rule:
  input:
    data=lambda wildcards: subset_task_lookup[(wildcards.sample, wildcards.subset)]["input_path"]
  output:
    umap=os.path.join(results_folder, "{sample}", "reclustering", "{subset}", "umap_recluster.png"),
    spatial=os.path.join(results_folder, "{sample}", "reclustering", "{subset}", "spatial_clusters.png"),
    markers=os.path.join(results_folder, "{sample}", "reclustering", "{subset}", "marker_genes.csv"),
    assignments=os.path.join(results_folder, "{sample}", "reclustering", "{subset}", "cluster_assignments.csv"),
    zarr=directory(os.path.join(results_folder, "{sample}", "reclustering", "{subset}", "{subset}.zarr"))
  threads: workflow_threads
  params:
    sample_id=lambda wildcards: wildcards.subset,
    output_dir=lambda wildcards: os.path.join(
      results_folder, wildcards.sample, "reclustering", wildcards.subset
    ),
    resolution=recluster_resolution,
    n_top_genes=recluster_n_top_genes,
    neighbors=recluster_neighbors,
    n_pcs=recluster_n_pcs,
    marker_method=recluster_marker_method,
    min_pct=recluster_min_pct,
    logfc_threshold=recluster_logfc_threshold
  shell:
      """
      OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
      python {spatialsnake_path}workflow/scripts/reclustering.py \
        --input {input.data} \
        --output_dir {params.output_dir} \
        --sample_id {params.sample_id} \
        --resolution {params.resolution} \
        --n_top_genes {params.n_top_genes} \
        --neighbors {params.neighbors} \
        --n_pcs {params.n_pcs} \
        --marker_method {params.marker_method} \
        --min_pct {params.min_pct} \
        --logfc_threshold {params.logfc_threshold} \
        --threads {threads}
      """
