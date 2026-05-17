sample_path_map = {sample: path for sample, path in zip(samples, downstream_file)}

rule reclustering_rule:
  input:
    data=lambda wildcards: sample_path_map[wildcards.sample]
  output:
    umap=os.path.join(results_folder, "{sample}", "reclustering", "umap_recluster.png"),
    spatial=os.path.join(results_folder, "{sample}", "reclustering", "spatial_clusters.png"),
    markers=os.path.join(results_folder, "{sample}", "reclustering", "marker_genes.csv"),
    assignments=os.path.join(results_folder, "{sample}", "reclustering", "cluster_assignments.csv"),
    zarr=directory(os.path.join(results_folder, "{sample}", "reclustering", "{sample}.zarr"))
  params:
    sample_id=lambda wildcards: wildcards.sample,
    resolution=recluster_resolution,
    n_top_genes=recluster_n_top_genes,
    neighbors=recluster_neighbors,
    n_pcs=recluster_n_pcs,
    marker_method=recluster_marker_method,
    min_pct=recluster_min_pct,
    logfc_threshold=recluster_logfc_threshold
  shell:
      """
      python {spatialsnake_path}workflow/scripts/reclustering.py \
        --input {input.data} \
        --output_dir {results_folder}/{wildcards.sample}/reclustering \
        --sample_id {params.sample_id} \
        --resolution {params.resolution} \
        --n_top_genes {params.n_top_genes} \
        --neighbors {params.neighbors} \
        --n_pcs {params.n_pcs} \
        --marker_method {params.marker_method} \
        --min_pct {params.min_pct} \
        --logfc_threshold {params.logfc_threshold}
      """
