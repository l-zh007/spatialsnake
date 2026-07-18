import os
import shlex


def _cpdb_optional_args(items):
  """Build safely quoted CLI arguments while omitting empty config values."""
  arguments = []
  for flag, value in items:
    value = "" if value is None else str(value).strip()
    if value:
      arguments.extend([flag, shlex.quote(value)])
  return " ".join(arguments)


cpdb_out_dir = os.path.join(results_folder, cpdb_sample_id, "cellPhoneDB_results") if channel=="single_analysis" else os.path.join(results_folder, "merge_data", "cellPhoneDB_results")
cpdb_heatmap_file = f"{cpdb_sample_id}_heatmap.png" if channel=="single_analysis" else "concentrate_heatmap.png"

cpdb_analysis_optional_args = _cpdb_optional_args([
  ("--output_name", output_name),
  ("--active_tf_path", active_tf_path),
  ("--microenvs_file_path", microenvs_file_path),
  ("--degs_file_path", degs_file_path),
])

cpdb_visual_optional_args = _cpdb_optional_args([
  ("--output_name", output_name),
  ("--cell_pairs", cell_pairs),
  ("--cell_type1", cell_type1),
  ("--cell_type2", cell_type2),
  ("--gene_family", gene_family),
  ("--cpdb_pathway", cpdb_pathway),
  ("--interaction_pairs", interaction_pairs),
  ("--cpdb_genes", cpdb_genes),
])

rule cellPhoneDB_rule:
  input:
    inputs=cellPhoneDB_input
  output:
    merge=os.path.join(cpdb_out_dir,"adata.h5ad")
  threads: workflow_threads
  params:
    sample_id = cpdb_sample_id,
    run_type = run_type,
    counts_data=counts_data,
    threshold= threshold,
    pvalue = pvalue,
    celltype_col=celltype_col,
    niche_col=niche_col,
    is_single_cell=is_single_cell,
    method=cpdb_method,
    species=species,
    iterations=iterations,
    optional_args=cpdb_analysis_optional_args
  shell:
      """
      OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
      python {spatialsnake_path:q}workflow/scripts/cellphoneDB.py \
        --input_dir {input.inputs:q} \
        --sample_id {params.sample_id:q} \
        --output_zarr_path {output.merge:q} \
        --type {params.run_type:q} \
        --counts_data {params.counts_data:q} \
        --celltype_col {params.celltype_col:q} \
        --niche_col {params.niche_col:q} \
        --is_single_cell {params.is_single_cell:q} \
        --method {params.method:q} \
        --species {params.species:q} \
        --iterations {params.iterations} \
        --threshold {params.threshold} \
        --threads {threads} \
        --pvalue {params.pvalue} \
        {params.optional_args}
      """
rule visualize_cellPhone:
  input:
    inputs=os.path.join(cpdb_out_dir,"adata.h5ad")
  output:
    merge=os.path.join(cpdb_out_dir,cpdb_heatmap_file)
  threads: workflow_threads
  params:
    sample_id = cpdb_sample_id,
    run_type = run_type,
    counts_data=counts_data,
    celltype=celltype,
    method=cpdb_method,
    pvalue=pvalue,
    optional_args=cpdb_visual_optional_args
  shell:
      """
      OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
      python {spatialsnake_path:q}workflow/scripts/cellphoneDB_visualize.py \
        --input_dir {input.inputs:q} \
        --sample_id {params.sample_id:q} \
        --output_zarr_path {output.merge:q} \
        --type {params.run_type:q} \
        --method {params.method:q} \
        --pvalue {params.pvalue} \
        --celltype {params.celltype:q} \
        {params.optional_args}
      """
