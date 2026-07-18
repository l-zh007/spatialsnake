COMPARE_OUTPUT_PATTERN = os.path.join(
    compare_result_root,
    compare_algorithm,
    "{target}",
    "{contrast}",
)


def get_compare_task(wildcards):
    key = (str(wildcards.target), str(wildcards.contrast))
    if key not in compare_task_lookup:
        raise ValueError(f"Unknown compare_gene task: {key}")
    return compare_task_lookup[key]


rule compare_gene_task:
    input:
        zarr=compare_input_zarr,
        sample_list=sample_list
    output:
        result_dir=directory(COMPARE_OUTPUT_PATTERN)
    threads:
        workflow_threads if compare_algorithm == "DESeq2" else 1
    params:
        celltype=lambda wildcards: get_compare_task(wildcards)["celltype"],
        comparison=lambda wildcards: get_compare_task(wildcards)["comparison"],
        reference=lambda wildcards: get_compare_task(wildcards)["reference"],
        run_type=run_type,
        algorithm=compare_algorithm,
        celltype_col=compare_celltype_col,
        sample_col=compare_sample_col,
        condition_col=compare_condition_col,
        count_layer=count_layer,
        min_replicates=min_replicates,
        min_cells_per_sample=min_cells_per_sample,
        min_total_counts_per_gene=min_total_counts_per_gene,
        cut_off_pvalue=config.get("cut_off_pvalue", 0.05),
        cut_off_logFC=config.get("cut_off_logFC", 0.5),
        de_top_n=de_top_n,
        gene_symbol_col=compare_gene_symbol_col,
        gene_id_type=compare_gene_id_type,
        species=species,
        go_ontology=compare_go_ontology,
        enrichment_top_n=compare_enrichment_top_n,
        python_script=os.path.join(spatialsnake_path, "workflow", "scripts", "DESeq2.py"),
        differential_r=os.path.join(spatialsnake_path, "workflow", "scripts", "differential_analysis.R"),
        enrichment_r=os.path.join(spatialsnake_path, "workflow", "scripts", "compare_enrichment.R"),
        rscript=config.get("rscript", "Rscript")
    shell:
        r"""
        OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} \
        MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
        python {params.python_script:q} \
          --input_dir {input.zarr:q} \
          --sample_list {input.sample_list:q} \
          --output_dir {output.result_dir:q} \
          --type {params.run_type:q} \
          --algorithm {params.algorithm:q} \
          --celltype {params.celltype:q} \
          --comparison {params.comparison:q} \
          --reference {params.reference:q} \
          --celltype_col {params.celltype_col:q} \
          --sample_col {params.sample_col:q} \
          --condition_col {params.condition_col:q} \
          --count_layer {params.count_layer:q} \
          --min_replicates {params.min_replicates} \
          --min_cells_per_sample {params.min_cells_per_sample} \
          --min_total_counts_per_gene {params.min_total_counts_per_gene} \
          --cut_off_pvalue {params.cut_off_pvalue} \
          --cut_off_logFC {params.cut_off_logFC} \
          --de_top_n {params.de_top_n} \
          --threads {threads} \
          --gene_symbol_col {params.gene_symbol_col:q} \
          --gene_id_type {params.gene_id_type:q} \
          --species {params.species:q} \
          --go_ontology {params.go_ontology:q} \
          --enrichment_top_n {params.enrichment_top_n} \
          --differential_r {params.differential_r:q} \
          --enrichment_r {params.enrichment_r:q} \
          --rscript {params.rscript:q}
        """
