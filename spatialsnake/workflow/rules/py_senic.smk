def split_pyscenic_rankings():
    rankings = config.get("feather_input", "")
    if isinstance(rankings, (list, tuple)):
        return [str(x).strip() for x in rankings if str(x).strip()]
    return [x.strip() for x in str(rankings).split(",") if x.strip()]


rule convert_to_loom:
    input:
        inputs = input_pysenic
    output:
        loom = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.loom")
    threads: workflow_threads
    params:
        sample_id = cpdb_sample_id,
        types = run_type
    shell:
      """
      python {spatialsnake_path}workflow/scripts/pysenic.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --loom {output.loom} \
        --types {params.types}
      """

rule pyscenic_grn:
    input:
        loom = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.loom"),
        tfs = config.get("tfs_input","")
    output:
        grn = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.grn.tsv")
    threads: workflow_threads
    params:
        gene_attr = config.get("gene_attr", "var_names"),
        cell_attr = config.get("cell_attr", "cell_id")
    resources:
        mem_mb=32000 
    shell:
        """
        python {spatialsnake_path}workflow/scripts/pyscenic_numpy_compat.py arboreto \
            --num_workers {threads} \
            --output {output.grn} \
            --method grnboost2 \
            --sparse \
            --gene_attribute {params.gene_attr} \
            --cell_id_attribute {params.cell_attr} \
            {input.loom} {input.tfs}
        """

rule pyscenic_ctx:
    input:
        grn = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.grn.tsv"),
        rankings = split_pyscenic_rankings(),
        motifs = config.get("motifs_input",""),
        loom = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.loom")
    output:
        regulons = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.regulons.csv")
    threads: workflow_threads
    params:
        gene_attr = config.get("gene_attr", "var_names"),
        cell_attr = config.get("cell_attr", "cell_id")
    shell:
        """
        python {spatialsnake_path}workflow/scripts/pyscenic_numpy_compat.py pyscenic ctx {input.grn} {input.rankings} \
            --annotations_fname {input.motifs} \
            --expression_mtx_fname {input.loom} \
            --output {output.regulons} \
            --num_workers {threads} \
            --mask_dropouts \
            --gene_attribute {params.gene_attr} \
            --cell_id_attribute {params.cell_attr}
        """
#            --mode "{params.mode}" \
rule pyscenic_aucell:
    input:
        regulons = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.regulons.csv")
    output:
        aucell = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.aucell.loom")
    threads: workflow_threads
    params:
        inputs = cellPhoneDB_input if cellPhoneDB_input else config.get("senic_input", ""),
        types = run_type,
        celltype = celltype_col,
        sample_id = cpdb_sample_id,
        top_regulons = config.get("pyscenic_top_regulons", 20),
        min_regulon_genes = config.get("pyscenic_min_regulon_genes", 10)
    shell:
        """
        python {spatialsnake_path}workflow/scripts/pysenic_visualize.py \
            --regulons {input.regulons} \
            --input_dir {params.inputs} \
            --sample_id {params.sample_id} \
            --celltype {params.celltype} \
            --outputs {output.aucell} \
            --num_workers {threads} \
            --types {params.types} \
            --top_regulons {params.top_regulons} \
            --min_regulon_genes {params.min_regulon_genes}
        """
