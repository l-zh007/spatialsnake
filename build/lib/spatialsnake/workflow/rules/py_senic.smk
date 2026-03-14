rule convert_to_loom:
    input:
        inputs = input_pysenic
    output:
        loom = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.loom")
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
    params:
        workers = config.get("senic_workers", 8),
        gene_attr = config.get("gene_attr", "var_names"),
        cell_attr = config.get("cell_attr", "cell_id")
    resources:
        mem_mb=32000 
    shell:
        """
        arboreto_with_multiprocessing.py \
            --num_workers {params.workers} \
            --output {output.grn} \
            --method genie3 \
            --sparse \
            --gene_attribute {params.gene_attr} \
            --cell_id_attribute {params.cell_attr} \
            {input.loom} {input.tfs}
        """

rule pyscenic_ctx:
    input:
        grn = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.grn.tsv"),
        rankings = config.get("feather_input",""),
        motifs = config.get("motifs_input",""),
        loom = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.loom")
    output:
        regulons = os.path.join(results_folder,"pysenic_results",f"{cpdb_sample_id}.regulons.csv")
    params:
        workers = config.get("senic_workers", 10),
        gene_attr = config.get("gene_attr", "var_names"),
        cell_attr = config.get("cell_attr", "cell_id")
    shell:
        """
        pyscenic ctx {input.grn} {input.rankings} \
            --annotations_fname {input.motifs} \
            --expression_mtx_fname {input.loom} \
            --output {output.regulons} \
            --num_workers {params.workers} \
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
    params:
        inputs = cellPhoneDB_input if cellPhoneDB_input else config.get("senic_input", ""),
        num_workers = config.get("senic_workers",8),
        types = run_type,
        celltype = celltype_col,
        sample_id = cpdb_sample_id
    shell:
        """
        python {spatialsnake_path}workflow/scripts/pysenic_visualize.py \
            --regulons {input.regulons} \
            --input_dir {params.inputs} \
            --sample_id {params.sample_id} \
            --celltype {params.celltype} \
            --outputs {output.aucell} \
            --num_workers {params.num_workers} \
            --types {params.types}
        """
