rule convert_to_loom:
    input:
        inputs = config["senic_input"]
    output:
        loom = os.path.join(results_folder,"pysenic_results",f"{sample_id}.loom")
    params:
        sample_id = sample_id
    shell:
      """
      python workflow/scripts/pysenic.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --loom {output.loom}
      """

rule pyscenic_grn:
    input:
        loom = os.path.join(results_folder,"pysenic_results",f"{sample_id}.loom"),
        tfs = config["tfs_input"]   #os.path.join("data",'hs_hgnc_tfs.txt')
    output:
        grn = os.path.join(results_folder,"pysenic_results",f"{sample_id}.grn.tsv")
    params:
        workers = config.get("num_workers", 32),
        gene_attr = config.get("gene_attribute", "var_names"),
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
        grn = os.path.join(results_folder,"pysenic_results",f"{sample_id}.grn.tsv"),
        rankings = config["feather_input"],  #os.path.join("data",'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'),
        motifs = config["motifs_input"],  #os.path.join("data",'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'),
        loom = os.path.join(results_folder,"pysenic_results",f"{sample_id}.loom")
    output:
        regulons = os.path.join(results_folder,"pysenic_results",f"{sample_id}.regulons.csv")
    params:
        workers = config.get("num_workers", 10),
        gene_attr = config.get("gene_attribute", "var_names"),
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
        loom =  os.path.join(results_folder,"pysenic_results",f"{sample_id}.loom"),
        regulons = os.path.join(results_folder,"pysenic_results",f"{sample_id}.regulons.csv")
    output:
        aucell = os.path.join(results_folder,"pysenic_results",f"{sample_id}.aucell.loom")
    params:
        workers = config.get("num_workers", 10),
        gene_attr = config.get("gene_attribute", "var_names"),
        cell_attr = config.get("cell_attr", "cell_id")
    shell:
        """
        pyscenic aucell {input.loom} {input.regulons} \
            --output {output.aucell} \
            --num_workers {params.workers} \
            --gene_attribute {params.gene_attr} \
            --cell_id_attribute {params.cell_attr}
        """
