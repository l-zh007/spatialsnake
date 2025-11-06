def annotion_help_output(run_type,file_name):
  if channel == 'single_analysis':
    if run_type=="visium_HD":
      return(os.path.join(results_folder, "{sample}_{bin}um", 'clustering',file_name))
    elif run_type=="visium" or run_type=="xenium" or run_type=="visium_segment":
      print("correct")
      return(os.path.join(results_folder, "{sample}",'clustering', file_name))
    elif run_type=="slide_seq":
      return(os.path.join(results_folder, "{sample}", 'clustering',file_name))
  if channel=="compare_analysis":
    return(os.path.join(results_folder,"merge_data",'clustering',file_name))


rule annotion_help:
  input:
    inputs=nomal_file(run_type,"clustering")
  output:
    merge=annotion_help_output(run_type,'marker_genes_pval.csv')
    # os.path.join(results_folder,"{sample}",'clustering','marker_genes_pval.csv') if channel == 'single_analysis' else os.path.join(results_folder,"merge_data",'clustering','marker_genes_pval.csv') 
  params:
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" and seg_filter==False else wildcards.sample,
    run_type = run_type,
    image_type = image_type,
    image_slice = image_slice,
    coord = "--coord {coord}" if image_slice==True else "",
    sample_cnt = len(samples) if channel=="comparision_analysis" else 1,
    markers_algorithm = markers_algorithm,
    shape_type=shape_type
  shell:
      """
      python {spatialsnake_path}workflow/scripts/cluster_visualize.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_zarr_path {output.merge} \
        --type {params.run_type} \
        --image_type {params.image_type} \
        --image_slice {params.image_slice} \
        --markers_algorithm {params.markers_algorithm} \
        --sample_cnt {params.sample_cnt} \
        --shape_type {params.shape_type} \
        {params.coord}
      """

rule enrich_rule:
  input:
    inputs = annotion_help_output(run_type,'marker_genes_pval.csv')
  output:
    merge = annotion_help_output(run_type,'kegg_data.csv')
    #os.path.join(results_folder,"{sample}",'clustering','kegg_data.csv') if channel == 'single_analysis' else os.path.join(results_folder,"merge_data",'clustering','kegg_data.csv')
  params:
    spacies = spacies,
    run_type = run_type,
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" and seg_filter==False else wildcards.sample,
  shell:
      """
      Rscript {spatialsnake_path}workflow/scripts/enrichment.R \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_path {output.merge} \
        --type {params.run_type} \
        --spacies {spacies}
      """

    


















