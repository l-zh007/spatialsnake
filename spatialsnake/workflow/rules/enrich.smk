def enrich_input(option):
  if option=="enrichment":
    if channel == 'single_analysis': 
      return(os.path.join(results_folder,"{sample}",'clustering','marker_csv',"{cluster}_marker.csv.csv"))
    if channel == "comparison"
      return(directory(os.path.join(results_folder,"merge_data",'clustering','marker_csv',"{cluster}_marker.csv.csv")))
  if option=="DEseq2":
    return(os.path.join(results_folder,"merge_data","DEseq2","marker_gene.csv"))
def enrich_output(option):
  if channel == 'single_analysis':
    os.path.join(results_folder,"{sample}",'clustering','enrich_output','{cluster}_kegg_plot.png') 
  else:
    return(os.path.join(results_folder,"merge_data",'clustering','marker_csv',"{cluster}_kegg_plot.png"))
  if option=="DEseq2":
    return(os.path.join(results_folder,"merge_data","DEseq2","{cluster}_kegg_plot.png"))
rule enrich_rule:
  input:
    inputs = enrich_input(option)
  output:
    merge = 
  params:
    spacies = spacies,
    GO_ont = GO_ont,
    type = run_type,
    sample_id = lambda wildcards: "concentrate" if channel=="compare_analysis" and seg_filter==False else wildcards.sample,
  shell:
      """
      R workflow/scripts/enrichment.R.py \
        --input_dir {input.inputs} \
        --sample_id {params.sample_id} \
        --output_path {output.merge} \
        --type {params.run_type} \
        --GO_ont {params.GO_ont} \
        --cluster {wildcards.cluster} \
        --spacies {spacies}
      """
