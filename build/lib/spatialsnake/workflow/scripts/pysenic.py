import os
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import loompy
import loompy as lp
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--loom', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
args = parser.parse_args()


output_dir = os.path.dirname(args.loom)
print(output_dir)
os.makedirs(output_dir,exist_ok=True)


if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]

print(adata.obs)
print(adata.var)
print("###############")
adata_T = adata.T
print(adata_T.obs)
print(adata_T.var)
adata.write_loom(args.loom)



print(adata_T.obs)


# with loompy.connect("./Non_Lession.loom") as ds:
#     print("Row attributes available:", ds)
#     print("Column attributes available:", ds.ca)


# pyscenic aucell \
# > Non_Lession.loom \
# > endo.reg.csv \
# > --output sample_endo.loom \
# > --num_workers 10 \
# > --gene_attribute var_names \
# > --cell_id_attribute spot_id
# 
# 
# pyscenic ctx scenic.pso.tsv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname Non_Lession.loom --mode "dask_multiprocessing" --output endo.reg.csv --num_workers 36 --mask_dropouts --gene_attribute var_names --cell_id_attribute spot_id
# 
# 
# pyscenic grn --num_workers 36 --output scenic.pso.tsv --method grnboost2 --gene_attribute var_names --cell_id_attribute spot_id Non_Lession.loom hs_hgnc_tfs.txt






# resources/hs_hgnc_tfs.txt下载转录因子列表：
# hs_hgnc_tfs.txt



# https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
# 从https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/下载reference数据库文件：
# 
# hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
# 
# 
# 
# 从https://resources.aertslab.org/cistarget/motif2tf/下载人的TF注释文件：
# 
# motifs-v9-nr.hgnc-m0.001-o0.0.tbl

# https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl





#pyscenic grn --num_workers 4 --output results/pysenic_results/scenic.pso.tsv --method grnboost2 --gene_attribute var_names --cell_id_attribute cell_id results/pysenic_results/concatenated_sdata.loom data/hs_hgnc_tfs.txt
#miniconda3/envs/spatial_env/lib/python3.10/site-packages/pyscenic


# 
# # 以人类基因组hg19版本为参考，转录起始位点（TSS）上下游5kb区域的基因与基序（motifs）的排名数据，数据整合7个物种信息用于评估基因与基序结合可能性，确定转录因子和靶基因间调控关系
# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather
# # 基序到转录因子的映射关系，通过分析转录因子结合位点基序，识别可能结合特定基序的转录因子。
# wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
# # 人类基因组中所有转录因子列表
# wget  https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
