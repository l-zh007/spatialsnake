###################preprocessing###########################################
import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "8"
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import gc
import matplotlib.pyplot as plt
import argparse
import bbknn
import anndata
import geosketch as sketch
import warnings
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')                   

parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')

parser.add_argument('--min_genes', type=int, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--min_cells', type=int, required=False,
                   help='Path for the output zarr file') 
parser.add_argument('--variable', type=str, required=True,
                   help='Path for the output zarr file')            ########  ?????????????/bool类型错误                    
parser.add_argument('--filter_dict',type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--NEIGHBORS', type=int, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--batch_method', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--mt_threshold', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--n_top_genes', type=int, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--n_comps', type=int, required=False,
                   help='Path for the output zarr file')   
parser.add_argument('--sketch', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_rate', type=float, required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()
type=args.type
dir_path=os.path.dirname(args.output_zarr_path)
sc.settings.n_jobs = 8


if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
print("init adata")
print(adata)

if args.filter_dict!=None:
  adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
  sc.pp.calculate_qc_metrics(
    adata,
    qc_vars="mt", 
    percent_top=(10, 20, 50, 150),
    inplace=True, 
    log1p=True)

if args.filter_dict!=None:
    print(args.filter_dict,args.sample_id,"filter the each sample with different params")
    samples = adata.obs["region"].unique()
    filtered_adatas = []
    print(args.filter_dict)
    filter_dict = json.loads(args.filter_dict)
    for sample in samples:
      sample_adata = adata[adata.obs["region"] == sample].copy()
      min_cells=filter_dict[sample][0]
      min_genes=filter_dict[sample][1]
      mt_threshold=filter_dict[sample][2]
      sc.pp.filter_genes(sample_adata,min_cells=min_cells)
      sc.pp.filter_cells(sample_adata,min_counts=min_genes)
      sample_adata=sample_adata[sample_adata.obs.pct_counts_mt<mt_threshold, :]
      filtered_adatas.append(sample_adata)
      del sample_adata
    adata = anndata.concat(filtered_adatas, join='inner', index_unique=None)
else:
    min_genes=args.min_genes
    min_cells=args.min_cells
    sc.pp.filter_genes(adata,min_cells=min_cells)
    sc.pp.filter_cells(adata,min_counts=min_genes)
    maxi=float(args.mt_threshold)
    adata=adata[adata.obs.pct_counts_mt<maxi, :]


sc.pl.violin(adata=adata, keys=["log1p_total_counts"], stripplot=False, inner="box",show=False, groupby="region",)
plt.title("Total UMI by Sample")
plt.axhline(y=4, color='r', linestyle='-')
plt.axhline(y=8, color='r', linestyle='-')
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}filtered_Total_UMI.png"),
  dpi=300,
  bbox_inches='tight')
plt.close()

sc.pl.violin(adata=adata, keys=["log1p_n_genes_by_counts"], groupby="region", stripplot=False, inner="box",show=False)
plt.title("Total Genes by Sample")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}filtered_Total_Genes.png"),
    dpi=300,
    bbox_inches='tight')
plt.close()

sc.pl.violin(adata=adata, keys=["log1p_total_counts_mt"], groupby="region", stripplot=False, inner="box",show=False)
plt.title("Mitochondrial Genes by Sample")
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}_Mitochondrial_Genes.png"),
  dpi=300,
  bbox_inches='tight')
plt.close()


sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}_scatter.png"),
    dpi=300,
    bbox_inches='tight')
plt.show()
plt.close()


sc.pp.normalize_total(adata, target_sum = None)
sc.pp.log1p(adata)





if args.variable=='True':
  sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,n_top_genes=args.n_top_genes)
  sc.pl.highly_variable_genes(adata)
  plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}_highly_variable.png"),
  dpi=300,
  bbox_inches='tight')
  plt.close()
  adata = adata[:, adata.var.highly_variable]
  print("select highly gene")


sc.tl.pca(adata,use_highly_variable=False)


print(adata)

# sampling 0.25
run_sketch = args.sketch=="True"
if run_sketch:
  SAMPLING_RATE = args.sample_rate
  sketched_adatas = {}
  value=adata.obs["group"].unique()
  print(value)
  for val in value:
    subset = adata.obs["group"] == val
    subset_adata = adata[subset, :]
    N = SAMPLING_RATE * subset_adata.shape[0]
    sketch_index = sketch.gs(X=subset_adata.obsm["X_pca"], N=int(N), seed=1)
    sketched_adatas[val] = subset_adata[sketch_index, :]
  sdata = sc.concat(sketched_adatas)
  sdata.uns=adata.uns
  sdata.varm["PCs"]=adata.varm["PCs"]
  adata.write(os.path.join(dir_path, f"sketch.h5ad"))
else:
  sdata=adata
  del adata
if args.batch_method=="harmony":
    sce.pp.harmony_integrate(sdata, key="region", basis="X_pca",max_iter_harmony=20)
elif args.batch_method=="BBkNN":
    bbknn.bbknn(sdata,batch_key="region")

# sc.pp.neighbors(sdata, n_neighbors=args.NEIGHBORS, n_pcs=20)
sc.pp.neighbors(sdata, n_neighbors=args.NEIGHBORS, use_rep="X_pca",metric="correlation",n_pcs=30)

sc.pl.pca_variance_ratio(sdata, log=True,n_pcs=20)
plt.title("pca_variance_ratio")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}pca_variance_ratio.png"),
    dpi=300,
    bbox_inches='tight')
plt.show()
plt.close()

print("final adata")
print(sdata)

if type!="slide_seq":
  sdata.obs['cell_id'] = sdata.obs['cell_id'].astype(str)
  sdata.obs['region'] = sdata.obs['region'].astype('category')
  concatenated_sdata[table]=sdata
  concatenated_sdata.write(os.path.join(args.output_zarr_path),overwrite=True)
else:
  sdata.write(os.path.join(args.output_zarr_path))
