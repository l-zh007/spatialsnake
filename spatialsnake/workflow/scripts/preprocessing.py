###################preprocessing###########################################
import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "8"
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import matplotlib.pyplot as plt
import argparse
import bbknn
import geosketch as sketch
import warnings
from spatialsnake.workflow.function.pca_choosing import select_pca_dimensions
from spatialsnake.workflow.function.stereoseq_selection import resolve_stereoseq_table_key
warnings.filterwarnings("ignore")

STEREOSEQ_TYPES = {"stereoseq", "StereoSeq", "Stereo-seq"}


def pick_table_key(sdata, run_type=None, input_spec=None):
  table_keys = list(getattr(sdata, "tables", {}).keys())
  if not table_keys:
    raise RuntimeError("No table found in SpatialData input.")
  if run_type in STEREOSEQ_TYPES:
    return resolve_stereoseq_table_key(input_spec, table_keys)
  return "table" if "table" in table_keys else table_keys[0]

def parse_bool(value):
  if isinstance(value, bool):
    return value
  if value is None:
    return False
  value_str = str(value).strip().lower()
  if value_str in ["true", "1", "yes", "y", "t"]:
    return True
  if value_str in ["false", "0", "no", "n", "f", "none", "null", ""]:
    return False
  raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")

parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')                   

parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--input_spec', type=str, required=False,
                   help='Stereo-seq input mode passed from sample.txt')

parser.add_argument('--min_genes', type=int, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--min_cells', type=int, required=False,
                   help='Path for the output zarr file') 
parser.add_argument('--variable', type=parse_bool, required=True,
                   help='Path for the output zarr file')            
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
parser.add_argument('--sketch', type=parse_bool, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_rate', type=float, required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()
type=args.type
dir_path=os.path.dirname(args.output_zarr_path)
sc.settings.n_jobs = 8
print(args.min_genes)
print(args.min_cells)

if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  table = pick_table_key(concatenated_sdata, run_type=type, input_spec=args.input_spec)
  adata = concatenated_sdata[table]
sample_cnt=adata.obs["sample"].astype(str).nunique() if "sample" in adata.obs else 1
# Calculate QC metrics for all cases
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"], 
    percent_top=(10, 20, 50),
    inplace=True, 
    log1p=True
)


if args.filter_dict!=None:
    print(args.filter_dict,args.sample_id,"filter the each sample with different params")
    samples = adata.obs["region"].unique()
    print(args.filter_dict)
    filter_dict = json.loads(args.filter_dict)
    kept_cell_ids = set()
    for sample in samples:
      sample_name = str(sample)
      if sample_name not in filter_dict:
        raise RuntimeError(f"Missing filter params for sample '{sample_name}' in filter_dict")
      sample_adata = adata[adata.obs["region"] == sample].copy()
      min_cells = filter_dict[sample_name][0]
      min_genes = filter_dict[sample_name][1]
      mt_threshold = filter_dict[sample_name][2]
      sc.pp.filter_genes(sample_adata, min_cells=min_cells)
      sc.pp.filter_cells(sample_adata, min_counts=min_genes)
      sample_adata = sample_adata[sample_adata.obs.pct_counts_mt < mt_threshold, :]
      kept_cell_ids.update(sample_adata.obs_names.astype(str).tolist())
      del sample_adata
    if len(kept_cell_ids) == 0:
      raise RuntimeError("No cells left after per-sample filtering.")
    # Keep original AnnData container metadata (uns/obsm/varm), including spatialdata_attrs.
    keep_mask = adata.obs_names.astype(str).isin(kept_cell_ids)
    adata = adata[keep_mask, :].copy()
else:
    min_genes=args.min_genes
    min_cells=args.min_cells
    print(f"filter with {min_genes},{min_cells}")
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


if "counts" not in adata.layers:
    adata.raw = adata.copy()
else:
    raw_adata = adata.copy()
    raw_adata.X = adata.layers['counts']
    adata.raw = raw_adata
    del raw_adata

sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)



use_hvg = False
if args.variable:
  sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,n_top_genes=args.n_top_genes)
  sc.pl.highly_variable_genes(adata)
  plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}_highly_variable.png"),
  dpi=300,
  bbox_inches='tight')
  plt.close()
  use_hvg = True
  print("select highly gene")

#sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata,use_highly_variable=use_hvg)


print(adata)

# sampling 0.25
run_sketch = args.sketch
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
print(args.batch_method)
if args.batch_method=='None' or sample_cnt==1:
    print("No Batch Effect")
elif args.batch_method=="harmony":
    print("using harmony_integrate")
    sce.pp.harmony_integrate(sdata, key="region", basis="X_pca",max_iter_harmony=20)
elif args.batch_method=="BBkNN":
    bbknn.bbknn(sdata,batch_key="region")
    print("using BBkNN integrate")





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


pca_selection=select_pca_dimensions(sdata)
print("******* recommand pcs :",pca_selection)


if type!="slide_seq":
  spatial_attrs = sdata.uns.get('spatialdata_attrs', {})
  instance_key = spatial_attrs.get('instance_key')
  sdata.obs[instance_key] = sdata.obs[instance_key].astype(str)
  sdata.obs['region'] = sdata.obs['region'].astype('category')
  concatenated_sdata[table]=sdata
  concatenated_sdata.write(os.path.join(args.output_zarr_path),overwrite=True)
else:
  sdata.write(os.path.join(args.output_zarr_path))
