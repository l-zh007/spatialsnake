import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")
import spatialdata as spd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from spatialdata.models import PointsModel
from spatialdata.transformations import get_transformation
from spatialsnake.workflow.function.plot import cluster_proportion, plot_spatial_clusters
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.stereoseq_selection import (
    extract_table_spatial_metadata,
    resolve_stereoseq_table_key,
    resolve_visual_target,
)
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step
import spatialdata_plot
import argparse
STEREOSEQ_TYPES = {"stereoseq", "StereoSeq", "Stereo-seq"}
logger = setup_logger("annotation_help")

parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--image_type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--coord', type=float, nargs=4, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_cnt', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--markers_algorithm', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--shape_type', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--image_slice', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--vis_mode', type=str, required=False,
                   help='Spatial visualization mode: auto, point, or shape')
parser.add_argument('--point_size', type=float, required=False, default=2.5,
                   help='Point size for point-based spatial visualization')
parser.add_argument('--input_spec', type=str, required=False,
                   help='Stereo-seq input mode passed from sample.txt')
parser.add_argument('--threads', type=int, default=4,
                   help='Threads allocated to marker analysis and visualization')
args = parser.parse_args()
type=args.type
if args.threads < 1:
    parser.error("--threads must be >= 1")
sc.settings.n_jobs = args.threads
dir_path=os.path.dirname(args.output_zarr_path)

def pick_table_key(sdata, run_type=None, input_spec=None):
    table_keys = list(getattr(sdata, "tables", {}).keys())
    if not table_keys:
        raise RuntimeError("No table found in SpatialData input.")
    if run_type in STEREOSEQ_TYPES:
        return resolve_stereoseq_table_key(input_spec, table_keys)
    return "table" if "table" in table_keys else table_keys[0]



def select_elements(items, keyword):
    if not keyword:
        return list(items)
    selected = [item for item in items if keyword in item]
    return selected if selected else list(items)

def sanitize_obs_index_name(adata, instance_key):
    index_name = adata.obs.index.name
    if not instance_key or instance_key not in adata.obs.columns or index_name != instance_key:
        return adata
    safe_index_name = "__obs_index__"
    suffix = 1
    while safe_index_name in adata.obs.columns:
        safe_index_name = f"__obs_index__{suffix}"
        suffix += 1
    adata.obs_names.name = safe_index_name
    return adata

def crop0(x,crs,x1,x2,y1,y2):
    return spd.bounding_box_query(
        x,
        min_coordinate=[x1, y1],
        max_coordinate=[x2, y2],
        axes=("x", "y"),
        target_coordinate_system=crs)

#args.input_dir="results/useful_results/data_splitting/ROI_Selection_1_cells_stats.zarr"
log_step(logger, 1, 5, "loading clustered SpatialData")
concatenated_sdata = spd.read_zarr(args.input_dir)
table = pick_table_key(concatenated_sdata, run_type=type, input_spec=args.input_spec)
adata = concatenated_sdata[table]
logger.info(f"Loaded {adata.n_obs} observations for marker and spatial visualization")
table_meta = extract_table_spatial_metadata(adata)
region_name = table_meta["region_name"]
instance_key = table_meta["instance_key"]
adata = sanitize_obs_index_name(adata, instance_key)
concatenated_sdata[table] = adata
log_step(logger, 2, 5, "rendering spatial cluster plots")
plot_spatial_clusters(
    concatenated_sdata, table, "clusters", dir_path, args.sample_id,
    image_type=args.image_type, vis_mode=args.vis_mode, point_size=args.point_size,
    coord=args.coord if str(args.image_slice).lower() == "true" else None,
)

log_step(logger, 3, 5, "saving cluster proportion plot")

cluster_proportion(adata, sample_col='region', 
                      cluster_col='clusters',
                      figsize=(12, 8),
                      palette='tab20',
                      sort_samples=None,
                      sort_clusters=None,
                      title="Clusters_proportion",
                      save_path=os.path.join(dir_path, f"Clusters_proportion.png"),
                      dpi=300)

log_step(logger, 4, 5, f"ranking marker genes with {args.markers_algorithm}")
markers_algorithm = args.markers_algorithm
sc.tl.rank_genes_groups(
    adata=adata,
    groupby="clusters",
    method=markers_algorithm,
    use_raw=False,
)
sc.pl.rank_genes_groups_dotplot(
    adata=adata, 
    groupby="clusters", 
    standard_scale="var", 
    n_genes=5,
    use_raw=False,
    show=False)
plt.savefig(os.path.join(dir_path, f"{args.sample_id}_marker_dotplot.png"),dpi=300,bbox_inches='tight')
plt.close()



log_step(logger, 5, 5, "exporting marker tables and cluster assignments")
export_cluster_csv(
    concatenated_sdata,
    args.type,
    dir_path,
    cell_id_col="cell_id",
    info_col="clusters",
    sample_col="sample",
    sample_id=args.sample_id
)

df_marker_genes = sc.get.rank_genes_groups_df(adata = adata,group = None,pval_cutoff=0.05)
output_filename = f"marker_genes_pval.csv"
output_path = os.path.join(dir_path,output_filename)
df_marker_genes.to_csv(output_path)



for group_value, sub_df in df_marker_genes.groupby("group"):
    sub_df_filtered = sub_df.drop(columns=["group"])
    output_filename = f"cluster_{group_value}.csv"
    output_path = os.path.join(dir_path,f'{group_value}')
    os.makedirs(os.path.join(output_path), exist_ok=True)
    sub_df_filtered.to_csv(os.path.join(output_path,output_filename), index=False)
logger.info("Annotation help module completed")
