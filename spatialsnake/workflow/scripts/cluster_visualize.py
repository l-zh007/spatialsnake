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
from spatialdata.transformations import Identity
from spatialsnake.workflow.function.plot import cluster_proportion
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

#args.input_dir="results/useful_results/ROI_Selection_1_cells_stats.zarr"
log_step(logger, 1, 5, "loading clustered SpatialData")
concatenated_sdata = spd.read_zarr(args.input_dir)
table = pick_table_key(concatenated_sdata, run_type=type, input_spec=args.input_spec)
adata = concatenated_sdata[table]
logger.info(f"Loaded {adata.n_obs} observations for marker and spatial visualization")
table_meta = extract_table_spatial_metadata(adata)
region_name = table_meta["region_name"]
instance_key = table_meta["instance_key"]
if region_name in concatenated_sdata.shapes and instance_key in adata.obs:
    adata.obs[instance_key] = adata.obs[instance_key].astype(str)
    shapes_df = concatenated_sdata.shapes[region_name]
    shapes_df.index = shapes_df.index.astype(str)
    common_index = shapes_df.index.intersection(adata.obs[instance_key])
    if len(common_index) > 0:
        if len(common_index) != len(shapes_df):
            shapes_df = shapes_df.loc[common_index]
        if len(common_index) != adata.n_obs:
            adata = adata[adata.obs[instance_key].isin(common_index)].copy()
        concatenated_sdata.shapes[region_name] = shapes_df
adata = sanitize_obs_index_name(adata, instance_key)
concatenated_sdata[table] = adata
sample_cnt=adata.obs["sample"].astype(str).nunique() if "sample" in adata.obs else 1
image_elements = list(concatenated_sdata.images.keys())
shape_elements = list(concatenated_sdata.shapes.keys())
point_elements = list(concatenated_sdata.points.keys())
valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
image_counts=len(image_elements)//sample_cnt
shape_count=len(shape_elements)//sample_cnt
shapes=[]
images=[]
systems=[]
print(args.image_type)
images = select_elements(image_elements, args.image_type)
systems = select_elements(valid_coord_systems, args.image_type)

logger.info(f"Found {len(images)} image(s), {len(shape_elements)} shape element(s), and {len(systems)} coordinate system(s)")
if shape_count>1 and type!="visium" and type != "xenium":
    shapes = select_elements(shape_elements, args.shape_type)
else:
    shapes=shape_elements

print(concatenated_sdata)
print(systems,images,shapes)

if isinstance(region_name, (list, tuple)):
  region_name = region_name[0] if region_name else None
active_vis_mode, visual_region_name = resolve_visual_target(concatenated_sdata, region_name, args.vis_mode)
point_visual_elements = []
if active_vis_mode == "shape":
  visual_elements = [visual_region_name] if visual_region_name in shape_elements else []
else:
  point_visual_elements = select_elements(point_elements, region_name)
  if not point_visual_elements:
    point_visual_elements = [visual_region_name] if visual_region_name in getattr(concatenated_sdata, "points", {}) else point_elements
  visual_elements = list(point_visual_elements)
  keep_ids = adata.obs_names.astype(str)
  point_transformations = {coord_system: Identity() for coord_system in valid_coord_systems} if valid_coord_systems else {"global": Identity()}
  for current_visual_region in point_visual_elements:
    source_points = concatenated_sdata.points[current_visual_region]
    points_df = source_points.compute() if hasattr(source_points, "compute") else source_points.copy()
    points_df = points_df.copy()
    points_df.index = points_df.index.astype(str)
    points_df = points_df.reindex(keep_ids).copy()
    points_df["clusters"] = adata.obs["clusters"].astype(str).reindex(points_df.index).fillna("unassigned").astype("category")
    points_df = points_df.dropna(subset=["x", "y"])
    concatenated_sdata.points[current_visual_region] = PointsModel.parse(
      points_df.sort_index(),
      coordinates={"x": "x", "y": "y"},
      transformations=point_transformations,
    )
total = min(len(visual_elements)*sample_cnt, len(images), max(1, len(systems)))
log_step(logger, 2, 5, f"rendering spatial cluster plots for {sample_cnt} sample(s)")
for i in range(sample_cnt):
  sys = systems[0] if sample_cnt < 2 else systems[i]
  title = images[i].replace("_hires_image", "")
  current_visual_region = visual_region_name
  if active_vis_mode == "point" and len(point_visual_elements) > 0:
    current_visual_region = point_visual_elements[0] if len(point_visual_elements) < 2 else point_visual_elements[min(i, len(point_visual_elements) - 1)]
  current_sdata = concatenated_sdata
  if str(args.image_slice).lower() == "true" and args.coord is not None:
    current_sdata = crop0(current_sdata, sys, args.coord[0], args.coord[1], args.coord[2], args.coord[3])
  fig, axes = plt.subplots(2, 1, figsize=(20, 13))
  axes = np.ravel(axes)
  current_sdata.pl.render_images(images[i]).pl.show(
    ax=axes[0],
    title="image",
    coordinate_systems=sys,
  )
  if active_vis_mode == "shape":
    current_sdata.pl.render_shapes(
      shapes[i],
      color="clusters",
      table_name=table,
      outline_alpha=0,
      fill_alpha=1,
    ).pl.show(
      ax=axes[1],
      coordinate_systems=sys,
      title=title,
    )
  else:
    current_sdata.pl.render_images(images[i]).pl.render_points(
      current_visual_region,
      color="clusters",
      size=args.point_size,
      alpha=0.8,
    ).pl.show(
      ax=axes[1],
      coordinate_systems=sys,
      title=title,
    )
  plt.savefig(
    os.path.join(dir_path, f"{images[i]}_Clusters.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
    
    
    
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
plt.savefig(os.path.join(dir_path, f"{args.sample_id}rank_genes_groups_dotplot.png"),dpi=300,bbox_inches='tight')
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
