import os
import spatialdata as spd
import scanpy as sc
import matplotlib.pyplot as plt
import argparse
import warnings
warnings.filterwarnings("ignore")
from spatialdata.models import Image2DModel, TableModel, ShapesModel
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
import re
import pandas as pd



parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--INPUT_FIlE', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_dir', type=str, required=True,
                   help='Path for the output zarr file')

parser.add_argument('--split_by', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--barcodes', type=str, default="",
                   help='Path for the output zarr file')
parser.add_argument('--max_x', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--min_x', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--max_y', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--min_y', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--shape_elements', type=str, required=True,
                   help='Path for the output zarr file')                                    
parser.add_argument('--roi_csv', type=str, default="",
                   help='Path for roi csv file')
args = parser.parse_args()

print(args.split_by)
print("@@@@@@@@@@@")
print(args.barcodes)
def crop0(x,min_x,max_x,min_y,max_y,src):
    return spd.bounding_box_query(
        x,
        min_coordinate=[min_x, min_y],
        max_coordinate=[max_x, max_y],
        axes=("x", "y"),
        target_coordinate_system=src)
    
    
def slice_by_cluster(adata,file_type,barcode,dir_path,barcode_values=None,concatenated_sdata=None):
  obs_values = adata.obs[barcode].astype(str)
  if barcode_values is None or len(barcode_values) == 0:
    values = obs_values.unique().tolist()
    for val in values:
      subset = obs_values == str(val)
      subset_adata = adata[subset, :]
      if file_type==".zarr":
        subset_adata.uns['spatialdata_attrs'] = {
      'region': list(subset_adata.obs["region"].unique()),
      'region_key': subset_adata.uns['spatialdata_attrs']["region_key"],
      'instance_key': subset_adata.uns['spatialdata_attrs']["instance_key"]}
        concatenated_sdata[table]=subset_adata
        safe_name = re.sub(r"[^\w\.-]+", "_", str(val))
        concatenated_sdata.write(os.path.join(dir_path,f"cluster_{safe_name}.zarr"),overwrite=True)
      else:
        safe_name = re.sub(r"[^\w\.-]+", "_", str(val))
        subset_adata.write(os.path.join(dir_path,f"cluster_{safe_name}.h5ad"))
    return
  else:
    values = [str(v) for v in barcode_values]
    missing = [v for v in values if v not in set(obs_values.unique())]
    if len(missing) > 0:
      raise ValueError(f"values not found in {barcode}: {','.join(missing)}")
  subset_mask = obs_values.isin([str(v) for v in values])
  subset_adata = adata[subset_mask, :]
  safe_name = re.sub(r"[^\w\.-]+", "_", "_".join(values))
  if file_type==".zarr":
    subset_adata.uns['spatialdata_attrs'] = {
  'region': list(subset_adata.obs["region"].unique()),
  'region_key': subset_adata.uns['spatialdata_attrs']["region_key"],
  'instance_key': subset_adata.uns['spatialdata_attrs']["instance_key"]}
    concatenated_sdata[table]=subset_adata
    concatenated_sdata.write(os.path.join(dir_path,f"{barcode}_selected_{safe_name}.zarr"),overwrite=True)
  else:
    subset_adata.write(os.path.join(dir_path,f"{barcode}_selected_{safe_name}.h5ad"))
def check_barcode_exit(adata,barcode):
  if barcode in adata.obs.columns:
    return True
  return False

def resolve_systems(concatenated_sdata):
  shape_elements = list(concatenated_sdata.shapes.keys())
  valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
  if len(valid_coord_systems)>1:
      systems = [cs for cs in valid_coord_systems if cs in shape_elements]
      if not systems:
        systems = valid_coord_systems
  else:
      systems = valid_coord_systems
  return systems

def slice_by_sample(file_type,barcode,concatenated_sdata):
  systems = resolve_systems(concatenated_sdata)
  if barcode in ["samples","sample","region"]:
    for i in range(len(systems)):
      subset_sample=concatenated_sdata.filter_by_coordinate_system(systems[i])
      subset_sample.write(os.path.join(dir_path,f"{systems[i]}.zarr"),overwrite=True)
    return
  if barcode == "group":
    obs = concatenated_sdata[table].obs
    if "group" not in obs.columns:
      raise ValueError("group column not found in table obs")
    if "region" not in obs.columns:
      raise ValueError("region column not found in table obs")
    region_col = "region"
    group_values = obs["group"].unique()
    for group_val in group_values:
      group_mask = obs["group"] == group_val
      group_regions = obs.loc[group_mask, region_col].unique().tolist()
      group_region_set = set(group_regions)
      systems_for_group = [cs for cs in systems if cs in group_region_set]
      if not systems_for_group:
        continue
      if len(systems_for_group) == 1:
        subset_group = concatenated_sdata.filter_by_coordinate_system(systems_for_group[0])
      else:
        sdatas = [concatenated_sdata.filter_by_coordinate_system(cs) for cs in systems_for_group]
        subset_group = spd.concatenate(sdatas, concatenate_tables=True)
      safe_name = re.sub(r"[^\w\.-]+", "_", str(group_val))
      subset_group.write(os.path.join(dir_path,f"group_{safe_name}.zarr"),overwrite=True)

def load_roi_csvs(roi_csv):
  if not roi_csv:
    raise ValueError("roi_csv is required for ROI splitting")
  if os.path.isdir(roi_csv):
    csv_paths = [os.path.join(roi_csv, f) for f in os.listdir(roi_csv) if f.endswith(".csv")]
    csv_paths.sort()
  else:
    csv_paths = [roi_csv]
  frames = []
  for csv_path in csv_paths:
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as handle:
      lines = handle.readlines()
    header_row = None
    for idx, line in enumerate(lines):
      lower = line.lower()
      if ("barcode" in lower or "cell id" in lower) and "," in line:
        header_row = idx
        break
    if header_row is None:
      df = pd.read_csv(csv_path, sep=",", engine="python", on_bad_lines="skip")
    else:
      df = pd.read_csv(csv_path, sep=",", engine="python", skiprows=header_row)
    lower_cols = {str(col).strip().lower(): col for col in df.columns}
    if "cell_id" not in df.columns:
      if "cell id" in lower_cols:
        df = df.rename(columns={lower_cols["cell id"]: "cell_id"})
      elif "cell_id" in lower_cols:
        df = df.rename(columns={lower_cols["cell_id"]: "cell_id"})
      elif "barcode" in lower_cols:
        df = df.rename(columns={lower_cols["barcode"]: "cell_id"})
    if "cell_id" not in df.columns:
      raise ValueError(f"cell_id column not found in {csv_path}")
    name_col = None
    for candidate in ["roi", "ROI", "region", "sample", "group"]:
      if candidate in df.columns:
        name_col = candidate
        break
    if name_col:
      roi_name = df[name_col].astype(str)
    else:
      roi_name = os.path.splitext(os.path.basename(csv_path))[0]
    df_out = pd.DataFrame({
      "cell_id": df["cell_id"].astype(str),
      "roi": roi_name if isinstance(roi_name, str) else roi_name.astype(str)
    })
    frames.append(df_out)
  if not frames:
    raise ValueError("no csv files found for ROI splitting")
  return pd.concat(frames, axis=0, ignore_index=True)

def slice_by_csv(adata,file_type,dir_path,roi_csv,concatenated_sdata=None):
  all_roi = load_roi_csvs(roi_csv)
  if "cell_id" in adata.obs.columns:
    obs_ids = adata.obs["cell_id"].astype(str)
  else:
    obs_ids = adata.obs_names.astype(str)
  for roi_name, roi_df in all_roi.groupby("roi"):
    roi_cell_ids = set(roi_df["cell_id"].astype(str).tolist())
    subset_mask = obs_ids.isin(roi_cell_ids)
    subset_adata = adata[subset_mask, :]
    if subset_adata.n_obs == 0:
      continue
    safe_name = re.sub(r"[^\w\.-]+", "_", str(roi_name))
    if file_type==".zarr":
      if "spatialdata_attrs" in subset_adata.uns:
        region_key = subset_adata.uns["spatialdata_attrs"].get("region_key", "region")
        instance_key = subset_adata.uns["spatialdata_attrs"].get("instance_key", "instance_id")
      else:
        region_key = "region"
        instance_key = "instance_id"
      if instance_key not in subset_adata.obs.columns:
        instance_key = "cell_id" if "cell_id" in subset_adata.obs.columns else subset_adata.obs_names.name or "cell_id"
      if instance_key in subset_adata.obs.columns:
        subset_adata.obs[instance_key] = subset_adata.obs[instance_key].astype(str)
      if region_key in subset_adata.obs.columns:
        region_vals = list(subset_adata.obs[region_key].unique())
      elif "region" in subset_adata.obs.columns:
        region_vals = list(subset_adata.obs["region"].unique())
      else:
        region_vals = []
      subset_adata.uns["spatialdata_attrs"] = {
        "region": region_vals,
        "region_key": region_key,
        "instance_key": instance_key
      }
      print(subset_adata)
      print("####")
      concatenated_sdata[table]=subset_adata
      concatenated_sdata.write(os.path.join(dir_path,f"ROI_{safe_name}.zarr"),overwrite=True)
    else:
      subset_adata.write(os.path.join(dir_path,f"ROI_{safe_name}.h5ad"))


dir_path=args.output_dir
os.makedirs(dir_path, exist_ok=True)
paths = str(args.INPUT_FIlE)
root, file_type = os.path.splitext(paths)
print("6666",args.INPUT_FIlE)
print(paths)
print(file_type)
if file_type!=".zarr":
  adata = sc.read_h5ad(paths)
else:
  concatenated_sdata = spd.read_zarr(paths)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
    print(adata.obs)


if check_barcode_exit(adata,args.split_by) and args.split_by in ["celltype","cluster","clusters","cell_type","celltypes","cell_types",'clusters','layer','layers']:
  barcode_values = [v.strip() for v in args.barcodes.split(",") if v.strip() != ""]
  if file_type==".zarr":
    slice_by_cluster(adata,file_type,args.split_by,dir_path,barcode_values=barcode_values,concatenated_sdata=concatenated_sdata)
  else:
    slice_by_cluster(adata,file_type,args.split_by,dir_path,barcode_values=barcode_values)
elif check_barcode_exit(adata,args.split_by) and args.split_by in ["samples","sample","region","group"]:
  if file_type==".zarr":
      slice_by_sample(file_type,args.split_by,concatenated_sdata)
  else:
      slice_by_cluster(adata,file_type,args.split_by,dir_path)
elif args.split_by in ["image","images"]:
  shape_elements = list(concatenated_sdata.shapes.keys())
  valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
  if args.shape_elements and args.shape_elements in shape_elements and args.shape_elements in valid_coord_systems:
      src=args.shape_elements
  else:
      src=valid_coord_systems[0]
  image_id=f"{args.min_x}_{args.max_x}_{args.min_y}_{args.max_y}"
  subset_image=crop0(concatenated_sdata,args.min_x,args.max_x,args.min_y,args.max_y,src)
  subset_image.write(os.path.join(dir_path,f"spatial{image_id}.zarr"),overwrite=True)
  axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
  subset_image.pl.render_images().pl.show(ax=axes[0], title="image",coordinate_systems=src)
  subset_image.pl.render_images().pl.render_shapes(color="clusters").pl.show(ax=axes[1],coordinate_systems=src)
  # subset_image.pl.render_shapes(color="clusters").pl.show()
  plt.savefig(
        os.path.join(dir_path, f"{image_id}_shape.png"),
        dpi=300,
        bbox_inches='tight')
  plt.close()
elif args.split_by in ["ROI","ROIs"]:
  if file_type==".zarr":
    print(args.roi_csv,"#####")
    slice_by_csv(adata,file_type,dir_path,args.roi_csv,concatenated_sdata=concatenated_sdata)
  else:
    slice_by_csv(adata,file_type,dir_path,args.roi_csv)
