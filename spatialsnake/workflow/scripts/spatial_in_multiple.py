import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")
import pandas as pd
import spatialdata as spd
import scanpy as sc
import gc
from spatialdata.models import TableModel
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
from spatialsnake.workflow.function.stereoseq_selection import resolve_stereoseq_table_key
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

logger = setup_logger("merge_integrate")
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_path', nargs='+', required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--group', nargs='+', required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', nargs='+', required=False,
                   help='Path for the output zarr file')
parser.add_argument('--input_spec', nargs='+', required=False,
                   help='Stereo-seq input mode list for compare analysis')
args = parser.parse_args()

type=args.type
group=args.group
sample=args.sample_id
input_specs=args.input_spec

log_step(logger, 1, 4, f"loading {len(sample)} integrated sample(s)")

def normalize_input_spec_list(input_specs, sample_count):
  if not input_specs:
    return [None] * sample_count
  specs = list(input_specs)
  if len(specs) == 1 and sample_count > 1:
    specs = specs * sample_count
  if len(specs) != sample_count:
    raise ValueError("Stereo-seq compare_analysis requires one input_spec per sample.")
  return specs


def QC_plot(type,sdata,zarr_name):
  dir_path=os.path.dirname(zarr_name)
  for table in sdata.tables.keys():
    adata = sdata[table]
  adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
  adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
  adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
  sc.pp.calculate_qc_metrics(
      adata, 
      qc_vars=["mt", "ribo", "hb"], 
      percent_top=(10, 20, 50, 150),
      inplace=True, 
      log1p=True)
  
  if type=="xenium":
    cprobes = (
      adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100)
    cwords = (adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100)

  if type=='xenium':
    image_num=4
  else:
    image_num=2
  
  fig, axs = plt.subplots(1, image_num, figsize=(15, 4))
  axs[0].set_title("Total transcripts per cell")
  sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0])

  axs[1].set_title("Unique transcripts per cell")
  sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1])
  if type=='xenium':
    axs[2].set_title("Area of segmented cells")
    sns.histplot(
      adata.obs["cell_area"],
      kde=False,
      ax=axs[2])

    axs[3].set_title("Nucleus ratio")
    sns.histplot(
      adata.obs["nucleus_area"] / adata.obs["cell_area"],
      kde=False,
      ax=axs[3])

  plt.savefig(
      os.path.join(dir_path, "total.png"),
      dpi=300,
      bbox_inches='tight')
  plt.close()
  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts"], 
    groupby="sample", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Total Counts by Sample (log1p)")
  # Ingestion is descriptive and platform-agnostic. Fixed count thresholds are
  # therefore not drawn here; filtering thresholds are selected and displayed
  # during preprocessing.
  plt.savefig(
    os.path.join(dir_path, "total_umi_by_sample.png"),
    dpi=300, 
    bbox_inches='tight')
  plt.close()



  sc.pl.violin(
    adata=adata, 
    keys=["log1p_n_genes_by_counts"], 
    groupby="sample", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Detected Genes by Sample (log1p)")
  plt.savefig(
    os.path.join(dir_path, "total_genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()

  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts_mt"], 
    groupby="sample", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Mitochondrial Counts by Sample (log1p)")
  plt.savefig(
    os.path.join(dir_path, "genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
  
  
  sc.pl.scatter(
    adata,
    "log1p_total_counts_mt",
    "log1p_n_genes_by_counts",
    color="pct_counts_mt",
    show=False,
  )
  plt.savefig(
    os.path.join(dir_path, "scatter.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
  instance_key = adata.uns["spatialdata_attrs"].get("instance_key")
  adata.obs[instance_key] = adata.obs[instance_key].astype(str)
  adata.obs['region'] = adata.obs['region'].astype('category')
  for table in sdata.tables.keys():
    sdata[table]=adata
  return sdata


def rename_elements(element_map, prefix):
  renamed = {}
  name_mapping = {}
  for original_name, value in element_map.items():
    new_name = f"{prefix}_{original_name}"
    renamed[new_name] = value
    name_mapping[str(original_name)] = new_name
  return renamed, name_mapping


def update_region_name(region_name, name_mapping):
  if region_name is None:
    return None
  region_text = str(region_name)
  return name_mapping.get(region_text, region_text)


def prepare_table_for_merge(table, sample_id, group_id, region_name, instance_key="cell_id"):
  table = table.copy()
  table.var_names_make_unique()
  table.obs_names = table.obs_names.astype(str)
  prefixed_ids = pd.Index([f"{sample_id}_{cell_id}" for cell_id in table.obs_names], dtype="object")
  table.obs_names = prefixed_ids
  table.obs[instance_key] = prefixed_ids.astype(str)
  table.obs["cell_id"] = prefixed_ids.astype(str)
  table.obs["sample"] = sample_id
  table.obs["group"] = group_id
  if region_name is None:
    region_name = str(sample_id)
  table.obs["region"] = region_name
  table.obs["region"] = table.obs["region"].astype("category")
  table.uns.pop("spatialdata_attrs", None)
  return table


def reindex_region_element(sdata, region_name, instance_ids):
  if region_name in getattr(sdata, "shapes", {}):
    sdata.shapes[region_name].index = instance_ids
  if region_name in getattr(sdata, "points", {}):
    sdata.points[region_name].index = instance_ids
  return sdata


def standardize_table_sample(
  sdata,
  sample_id,
  group_id,
  source_table_key,
  region_name=None,
  table_key="table",
  instance_key="cell_id",
):
  table = sdata[source_table_key]
  spatial_attrs = table.uns.get("spatialdata_attrs", {}) if hasattr(table, "uns") else {}
  spatial_region = spatial_attrs.get("region")
  if isinstance(spatial_region, (list, tuple)):
    spatial_region = spatial_region[0] if spatial_region else None
  final_region = region_name if region_name is not None else spatial_region
  table = prepare_table_for_merge(table, sample_id, group_id, final_region, instance_key=instance_key)
  final_region = str(table.obs["region"].iloc[0])
  sdata = reindex_region_element(sdata, final_region, table.obs[instance_key])
  sdata.tables = {
    table_key: TableModel.parse(
      table,
      region=final_region,
      region_key="region",
      instance_key=instance_key,
    )
  }
  return sdata


def prepare_xenium_like_sample(sdata, sample_id, group_id):
  renamed_images, _ = rename_elements(getattr(sdata, "images", {}), sample_id)
  renamed_shapes, shape_mapping = rename_elements(getattr(sdata, "shapes", {}), sample_id)
  renamed_labels, _ = rename_elements(getattr(sdata, "labels", {}), sample_id)
  renamed_points, point_mapping = rename_elements(getattr(sdata, "points", {}), sample_id)
  sdata.images = renamed_images
  sdata.shapes = renamed_shapes
  sdata.labels = renamed_labels
  sdata.points = renamed_points

  table_key = next(iter(sdata.tables.keys()))
  spatial_attrs = sdata[table_key].uns.get("spatialdata_attrs", {}) if hasattr(sdata[table_key], "uns") else {}
  region_name = spatial_attrs.get("region")
  if isinstance(region_name, (list, tuple)):
    region_name = region_name[0] if region_name else None
  region_name = update_region_name(region_name, shape_mapping) or update_region_name(region_name, point_mapping)
  if region_name is None and shape_mapping:
    region_name = next(iter(shape_mapping.values()))
  if region_name is None and point_mapping:
    region_name = next(iter(point_mapping.values()))
  return standardize_table_sample(
    sdata,
    sample_id=sample_id,
    group_id=group_id,
    source_table_key=table_key,
    region_name=region_name,
    table_key="table",
    instance_key="cell_id",
  )


sdatas = []
if type=='visium_segment':
  for i in range(len(sample)):
    sdata=spd.read_zarr(args.input_path[i])
    for table in sdata.tables.values():
      table.var_names_make_unique()
      table.obs['cell_id'] = table.obs.index
      table.obs["sample"]=sample[i]
      table.obs["group"]=group[i]
    sdatas.append(sdata)
    del sdata,table
elif type=='visium_HD':
  for i in range(len(sample)):
    sdata=spd.read_zarr(args.input_path[i])
    TABLE_KEY = 'segmentation_counts'
    for shapes in sdata.shapes.keys():
      shape_key=shapes
    for table in sdata.tables.values():
        # table.obs['cell_id'] = table.obs['location_id']
        # table.var_names_make_unique()
        table.obs['cell_id'] = table.obs.index
        table.obs["sample"]=sample[i]
        table.obs["group"]=group[i]
        table.obs['region']=shape_key
        sdata.shapes[shape_key].index=table.obs.index
    del table.uns['spatialdata_attrs']
    sdata.tables={
            TABLE_KEY: TableModel.parse(
                table,
                region=shape_key, # Link table to shapes element
                region_key='region', # Column in adata.obs indicating region name
                instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
            )
        }
    sdatas.append(sdata)
elif type=="visium":
  for i in range(len(sample)):
      sdata=spd.read_zarr(args.input_path[i])
      SHAPES_KEY = sample[i]
      TABLE_KEY = 'segmentation_counts'
      for table in sdata.tables.values():
          table.obs["sample"] = sample[i]
          table.obs["group"]=group[i]
          table.obs['cell_id'] = table.obs.index
          sdata.shapes[sample[i]].index=table.obs['cell_id']
      del table.uns['spatialdata_attrs']
      sdata.tables={
              TABLE_KEY: TableModel.parse(
                  table,
                  region=SHAPES_KEY, # Link table to shapes element
                  region_key='region', # Column in adata.obs indicating region name
                  instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
              )
          }
      sdatas.append(sdata)
elif type=="xenium":
  for i in range(len(sample)):
    sdata=spd.read_zarr(args.input_path[i])
    new_images = {}
    for img_name in sdata.images.keys():
        new_name = f"{sample[i]}_{img_name}"
        new_images[new_name] = sdata.images[img_name]
    sdata.images = new_images
    new_images = {}
    for shapes_name in sdata.shapes.keys():
        new_name = f"{sample[i]}_{shapes_name}"
        new_images[new_name] = sdata.shapes[shapes_name]
    sdata.shapes = new_images
    new_images = {}
    for labels_name in sdata.labels.keys():
        new_name = f"{sample[i]}_{labels_name}"
        new_images[new_name] = sdata.labels[labels_name]
    sdata.labels = new_images

    new_images = {}
    for points_name in sdata.points.keys():
        new_name = f"{sample[i]}_{points_name}"
        new_images[new_name] = sdata.points[points_name]
    sdata.points = new_images

    new_images = {}
    SHAPES_KEY = f"{sample[i]}_{shapes_name}"
    TABLE_KEY = 'segmentation_counts'
    for table in sdata.tables.values():
        table.var_names_make_unique()
        table.obs["sample"] = sample[i]
        table.obs["group"]=group[i]
        table.obs['cell_id'] = table.obs.index.astype(str)
        table.obs['region'] = SHAPES_KEY
        table.obs['region'] = table.obs['region'].astype('category')
    if SHAPES_KEY in sdata.shapes:
        sdata.shapes[SHAPES_KEY].index = table.obs['cell_id']
    del table.uns['spatialdata_attrs']
    sdata.tables={
            TABLE_KEY: TableModel.parse(
                table,
                region=SHAPES_KEY,
                region_key='region',
                instance_key='cell_id'
            )
        }
    sdatas.append(sdata)
elif type in ["Merfish", "merscope", "MERFISH"]:
  for i in range(len(sample)):
    sdata = spd.read_zarr(args.input_path[i])
    for table in sdata.tables.values():
      table.obs["sample"] = sample[i]
      table.obs["group"] = group[i]
      table.obs['region'] = table.obs['region'].astype('category')
      table.obs.index.name = None
    sdatas.append(sdata)
elif type in ["stereoseq", "StereoSeq", "Stereo-seq"]:
  for i in range(len(sample)):
    sdata = spd.read_zarr(args.input_path[i])
    logger.info(f"Loaded Stereo-seq sample {sample[i]} with {len(sdata.tables)} table(s)")
    new_images = {}
    for img_name in sdata.images.keys():
        new_name = f"{sample[i]}_{img_name}"
        new_images[new_name] = sdata.images[img_name]
    sdata.images = new_images

    new_images = {}
    for points_name in sdata.points.keys():
        new_name = f"{sample[i]}_{points_name}"
        new_images[new_name] = sdata.points[points_name]
    sdata.points = new_images
    SHAPES_KEY = f"{sample[i]}_{points_name}"
    new_images = {}
    if sdata.shapes.keys():
      for shapes_name in sdata.shapes.keys():
          new_name = f"{sample[i]}_{shapes_name}"
          new_images[new_name] = sdata.shapes[shapes_name]
      sdata.shapes = new_images
      SHAPES_KEY = f"{sample[i]}_{shapes_name}"
    TABLE_KEY = 'segmentation_counts'
    for table in sdata.tables.values():
        instance_key = table.uns['spatialdata_attrs'].get("instance_key")
        table.var_names_make_unique()
        table.obs["sample"] = sample[i]
        table.obs["group"]=group[i]
        table.obs['region'] = SHAPES_KEY
        table.obs['region'] = table.obs['region'].astype('category')
    del table.uns['spatialdata_attrs']
    sdata.tables={
            TABLE_KEY: TableModel.parse(
                table,
                region=SHAPES_KEY,
                region_key='region',
                instance_key=instance_key
            )
        }
    sdatas.append(sdata)
log_step(logger, 2, 4, "concatenating SpatialData objects")
concatenated_sdata = spd.concatenate(sdatas, concatenate_tables=True)
logger.info(f"Merged object contains {len(concatenated_sdata.tables)} table(s)")
log_step(logger, 3, 4, "calculating merged quality-control plots")
concatenated_sdata=QC_plot(type,concatenated_sdata,args.output_zarr_path)
log_step(logger, 4, 4, f"saving merged data to {args.output_zarr_path}")
concatenated_sdata.write(args.output_zarr_path, overwrite=True)
del concatenated_sdata, sdatas
gc.collect()
logger.info("Merge integrate module completed")
