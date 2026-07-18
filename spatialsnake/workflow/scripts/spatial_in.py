import os
import spatialdata as spd
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata_io
import json
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
import argparse
import seaborn as sns
from spatialsnake.workflow.function.merfish_utils import (
    align_merfish_image_to_shape_space,
    find_merfish_transform_csv,
)
from spatialsnake.workflow.function.stereoseq_spec import parse_stereoseq_input_spec
from spatialsnake.workflow.function.stereoseq_v8 import stereoseq_v8
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

logger = setup_logger("integrate")

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
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--count_file', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--channel', type=str, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--threads', type=int, default=4,
                   help='Threads allocated to this integration job')
########## args of visium HD
parser.add_argument("--bin_size", type=str, required=False, help="bin")
parser.add_argument("--input_spec", type=str, required=False, help="Stereo-seq input mode")

####### args of xenium
parser.add_argument("--cells_boundaries", type=parse_bool, required=False, help="bin")
parser.add_argument("--nucleus_boundaries", type=parse_bool, required=False, help="bin")
parser.add_argument("--nucleus_labels", type=parse_bool, required=False, help="bin")
parser.add_argument("--morphology_mip", type=parse_bool, required=False, help="bin")

##########args of segment visium
parser.add_argument("--scale_factors", type=str, required=False, help="bin")
parser.add_argument("--image", type=str, required=False, help="bin")
parser.add_argument("--geojson", type=str, required=False, help="bin")


####### args of merscope/merfish
parser.add_argument("--merscope_z_layers", type=int, required=False, help="Optional z layers for spatialdata_io.merscope")
parser.add_argument("--merscope_region_name", type=str, required=False, help="Optional region name for spatialdata_io.merscope")
parser.add_argument("--merscope_transcripts", type=parse_bool, required=False, help="Load transcripts in spatialdata_io.merscope")
parser.add_argument("--merscope_cells_boundaries", type=parse_bool, required=False, help="Load cell boundaries in spatialdata_io.merscope")
parser.add_argument("--merscope_cells_table", type=parse_bool, required=False, help="Load cells table in spatialdata_io.merscope")
parser.add_argument("--merscope_mosaic_images", type=parse_bool, required=False, help="Load mosaic images in spatialdata_io.merscope")

args = parser.parse_args()


count_file= args.count_file
type=args.type
log_step(logger, 1, 3, f"preparing {type} input for sample {args.sample_id}")
_input_norm = os.path.normpath(args.input_dir)
if type == "visium_HD":
  parts = _input_norm.split(os.path.sep)
  if "binned_outputs" in parts:
    idx = parts.index("binned_outputs")
    real_dir = os.path.join(*parts[:idx]) if idx > 0 else os.path.dirname(_input_norm)
  else:
    real_dir = os.path.dirname(_input_norm)
else:
  real_dir = _input_norm if os.path.isdir(_input_norm) else os.path.dirname(_input_norm)


def QC_plot(type,sdata,zarr_name):
  dir_path=os.path.dirname(zarr_name)
  for table in sdata.tables.keys():
    adata = sdata[table]
  adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
  sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars="mt",
    percent_top=(10, 20, 50),
    inplace=True, 
    log1p=True)
  
  if type=="xenium":
    cprobes = (
      adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100)
    cwords = (adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100)
    logger.info(f"Negative DNA probe count percent: {cprobes:.3f}")
    logger.info(f"Negative decoding count percent: {cwords:.3f}")

  if type=='xenium':
    image_num=4
  else:
    image_num=2
  _, axs = plt.subplots(1, image_num, figsize=(15, 4))
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
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Total Counts by Region (log1p)")
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
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Detected Genes by Region (log1p)")
  plt.savefig(
    os.path.join(dir_path, "total_genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()

  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts_mt"], 
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Mitochondrial Counts by Region (log1p)")
  plt.savefig(
    os.path.join(dir_path, "genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
  
  
  sc.pl.scatter(
    adata,
    "total_counts",
    "n_genes_by_counts",
    color="pct_counts_mt",
    show=False,
  )
  plt.savefig(
    os.path.join(dir_path, "scatter.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
  sdata[table]=adata
  return sdata


def align_shapes_with_table_cell_id(sdata, sample_id=None, table_key="table"):
  shape_keys = list(getattr(sdata, "shapes", {}).keys())
  table_items = list(getattr(sdata, "tables", {}).items())
  if not shape_keys or not table_items:
    return sdata

  shape_key = shape_keys[0]
  _, table = table_items[0]
  table.var_names_make_unique()
  table.obs_names = table.obs_names.astype(str)
  table.obs["cell_id"] = table.obs_names
  if sample_id is not None:
    table.obs["sample"] = sample_id
    table.obs["group"] = sample_id
  table.obs["region"] = shape_key
  table.obs["region"] = table.obs["region"].astype("category")
  sdata.shapes[shape_key].index = table.obs["cell_id"]
  table.uns.pop("spatialdata_attrs", None)
  sdata.tables = {
    table_key: TableModel.parse(
      table,
      region=shape_key,
      region_key="region",
      instance_key="cell_id",
    )
  }
  return sdata












def create_zarr_bin(path_to_inputs,sample_id,zarr_name,filtered_counts_file,bin_size):
    logger.info(f"Reading Visium HD data from {path_to_inputs}")
    sdata = spatialdata_io.visium_hd(path_to_inputs, dataset_id=sample_id, filtered_counts_file=filtered_counts_file, bin_size=bin_size)
    sdata = align_shapes_with_table_cell_id(
        sdata,
        sample_id=sample_id,
        table_key="table",
    )
    sdata=QC_plot(type,sdata,zarr_name)
    logger.info(f"Writing Visium HD SpatialData to {zarr_name}")
    sdata.write(zarr_name, overwrite=True)

def creat_zarr_visium(path_to_inputs,sample_id,zarr_name,h5_name):
    logger.info(f"Reading Visium data from {path_to_inputs}")
    sdata=spatialdata_io.visium(path_to_inputs,dataset_id=sample_id,counts_file=h5_name)
    sdata=QC_plot(type,sdata,zarr_name)
    SHAPES_KEY = sample_id
    TABLE_KEY = 'table'
    for table in sdata.tables.values():
          table.obs["sample"] = sample_id
          table.obs["group"] = sample_id
          table.obs['cell_id'] = table.obs.index
          sdata.shapes[sample_id].index=table.obs['cell_id']
    del table.uns['spatialdata_attrs']
    sdata.tables={
              TABLE_KEY: TableModel.parse(
                  table,
                  region=SHAPES_KEY, # Link table to shapes element
                  region_key='region', # Column in adata.obs indicating region name
                  instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
              )
          }
    logger.info(f"Prepared Visium table with {sdata['table'].n_obs} observations")
    sdata.write(zarr_name, overwrite=True)

def create_zarr(count_matrix_path,image_path,scale_factors_path,geojson_path,sample_name,zarr_name):
    logger.info(f"Reading segmented Visium count matrix from {count_matrix_path}")
    COUNT_MATRIX_PATH = count_matrix_path
    IMAGE_PATH = image_path
    SCALE_FACTORS_PATH = scale_factors_path
    GEOJSON_PATH = geojson_path

    # Load AnnData
    adata = sc.read_10x_h5(COUNT_MATRIX_PATH)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.obs.index = sample_name +"_" + adata.obs.index.astype(str)

    # Load and preprocess image data
    image_data = np.array(Image.open(IMAGE_PATH))
    if image_data.ndim == 2:
        image_data = image_data[np.newaxis, :, :] # Add channel dimension for grayscale
    elif image_data.ndim == 3:
        image_data = np.transpose(image_data, (2, 0, 1)) # (H, W, C) -> (C, H, W) for spatialdata

    # Load scale factors
    with open(SCALE_FACTORS_PATH, 'r') as f:
        scale_data = json.load(f)

    # Load GeoJSON data
    with open(GEOJSON_PATH, 'r') as f:
        geojson_data = json.load(f)

    # Define coordinate systems:
    # `downscale_to_hires`: The coordinate system where shapes are located, scaled relative to the hires resolution.

    hires_scale = scale_data['tissue_hires_scalef']

    # Transformation for shapes (from pixel to downscale_to_hires)
    shapes_transformations = {
       "downscale_to_hires": Scale(np.array([hires_scale, hires_scale]), axes=("x", "y")) # if the high-resolution microscope image is being used and Identity() transform would be performed.
    }

    # Transformation for the 'hires_tissue_image' (it's already in the 'downscale_to_hires' space visually)
    image_transformations = {
        "downscale_to_hires": Identity()
    }

    # Process Cell Segmentation (GeoJSON) and Integrate with AnnData

    # Create a mapping from adata.obs.index to geojson features
    geojson_features_map = {
        f"{sample_name}_cellid_{feature['properties']['cell_id']:09d}-1": feature
        for feature in geojson_data['features']
    }

    # Prepare data for GeoDataFrame and update adata.obs
    geometries = []
    cell_ids_ordered = []

    for obs_index_str in adata.obs.index:
        feature = geojson_features_map.get(obs_index_str)
        if feature:
            # Create shapely Polygon from coordinates
            polygon_coords = np.array(feature['geometry']['coordinates'][0])
            geometries.append(Polygon(polygon_coords))
            cell_ids_ordered.append(obs_index_str)
        else:
            geometries.append(None) # Or a suitable placeholder
            cell_ids_ordered.append(obs_index_str)

    # Remove None entries if any (or handle them upstream)
    valid_indices = [i for i, geom in enumerate(geometries) if geom is not None]
    geometries = [geometries[i] for i in valid_indices]
    cell_ids_ordered = [cell_ids_ordered[i] for i in valid_indices]


    # Create GeoDataFrame for shapes
    shapes_gdf = gpd.GeoDataFrame({
        'cell_id': cell_ids_ordered,
        'geometry': geometries
    }, index=cell_ids_ordered)
    # Update adata.obs with cluster information and spatial identifiers
    adata.obs['cell_id'] = adata.obs.index
    adata.obs['region'] = sample_name + '_cell_boundaries'
    adata.obs['region'] = adata.obs['region'].astype('category')
    adata.obs['group'] = sample_name
    adata = adata[shapes_gdf.index].copy() # Filter adata to match shapes_gdf

    # Define names for SpatialData elements
    IMAGE_KEY =  sample_name + '_hires_tissue_image'
    TABLE_KEY =  'table'
    SHAPES_KEY = sample_name + '_cell_boundaries'

    # Create SpatialData elements directly
    sdata = spd.SpatialData(
        images={
            IMAGE_KEY: Image2DModel.parse(image_data, transformations=image_transformations)
        },
        tables={
            TABLE_KEY: TableModel.parse(
                adata,
                region=SHAPES_KEY, # Link table to shapes element
                region_key='region', # Column in adata.obs indicating region name
                instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
            )
        },
        shapes={
            SHAPES_KEY: ShapesModel.parse(shapes_gdf, transformations=shapes_transformations)
        }
    )
    if args.channel=="single_analysis":
      sdata=QC_plot(type,sdata,zarr_name)
    sdata.write(zarr_name, overwrite=True)

def create_zarr_xenium(path_to_inputs,zarr_name,cells_boundaries,nucleus_boundaries,nucleus_labels,morphology_mip):
    logger.info(f"Reading Xenium data from {path_to_inputs}")
    sdata = spatialdata_io.xenium(path_to_inputs,cells_boundaries=cells_boundaries, nucleus_boundaries=nucleus_boundaries,nucleus_labels=nucleus_labels,morphology_mip=morphology_mip,n_jobs=args.threads,cells_as_circles=True,morphology_focus = True)
    new_images = {}
    for img_name in sdata.images.keys():
        new_name = f"{args.sample_id}_{img_name}"
        new_images[new_name] = sdata.images[img_name]
    sdata.images = new_images
    new_images = {}
    for shapes_name in sdata.shapes.keys():
        new_name = f"{args.sample_id}_{shapes_name}"
        new_images[new_name] = sdata.shapes[shapes_name]
    sdata.shapes = new_images
    new_images = {}
    for labels_name in sdata.labels.keys():
        new_name = f"{args.sample_id}_{labels_name}"
        new_images[new_name] = sdata.labels[labels_name]
    sdata.labels = new_images

    new_images = {}
    for points_name in sdata.points.keys():
        new_name = f"{args.sample_id}_{points_name}"
        new_images[new_name] = sdata.points[points_name]
    sdata.points = new_images

    new_images = {}
    SHAPES_KEY = f"{args.sample_id}"+'_cell_circles'
    TABLE_KEY = 'table'
    for table in sdata.tables.values():
        table.var_names_make_unique()
        table.obs['cell_id'] = table.obs['cell_id'].astype(str)
        table.obs["sample"] = args.sample_id
        table.obs['region'] = SHAPES_KEY
        table.obs['region'] = table.obs['region'].astype('category')
        logger.info(f"Prepared Xenium table with {table.n_obs} observations")
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
    sdata=QC_plot(type,sdata,zarr_name)
    sdata.write(zarr_name,overwrite=True)
    
def _create_mock_spatialdata(sample_id: str, platform: str):
    """
    Create a tiny SpatialData object for smoke testing imaging platform branches.
    This is used only when MOCK_DATASET.txt exists under --input_dir.
    """
    import numpy as np
    import pandas as pd
    import anndata as ad
    from spatialdata import SpatialData
    from spatialdata.models import TableModel, Image2DModel

    # ----- tiny counts -----
    genes = [f"Gene{i}" for i in range(5)]
    cells = [f"{platform}_cell{i}" for i in range(5)]
    X = np.arange(25, dtype=np.float32).reshape(5, 5)

    region_value = str(sample_id)

    obs = pd.DataFrame(index=cells)
    # region ����ͬʱ��1) obs ���� region �� 2) TableModel.parse ��ʽ�� region=...
    obs["region"] = pd.Categorical([region_value] * len(cells))
    obs["cell_id"] = obs.index
    obs["group"] = region_value

    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X, obs=obs, var=var)

    # ----- coordinates -----
    coords = np.stack([np.arange(len(cells)), np.arange(len(cells))], axis=1).astype(np.float32)
    adata.obsm["spatial"] = coords

    table = TableModel.parse(
        adata,
        region=region_value,
        region_key="region",
        instance_key="cell_id",
    )

    # ----- tiny image (10x10) -----
    img = (np.arange(100).reshape(10, 10) % 255).astype(np.uint8)
    img = img[None, :, :]  # ���� channel ά�� -> (1, 10, 10)
    img = Image2DModel.parse(img, dims=("c", "y", "x"))

    sdata = SpatialData(tables={"table": table}, images={"mock_image": img})
    return sdata



def create_zarr_merscope(path_to_inputs, zarr_name):
    """Read Vizgen MERSCOPE / MERFISH output directory via spatialdata-io."""
    logger.info(f"Reading MERSCOPE/MERFISH data from {path_to_inputs}")
    # Mock mode
    if os.path.exists(os.path.join(path_to_inputs, 'MOCK_DATASET.txt')):
        sdata = _create_mock_spatialdata(args.sample_id, 'merscope')
    else:
        sdata = spatialdata_io.merscope(
            path_to_inputs,
            z_layers=args.merscope_z_layers,
            transcripts=args.merscope_transcripts,
            cells_boundaries=args.merscope_cells_boundaries,
            cells_table=args.merscope_cells_table,
            mosaic_images=args.merscope_mosaic_images,
        )

    for _, table in sdata.tables.items():
        table.obs['cell_id'] = table.obs.index
        if 'sample' not in table.obs.columns:
            table.obs['sample'] = args.sample_id
        else:
            table.obs['sample'] = table.obs['sample'].astype(str)
    logger.info(f"Loaded MERSCOPE/MERFISH object with {len(sdata.tables)} table(s)")
    transform_csv = find_merfish_transform_csv(path_to_inputs)
    if transform_csv and os.path.isfile(transform_csv):
        sdata = align_merfish_image_to_shape_space(
            sdata,
            transform_csv=transform_csv,
            coordinate_system="global",
            register_points=True,
        )
    else:
        logger.warning(f"MERFISH transform csv not found, skip alignment: {transform_csv}")

    if args.channel == "single_analysis":
        sdata = QC_plot(type, sdata, zarr_name)
    sdata.write(zarr_name, overwrite=True)

def create_zarr_stereoseq(path_to_inputs, zarr_name, bin_size=None):
  logger.info(f"Reading Stereo-seq data from {path_to_inputs}")
  input_specs = parse_stereoseq_input_spec(args.input_spec if args.input_spec not in [None, ""] else bin_size)
  sdata = stereoseq_v8(
      path_to_inputs,
      bin_sizes=input_specs,
      load_analysis=False,
  )

  for _, table in sdata.tables.items():
    spatial_attrs = table.uns.get('spatialdata_attrs', {}) if hasattr(table, 'uns') else {}
    instance_key = spatial_attrs.get('instance_key')
    region_key = spatial_attrs.get('region_key', 'region')
    region_name = spatial_attrs.get('region', args.sample_id)
    if instance_key == 'cell_id' and 'cell_id' not in table.obs.columns:
      table.obs['cell_id'] = table.obs.index
    if 'cell_id' not in table.obs.columns and table.obs.index.is_unique:
      table.obs['cell_id'] = table.obs.index.astype(str)
    if region_key not in table.obs.columns:
      table.obs[region_key] = region_name
    if 'sample' not in table.obs.columns:
      table.obs['sample'] = args.sample_id
    table.obs['group'] = args.sample_id
  if args.channel == "single_analysis":
    sdata = QC_plot(type, sdata, zarr_name)
  sdata.write(zarr_name, overwrite=True)

if type=='visium_segment':
        log_step(logger, 2, 3, "reading segmented Visium files")
        # ``--input_dir`` is resolved from sample.txt by the Snakemake rule and
        # may point either to segmented_outputs or directly to the count file.
        # All companion files must therefore be resolved relative to that user
        # path rather than to a hard-coded data/{sample_id} directory.
        segmented_dir = real_dir
        count_file = _input_norm if os.path.isfile(_input_norm) else os.path.join(segmented_dir, count_file)
        image_file=os.path.join(segmented_dir, 'spatial', args.image)
        scale_factors_file=os.path.join(segmented_dir, 'spatial', args.scale_factors)
        geojson_file=os.path.join(segmented_dir, args.geojson)
        required_paths = [count_file, image_file, scale_factors_file, geojson_file]
        missing_paths = [path for path in required_paths if not os.path.exists(path)]
        if missing_paths:
          raise FileNotFoundError("Missing segmented Visium input(s): " + ", ".join(missing_paths))
        create_zarr(count_matrix_path=count_file,
                image_path=image_file,
                scale_factors_path=scale_factors_file,
                geojson_path=geojson_file,
                sample_name=args.sample_id,
                zarr_name=args.output_zarr_path)
elif type=='visium_HD':
      log_step(logger, 2, 3, "reading Visium HD files")
      is_filtered = True if count_file == "filtered_feature_bc_matrix.h5" else False
      create_zarr_bin(path_to_inputs=real_dir,
                sample_id=args.sample_id,
                zarr_name=args.output_zarr_path,
                filtered_counts_file=is_filtered,
                bin_size=args.bin_size)
elif type=="visium":
      log_step(logger, 2, 3, "reading Visium files")
      creat_zarr_visium(
                path_to_inputs=real_dir,
                sample_id=args.sample_id,
                zarr_name=args.output_zarr_path,
                h5_name=count_file)
elif type=="xenium":
      log_step(logger, 2, 3, "reading Xenium files")
      create_zarr_xenium(path_to_inputs=real_dir,zarr_name=args.output_zarr_path,cells_boundaries=args.cells_boundaries,nucleus_boundaries=args.nucleus_boundaries,nucleus_labels=args.nucleus_labels,morphology_mip=args.morphology_mip)
elif type in ["Merfish", "merscope", "MERFISH"]:
      log_step(logger, 2, 3, "reading MERSCOPE/MERFISH files")
      create_zarr_merscope(path_to_inputs=real_dir, zarr_name=args.output_zarr_path)

elif type in ["stereoseq", "StereoSeq", "Stereo-seq"]:
      log_step(logger, 2, 3, "reading Stereo-seq files")
      create_zarr_stereoseq(path_to_inputs=real_dir, zarr_name=args.output_zarr_path, bin_size=args.bin_size)

log_step(logger, 3, 3, f"integrated data saved to {args.output_zarr_path}")
logger.info("Integrate module completed")



      
