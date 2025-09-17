import pandas as pd
import glob
import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "1"
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

import json
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils
compare=False



display(HTML(db_releases_utils.get_remote_database_versions_html()['db_releases_html_table']))
cpdb_version ='v5.0.0'
cpdb_target_dir = os.path.join('./cellphonedb_v500_NatProtocol/', cpdb_version)
from cellphonedb.utils import db_utils
db_utils.download_database(cpdb_target_dir, cpdb_version)


if compare:
  concatenated_sdata = spd.read_zarr("concatenated_sdata")
  unique_samples = adata.obs['region'].unique()
  for sample in unique_samples:
      adata_sample = adata[adata.obs['region'] == sample, :].copy()
      h5ad_path = os.path.join(output_dir, f"{region}.h5ad")
      adata_sample.write_h5ad(h5ad_path)
      print(f"已保存 AnnData 子集：{h5ad_path}")
    
      df_extract = adata_sample.obs[['cell_id', 'grouped_clusters']].copy()
      txt_path = os.path.join(output_dir, f"{sample}_cellid_groupedclusters.txt")
      df_extract.to_csv(txt_path, sep="\t", index=False)
      print(f"已保存两列数据：{txt_path}")
  exit()



samples = {"Colon_Cancer_P1":["data/Visium_HD_Human_Colon_Cancer_P1",8,"Colon_Cancer_P1.zarr"],
            "Colon_Cancer_P2":["data/Visium_HD_Human_Colon_Cancer_P2",8,"Colon_Cancer_P2.zarr"]}
for key, inputs in samples.items():
if type=="slide_seq":
  print("secceed")
else:
  for key, inputs in samples.items():
    concatenated_sdata = spd.read_zarr(os.path.join(pre_dir,inputs[-1]))
    for table in concatenated_sdata.tables.keys():
      print(table)
      adata = concatenated_sdata[table]
      file_name=inputs[-1].split(".")[0]
      adata.write(f"{file_name}.h5ad")
      if type!="visium":
          df_extract = adata_sample.obs[['cell_id', 'cell_type']].copy()
      else:
          df_extract = adata_sample.obs[['spot_id', 'cell_type']].copy()
      txt_path = os.path.join(output_dir, f"{key}_cellid_cell_type.txt")
      df_extract.to_csv(txt_path, sep="\t", index=False)
      print(f"已保存两列数据：{txt_path}")











