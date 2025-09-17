import os
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

output_dir = "DEG_results"
os.makedirs(output_dir, exist_ok=True)


for key, inputs in samples.items():
  if type=="slide_seq":
    print(os.path.join(pre_dir,inputs[-1]))
    adata = sc.read_h5ad(os.path.join(pre_dir,inputs[-1]))
  else:
    concatenated_sdata = spd.read_zarr(os.path.join(pre_dir,inputs[-1]))
    for table in concatenated_sdata.tables.keys():
      print(table)
      adata = concatenated_sdata[table]

  sc.tl.rank_genes_groups(adata = adata, groupby="clusters", method="wilcoxon")



  fig = sc.pl.rank_genes_groups_dotplot(
    adata=adata, 
    groupby="clusters", 
    standard_scale="var", 
    n_genes=5,
    show=False)
  plt.savefig(
    os.path.join(output_dir, f"{key}rank_genes_groups_dotplot.png"),
    dpi=300,
    bbox_inches='tight')
  df_marker_genes = sc.get.rank_genes_groups_df(adata = adata,group = None,pval_cutoff=0.05)
  output_filename = f"marker_genes_pval.csv"
  output_path = os.path.join(output_dir, output_filename)
  df_marker_genes.to_csv(output_path)



  for group_value, sub_df in df_marker_genes.groupby("clusters"):
      sub_df_filtered = sub_df.drop(columns=["clusters"])
      output_filename = f"cluster_{group_value}.csv"
      output_path = os.path.join(output_dir, output_filename)
      sub_df_filtered.to_csv(output_path, index=False)




  gene_name = ["AREG", "MET"]
  for name in gene_name:
    sdata.pl.render_images("morphology_focus").pl.render_shapes(
        "cell_circles",
        color=name,
        table_name="table",
        use_raw=False,
    ).pl.show(
        title=f"{name} expression over Morphology image",
        coordinate_systems="global",
        figsize=(10, 5),
    )







