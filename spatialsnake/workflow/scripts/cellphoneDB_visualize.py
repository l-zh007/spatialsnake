import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import anndata as ad
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--output_name', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--cell_type1', type=str, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--cell_type2', type=str, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--gene_family', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--celltype', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--degs_file_path', type=str, required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()

adata = ad.read_h5ad(args.input_dir)
output_dir=os.path.dirname(args.output_zarr_path)
if os.path.isfile(args.degs_file_path):
  means = pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"degs_analysis_means_{args.output_name}.txt"), sep="\t")
  pvals = pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"degs_analysis_pvalues_{args.output_name}.txt"), sep="\t")
  decon = pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"degs_analysis_deconvoluted_{args.output_name}.txt"), sep="\t")
  interactions=pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"degs_analysis_relevant_interactions_{args.output_name}.txt"), sep="\t")
else:
  means = pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"statistical_analysis_means_{args.output_name}.txt"), sep="\t")
  pvals = pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"statistical_analysis_pvalues_{args.output_name}.txt"), sep="\t")
  decon = pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"statistical_analysis_deconvoluted_{args.output_name}.txt"), sep="\t")
  interactions=pd.read_csv(os.path.join(output_dir,f"cellphonedb_output",f"statistical_analysis_relevant_interactions_{args.output_name}.txt"), sep="\t")



p1=kpy.plot_cpdb_heatmap(pvals=pvals, figsize=(5, 5), title="Sum of significant interactions")
print(type(p1))


p1.savefig(
    os.path.join(output_dir,f"{args.sample_id}_heatmap.png"),
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    pad_inches=0.1
)


p3=kpy.plot_cpdb(
    adata = adata,
    cell_type1 = args.cell_type1,
    cell_type2 = args.cell_type2, 
    means = means,
    pvals = pvals,
    celltype_key = args.celltype,
    figsize = (12,40),
    max_size = 5,
    highlight_size = 0.75,
    degs_analysis = False,
    standard_scale = False,
    interaction_scores = interactions,
    scale_alpha_by_interaction_scores=True)
p3.save(os.path.join(output_dir,f"{args.sample_id}_dot_plot.png"),dpi=300)




p4=kpy.plot_cpdb(
    adata=adata,
    cell_type1=args.cell_type1,
    cell_type2=args.cell_type2,
    means=means,
    pvals=pvals,
    celltype_key=args.celltype,
    gene_family=args.gene_family,
    highlight_size=1,
    figsize=(20, 8),
)
p4.save(os.path.join(output_dir,f"{args.sample_id}_dot_family_plot.png"),dpi=300)








# p=kpy.plot_cpdb_chord(
#     adata=adata,
#     cell_type1="Colorecta",
#     cell_type2=".",
#     means=means,
#     pvals=pvals,
#     deconvoluted=decon,
#     celltype_key="celltype",
#     interaction=["PTPRC", "CD40", "CLEC2D"],
#     link_kwargs={"direction": 1, "allow_twist": True, "r1": 90, "r2": 90},
#     sector_text_kwargs={"color": "black", "size": 12, "r": 105, "adjust_rotation": True},
#     legend_kwargs={"loc": "center", "bbox_to_anchor": (1, 1), "fontsize": 8},
#     link_offset=1,
#     same_producer_colors=True)
# p.savefig(os.path.join(output_dir,f"{args.sample_id}_chord_plot.png"),dpi=300)





