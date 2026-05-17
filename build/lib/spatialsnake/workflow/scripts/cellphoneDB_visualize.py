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
parser.add_argument('--interaction_pairs', type=str, required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()

adata = ad.read_h5ad(args.input_dir)
output_dir=os.path.dirname(args.output_zarr_path)
base_dir = os.path.join(output_dir, "cellphonedb_output")
prefix = "degs_analysis" if args.degs_file_path and os.path.isfile(args.degs_file_path) else "statistical_analysis"

def pick_file(candidates, required):
  for name in candidates:
    path = os.path.join(base_dir, name)
    if os.path.isfile(path):
      return path
  if required:
    raise FileNotFoundError(f"Missing required file. Tried: {candidates}")
  return None

output_name = args.output_name
means_path = pick_file([f"{prefix}_means_{output_name}.txt"], True)
pvals_path = pick_file([f"{prefix}_pvalues_{output_name}.txt"], True)
decon_path = pick_file([f"{prefix}_deconvoluted_{output_name}.txt", f"{prefix}_deconvoluted_percents_{output_name}.txt"], False)
interactions_path = pick_file([f"{prefix}_interaction_scores_{output_name}.txt", f"{prefix}_relevant_interactions_{output_name}.txt"], False)

means = pd.read_csv(means_path, sep="\t")
pvals = pd.read_csv(pvals_path, sep="\t")
decon = pd.read_csv(decon_path, sep="\t") if decon_path else None
interactions = pd.read_csv(interactions_path, sep="\t") if interactions_path else pd.DataFrame()

def split_values(value):
  if value is None:
    return []
  return [v.strip() for v in str(value).split(",") if v.strip()]

def resolve_celltype_key(adata, user_key):
  if user_key and user_key in adata.obs.columns:
    return user_key
  for candidate in ["celltype", "cell_type", "celltypes", "cell_type_annotation"]:
    if candidate in adata.obs.columns:
      return candidate
  raise ValueError("celltype column not found in adata.obs")

def resolve_celltypes(adata, celltype_key, cell_type1, cell_type2):
  value_counts = adata.obs[celltype_key].value_counts()
  all_types = value_counts.index.astype(str).tolist()
  if len(all_types) == 0:
    raise ValueError("no cell types found in adata.obs")
  ct1 = cell_type1 if cell_type1 in all_types else ""
  ct2 = cell_type2 if cell_type2 in all_types else ""
  if not ct1 and not ct2:
    ct1 = all_types[0]
    ct2 = all_types[1] if len(all_types) > 1 else all_types[0]
  elif ct1 and not ct2:
    ct2 = next((t for t in all_types if t != ct1), ct1)
  elif ct2 and not ct1:
    ct1 = next((t for t in all_types if t != ct2), ct2)
  return ct1, ct2

def resolve_interaction_pairs(interactions, pvals, user_pairs, max_items=8):
  pairs = split_values(user_pairs)
  if pairs:
    return pairs[:max_items]
  if "interacting_pair" in interactions.columns:
    candidates = interactions["interacting_pair"].dropna().astype(str).unique().tolist()
    if candidates:
      return candidates[:max_items]
  if "interaction" in interactions.columns:
    candidates = interactions["interaction"].dropna().astype(str).unique().tolist()
    if candidates:
      return candidates[:max_items]
  if hasattr(pvals, "index"):
    candidates = [str(v) for v in pvals.index.tolist() if str(v)]
    if candidates:
      return candidates[:max_items]
  return []

def resolve_gene_family(user_gene_family):
  allowed = {"chemokines", "th1", "th2", "th17", "treg", "costimulatory", "coinhibitory"}
  if user_gene_family and str(user_gene_family) in allowed:
    return str(user_gene_family)
  return "chemokines"

celltype_key = resolve_celltype_key(adata, args.celltype)
cell_type1, cell_type2 = resolve_celltypes(adata, celltype_key, args.cell_type1, args.cell_type2)
interaction_pairs = resolve_interaction_pairs(interactions, pvals, args.interaction_pairs)
selected_gene_family = resolve_gene_family(args.gene_family)



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
    cell_type1 = cell_type1,
    cell_type2 = cell_type2, 
    means = means,
    pvals = pvals,
    celltype_key = celltype_key,
    figsize = (10,10),
    max_size = 5,
    highlight_size = 0.75,
    degs_analysis = False,
    standard_scale = False,
    interaction_scores = interactions if interactions_path else None,
    scale_alpha_by_interaction_scores=bool(interactions_path))
p3.save(os.path.join(output_dir,f"{args.sample_id}_dot_plot.png"),dpi=300,limitsize=False)

# chemokines', 'th1', 'th2', 'th17', 'treg', 'costimulatory', 'coinhibitory
#selected_gene_family = 'costimulatory'
print(selected_gene_family)
print("@@@@")
p4=kpy.plot_cpdb(
    adata=adata,
    cell_type1=cell_type1,
    cell_type2=cell_type2,
    means=means,
    pvals=pvals,
    celltype_key=celltype_key,
    gene_family=selected_gene_family,
    highlight_size=0.9,
    figsize=(12, 6),
)
p4.save(os.path.join(output_dir,f"{args.sample_id}_dot_family_plot.png"),dpi=300,limitsize=False)

if interaction_pairs and decon is not None:
  chord_fig = kpy.plot_cpdb_chord(
      adata=adata,
      cell_type1=cell_type1,
      cell_type2=cell_type2,
      means=means,
      pvals=pvals,
      deconvoluted=decon,
      celltype_key=celltype_key,
      interaction=interaction_pairs,
      link_kwargs={"direction": 1, "allow_twist": True, "r1": 95, "r2": 90},
      sector_text_kwargs={"color": "black", "size": 12, "r": 105, "adjust_rotation": True},
      legend_kwargs={"loc": "center", "bbox_to_anchor": (1, 1), "fontsize": 8},
      link_offset=1,
  )
  chord_fig.savefig(os.path.join(output_dir,f"{args.sample_id}_chord_plot.png"),dpi=300)






