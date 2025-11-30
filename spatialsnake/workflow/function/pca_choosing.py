import numpy as np
import scanpy as sc
import numpy as np

def select_pca_dimensions(adata):
    """
    revise from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
    the best n_pcs to choose
    """
    if 'pca' not in adata.uns or 'variance' not in adata.uns['pca']:
        raise ValueError("please run the [sc.tl.pca(adata)]")
    variance = adata.uns['pca']['variance']
    stdev = np.sqrt(variance)
    pct = stdev / stdev.sum() * 100
    cum = np.cumsum(pct)
    criteria_co1 = (cum > 90) & (pct < 5)
    co1_indices = np.where(criteria_co1)[0]
    
    if len(co1_indices) > 0:
        co1 = co1_indices[0] + 1
    else:
        co1 = len(pct)
    diffs = pct[:-1] - pct[1:]
    co2_indices = np.where(diffs > 0.1)[0]
    
    if len(co2_indices) > 0:
        last_idx = co2_indices.max()
        co2 = last_idx + 2
    else:
        co2 = len(pct)
    n_pcs = min(co1, co2)
    print(f"Cutoff 1 (>90% & PC<5%): PC {co1}")
    print(f"Cutoff 2 (elbow): PC {co2}")
    print(f"recommand: {n_pcs}")
    return n_pcs
