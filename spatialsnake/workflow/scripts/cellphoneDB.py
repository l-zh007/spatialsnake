import argparse
import os
import urllib.request

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import pandas as pd
import scanpy as sc
import spatialdata as spd
from cellphonedb.src.core.methods import (
    cpdb_degs_analysis_method,
    cpdb_statistical_analysis_method,
)
from cellphonedb.utils import db_utils

from spatialsnake.workflow.function.logging_utils import log_step, setup_logger


logger = setup_logger("cellphonedb")
parser = argparse.ArgumentParser(description="Run CellPhoneDB on an annotated h5ad or SpatialData zarr object")
parser.add_argument("--input_dir", required=True, help="Annotated h5ad or SpatialData zarr input")
parser.add_argument("--output_zarr_path", required=True, help="Intermediate h5ad used by CellPhoneDB")
parser.add_argument("--sample_id", required=True, help="Sample identifier")
parser.add_argument("--type", required=False, default="", help="Input technology (kept for workflow compatibility)")
parser.add_argument("--counts_data", default="hgnc_symbol", help="CellPhoneDB gene identifier type")
parser.add_argument("--iterations", type=int, default=1000, help="Permutation iterations")
parser.add_argument("--threads", type=int, default=8, help="CellPhoneDB worker threads")
parser.add_argument("--pvalue", type=float, default=0.05, help="Statistical significance threshold")
parser.add_argument("--output_name", default="", help="Suffix for CellPhoneDB result tables")
parser.add_argument("--microenvs_file_path", default="", help="Optional CellPhoneDB microenvironment file")
parser.add_argument("--active_tf_path", default="", help="Optional active transcription-factor file")
parser.add_argument("--threshold", type=float, default=0.1, help="Minimum expressing-cell proportion")
parser.add_argument("--degs_file_path", default="", help="Optional DEG file for CellPhoneDB DEG mode")
parser.add_argument("--celltype_col", required=True, help="Cell-type annotation column in obs")
parser.add_argument("--niche_col", default="spatial_cluster", help="Spatial niche column in obs")
parser.add_argument("--is_single_cell", default="False", help="Whether the input is single-cell rather than spatial")
parser.add_argument("--method", default="statistical", choices=["statistical", "degs"], help="CellPhoneDB method")
parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Input species")
args = parser.parse_args()


def parse_bool(value):
    if isinstance(value, bool):
        return value
    value = "" if value is None else str(value).strip().lower()
    if value in {"true", "1", "yes", "y", "t"}:
        return True
    if value in {"false", "0", "no", "n", "f", ""}:
        return False
    raise ValueError(f"Invalid boolean value: {value}")


def ensure_cell_id_column(adata):
    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = adata.obs_names.astype(str)
    else:
        adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    if adata.obs["cell_id"].str.strip().eq("").any() or adata.obs["cell_id"].duplicated().any():
        raise ValueError("cell_id must contain unique, non-empty values")
    return adata


def load_mouse_human_orthologs(cache_dir):
    """Download and parse the MGI one-to-one mouse-human protein-coding map."""
    os.makedirs(cache_dir, exist_ok=True)
    mapping_path = os.path.join(cache_dir, "HOM_ProteinCoding.rpt")
    if not os.path.isfile(mapping_path):
        logger.info("Downloading the MGI one-to-one mouse-human ortholog table")
        temporary_path = f"{mapping_path}.tmp"
        urllib.request.urlretrieve(
            "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt",
            temporary_path,
        )
        os.replace(temporary_path, mapping_path)

    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    mouse_col = next((column for column in mapping.columns if "mouse" in column.lower() and "symbol" in column.lower()), None)
    human_col = next((column for column in mapping.columns if "human" in column.lower() and "symbol" in column.lower()), None)
    if mouse_col is None or human_col is None:
        columns = [
            "MGI Marker Accession ID",
            "Mouse Gene Symbol",
            "Mouse NCBI Gene ID",
            "HGNC ID",
            "Human Gene Symbol",
            "Human NCBI Gene ID",
        ]
        mapping = pd.read_csv(mapping_path, sep="\t", dtype=str, header=None, names=columns)
        mouse_col, human_col = "Mouse Gene Symbol", "Human Gene Symbol"

    mapping = mapping[[mouse_col, human_col]].dropna().drop_duplicates()
    mapping.columns = ["mouse_symbol", "human_symbol"]
    mapping = mapping.drop_duplicates("mouse_symbol", keep=False).drop_duplicates("human_symbol", keep=False)
    return dict(zip(mapping["mouse_symbol"], mapping["human_symbol"]))


def convert_mouse_anndata(adata, orthologs):
    original_symbols = pd.Index(adata.var_names.astype(str))
    human_symbols = original_symbols.map(orthologs)
    keep = human_symbols.notna()
    converted = adata[:, keep].copy()
    converted.var["mouse_gene_symbol"] = original_symbols[keep].to_numpy()
    converted.var_names = pd.Index(human_symbols[keep].astype(str), name=adata.var_names.name)
    if converted.n_vars == 0:
        raise ValueError("No mouse gene symbols could be mapped to one-to-one human orthologs")
    logger.info("Mouse-to-human ortholog conversion retained %d of %d genes", converted.n_vars, adata.n_vars)
    return converted


def convert_mouse_gene_file(input_path, orthologs, output_path, label):
    if not input_path or not os.path.isfile(input_path):
        return None
    table = pd.read_csv(input_path, sep=None, engine="python", dtype=str)
    if table.shape[1] < 2:
        raise ValueError(f"{label} file must contain at least two columns")
    gene_column = table.columns[1]
    table[gene_column] = table[gene_column].map(orthologs)
    table = table.dropna(subset=[gene_column])
    table.to_csv(output_path, sep="\t", index=False)
    logger.info("Converted %d %s records to human orthologs", len(table), label)
    return output_path


output_dir = os.path.dirname(args.output_zarr_path)
os.makedirs(output_dir, exist_ok=True)
output_name = args.output_name.strip() or args.sample_id
is_single_cell = parse_bool(args.is_single_cell)
threads = max(1, args.threads)

log_step(logger, 1, 5, "preparing CellPhoneDB database")
cpdb_version = "v5.0.0"
cpdb_target_dir = os.path.join(output_dir, "cellphonedb_v500_NatProtocol", cpdb_version)
cpdb_file_path = os.path.join(cpdb_target_dir, "cellphonedb.zip")
os.makedirs(cpdb_target_dir, exist_ok=True)
if not os.path.isfile(cpdb_file_path):
    db_utils.download_database(cpdb_target_dir, cpdb_version)

log_step(logger, 2, 5, "loading expression data")
input_ext = os.path.splitext(args.input_dir)[1].lower()
if input_ext in {".h5ad", ".h5"}:
    adata = sc.read_h5ad(args.input_dir)
else:
    spatial_object = spd.read_zarr(args.input_dir)
    table_names = list(spatial_object.tables)
    if not table_names:
        raise ValueError("No table was found in the SpatialData object")
    adata = spatial_object.tables[table_names[0]]

adata = ensure_cell_id_column(adata)
if args.celltype_col not in adata.obs.columns:
    raise ValueError(f"celltype_col not found in obs: {args.celltype_col}")

orthologs = None
counts_data = args.counts_data
degs_file_path = args.degs_file_path if args.degs_file_path and os.path.isfile(args.degs_file_path) else None
active_tf_path = args.active_tf_path if args.active_tf_path and os.path.isfile(args.active_tf_path) else None
if args.species == "mouse":
    orthologs = load_mouse_human_orthologs(cpdb_target_dir)
    adata = convert_mouse_anndata(adata, orthologs)
    counts_data = "hgnc_symbol"
    degs_file_path = convert_mouse_gene_file(
        degs_file_path,
        orthologs,
        os.path.join(output_dir, f"{args.sample_id}_degs_human_orthologs.txt"),
        "DEG",
    )
    active_tf_path = convert_mouse_gene_file(
        active_tf_path,
        orthologs,
        os.path.join(output_dir, f"{args.sample_id}_active_tfs_human_orthologs.txt"),
        "active TF",
    )
    logger.info("Running CellPhoneDB with human orthologs projected from mouse genes")

microenvs_file_path = None
generated_microenvs_path = os.path.join(output_dir, f"{args.sample_id}_microenvs.txt")
if not is_single_cell:
    if args.niche_col in adata.obs.columns:
        microenvs = adata.obs[[args.celltype_col, args.niche_col]].dropna().drop_duplicates().copy()
        microenvs.columns = ["cell_type", "microenvironment"]
        microenvs = microenvs.astype(str)
        microenvs = microenvs[
            microenvs["cell_type"].str.strip().ne("")
            & microenvs["microenvironment"].str.strip().ne("")
        ]
        microenvs.to_csv(generated_microenvs_path, sep="\t", index=False)
        microenvs_file_path = generated_microenvs_path
        logger.info("Using obs[%s] as the spatial microenvironment", args.niche_col)
    elif args.microenvs_file_path and os.path.isfile(args.microenvs_file_path):
        microenvs_file_path = os.path.abspath(args.microenvs_file_path)
        logger.info("Using the user-provided microenvironment file")
    else:
        logger.info(
            "obs[%s] was not found and no microenvironment file was provided; running without spatial restriction",
            args.niche_col,
        )
else:
    logger.info("Single-cell mode selected; running without automatic spatial niche restriction")

adata.uns["cellphonedb"] = {
    "method": args.method,
    "output_name": output_name,
    "species": args.species,
    "celltype_col": args.celltype_col,
    "niche_col": args.niche_col,
    "microenvs_file_path": microenvs_file_path or "",
    "pvalue": float(args.pvalue),
}

log_step(logger, 3, 5, "writing CellPhoneDB input files")
adata.write(args.output_zarr_path)
metadata = adata.obs[["cell_id", args.celltype_col]].copy()
metadata.columns = ["cell_id", "cell_type"]
metadata_path = os.path.join(output_dir, f"{args.sample_id}_cellid_cell_type.txt")
metadata.to_csv(metadata_path, sep="\t", index=False)

out_path = os.path.join(output_dir, "cellphonedb_output")
os.makedirs(out_path, exist_ok=True)
common_kwargs = {
    "cpdb_file_path": cpdb_file_path,
    "meta_file_path": metadata_path,
    "counts_file_path": args.output_zarr_path,
    "counts_data": counts_data,
    "threshold": args.threshold,
    "threads": threads,
    "score_interactions": True,
    "output_path": out_path,
    "output_suffix": output_name,
}
if active_tf_path:
    common_kwargs["active_tfs_file_path"] = active_tf_path
if microenvs_file_path:
    common_kwargs["microenvs_file_path"] = microenvs_file_path

log_step(logger, 4, 5, f"running CellPhoneDB {args.method} analysis")
if args.method == "degs":
    if not degs_file_path:
        raise ValueError("degs_file_path is required when method is degs")
    cpdb_results = cpdb_degs_analysis_method.call(degs_file_path=degs_file_path, **common_kwargs)
else:
    cpdb_results = cpdb_statistical_analysis_method.call(
        iterations=args.iterations,
        pvalue=args.pvalue,
        **common_kwargs,
    )

log_step(logger, 5, 5, f"CellPhoneDB outputs saved to {out_path}")
logger.info("CellPhoneDB module completed")
