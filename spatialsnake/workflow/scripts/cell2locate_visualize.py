import copy
import os
import argparse
import numpy as np
import pandas as pd
import spatialdata as spd
from spatialsnake.workflow.function.transform import zarr_to_h5ad
from spatialsnake.workflow.function.plot import plot_cell2location_dotplot
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

logger = setup_logger("cell2location_visualize")

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", required=True)
parser.add_argument("--sample_id", required=True)
parser.add_argument("--output_zarr_path", required=True)
parser.add_argument("--type", required=True)
parser.add_argument("--image_type", required=True)
parser.add_argument("--image_slice", required=False)
parser.add_argument("--sample_cnt", type=int, required=True)
parser.add_argument("--coord", type=int, nargs=4, required=False)
parser.add_argument("--input_st", required=False)
parser.add_argument("--shape_type", required=False)
parser.add_argument("--celltype_col", required=False, default="celltype")
parser.add_argument("--batch_key_st", required=False, default="sample")
parser.add_argument("--microenvironment_output", required=True)
parser.add_argument("--microenvironment_threshold", type=float, default=0.10)
parser.add_argument("--dotplot_max_cell_types", type=int, default=30)
parser.add_argument("--dotplot_enrichment_clip", type=float, default=2.5)
args = parser.parse_args()
if not np.isfinite(args.microenvironment_threshold) or not 0 < args.microenvironment_threshold <= 1:
    parser.error("--microenvironment_threshold must be in the interval (0, 1]")
if args.dotplot_max_cell_types < 0:
    parser.error("--dotplot_max_cell_types must be >= 0; use 0 to display all cell types")
if not np.isfinite(args.dotplot_enrichment_clip) or args.dotplot_enrichment_clip <= 0:
    parser.error("--dotplot_enrichment_clip must be a positive finite number")


def prepare_scanpy_adata(zarr_path, table_adata, output_dir):
    converted = zarr_to_h5ad(
        [zarr_path],
        os.path.join(output_dir, "convert_data.h5ad"),
        save_image=True
    )
    if converted is None:
        return None
    converted.obs_names = converted.obs_names.astype(str)
    table_copy = table_adata.copy()
    table_copy.obs_names = table_copy.obs_names.astype(str)
    if not converted.obs_names.is_unique or not table_copy.obs_names.is_unique:
        raise ValueError("Cannot prepare co-location input with duplicated observation IDs")
    if set(converted.obs_names) != set(table_copy.obs_names):
        missing = table_copy.obs_names[~table_copy.obs_names.isin(converted.obs_names)].tolist()[:10]
        extra = converted.obs_names[~converted.obs_names.isin(table_copy.obs_names)].tolist()[:10]
        raise ValueError(
            "Converted AnnData observations do not match the cell2location table. "
            f"Missing examples: {missing}; extra examples: {extra}"
        )
    table_copy = table_copy[converted.obs_names].copy()
    for column in table_copy.obs.columns:
        converted.obs[column] = table_copy.obs[column].to_numpy()
    for key in ("q05_cell_abundance_w_sf", "means_cell_abundance_w_sf", "q95_cell_abundance_w_sf", "stds_cell_abundance_w_sf"):
        if key in table_copy.obsm:
            converted.obsm[key] = table_copy.obsm[key].copy()
    if "mod" in table_copy.uns:
        converted.uns["mod"] = copy.deepcopy(table_copy.uns["mod"])
    if args.batch_key_st not in converted.obs:
        raise ValueError(
            f"Co-location sample column '{args.batch_key_st}' was not found in the converted AnnData"
        )
    batch_values = converted.obs[args.batch_key_st]
    if batch_values.isna().any() or batch_values.astype(str).str.strip().eq("").any():
        raise ValueError(
            f"Co-location sample column '{args.batch_key_st}' contains missing or blank values"
        )
    return converted


def transfer_colocation_results(source, target):
    source_names = pd.Index(source.obs_names.astype(str))
    target_names = pd.Index(target.obs_names.astype(str))
    if not source_names.is_unique or not target_names.is_unique:
        raise ValueError("Cannot transfer co-location results with duplicated observation IDs")
    if set(source_names) != set(target_names):
        missing = target_names[~target_names.isin(source_names)].tolist()[:10]
        extra = source_names[~source_names.isin(target_names)].tolist()[:10]
        raise ValueError(
            "Co-location result observations do not match the original spatial table. "
            f"Missing examples: {missing}; extra examples: {extra}"
        )

    source_obs = source.obs.copy()
    source_obs.index = source_names
    n_fact12_slot = source.uns.get("mod_coloc_n_fact12")
    if not isinstance(n_fact12_slot, dict) or "fact_names" not in n_fact12_slot:
        raise KeyError("Expected n_fact=12 factor names in uns['mod_coloc_n_fact12']")
    factor_names = [str(value) for value in n_fact12_slot["fact_names"]]
    if len(factor_names) != 12 or len(set(factor_names)) != len(factor_names):
        raise ValueError("uns['mod_coloc_n_fact12']['fact_names'] must contain 12 unique names")
    factor_columns = [f"mean_nUMI_factors{factor_name}" for factor_name in factor_names]
    missing_columns = [column for column in factor_columns if column not in source_obs.columns]
    if missing_columns:
        raise KeyError(
            "Expected continuous n_fact=12 co-location scores missing from adata_coloc.obs: "
            f"{missing_columns}"
        )
    stale_columns = [
        column
        for column in target.obs.columns
        if str(column).startswith("mean_nUMI_factorsfact_") and column not in factor_columns
    ]
    if stale_columns:
        target.obs.drop(columns=stale_columns, inplace=True)
        logger.warning("Removed %d stale co-location factor score column(s)", len(stale_columns))
    for column in factor_columns:
        target.obs[column] = source_obs.loc[target_names, column].to_numpy()

    for key in ("mod_coloc_n_fact11", "mod_coloc_n_fact12"):
        if key not in source.uns:
            raise KeyError(f"Expected co-location model output missing from uns['{key}']")
        target.uns[key] = copy.deepcopy(source.uns[key])
    logger.info(f"Transferred {len(factor_columns)} continuous n_fact=12 co-location score columns")
    return target


def get_cell_type_fractions(res_dict, adata_coloc, n_fact=12):
    model_result = res_dict.get(f"n_fact{n_fact}", {})
    model = model_result.get("mod") if isinstance(model_result, dict) else None
    fractions = getattr(model, "cell_type_fractions", None)
    if isinstance(fractions, pd.DataFrame):
        return fractions.copy()

    slot_name = f"mod_coloc_n_fact{n_fact}"
    if slot_name not in adata_coloc.uns:
        raise KeyError(f"Expected co-location model output missing from uns['{slot_name}']")
    slot = adata_coloc.uns[slot_name]
    try:
        loadings = np.asarray(slot["post_sample_means"]["cell_type_factors"], dtype=float)
        cell_types = pd.Index(slot["var_names"].astype(str) if hasattr(slot["var_names"], "astype") else slot["var_names"])
        factor_names = pd.Index(slot["fact_names"].astype(str) if hasattr(slot["fact_names"], "astype") else slot["fact_names"])
    except (KeyError, TypeError) as error:
        raise KeyError(
            f"Cannot reconstruct cell_type_fractions from uns['{slot_name}']"
        ) from error
    if loadings.shape != (len(cell_types), len(factor_names)):
        raise ValueError(
            f"Invalid cell_type_factors shape {loadings.shape}; expected "
            f"({len(cell_types)}, {len(factor_names)})"
        )
    row_sums = loadings.sum(axis=1)
    if np.any(row_sums <= 0):
        raise ValueError("Cannot calculate cell_type_fractions from non-positive cell-type loadings")
    return pd.DataFrame(loadings / row_sums[:, None], index=cell_types, columns=factor_names)


def write_microenvironments(res_dict, adata_coloc, output_path, threshold, n_fact=12):
    if not np.isfinite(threshold) or threshold <= 0 or threshold > 1:
        raise ValueError("cell2location_microenvironment_threshold must be in the interval (0, 1]")

    fractions = get_cell_type_fractions(res_dict, adata_coloc, n_fact=n_fact)
    fractions.index = fractions.index.astype(str)
    fractions.columns = fractions.columns.astype(str)
    if not fractions.index.is_unique or not fractions.columns.is_unique:
        raise ValueError("cell_type_fractions contains duplicated cell-type or factor names")
    values = fractions.to_numpy(dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("cell_type_fractions must contain finite, non-negative values")

    relationships = []
    for factor_name in fractions.columns:
        members = fractions.index[fractions[factor_name].to_numpy() >= threshold]
        if len(members) == 0:
            logger.warning(
                f"Co-location factor '{factor_name}' has no cell types at loading threshold {threshold:g}"
            )
            continue
        relationships.extend(
            {"cell_type": cell_type, "microenvironment": factor_name}
            for cell_type in members
        )

    output_table = pd.DataFrame(
        relationships,
        columns=["cell_type", "microenvironment"],
    ).drop_duplicates()
    if output_table.empty:
        raise ValueError(
            f"No cell types passed the n_fact={n_fact} microenvironment threshold {threshold:g}"
        )
    assigned_cell_types = set(output_table["cell_type"])
    omitted_cell_types = [cell_type for cell_type in fractions.index if cell_type not in assigned_cell_types]
    if omitted_cell_types:
        logger.warning(
            "Microenvironment threshold omitted %d of %d reference cell types; omitted types will not "
            "participate in CellPhoneDB pairs when this file is used. Examples: %s",
            len(omitted_cell_types),
            len(fractions.index),
            omitted_cell_types[:10],
        )
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    output_table.to_csv(output_path, sep="\t", index=False)
    logger.info(
        "Saved %d cell type-factor relationship(s) covering %d/%d reference cell types to %s",
        len(output_table),
        len(assigned_cell_types),
        len(fractions.index),
        output_path,
    )
    return output_table



output_dir = os.path.dirname(args.output_zarr_path)
figure_dir = os.path.join(output_dir, "figure")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(figure_dir, exist_ok=True)

log_step(logger, 1, 6, "loading cell2location SpatialData output")
concatenated_sdata = spd.read_zarr(args.input_dir)
table_key = next(iter(concatenated_sdata.tables.keys()))
adata_vis = concatenated_sdata[table_key]
log_step(logger, 2, 6, "saving relative regional cell-abundance dotplot")
if (
    args.celltype_col in adata_vis.obs.columns
    and "q05_cell_abundance_w_sf" in adata_vis.obsm
):
    plot_cell2location_dotplot(
        adata_vis,
        unsupervised_cluster_col=args.celltype_col,
        sample_col=args.batch_key_st,
        abundance_key="q05_cell_abundance_w_sf",
        max_cell_types=args.dotplot_max_cell_types,
        enrichment_clip=args.dotplot_enrichment_clip,
        save_path=os.path.join(figure_dir, "cell2location_relative_abundance_dotplot.pdf"),
        source_data_path=os.path.join(
            figure_dir,
            "cell2location_relative_abundance_dotplot_source.tsv",
        ),
    )
else:
    logger.warning(
        f"Skip cell2location regional-abundance dotplot because '{args.celltype_col}' was not found in obs, "
        "or 'q05_cell_abundance_w_sf' was missing in obsm."
    )

log_step(logger, 3, 6, "preparing Scanpy-compatible AnnData for co-location")
adata_for_colocation = prepare_scanpy_adata(args.input_dir, adata_vis, output_dir)
if adata_for_colocation is None:
    raise ValueError("Failed to prepare a Scanpy-compatible AnnData object for cell2location co-location.")

from cell2location import run_colocation

log_step(logger, 4, 6, "running cell2location co-location analysis")
res_dict, adata_coloc = run_colocation(
    adata_for_colocation,
    model_name="CoLocatedGroupsSklearnNMF",
    train_args={
        "n_fact": np.arange(11, 13),
        "sample_name_col": args.batch_key_st,
        "n_restarts": 16
    },
    model_kwargs={"alpha": 0.01, "init": "random", "nmf_kwd_args": {"tol": 0.000001}},
    export_args={
        "path": f"{output_dir}/CoLocatedComb/",
        "plot_histology": True,
        "scanpy_alpha_img": 1.0,
    }
)

log_step(logger, 5, 6, "saving co-location results and CellPhoneDB microenvironments")
microenvironment_table = write_microenvironments(
    res_dict,
    adata_coloc,
    args.microenvironment_output,
    args.microenvironment_threshold,
    n_fact=12,
)
if args.celltype_col in adata_vis.obs.columns:
    metadata_labels = set(
        adata_vis.obs[args.celltype_col].dropna().astype(str).str.strip()
    )
    incompatible_labels = sorted(set(microenvironment_table["cell_type"]) - metadata_labels)
    if incompatible_labels:
        logger.warning(
            "The microenvironment TSV uses reference cell-type labels, but %d label(s) are absent "
            "from obs['%s']. Do not pass this TSV to CellPhoneDB unless its metadata uses the same "
            "reference labels. Examples: %s",
            len(incompatible_labels),
            args.celltype_col,
            incompatible_labels[:10],
        )
adata_vis = transfer_colocation_results(adata_coloc, adata_vis)
concatenated_sdata[table_key] = adata_vis
concatenated_sdata.write(args.output_zarr_path)

log_step(logger, 6, 6, "cell2location visualization outputs saved")
logger.info("Cell2location visualization module completed")
