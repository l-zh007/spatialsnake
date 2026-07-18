"""Replicate-aware pseudobulk differential expression for compare_analysis.

One invocation analyses exactly one cell-type/region label and one directed
comparison.  Snakemake expands these invocations so that targets can run in
parallel.  Plotting and enrichment are delegated to R; this module never opens
a graphics device and does not retain pseudobulk helper files.
"""

from __future__ import annotations

import argparse
import inspect
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import spatialdata as spd
from scipy import sparse

from spatialsnake.workflow.function.get_sample import read_compare_sample_table
from spatialsnake.workflow.function.logging_utils import log_step, setup_logger


logger = setup_logger("compare_gene")

FINAL_COLUMNS = [
    "gene",
    "gene_symbol",
    "base_mean",
    "log2FoldChange",
    "statistic",
    "pvalue",
    "padj",
    "comparison",
    "reference",
    "celltype",
    "algorithm",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one replicate-aware pseudobulk cell-type/region contrast."
    )
    parser.add_argument("--input_dir", required=True, help="Annotated merged SpatialData Zarr")
    parser.add_argument(
        "--sample_list",
        required=True,
        help="Three-column compare sample table (space or tab separated)",
    )
    parser.add_argument("--output_dir", required=True, help="Final cell-type/contrast directory")
    parser.add_argument("--type", required=True, help="Spatial platform name, used for result semantics")
    parser.add_argument("--algorithm", required=True, help="DESeq2 or edgeR")
    parser.add_argument("--celltype", required=True, help="Exact target label")
    parser.add_argument("--comparison", required=True, help="Numerator condition")
    parser.add_argument("--reference", required=True, help="Reference condition")
    parser.add_argument("--celltype_col", default="celltype")
    parser.add_argument("--sample_col", default="sample")
    parser.add_argument("--condition_col", default="group")
    parser.add_argument("--count_layer", default="counts")
    parser.add_argument("--min_replicates", type=int, default=2)
    parser.add_argument("--min_cells_per_sample", type=int, default=30)
    parser.add_argument("--min_total_counts_per_gene", type=int, default=10)
    parser.add_argument("--cut_off_pvalue", type=float, default=0.05)
    parser.add_argument("--cut_off_logFC", type=float, default=0.5)
    parser.add_argument("--de_top_n", type=int, default=20)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--gene_symbol_col",
        default="",
        nargs="?",
        const="",
        help="Optional var column containing gene symbols",
    )
    parser.add_argument("--gene_id_type", default="auto")
    parser.add_argument("--species", default="human")
    parser.add_argument("--go_ontology", default="BP")
    parser.add_argument("--enrichment_top_n", type=int, default=10)
    parser.add_argument("--differential_r", required=True)
    parser.add_argument("--enrichment_r", required=True)
    parser.add_argument("--rscript", default="Rscript")
    return parser.parse_args()


def normalize_algorithm(value: str) -> str:
    text = str(value).strip().lower()
    if text in {"deseq2", "deseq"}:
        return "DESeq2"
    if text == "edger":
        return "edgeR"
    raise ValueError("compare_algorithm must be DESeq2 or edgeR")


def read_sample_table(path: str) -> pd.DataFrame:
    try:
        rows = read_compare_sample_table(path)
    except (OSError, RuntimeError, ValueError) as exc:
        raise ValueError(str(exc)) from exc
    table = pd.DataFrame(rows, columns=["sample_id", "input_path", "group"])
    if table["group"].nunique() < 2:
        raise ValueError("compare_analysis requires at least two biological groups")
    return table


def get_primary_table(sdata):
    if "table" in sdata.tables:
        return sdata.tables["table"]
    table_names = list(sdata.tables.keys())
    if not table_names:
        raise ValueError("Input SpatialData contains no table")
    return sdata.tables[table_names[0]]


def validate_count_matrix(matrix, row_indices: np.ndarray | None = None) -> None:
    if len(matrix.shape) != 2:
        raise ValueError("Raw counts must be a two-dimensional matrix")
    selected = matrix if row_indices is None else matrix[row_indices, :]
    if sparse.issparse(selected):
        values = np.asarray(selected.data)
        if values.size and (
            not np.isfinite(values).all()
            or (values < 0).any()
            or not np.allclose(values, np.rint(values), atol=1e-6)
        ):
            raise ValueError("Raw counts contain non-finite, negative, or non-integer values")
        return
    for start in range(0, selected.shape[0], 1024):
        values = np.asarray(selected[start : start + 1024, :])
        if values.size and (
            not np.isfinite(values).all()
            or (values < 0).any()
            or not np.allclose(values, np.rint(values), atol=1e-6)
        ):
            raise ValueError("Raw counts contain non-finite, negative, or non-integer values")


def select_count_matrix(
    adata,
    requested_layer: str,
    row_indices: np.ndarray | None = None,
):
    requested = str(requested_layer or "").strip()
    candidates: list[tuple[str, Any, Any, Any]] = []
    normalized_request = requested.lower()
    standard_requests = {"", "auto", "counts", "raw_counts", "none", "null"}
    strict_request = normalized_request not in standard_requests
    if requested and normalized_request not in {"auto", "none", "null"}:
        if requested in adata.layers:
            candidates.append(
                (f"layers[{requested!r}]", adata.layers[requested], adata.var_names, adata.var)
            )
        elif strict_request:
            raise ValueError(f"Requested count_layer {requested!r} is absent from the SpatialData table")
        else:
            logger.warning(
                "Preferred count layer %r is absent; checking standard raw-count sources",
                requested,
            )

    existing_sources = {source for source, _, _, _ in candidates}
    for layer_name in ("counts", "raw_counts"):
        source = f"layers[{layer_name!r}]"
        if layer_name in adata.layers and source not in existing_sources:
            candidates.append((source, adata.layers[layer_name], adata.var_names, adata.var))
    if adata.raw is not None and adata.raw.n_obs == adata.n_obs:
        candidates.append(("raw.X", adata.raw.X, adata.raw.var_names, adata.raw.var))
    candidates.append(("X", adata.X, adata.var_names, adata.var))

    errors = []
    for source, matrix, var_names, var in candidates:
        try:
            validate_count_matrix(matrix, row_indices=row_indices)
            names = pd.Index(var_names.astype(str))
            if names.has_duplicates:
                duplicated = names[names.duplicated()].unique().tolist()[:10]
                raise ValueError(f"duplicate gene IDs: {', '.join(duplicated)}")
            return matrix, names.to_numpy(), var.copy(), source
        except ValueError as exc:
            errors.append(f"{source}: {exc}")
            if strict_request:
                break
    detail = "; ".join(errors) if errors else "no candidate count matrix was found"
    raise ValueError(f"compare_analysis requires non-negative integer raw counts ({detail})")


def validate_obs_metadata(obs: pd.DataFrame, sample_table: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    for column in (args.sample_col, args.celltype_col):
        if column not in obs.columns:
            raise ValueError(f"Required SpatialData obs column is missing: {column}")
    result = obs.copy()
    result[args.sample_col] = result[args.sample_col].astype(str)
    result[args.celltype_col] = result[args.celltype_col].astype(str)
    sample_ids = set(sample_table["sample_id"])
    observed_ids = set(result[args.sample_col])
    missing = sorted(sample_ids - observed_ids)
    extra = sorted(observed_ids - sample_ids)
    if missing:
        raise ValueError(f"Samples in sample.txt are absent from obs[{args.sample_col!r}]: {', '.join(missing)}")
    if extra:
        raise ValueError(
            f"obs[{args.sample_col!r}] contains samples absent from sample.txt: {', '.join(extra[:20])}"
        )

    if args.condition_col in result.columns:
        condition = result[[args.sample_col, args.condition_col]].copy()
        condition[args.condition_col] = condition[args.condition_col].astype(str).str.strip()
        condition = condition[~condition[args.condition_col].isin({"", "nan", "None", "NA"})]
        if not condition.empty:
            counts = condition.groupby(args.sample_col)[args.condition_col].nunique()
            ambiguous = counts[counts > 1].index.tolist()
            obs_map = condition.drop_duplicates(args.sample_col).set_index(args.sample_col)[args.condition_col]
            expected_map = sample_table.set_index("sample_id")["group"]
            conflicts = [
                sample
                for sample in expected_map.index
                if sample in obs_map.index and str(obs_map[sample]) != str(expected_map[sample])
            ]
            if ambiguous or conflicts:
                details = []
                if ambiguous:
                    details.append("multiple obs labels: " + ", ".join(ambiguous[:10]))
                if conflicts:
                    details.append("different obs/sample.txt labels: " + ", ".join(conflicts[:20]))
                logger.warning(
                    "Ignoring obs[%r] group metadata (%s); sample.txt is the authoritative group source",
                    args.condition_col,
                    "; ".join(details),
                )
    return result


def aggregate_by_sample(matrix, row_indices: np.ndarray, samples: np.ndarray, sample_order: list[str], genes) -> pd.DataFrame:
    subset = matrix[row_indices, :]
    if sparse.issparse(subset):
        subset = subset.tocsr().astype(np.int64)
    else:
        subset = sparse.csr_matrix(np.asarray(subset, dtype=np.int64))
    codes = pd.Categorical(samples, categories=sample_order).codes
    if (codes < 0).any():
        raise ValueError("Target observations contain sample IDs absent from sample.txt")
    selector = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.int64), (codes, np.arange(len(codes)))),
        shape=(len(sample_order), len(codes)),
    )
    aggregated = selector @ subset
    dense = aggregated.toarray()
    if not np.isfinite(dense).all() or (dense < 0).any() or not np.equal(dense, np.rint(dense)).all():
        raise ValueError("Pseudobulk aggregation produced invalid counts")
    return pd.DataFrame(dense.astype(np.int64), index=sample_order, columns=genes)


def build_sample_metadata(sample_table: pd.DataFrame) -> pd.DataFrame:
    """Use each sample_id as one independent biological replicate."""
    metadata = sample_table[["sample_id", "group"]].rename(columns={"group": "condition"})
    metadata = metadata.rename_axis(None).set_index("sample_id", drop=False)
    metadata["replicate"] = metadata["sample_id"]
    return metadata


def resolve_gene_symbols(
    var_names: np.ndarray,
    var: pd.DataFrame,
    column: str,
    gene_id_type: str,
) -> np.ndarray:
    requested = str(column or "").strip()
    if requested:
        if requested not in var.columns:
            raise ValueError(f"compare_gene_symbol_col {requested!r} is absent from SpatialData var")
        symbols = var[requested].astype(str).str.strip().to_numpy()
        missing = np.isin(symbols, ["", "nan", "None", "NA"])
        if missing.any():
            raise ValueError(
                f"compare_gene_symbol_col {requested!r} contains {int(missing.sum())} empty values"
            )
        return symbols

    genes = pd.Index(var_names.astype(str))
    identifier_type = str(gene_id_type or "auto").strip().upper()
    if identifier_type == "SYMBOL":
        return genes.to_numpy()
    if identifier_type in {"ENSEMBL", "ENTREZID"}:
        return np.full(len(genes), "", dtype=object)

    # In auto mode, retain var_names as symbols only when they do not look like
    # a predominantly ENSEMBL or ENTREZ identifier set.  This keeps the public
    # gene_symbol column useful for common SYMBOL-indexed objects while allowing
    # the enrichment script to detect and map non-SYMBOL gene IDs correctly.
    ensembl_fraction = genes.str.match(r"^ENS[A-Z0-9]*G[0-9]+(?:\.[0-9]+)?$").mean()
    entrez_fraction = genes.str.match(r"^[0-9]+$").mean()
    if max(float(ensembl_fraction), float(entrez_fraction)) >= 0.8:
        return np.full(len(genes), "", dtype=object)
    return genes.to_numpy()


def build_pydeseq2_inference(threads: int, dds_class, stats_class):
    dds_kwargs: dict[str, Any] = {}
    stats_kwargs: dict[str, Any] = {}
    if "inference" in inspect.signature(dds_class).parameters:
        from pydeseq2.default_inference import DefaultInference

        inference = DefaultInference(n_cpus=threads)
        dds_kwargs["inference"] = inference
        if "inference" in inspect.signature(stats_class).parameters:
            stats_kwargs["inference"] = inference
    elif "n_cpus" in inspect.signature(dds_class).parameters:
        dds_kwargs["n_cpus"] = threads
        if "n_cpus" in inspect.signature(stats_class).parameters:
            stats_kwargs["n_cpus"] = threads
    return dds_kwargs, stats_kwargs


def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: np.ndarray,
    gene_symbols: np.ndarray,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("PyDESeq2 0.5.2 is required for compare_algorithm=DESeq2") from exc

    dds_kwargs, stats_kwargs = build_pydeseq2_inference(args.threads, DeseqDataSet, DeseqStats)
    model_metadata = metadata[["condition"]].copy()
    model_metadata["condition"] = pd.Categorical(
        model_metadata["condition"].astype(str),
        categories=[args.reference, args.comparison],
    )
    dds = DeseqDataSet(
        counts=counts.astype(np.int64),
        metadata=model_metadata,
        design="~condition",
        refit_cooks=True,
        **dds_kwargs,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["condition", args.comparison, args.reference],
        alpha=args.cut_off_pvalue,
        **stats_kwargs,
    )
    stats.summary()
    raw = stats.results_df.reindex(genes).copy()
    result = pd.DataFrame(
        {
            "gene": genes,
            "gene_symbol": gene_symbols,
            "base_mean": pd.to_numeric(raw["baseMean"], errors="coerce").to_numpy(),
            "log2FoldChange": pd.to_numeric(raw["log2FoldChange"], errors="coerce").to_numpy(),
            "statistic": pd.to_numeric(raw["stat"], errors="coerce").to_numpy(),
            "pvalue": pd.to_numeric(raw["pvalue"], errors="coerce").to_numpy(),
            "padj": pd.to_numeric(raw["padj"], errors="coerce").to_numpy(),
            "comparison": args.comparison,
            "reference": args.reference,
            "celltype": args.celltype,
            "algorithm": "DESeq2",
        }
    )
    try:
        dds.vst(use_design=False)
        normalized_values = np.asarray(dds.layers["vst_counts"], dtype=float)
        logger.info("Using PyDESeq2 variance-stabilized expression for the heatmap")
    except Exception as exc:  # noqa: BLE001 - retain a valid normalized transform
        if "normed_counts" in dds.layers:
            linear_normalized = np.asarray(dds.layers["normed_counts"], dtype=float)
        else:
            size_factors = np.asarray(dds.obs["size_factors"], dtype=float)
            linear_normalized = counts.to_numpy(dtype=float) / size_factors[:, None]
        normalized_values = np.log2(linear_normalized + 1.0)
        logger.warning(
            "PyDESeq2 VST was unavailable (%s); using log2(normalized counts + 1) for the heatmap",
            exc,
        )
    normalized = pd.DataFrame(normalized_values, index=counts.index, columns=genes)
    return result[FINAL_COLUMNS], normalized


def write_tsv(df: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression=compression, na_rep="NA")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, na_rep="NA")


def run_command(command: list[str], label: str) -> None:
    logger.info("Running %s", label)
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed with exit code {completed.returncode}")


def clean_known_outputs(output_dir: Path) -> None:
    for name in (
        "differential_expression.csv",
        "differential_expression.tsv.gz",
        "volcano.pdf",
        "heatmap.pdf",
        "GO_enrichment.csv",
        "GO_enrichment.pdf",
        "GO_BP_enrichment.pdf",
        "GO_CC_enrichment.pdf",
        "GO_MF_enrichment.pdf",
        "KEGG_enrichment.csv",
        "KEGG_enrichment.pdf",
        # Remove the previous combined outputs when reusing a result directory.
        "enrichment.tsv",
        "enrichment.pdf",
        "significant_genes.tsv",
        "GO_enrichment.tsv",
        "KEGG_enrichment.tsv",
    ):
        path = output_dir / name
        if path.exists() or path.is_symlink():
            path.unlink()


def go_enrichment_plot_paths(output_dir: Path, ontology: str) -> list[Path]:
    selected = str(ontology or "BP").strip().upper()
    if selected == "ALL":
        ontologies = ["BP", "CC", "MF"]
    elif selected in {"BP", "CC", "MF"}:
        ontologies = [selected]
    else:
        raise ValueError("compare_go_ontology must be BP, CC, MF, or ALL")
    return [output_dir / f"GO_{value}_enrichment.pdf" for value in ontologies]


def main() -> None:
    args = parse_args()
    algorithm = normalize_algorithm(args.algorithm)
    if args.threads < 1:
        raise ValueError("threads must be >= 1")
    if args.min_replicates < 2:
        raise ValueError("min_replicates must be >= 2")
    if args.comparison == args.reference:
        raise ValueError("comparison and reference must be different groups")
    if not 0 < args.cut_off_pvalue <= 1:
        raise ValueError("cut_off_pvalue must be in (0, 1]")
    if args.cut_off_logFC < 0:
        raise ValueError("cut_off_logFC must be >= 0")

    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(args.threads)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    clean_known_outputs(output_dir)
    hidden_files = [
        output_dir / ".pseudobulk_counts.tsv.gz",
        output_dir / ".sample_metadata.tsv",
        output_dir / ".normalized_expression.tsv.gz",
        output_dir / ".gene_symbols.tsv",
    ]

    try:
        log_step(logger, 1, 6, "validating sample metadata and SpatialData input")
        sample_table = read_sample_table(args.sample_list)
        available_groups = set(sample_table["group"])
        for group in (args.comparison, args.reference):
            if group not in available_groups:
                raise ValueError(f"Contrast group {group!r} is absent from sample.txt")

        sdata = spd.read_zarr(args.input_dir)
        adata = get_primary_table(sdata)
        obs = validate_obs_metadata(adata.obs, sample_table, args)
        target_mask = obs[args.celltype_col].astype(str).to_numpy() == str(args.celltype)
        if not target_mask.any():
            raise ValueError(
                f"celltype {args.celltype!r} is absent from obs[{args.celltype_col!r}]"
            )
        row_indices = np.flatnonzero(target_mask)
        matrix, genes, var, count_source = select_count_matrix(
            adata,
            args.count_layer,
            row_indices=row_indices,
        )
        gene_symbols = resolve_gene_symbols(
            genes,
            var,
            args.gene_symbol_col,
            args.gene_id_type,
        )
        logger.info("Using validated raw counts from %s", count_source)

        log_step(logger, 2, 6, f"aggregating target label {args.celltype!r} by sample")
        sample_values = obs.loc[target_mask, args.sample_col].astype(str).to_numpy()
        sample_order = sample_table["sample_id"].tolist()
        counts = aggregate_by_sample(matrix, row_indices, sample_values, sample_order, genes)
        units = pd.Series(sample_values).value_counts().reindex(sample_order, fill_value=0)
        metadata = build_sample_metadata(sample_table)

        selected = metadata["condition"].isin([args.comparison, args.reference])
        selected &= units.reindex(metadata.index).to_numpy() >= args.min_cells_per_sample
        metadata = metadata.loc[selected].copy()
        counts = counts.loc[metadata.index].copy()
        if metadata.empty:
            raise ValueError("No samples remain after target and min_cells_per_sample filtering")

        replicate_counts = metadata["condition"].value_counts()
        insufficient = {
            group: int(replicate_counts.get(group, 0))
            for group in (args.comparison, args.reference)
            if int(replicate_counts.get(group, 0)) < args.min_replicates
        }
        if insufficient:
            details = ", ".join(f"{key}={value}" for key, value in insufficient.items())
            raise ValueError(
                f"Each group requires at least {args.min_replicates} biological replicates after filtering: {details}"
            )

        gene_keep = counts.sum(axis=0).to_numpy() >= args.min_total_counts_per_gene
        gene_keep &= (counts > 0).sum(axis=0).to_numpy() >= args.min_replicates
        if int(gene_keep.sum()) < 2:
            raise ValueError("Fewer than two genes remain after pseudobulk expression filtering")
        counts = counts.loc[:, gene_keep].copy()
        genes = genes[gene_keep]
        gene_symbols = gene_symbols[gene_keep]

        log_step(logger, 3, 6, f"running {algorithm} for {args.comparison} versus {args.reference}")
        counts_file, metadata_file, normalized_file, gene_symbols_file = hidden_files
        counts_export = counts.copy()
        counts_export.insert(0, "replicate", counts_export.index.astype(str))
        write_tsv(counts_export, counts_file, compression="gzip")
        metadata_export = metadata.copy()
        metadata_export.to_csv(metadata_file, sep="\t", index=False)
        write_tsv(
            pd.DataFrame({"gene": genes.astype(str), "gene_symbol": gene_symbols.astype(str)}),
            gene_symbols_file,
        )

        result_file = output_dir / "differential_expression.csv"
        if algorithm == "DESeq2":
            result, normalized = run_deseq2(
                counts,
                metadata,
                genes,
                gene_symbols,
                args,
            )
            write_csv(result, result_file)
            normalized_export = normalized.copy()
            normalized_export.insert(0, "replicate", normalized_export.index.astype(str))
            write_tsv(normalized_export, normalized_file, compression="gzip")

        r_command = [
            args.rscript,
            args.differential_r,
            "--algorithm",
            algorithm,
            "--counts",
            str(counts_file),
            "--metadata",
            str(metadata_file),
            "--output_dir",
            str(output_dir),
            "--comparison",
            args.comparison,
            "--reference",
            args.reference,
            "--celltype",
            args.celltype,
            "--p_cutoff",
            str(args.cut_off_pvalue),
            "--logfc_cutoff",
            str(args.cut_off_logFC),
            "--top_n",
            str(args.de_top_n),
            "--result",
            str(result_file),
            "--normalized",
            str(normalized_file),
            "--gene_symbols",
            str(gene_symbols_file),
        ]
        run_command(r_command, f"{algorithm} R statistics/visualization")
        differential_outputs = [
            result_file,
            output_dir / "volcano.pdf",
            output_dir / "heatmap.pdf",
        ]
        missing_differential = [path.name for path in differential_outputs if not path.exists()]
        if missing_differential:
            raise RuntimeError(
                "Differential analysis did not create required output(s): "
                + ", ".join(missing_differential)
            )

        log_step(logger, 4, 6, "running separate GO/KEGG enrichment from the final DEG table")
        go_plot_paths = go_enrichment_plot_paths(output_dir, args.go_ontology)
        enrichment_command = [
            args.rscript,
            args.enrichment_r,
            "--input",
            str(result_file),
            "--output_go_table",
            str(output_dir / "GO_enrichment.csv"),
            "--output_kegg_table",
            str(output_dir / "KEGG_enrichment.csv"),
            "--output_go_plot",
            str(output_dir / "GO_enrichment.pdf"),
            "--output_kegg_plot",
            str(output_dir / "KEGG_enrichment.pdf"),
            "--species",
            args.species,
            "--gene_id_type",
            args.gene_id_type,
            "--padj_cutoff",
            str(args.cut_off_pvalue),
            "--logfc_cutoff",
            str(args.cut_off_logFC),
            "--go_ontology",
            args.go_ontology,
            "--top_n",
            str(args.enrichment_top_n),
        ]
        run_command(enrichment_command, "GO/KEGG enrichment")
        enrichment_outputs = [
            output_dir / "GO_enrichment.csv",
            *go_plot_paths,
            output_dir / "KEGG_enrichment.csv",
            output_dir / "KEGG_enrichment.pdf",
        ]
        missing_enrichment = [path.name for path in enrichment_outputs if not path.exists()]
        if missing_enrichment:
            raise RuntimeError(
                "Enrichment did not create required output(s): "
                + ", ".join(missing_enrichment)
            )

        log_step(logger, 5, 6, "checking final result contract")
        for path in hidden_files:
            try:
                path.unlink(missing_ok=True)
            except OSError as exc:
                raise RuntimeError(f"Unable to remove internal helper file {path}: {exc}") from exc
        allowed_outputs = {
            "differential_expression.csv",
            "volcano.pdf",
            "heatmap.pdf",
            "GO_enrichment.csv",
            *(path.name for path in go_plot_paths),
            "KEGG_enrichment.csv",
            "KEGG_enrichment.pdf",
            ".snakemake_timestamp",
        }
        unexpected = [
            path.name
            for path in output_dir.iterdir()
            if path.name not in allowed_outputs
        ]
        if unexpected:
            raise RuntimeError("Unexpected compare result files were created: " + ", ".join(unexpected))
        log_step(logger, 6, 6, "compare_analysis contrast completed")
    finally:
        for path in hidden_files:
            try:
                path.unlink(missing_ok=True)
            except OSError as exc:
                logger.warning("Unable to remove internal helper file %s: %s", path, exc)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        logger.error(str(error))
        sys.exit(1)
