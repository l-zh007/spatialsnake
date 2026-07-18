import os
from typing import Any, Dict, Iterable, List, Optional

from spatialsnake.workflow.function.get_sample import read_compare_sample_table


class ConfigValidationError(ValueError):
    """Raised when CLI-level workflow configuration is invalid."""


NONE_VALUES = {"", "none", "null", "na", "nan"}
TRUE_VALUES = {"true", "1", "yes", "y", "t"}
FALSE_VALUES = {"false", "0", "no", "n", "f"}


CHOICE_RULES = {
    "option": {
        "choices": [
            "integrate",
            "preprocess",
            "clustering",
            "annotation_help",
            "annotation",
            "compare_stage",
            "advance_analysis",
            "reclustering",
            "splitting",
            "merge",
            "transform",
            "all",
        ],
        "case_sensitive": True,
    },
    "channel": {"choices": ["single_analysis", "compare_analysis"], "case_sensitive": True},
    "run_type": {
        "choices": ["visium", "visium_HD", "visium_segment", "xenium", "Merfish", "stereoseq"],
        "case_sensitive": True,
    },
    "batch_method": {"choices": ["harmony", "BBKNN"], "allow_none": True, "case_sensitive": False},
    "cluster_algorithm": {"choices": ["leiden", "louvain", "Kmeans", "kmeans"], "case_sensitive": False},
    "anno_algorithm": {"choices": ["manual", "reannotation", "cell2Location", "RCTD"], "case_sensitive": False},
    "RCTD_mode": {"choices": ["full", "doublet"], "case_sensitive": False},
    "device": {"choices": ["cpu", "cuda"], "case_sensitive": False},
    "compare_algorithm": {"choices": ["DESeq2", "edgeR"], "case_sensitive": False},
    "compare_gene_id_type": {
        "choices": ["auto", "SYMBOL", "ENSEMBL", "ENTREZID"],
        "case_sensitive": False,
    },
    "compare_go_ontology": {
        "choices": ["BP", "CC", "MF", "ALL"],
        "case_sensitive": False,
    },
    "species": {"choices": ["human", "mouse"], "case_sensitive": False},
    "runpipe": {
        "choices": ["cellPhoneDB", "pysenic", "liana", "cellcharter", "banksy", "cellchat", "compare_gene"],
        "case_sensitive": False,
    },
    "liana_method": {
        "choices": [
            "cellphonedb",
            "connectome",
            "logfc",
            "natmi",
            "singlecellsignalr",
            "rank_aggregate",
            "geometric_mean",
            "scseqcomm",
            "cellchat",
        ],
        "case_sensitive": False,
    },
    "cellchat_species": {"choices": ["human", "mouse"], "case_sensitive": False},
    "cellchat_db_subset": {
        "choices": ["all_interactions", "secreted_only", "secreted_ecm"],
        "case_sensitive": False,
    },
    "cellchat_compare_gene_plot_type": {
        "choices": ["dot", "violin", "bar"],
        "case_sensitive": False,
    },
    "transform_from": {"choices": ["zarr", "h5ad"], "case_sensitive": False},
    "transform_to": {"choices": ["zarr", "h5ad", "seurat"], "case_sensitive": False},
    "seurat_matrix": {"choices": ["auto", "raw", "X"], "case_sensitive": False},
    # ``split_by`` intentionally accepts any non-empty obs metadata column in
    # addition to the case-insensitive special modes handled by splitting.py.
    "subset_mode": {"choices": ["table", "spatial"], "case_sensitive": False},
    "annotation_format": {
        "choices": ["auto", "loupe", "xenium", "both", "none"],
        "case_sensitive": False,
    },
    "merge_by": {"choices": ["sample", "annotation"], "case_sensitive": False},
    "feature_join": {"choices": ["inner", "outer"], "case_sensitive": False},
    "label_policy": {"choices": ["preserve", "prefix", "offset"], "case_sensitive": False},
    "conflict_policy": {"choices": ["error", "first", "last"], "case_sensitive": False},
    "existing_policy": {
        "choices": ["overwrite_matched", "fill_missing", "error"],
        "case_sensitive": False,
    },
}


INT_RULES = {
    "jobs": {"min": 1},
    "threads": {"min": 1},
    "workers": {"min": 1},
    "max_cores": {"min": 1},
    "senic_workers": {"min": 1},
    "min_counts": {"min": 1},
    "min_cells": {"min": 1},
    "n_top_genes": {"min": 1},
    "n_comps": {"min": 1},
    "n_clusters": {"min": 2},
    "NEIGHBORS": {"min": 1},
    "pcs": {"min": 1},
    "recluster_neighbors": {"min": 1},
    "recluster_n_pcs": {"min": 1},
    "recluster_n_top_genes": {"min": 1},
    "min_replicates": {"min": 2},
    "min_cells_per_sample": {"min": 1},
    "min_total_counts_per_gene": {"min": 1},
    "de_top_n": {"min": 1},
    "compare_enrichment_top_n": {"min": 1},
    "iterations": {"min": 1},
    "liana_min_cells": {"min": 1},
    "liana_top_n": {"min": 1},
    "cellchat_min_cells": {"min": 1},
    "cellchat_workers": {"min": 1},
    "cellchat_interaction_length": {"min": 1},
    "cellchat_top_pathways": {"min": 1},
    "cellchat_top_cell_pairs": {"min": 1},
    "cellchat_bubble_top_lr": {"min": 0},
    "cellchat_compare_top_cell_pairs": {"min": 1},
    "cellchat_compare_top_pathways": {"min": 1},
    "cellchat_compare_top_lr": {"min": 1},
    "k_geom": {"min": 1},
    "rctd_dotplot_max_cell_types": {"min": 0},
}


FLOAT_RULES = {
    "mt_threshold": {"min": 0.0, "max": 100.0},
    "sample_rate": {"exclusive_min": 0.0, "max": 1.0},
    "resolution": {"exclusive_min": 0.0},
    "MIN_DIST": {"min": 0.0},
    "SPREAD": {"exclusive_min": 0.0},
    "recluster_resolution": {"exclusive_min": 0.0},
    "recluster_min_pct": {"min": 0.0, "max": 1.0},
    "recluster_logfc_threshold": {"min": 0.0},
    "point_size": {"exclusive_min": 0.0},
    "cut_off_pvalue": {"exclusive_min": 0.0, "max": 1.0},
    "cut_off_logFC": {"min": 0.0},
    "threshold": {"min": 0.0, "max": 1.0},
    "pvalue": {"exclusive_min": 0.0, "max": 1.0},
    "liana_expr_prop": {"min": 0.0, "max": 1.0},
    "liana_pvalue": {"exclusive_min": 0.0, "max": 1.0},
    "cellchat_trim": {"min": 0.0, "max": 1.0},
    "cellchat_spot_size": {"exclusive_min": 0.0},
    "cellchat_future_max_size_gb": {"min": 0.0},
    "max_m": {"min": 0.0},
    "rctd_dotplot_enrichment_clip": {"exclusive_min": 0.0},
    "min_match_rate": {"min": 0.0, "max": 1.0},
    "memory_limit_gb": {"min": 0.0},
}


NON_EMPTY_RULES = [
    "compare_sample_col",
    "compare_condition_col",
    "compare_celltype_col",
    "compare_input_zarr",
    "count_layer",
    "cell_focus",
    "celltype_col",
    "condition_col",
    "sample_col",
    "cellcharter_col",
]


def validate_workflow_config(config: Dict[str, Any]) -> None:
    """Validate key workflow config values before launching Snakemake."""
    errors: List[str] = []

    _validate_sample_list(config, errors)
    _validate_compare_sample_table(config, errors)
    _validate_compare_cellchat_table(config, errors)
    _validate_choice_rules(config, errors)
    _validate_int_rules(config, errors)
    _validate_float_rules(config, errors)
    _validate_non_empty_rules(config, errors)
    _validate_compare_contrasts(config, errors)
    _validate_sample_parameters(config, errors)
    _validate_image_slice(config, errors)
    _validate_removed_rctd_zarr_input(config, errors)

    if errors:
        joined = "\n".join(f"- {message}" for message in errors)
        raise ConfigValidationError(f"Invalid Spatialsnake configuration:\n{joined}")


def validate_splitting_config(config: Dict[str, Any]) -> None:
    """Validate useful-tool splitting options without workflow-only fields."""
    errors: List[str] = []
    split_by = config.get("split_by")
    if _is_none_like(split_by):
        errors.append("split_by must name an obs metadata column or a splitting mode.")

    for name in ("subset_mode", "annotation_format"):
        rule = CHOICE_RULES[name]
        value = config.get(name)
        choices = {str(choice).lower() for choice in rule["choices"]}
        normalized = str(value).strip().lower() if value is not None else ""
        # ``none`` is a real annotation output mode, not a missing value.
        if _is_none_like(value) and normalized not in choices:
            errors.append(f"{name} must not be empty. Allowed values: {_format_choices(rule['choices'])}.")
            continue
        if normalized not in choices:
            errors.append(f"{name} must be one of {_format_choices(rule['choices'])}; got {value!r}.")

    if errors:
        joined = "\n".join(f"- {message}" for message in errors)
        raise ConfigValidationError(f"Invalid Spatialsnake splitting configuration:\n{joined}")


def validate_merge_config(config: Dict[str, Any]) -> None:
    """Validate the direct-Python merge utility without workflow-only fields."""
    errors: List[str] = []
    for name in (
        "merge_by",
        "feature_join",
        "label_policy",
        "conflict_policy",
        "existing_policy",
    ):
        rule = CHOICE_RULES[name]
        value = config.get(name)
        normalized = str(value).strip().lower() if value is not None else ""
        choices = {str(choice).lower() for choice in rule["choices"]}
        if normalized not in choices:
            errors.append(
                f"{name} must be one of {_format_choices(rule['choices'])}; got {value!r}."
            )

    match_rate = _parse_float("min_match_rate", config.get("min_match_rate", 0.95), errors)
    if match_rate is not None and not 0.0 <= match_rate <= 1.0:
        errors.append(f"min_match_rate must be between 0 and 1; got {match_rate!r}.")

    for name in ("output_dir", "output_name", "sample_col", "target_col"):
        if _is_none_like(config.get(name)):
            errors.append(f"{name} must not be empty.")

    label_policy = str(config.get("label_policy", "preserve")).strip().lower()
    if label_policy != "preserve" and _is_none_like(config.get("label_cols")):
        errors.append(f"label_policy={label_policy!r} requires label_cols.")
    if str(config.get("merge_by", "")).strip().lower() == "sample" and not _is_none_like(
        config.get("annotation_csv")
    ):
        errors.append("annotation_csv is only valid when merge_by=annotation.")

    if errors:
        joined = "\n".join(f"- {message}" for message in errors)
        raise ConfigValidationError(f"Invalid Spatialsnake merge configuration:\n{joined}")


def validate_transform_config(config: Dict[str, Any]) -> None:
    """Validate the direct-Python transform utility and supported directions."""
    errors: List[str] = []
    source = str(config.get("transform_from", "")).strip().lower()
    target = str(config.get("transform_to", "")).strip().lower()
    allowed = {
        ("zarr", "h5ad"),
        ("zarr", "seurat"),
        ("h5ad", "zarr"),
        ("h5ad", "seurat"),
    }
    if (source, target) not in allowed:
        errors.append(
            "transform_from/transform_to must define one of zarr->h5ad, "
            "zarr->seurat, h5ad->zarr, or h5ad->seurat."
        )

    data_type = str(config.get("type", "auto")).strip().lower()
    if data_type not in {"auto", "st", "sc"}:
        errors.append(f"type must be auto, st, or sc; got {config.get('type')!r}.")

    matrix = str(config.get("seurat_matrix", "auto")).strip().lower()
    if matrix not in {"auto", "raw", "x"}:
        errors.append(
            f"seurat_matrix must be auto, raw, or X; got {config.get('seurat_matrix')!r}."
        )

    memory = _parse_float("memory_limit_gb", config.get("memory_limit_gb", 0), errors)
    if memory is not None and memory < 0:
        errors.append(f"memory_limit_gb must be >= 0; got {memory!r}.")

    for name in ("save_image", "keep_intermediate"):
        value = config.get(name, False)
        if not isinstance(value, bool) and str(value).strip().lower() not in TRUE_VALUES | FALSE_VALUES:
            errors.append(f"{name} must be a boolean; got {value!r}.")

    if _is_none_like(config.get("output_dir")):
        errors.append("output_dir must not be empty.")

    if errors:
        joined = "\n".join(f"- {message}" for message in errors)
        raise ConfigValidationError(f"Invalid Spatialsnake transform configuration:\n{joined}")


def _validate_sample_list(config: Dict[str, Any], errors: List[str]) -> None:
    sample_list = config.get("sample_list")
    if _is_none_like(sample_list):
        errors.append("sample_list is required for workflow commands.")
        return
    if not os.path.isfile(str(sample_list)):
        errors.append(f"sample_list file does not exist: {sample_list}")


def _validate_compare_sample_table(config: Dict[str, Any], errors: List[str]) -> None:
    """Validate the shared compare_gene sample contract before Snakemake starts."""
    if str(config.get("option", "")).strip() != "compare_stage":
        return
    if str(config.get("channel", "")).strip() != "compare_analysis":
        return
    if str(config.get("runpipe", "compare_gene")).strip().lower() != "compare_gene":
        return
    path = config.get("sample_list")
    if _is_none_like(path) or not os.path.isfile(str(path)):
        return
    try:
        rows = read_compare_sample_table(str(path))
    except (OSError, RuntimeError, ValueError) as exc:
        errors.append(str(exc))
        return
    groups = [row["group"] for row in rows]
    if len({group for group in groups if group}) < 2:
        errors.append("compare_gene sample.txt requires at least two biological groups.")


def _validate_compare_cellchat_table(config: Dict[str, Any], errors: List[str]) -> None:
    """CellChat compare-stage accepts exactly two condition-level RDS objects."""
    if str(config.get("option", "")).strip() != "compare_stage":
        return
    if str(config.get("channel", "")).strip() != "compare_analysis":
        return
    if str(config.get("runpipe", "")).strip().lower() != "cellchat":
        return
    path = config.get("sample_list")
    if _is_none_like(path) or not os.path.isfile(str(path)):
        return
    try:
        with open(str(path), encoding="utf-8") as handle:
            rows = [line.split() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    except OSError as exc:
        errors.append(f"Unable to read CellChat sample table: {exc}")
        return
    if rows and rows[0] and rows[0][0].strip().lower() == "sample_id":
        rows = rows[1:]
    if len(rows) != 2:
        errors.append(
            "compare-stage CellChat requires exactly two condition-level RDS rows; "
            "use advance_analysis for one condition or explicit pairwise runs for more conditions."
        )
        return
    malformed = [index + 1 for index, row in enumerate(rows) if len(row) < 2]
    if malformed:
        errors.append(f"CellChat sample table rows require sample_id and input_path columns: {malformed}.")


def _validate_choice_rules(config: Dict[str, Any], errors: List[str]) -> None:
    for name, rule in CHOICE_RULES.items():
        if name not in config:
            continue
        value = config.get(name)
        normalized_choices = {str(choice).lower() for choice in rule["choices"]}
        normalized_value = str(value).strip().lower() if value is not None else ""
        if _is_none_like(value) and normalized_value not in normalized_choices:
            if rule.get("allow_none", False):
                continue
            errors.append(f"{name} must not be empty. Allowed values: {_format_choices(rule['choices'])}.")
            continue
        choices = rule["choices"]
        if rule.get("case_sensitive", False):
            valid = str(value) in choices
        else:
            valid = normalized_value in normalized_choices
        if not valid:
            errors.append(f"{name} must be one of {_format_choices(choices)}; got {value!r}.")


def _validate_int_rules(config: Dict[str, Any], errors: List[str]) -> None:
    for name, rule in INT_RULES.items():
        if name not in config or _is_none_like(config.get(name)):
            continue
        value = _parse_int(name, config.get(name), errors)
        if value is None:
            continue
        if "min" in rule and value < rule["min"]:
            errors.append(f"{name} must be >= {rule['min']}; got {config.get(name)!r}.")
        if "max" in rule and value > rule["max"]:
            errors.append(f"{name} must be <= {rule['max']}; got {config.get(name)!r}.")


def _validate_float_rules(config: Dict[str, Any], errors: List[str]) -> None:
    for name, rule in FLOAT_RULES.items():
        if name not in config or _is_none_like(config.get(name)):
            continue
        value = _parse_float(name, config.get(name), errors)
        if value is None:
            continue
        if "min" in rule and value < rule["min"]:
            errors.append(f"{name} must be >= {rule['min']}; got {config.get(name)!r}.")
        if "exclusive_min" in rule and value <= rule["exclusive_min"]:
            errors.append(f"{name} must be > {rule['exclusive_min']}; got {config.get(name)!r}.")
        if "max" in rule and value > rule["max"]:
            errors.append(f"{name} must be <= {rule['max']}; got {config.get(name)!r}.")


def _validate_non_empty_rules(config: Dict[str, Any], errors: List[str]) -> None:
    for name in NON_EMPTY_RULES:
        if name in config and _is_none_like(config.get(name)):
            errors.append(f"{name} must not be empty.")


def _validate_compare_contrasts(config: Dict[str, Any], errors: List[str]) -> None:
    """Validate planned contrasts while accepting YAML mappings or CLI text."""
    if "compare_contrasts" not in config:
        return
    value = config.get("compare_contrasts")
    if value is None or value == [] or _is_none_like(value):
        return

    contrasts: List[Any]
    if isinstance(value, str):
        contrasts = [item.strip() for item in value.split(",") if item.strip()]
    elif isinstance(value, list):
        contrasts = value
    else:
        errors.append(
            "compare_contrasts must be a YAML list of comparison/reference mappings "
            "or comma-separated comparison:reference text."
        )
        return

    if not contrasts:
        return
    seen = set()
    for index, contrast in enumerate(contrasts, start=1):
        if isinstance(contrast, str):
            parts = [part.strip() for part in contrast.split(":")]
            if len(parts) != 2 or not all(parts):
                errors.append(
                    f"compare_contrasts item {index} must use comparison:reference; "
                    f"got {contrast!r}."
                )
                continue
            comparison, reference = parts
        elif isinstance(contrast, dict):
            comparison = str(contrast.get("comparison", "")).strip()
            reference = str(contrast.get("reference", "")).strip()
            if set(contrast) != {"comparison", "reference"} or not comparison or not reference:
                errors.append(
                    f"compare_contrasts item {index} must contain only non-empty "
                    "comparison and reference keys."
                )
                continue
        else:
            errors.append(
                f"compare_contrasts item {index} must be a mapping or "
                f"comparison:reference text; got {contrast!r}."
            )
            continue

        if comparison == reference:
            errors.append(
                f"compare_contrasts item {index} compares {comparison!r} with itself."
            )
        pair = (comparison, reference)
        if pair in seen:
            errors.append(f"compare_contrasts contains duplicate contrast {comparison}:{reference}.")
        seen.add(pair)


def _validate_sample_parameters(config: Dict[str, Any], errors: List[str]) -> None:
    """Validate YAML-only per-sample metadata and parameter overrides."""
    if "sample_parameters" not in config:
        return
    value = config.get("sample_parameters")
    if value is None:
        return
    if not isinstance(value, dict):
        errors.append("sample_parameters must be a YAML mapping keyed by sample_id.")
        return

    sample_ids = _read_sample_ids(config.get("sample_list"))
    reserved_keys = {"sample_id", "input_path", "group", "condition"}
    for sample_id, parameters in value.items():
        if _is_none_like(sample_id):
            errors.append("sample_parameters contains an empty sample_id key.")
            continue
        if sample_ids is not None and str(sample_id) not in sample_ids:
            errors.append(
                f"sample_parameters contains unknown sample_id {sample_id!r}; "
                "keys must match sample.txt."
            )
        if not isinstance(parameters, dict):
            errors.append(
                f"sample_parameters[{sample_id!r}] must be a mapping of parameter names to values."
            )
            continue

        invalid_keys = reserved_keys.intersection(parameters)
        if invalid_keys:
            errors.append(
                f"sample_parameters[{sample_id!r}] must not override "
                f"{', '.join(sorted(invalid_keys))}; keep these fields in sample.txt."
            )
        for field in ["replicate", "subject", "batch", "input_spec"]:
            if field in parameters and _is_none_like(parameters.get(field)):
                errors.append(
                    f"sample_parameters[{sample_id!r}][{field!r}] must not be empty."
                )
        for field in ["min_cells", "min_counts", "bin_size"]:
            if field not in parameters:
                continue
            parsed = _parse_int(
                f"sample_parameters[{sample_id!r}][{field!r}]",
                parameters.get(field),
                errors,
            )
            if parsed is not None and parsed < 1:
                errors.append(
                    f"sample_parameters[{sample_id!r}][{field!r}] must be >= 1; "
                    f"got {parameters.get(field)!r}."
                )
        if "mt_threshold" in parameters:
            parsed = _parse_float(
                f"sample_parameters[{sample_id!r}]['mt_threshold']",
                parameters.get("mt_threshold"),
                errors,
            )
            if parsed is not None and not 0.0 <= parsed <= 100.0:
                errors.append(
                    f"sample_parameters[{sample_id!r}]['mt_threshold'] must be "
                    f"between 0 and 100; got {parameters.get('mt_threshold')!r}."
                )


def _read_sample_ids(sample_list: Any) -> Optional[set]:
    """Read first-column sample IDs for cross-checking per-sample YAML keys."""
    if _is_none_like(sample_list) or not os.path.isfile(str(sample_list)):
        return None
    sample_ids = set()
    first_data_row = True
    with open(str(sample_list), encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if not fields:
                continue
            sample_id = fields[0].strip()
            if first_data_row and sample_id.lower() in {"sample", "sample_id"}:
                first_data_row = False
                continue
            first_data_row = False
            sample_ids.add(sample_id)
    return sample_ids


def _validate_image_slice(config: Dict[str, Any], errors: List[str]) -> None:
    if not _parse_bool(config.get("image_slice", False)):
        return
    coords = {}
    for name in ["x1", "x2", "y1", "y2"]:
        coords[name] = _parse_float(name, config.get(name), errors)
    if any(value is None for value in coords.values()):
        return
    if coords["x2"] <= coords["x1"]:
        errors.append(f"image_slice=True requires x2 > x1; got x1={coords['x1']}, x2={coords['x2']}.")
    if coords["y2"] <= coords["y1"]:
        errors.append(f"image_slice=True requires y2 > y1; got y1={coords['y1']}, y2={coords['y2']}.")


def _validate_removed_rctd_zarr_input(
    config: Dict[str, Any], errors: List[str]
) -> None:
    """Reject the removed duplicate RCTD spatial-input parameter."""
    if "zarr_input" not in config:
        return
    if str(config.get("option", "")).strip().lower() != "annotation":
        return
    if str(config.get("anno_algorithm", "")).strip().lower() != "rctd":
        return
    errors.append(
        "zarr_input was removed from RCTD. Move the SpatialData Zarr path to "
        "the second column of sample.txt and use: "
        "sample_id spatial_zarr sc_reference."
    )


def _parse_int(name: str, value: Any, errors: List[str]) -> Optional[int]:
    if isinstance(value, bool):
        errors.append(f"{name} must be an integer; got {value!r}.")
        return None
    try:
        text = str(value).strip()
        if not text or not text.lstrip("+-").isdigit():
            raise ValueError
        parsed = int(text)
    except (TypeError, ValueError):
        errors.append(f"{name} must be an integer; got {value!r}.")
        return None
    return parsed


def _parse_float(name: str, value: Any, errors: List[str]) -> Optional[float]:
    if isinstance(value, bool):
        errors.append(f"{name} must be numeric; got {value!r}.")
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        errors.append(f"{name} must be numeric; got {value!r}.")
        return None
    return parsed


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if _is_none_like(value):
        return False
    normalized = str(value).strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    return bool(value)


def _is_none_like(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().lower() in NONE_VALUES
    return False


def _format_choices(choices: Iterable[Any]) -> str:
    return ", ".join(str(choice) for choice in choices)
