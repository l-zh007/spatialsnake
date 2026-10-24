from __future__ import annotations


def parse_stereoseq_input_spec(input_spec):
    if input_spec in [None, "", False, "None", "False", "false", "NULL", "null"]:
        return None

    if isinstance(input_spec, (list, tuple)):
        items = list(input_spec)
    else:
        items = [part.strip() for part in str(input_spec).split(",") if part.strip()]

    parsed = []

    for item in items:
        key = str(item).strip().lower()
        if key in {"cellbin", "cell_bin", "cell"}:
            parsed.append("cellbin")
            continue
        if key in {"adjusted_cellbin", "adjusted.cellbin", "adjusted-cellbin", "adjusted"}:
            parsed.append("adjusted_cellbin")
            continue

        try:
            value = float(item)
            if not value.is_integer() or value <= 0:
                raise ValueError("bin size must be a positive integer")
            parsed.append(int(value))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid Stereo-seq input_spec '{item}'. Use cellbin, adjusted_cellbin, "
                "or one bin size such as 50."
            ) from exc

    if len(parsed) != 1:
        raise ValueError("Select one Stereo-seq bin size, cellbin, or adjusted_cellbin per analysis; run other resolutions separately.")
    return parsed
