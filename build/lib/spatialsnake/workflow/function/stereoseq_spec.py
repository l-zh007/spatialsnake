from __future__ import annotations


def parse_stereoseq_input_spec(input_spec):
    if input_spec in [None, "", False, "None", "False", "false", "NULL", "null"]:
        return None

    if isinstance(input_spec, (list, tuple)):
        items = list(input_spec)
    else:
        items = [part.strip() for part in str(input_spec).split(",") if part.strip()]

    parsed = []
    seen_modes = set()

    for item in items:
        key = str(item).strip().lower()
        if key in {"cellbin", "cell_bin", "cell"}:
            parsed.append("cellbin")
            seen_modes.add("cellbin")
            continue
        if key in {"adjusted_cellbin", "adjusted.cellbin", "adjusted-cellbin", "adjusted"}:
            parsed.append("adjusted_cellbin")
            seen_modes.add("cellbin")
            continue

        try:
            parsed.append(int(float(item)))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid Stereo-seq input_spec '{item}'. Use cellbin, adjusted_cellbin, "
                "or one/multiple comma-separated bin sizes such as 50 or 50,150."
            ) from exc
        seen_modes.add("bin")

    if len(seen_modes) > 1:
        raise ValueError(
            "Stereo-seq input_spec must use a single mode only: cellbin, adjusted_cellbin, "
            "or comma-separated bin sizes."
        )

    return parsed
