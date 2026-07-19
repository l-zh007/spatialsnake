#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import conda_pack
import conda_pack.core


def malformed_records(prefix: Path) -> set[str]:
    records: set[str] = set()
    for record in (prefix / "conda-meta").glob("*.json"):
        with record.open(encoding="utf-8") as handle:
            metadata = json.load(handle)
        if "link" not in metadata or not isinstance(metadata.get("files"), list):
            records.add(record.name)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run conda-pack while treating malformed conda-meta records as unmanaged files."
    )
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--compress-level", type=int, default=1)
    args = parser.parse_args()

    prefix = Path(args.prefix).resolve()
    conda_meta = (prefix / "conda-meta").resolve()
    skipped = malformed_records(prefix)

    original_listdir = conda_pack.core.os.listdir

    def filtered_listdir(path: str | bytes | os.PathLike[str]) -> list[str]:
        entries = original_listdir(path)
        if Path(path).resolve() == conda_meta:
            return [entry for entry in entries if entry not in skipped]
        return entries

    if skipped:
        print(
            "Treating malformed conda metadata as unmanaged files: "
            + ", ".join(sorted(skipped))
        )

    conda_pack.core.os.listdir = filtered_listdir
    try:
        conda_pack.pack(
            prefix=str(prefix),
            output=args.output,
            compress_level=args.compress_level,
            force=True,
            ignore_missing_files=True,
        )
    finally:
        conda_pack.core.os.listdir = original_listdir

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
