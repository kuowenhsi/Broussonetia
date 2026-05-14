#!/usr/bin/env python3
"""Relabel IQ-TREE Newick tip labels with Broussonetia sample metadata."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


TIP_RE = re.compile(r"(?<=[(,])([^():,;]+)(?=:)")


def quote_newick_label(label: str) -> str:
    return "'" + label.replace("'", "''") + "'"


def clean_geo_name(value: str) -> str:
    return value.replace("\\,", ",").strip()


def build_label(row: dict[str, str]) -> str:
    organism = row["Organism"].strip()
    geo = clean_geo_name(row["geo_loc_name"])
    haplotype = row["haplotype"].strip()
    run = row["Run"].strip()
    return f"{organism} | {geo} | {haplotype} | {run}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Duplicate an IQ-TREE Newick treefile with metadata-rich tip labels."
    )
    parser.add_argument("--tree", required=True, type=Path, help="Input Newick treefile.")
    parser.add_argument("--metadata", required=True, type=Path, help="Sample metadata CSV.")
    parser.add_argument("--out-tree", required=True, type=Path, help="Relabeled output treefile.")
    parser.add_argument("--out-key", required=True, type=Path, help="Old/new label key CSV.")
    args = parser.parse_args()

    metadata: dict[str, dict[str, str]] = {}
    with args.metadata.open(newline="") as handle:
        for row in csv.DictReader(handle):
            metadata[row["Run"]] = row

    tree = args.tree.read_text()
    seen: set[str] = set()
    replacements: list[tuple[str, str, dict[str, str]]] = []

    def replace_tip(match: re.Match[str]) -> str:
        old_label = match.group(1)
        if old_label not in metadata:
            raise SystemExit(f"No metadata row found for tree tip: {old_label}")
        new_label = build_label(metadata[old_label])
        if new_label in seen:
            raise SystemExit(f"New label is not unique: {new_label}")
        seen.add(new_label)
        replacements.append((old_label, new_label, metadata[old_label]))
        return quote_newick_label(new_label)

    relabeled_tree = TIP_RE.sub(replace_tip, tree)

    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    args.out_tree.write_text(relabeled_tree)

    with args.out_key.open("w", newline="") as handle:
        fieldnames = [
            "old_label",
            "new_label",
            "Organism",
            "geo_loc_name",
            "haplotype",
            "isolate",
            "Sample Name",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for old_label, new_label, row in replacements:
            writer.writerow(
                {
                    "old_label": old_label,
                    "new_label": new_label,
                    "Organism": row["Organism"],
                    "geo_loc_name": clean_geo_name(row["geo_loc_name"]),
                    "haplotype": row["haplotype"],
                    "isolate": row["isolate"],
                    "Sample Name": row["Sample Name"],
                }
            )

    print(f"Wrote relabeled tree: {args.out_tree}")
    print(f"Wrote label key: {args.out_key}")
    print(f"Relabeled tips: {len(replacements)}")


if __name__ == "__main__":
    main()
