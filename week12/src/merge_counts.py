#!/usr/bin/env python3
"""Merge direct and model-based count tables into named sample matrices."""

import argparse
import csv
from pathlib import Path


def read_map(path):
    tx_gene = {}
    gene_name = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tx_gene[row["transcript_id"]] = row["gene_id"]
            if row["gene_id"] in gene_name and gene_name[row["gene_id"]] != row["gene_name"]:
                raise ValueError(f"Conflicting name for {row['gene_id']}")
            gene_name[row["gene_id"]] = row["gene_name"]
    return tx_gene, gene_name


def read_counts(path, method, level, star_column):
    counts = {}
    with path.open() as handle:
        if method == "star":
            for row in csv.reader(handle, delimiter="\t"):
                if row[0].startswith("N_"):
                    continue
                counts[row[0]] = row[star_column - 1]
        else:
            rows = (line for line in handle if not line.startswith("#"))
            for row in csv.DictReader(rows, delimiter="," if method == "prepde" else "\t"):
                key = ("Geneid" if method == "featurecounts" else
                       "gene_id" if level == "gene" else "transcript_id")
                value = (list(row.values())[-1] if method == "featurecounts" else
                         row["expected_count"] if method == "rsem" else
                         row[list(row)[1]])
                if row[key] in counts:
                    raise ValueError(f"Duplicate {key} in {path}: {row[key]}")
                counts[row[key]] = value
    return counts


def write_matrix(output, ids, per_sample, samples, labels, level):
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([f"{level}_id", "gene_name", *samples])
        for key in sorted(ids):
            writer.writerow([key, labels.get(key, ""), *[per_sample[s].get(key, 0) for s in samples]])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("method", choices=["star", "featurecounts", "rsem", "prepde"])
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--level", choices=["gene", "transcript"], default="gene")
    parser.add_argument("--star-column", type=int, default=4)
    args = parser.parse_args()
    with args.design.open() as handle:
        samples = [row["sample"] for row in csv.DictReader(handle)]
    if len(samples) != len(set(samples)) or not samples:
        raise ValueError("Design must contain unique samples")
    tx_gene, gene_name = read_map(args.mapping)
    labels = ({t: gene_name[g] for t, g in tx_gene.items()}
              if args.level == "transcript" else gene_name)
    if args.method == "star":
        paths = {s: args.input / s / "ReadsPerGene.out.tab" for s in samples}
    elif args.method == "featurecounts":
        paths = {s: args.input / f"{s}.txt" for s in samples}
    elif args.method == "rsem":
        suffix = "genes.results" if args.level == "gene" else "isoforms.results"
        paths = {s: args.input / f"{s}.{suffix}" for s in samples}
    else:
        paths = {s: args.input for s in samples}
    if args.method == "prepde":
        # prepDE's gene CSV can misgroup transcripts absent from its first sample.
        # Sum its transcript counts using the checked GTF transcript-to-gene map.
        with args.input.open() as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ["transcript_id", *samples]:
                raise ValueError(f"Unexpected prepDE columns: {reader.fieldnames}")
            data = list(reader)
        if len(data) != len({row["transcript_id"] for row in data}):
            raise ValueError("Duplicate transcript IDs in prepDE matrix")
        unknown_tx = {row["transcript_id"] for row in data} - set(tx_gene)
        if unknown_tx:
            raise ValueError(f"Unknown prepDE transcripts: {sorted(unknown_tx)[:3]}")
        if args.level == "transcript":
            per_sample = {s: {row["transcript_id"]: row[s] or 0 for row in data}
                          for s in samples}
        else:
            per_sample = {s: {} for s in samples}
            for row in data:
                gene = tx_gene[row["transcript_id"]]
                for sample in samples:
                    counts = per_sample[sample]
                    counts[gene] = counts.get(gene, 0) + int(row[sample] or 0)
    else:
        for path in paths.values():
            if not path.is_file():
                raise FileNotFoundError(path)
        per_sample = {s: read_counts(paths[s], args.method, args.level,
                                     args.star_column) for s in samples}
    ids = set().union(*(set(x) for x in per_sample.values()))
    unknown = ids - set(labels)
    if unknown:
        raise ValueError(f"{len(unknown)} output IDs absent from annotation, e.g. {sorted(unknown)[:3]}")
    for sample in samples:
        for value in per_sample[sample].values():
            if float(value) < 0:
                raise ValueError("Negative count")
    write_matrix(args.output, ids, per_sample, samples, labels, args.level)
    print(f"{args.method}: {len(ids)} {args.level} rows across {len(samples)} samples")


if __name__ == "__main__":
    main()
