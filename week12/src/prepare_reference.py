#!/usr/bin/env python3
"""Check the course reference and create transcript-to-gene and Salmon decoy files."""

import argparse
import csv
import re
import shutil
from pathlib import Path


def attribute(text, name):
    match = re.search(rf'(?:^|;\s*){name} "([^"]+)"', text)
    return match.group(1) if match else ""


def fasta_records(path):
    name, parts = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--transcripts", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    tx, locations, exons = {}, {}, {}
    with args.gtf.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in ("transcript", "exon"):
                continue
            tid = attribute(fields[8], "transcript_id")
            if fields[2] == "exon":
                exons.setdefault(tid, []).append((int(fields[3]), int(fields[4])))
                continue
            gid = attribute(fields[8], "gene_id")
            name = attribute(fields[8], "gene_name")
            if not tid or not gid or tid in tx:
                raise ValueError(f"Missing or duplicate GTF transcript ID: {tid}")
            tx[tid] = (gid, name)
            locations[tid] = (fields[0], fields[6])

    transcripts = dict(fasta_records(args.transcripts))
    if len(transcripts) != sum(line.startswith(">") for line in args.transcripts.open()):
        raise ValueError("Duplicate transcript FASTA IDs")
    if transcripts.keys() != tx.keys() or exons.keys() != tx.keys():
        raise ValueError(
            f"Transcript IDs differ: {len(tx.keys() - transcripts.keys())} only in GTF, "
            f"{len(transcripts.keys() - tx.keys())} only in FASTA"
        )

    genome = dict(fasta_records(args.genome))
    decoys = list(genome)
    if not decoys or transcripts.keys() & genome.keys():
        raise ValueError("Genome FASTA missing contigs or IDs overlap transcripts")
    complement = str.maketrans("ACGTN", "TGCAN")
    for tid, pieces in exons.items():
        chrom, strand = locations[tid]
        sequence = "".join(genome[chrom][start - 1:end] for start, end in sorted(pieces))
        if strand == "-":
            sequence = sequence.translate(complement)[::-1]
        if sequence != transcripts[tid]:
            raise ValueError(f"Transcript FASTA disagrees with genome/GTF: {tid}")

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "tx2gene.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id", "gene_name"])
        for tid in sorted(tx):
            writer.writerow([tid, *tx[tid]])
    with (args.out / "decoys.txt").open("w") as handle:
        handle.write("\n".join(decoys) + "\n")
    with (args.out / "salmon_gentrome.fa").open("wb") as output:
        for path in (args.transcripts, args.genome):
            with path.open("rb") as source:
                shutil.copyfileobj(source, output)
            output.write(b"\n")
    print(f"Validated {len(transcripts)} transcript sequences, {len(set(tx.values()))} gene/name pairs, {len(decoys)} decoys")


if __name__ == "__main__":
    main()
