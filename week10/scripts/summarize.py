#!/usr/bin/env python3

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def run(command, output):
    with output.open("w") as handle:
        subprocess.run(command, check=True, stdout=handle)


def main(design_path, peak_dir, result_dir, suffix):
    peak_dir = Path(peak_dir)
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)

    with open(design_path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    samples = [row for row in rows if row["Assay"] != "Input"]
    groups = defaultdict(list)
    for row in samples:
        peak = peak_dir / f'{row["Sample"]}{suffix}_peaks.narrowPeak'
        if not peak.is_file():
            raise FileNotFoundError(peak)
        groups[(row["Assay"], row["Group"])].append(row)

    with (result_dir / f"peak_counts{suffix}.tsv").open("w") as handle:
        handle.write("sample\tassay\tgroup\treplicate\tpeaks\n")
        for row in samples:
            peak = peak_dir / f'{row["Sample"]}{suffix}_peaks.narrowPeak'
            count = sum(1 for line in peak.open() if line.strip())
            handle.write(
                f'{row["Sample"]}\t{row["Assay"]}\t{row["Group"]}\t'
                f'{row["Replicate"]}\t{count}\n'
            )

    consensus = {}
    with (result_dir / f"replicate_overlap{suffix}.tsv").open("w") as handle:
        handle.write(
            "assay\tgroup\treplicate_1\treplicate_2\tintersection_bp\t"
            "union_bp\tjaccard\toverlapping_intervals\n"
        )
        for (assay, group), members in sorted(groups.items()):
            members.sort(key=lambda row: int(row["Replicate"]))
            if len(members) != 2:
                raise ValueError(f"{assay} {group} requires exactly two replicates")
            peak1 = peak_dir / f'{members[0]["Sample"]}{suffix}_peaks.narrowPeak'
            peak2 = peak_dir / f'{members[1]["Sample"]}{suffix}_peaks.narrowPeak'
            result = subprocess.run(
                ["bedtools", "jaccard", "-a", peak1, "-b", peak2],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip().splitlines()[-1].split("\t")
            handle.write(
                f"{assay}\t{group}\t{members[0]['Sample']}\t"
                f"{members[1]['Sample']}\t" + "\t".join(result) + "\n"
            )

            label = f"{assay.lower()}_{group.lower()}"
            consensus_path = result_dir / f"{label}_consensus{suffix}.bed"
            run(["bedtools", "intersect", "-u", "-a", peak1, "-b", peak2], consensus_path)
            consensus[(assay, group)] = consensus_path

            rep1_only = result_dir / f"{label}_rep1_only{suffix}.bed"
            run(["bedtools", "intersect", "-v", "-a", peak1, "-b", peak2], rep1_only)

    shared = result_dir / f"shared_il6_atac_pstat3{suffix}.bed"
    run(
        [
            "bedtools",
            "intersect",
            "-u",
            "-a",
            consensus[("pSTAT3", "IL6")],
            "-b",
            consensus[("ATAC", "IL6")],
        ],
        shared,
    )

    pstat3_control = [
        peak_dir / f'{row["Sample"]}{suffix}_peaks.narrowPeak'
        for row in samples
        if row["Assay"] == "pSTAT3" and row["Group"] == "Control"
    ]
    atac_control = [
        peak_dir / f'{row["Sample"]}{suffix}_peaks.narrowPeak'
        for row in samples
        if row["Assay"] == "ATAC" and row["Group"] == "Control"
    ]

    gained = result_dir / f"pstat3_il6_gained{suffix}.bed"
    run(
        ["bedtools", "intersect", "-v", "-a", consensus[("pSTAT3", "IL6")],
         "-b", *pstat3_control],
        gained,
    )

    gained_open = result_dir / f"pstat3_il6_gained_atac_open{suffix}.bed"
    run(
        ["bedtools", "intersect", "-u", "-a", gained,
         "-b", consensus[("ATAC", "IL6")]],
        gained_open,
    )

    preexisting = result_dir / f"pstat3_il6_preexisting_open{suffix}.bed"
    run(
        ["bedtools", "intersect", "-u", "-a", gained_open,
         "-b", consensus[("ATAC", "Control")]],
        preexisting,
    )

    candidate_new = result_dir / f"pstat3_il6_candidate_new_open{suffix}.bed"
    run(
        ["bedtools", "intersect", "-v", "-a", gained_open,
         "-b", *atac_control],
        candidate_new,
    )

    ambiguous = result_dir / f"pstat3_il6_atac_control_ambiguous{suffix}.bed"
    run(
        ["bedtools", "intersect", "-v", "-a", gained_open,
         "-b", preexisting, candidate_new],
        ambiguous,
    )

    def line_count(path):
        return sum(1 for line in path.open() if line.strip())

    with (result_dir / f"accessibility_classes{suffix}.tsv").open("w") as handle:
        handle.write("class\tpstat3_peaks\n")
        handle.write(f"IL6_pSTAT3_consensus\t{line_count(consensus[('pSTAT3', 'IL6')])}\n")
        handle.write(f"IL6_gained_pSTAT3\t{line_count(gained)}\n")
        handle.write(f"IL6_gained_and_ATAC_open\t{line_count(gained_open)}\n")
        handle.write(f"preexisting_accessible\t{line_count(preexisting)}\n")
        handle.write(f"candidate_newly_accessible\t{line_count(candidate_new)}\n")
        handle.write(f"control_ATAC_ambiguous\t{line_count(ambiguous)}\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("usage: summarize.py DESIGN PEAK_DIR RESULT_DIR SUFFIX")
    main(*sys.argv[1:])

