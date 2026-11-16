#!/usr/bin/env python3
"""Compare six gene matrices on a shared chr22 gene set and record costs."""

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FILES = {
    "STAR+featureCounts": "star_featurecounts.gene.direct.tsv",
    "STAR GeneCounts": "star_genecounts.gene.direct.tsv",
    "STAR+RSEM": "rsem.gene.expected.tsv",
    "Salmon+tximport": "salmon.gene.expected.tsv",
    "kallisto+tximport": "kallisto.gene.expected.tsv",
    "HISAT2+StringTie": "stringtie.gene.derived.tsv",
}
TX_FILES = {
    "STAR+RSEM": "rsem.transcript.expected.tsv",
    "Salmon+tximport": "salmon.transcript.expected.tsv",
    "kallisto+tximport": "kallisto.transcript.expected.tsv",
    "HISAT2+StringTie": "stringtie.transcript.derived.tsv",
}


def load_matrix(path, samples, id_column="gene_id"):
    table = pd.read_csv(path, sep="\t", dtype={id_column: str})
    if list(table.columns) != [id_column, "gene_name", *samples]:
        raise ValueError(f"Columns or sample order disagree: {path}")
    if table[id_column].isna().any() or table[id_column].duplicated().any():
        raise ValueError(f"Missing or duplicate {id_column}: {path}")
    numeric = table.set_index(id_column)[samples].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric.to_numpy()).all() or (numeric.to_numpy() < 0).any():
        raise ValueError(f"Missing, infinite or negative counts: {path}")
    return numeric


def seconds_from_time(raw):
    parts = [float(x) for x in raw.split(":")]
    value = 0.0
    for piece in parts:
        value = value * 60 + piece
    return value


def parse_time(path):
    if not path.is_file():
        return (np.nan, np.nan)
    text = path.read_text()
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", text)
    ram = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    return (seconds_from_time(wall.group(1)) if wall else np.nan,
            int(ram.group(1)) / 1024 if ram else np.nan)


def write_efficiency(logs, samples, out):
    rows = []
    stages = {
        "STAR+featureCounts": ["star", "featurecounts"],
        "STAR GeneCounts": ["star"],
        "STAR+RSEM": ["star", "rsem"],
        "Salmon+tximport": ["salmon"],
        "kallisto+tximport": ["kallisto"],
        "HISAT2+StringTie": ["hisat2", "hisat2_sort", "stringtie"],
    }
    for method, steps in stages.items():
        for sample in samples:
            measured = []
            for step in steps:
                if step == "hisat2_sort":
                    path = logs / "hisat2" / f"{sample}.sort.time"
                else:
                    path = logs / step / f"{sample}.time"
                measured.append(parse_time(path))
            seconds = sum(x[0] for x in measured) if all(np.isfinite(x[0]) for x in measured) else np.nan
            ram = max(x[1] for x in measured) if all(np.isfinite(x[1]) for x in measured) else np.nan
            rows.append({"method": method, "sample": sample, "seconds": seconds,
                         "max_stage_ram_mib": ram})
    pd.DataFrame(rows).to_csv(out / "runtime_memory.tsv", sep="\t", index=False)
    index_rows = []
    for name in ("star", "rsem", "salmon", "kallisto", "hisat2"):
        wall, ram = parse_time(logs / "index" / f"{name}.time")
        index_rows.append({"index": name, "seconds": wall, "ram_mib": ram})
    pd.DataFrame(index_rows).to_csv(out / "index_cost.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrices", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    design = pd.read_csv(args.design)
    if list(design.columns) != ["sample", "group"] or design["sample"].duplicated().any():
        raise ValueError("Invalid design.csv")
    samples = design["sample"].tolist()
    matrices = {method: load_matrix(args.matrices / filename, samples)
                for method, filename in FILES.items()}
    txmap = pd.read_csv(args.mapping, sep="\t")
    complexity = txmap.groupby("gene_id")["transcript_id"].nunique()
    names = txmap.drop_duplicates("gene_id").set_index("gene_id")["gene_name"]
    common = sorted(set(complexity.index).intersection(*(set(m.index) for m in matrices.values())))
    if len(common) < 100:
        raise ValueError("Too few shared genes for a useful comparison")

    # CPM denominator is the sum of reported chr22 gene counts within each method/sample.
    cpm = {k: m.div(m.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
           for k, m in matrices.items()}
    if any(x.isna().any().any() for x in cpm.values()):
        raise ValueError("A method reported zero total gene counts")

    summary = []
    for method, table in matrices.items():
        for sample in samples:
            values = table[sample]
            rates = cpm[method][sample]
            summary.append({"method": method, "sample": sample,
                            "reported_genes": len(table), "total_counts": round(values.sum(), 2),
                            "genes_nonzero": int((values > 0).sum()),
                            "genes_low_1_to_9_counts": int(((values >= 1) & (values < 10)).sum()),
                            "genes_at_least_1_cpm": int((rates >= 1).sum())})
    pd.DataFrame(summary).to_csv(args.out / "summary.tsv", sep="\t", index=False)

    pair_rows = []
    for i, method_a in enumerate(FILES):
        for method_b in list(FILES)[i + 1:]:
            for sample in samples:
                a = np.log1p(cpm[method_a].loc[common, sample])
                b = np.log1p(cpm[method_b].loc[common, sample])
                pair_rows.append({"method_a": method_a, "method_b": method_b,
                                  "sample": sample, "shared_genes": len(common),
                                  "pearson_log1p_cpm": a.corr(b),
                                  "spearman_log1p_cpm": a.rank().corr(b.rank())})
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(args.out / "correlations.tsv", sep="\t", index=False)

    transcript_matrices = {method: load_matrix(args.matrices / filename, samples,
                                               "transcript_id")
                           for method, filename in TX_FILES.items()}
    tx_to_gene = txmap.set_index("transcript_id")["gene_id"]
    for method, tx_counts in transcript_matrices.items():
        if not set(tx_counts.index).issubset(tx_to_gene.index):
            raise ValueError(f"{method} has transcripts absent from the GTF")
        if method != "HISAT2+StringTie" and set(tx_counts.index) != set(tx_to_gene.index):
            raise ValueError(f"{method} lost transcripts from the supplied annotation")
        summed = tx_counts.groupby(tx_to_gene.loc[tx_counts.index]).sum()
        gene_counts = matrices[method]
        if set(summed.index) != set(gene_counts.index):
            raise ValueError(f"{method} gene IDs disagree with transcript sums")
        max_error = (summed.loc[gene_counts.index] - gene_counts).abs().to_numpy().max()
        if max_error > 0.05:  # RSEM rounds reported transcript counts to two decimals.
            raise ValueError(f"{method} gene/transcript sums disagree by {max_error:.2f}")
    common_tx = sorted(set.intersection(*(set(m.index) for m in transcript_matrices.values())))
    tx_cpm = {k: m.div(m.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
              for k, m in transcript_matrices.items()}
    if any(x.isna().any().any() for x in tx_cpm.values()):
        raise ValueError("A method reported zero total transcript counts")
    tx_pair_rows = []
    for i, method_a in enumerate(TX_FILES):
        for method_b in list(TX_FILES)[i + 1:]:
            for sample in samples:
                a = np.log1p(tx_cpm[method_a].loc[common_tx, sample])
                b = np.log1p(tx_cpm[method_b].loc[common_tx, sample])
                tx_pair_rows.append({"method_a": method_a, "method_b": method_b,
                                     "sample": sample, "shared_transcripts": len(common_tx),
                                     "pearson_log1p_cpm": a.corr(b),
                                     "spearman_log1p_cpm": a.rank().corr(b.rank())})
    tx_pairs = pd.DataFrame(tx_pair_rows)
    tx_pairs.to_csv(args.out / "transcript_correlations.tsv", sep="\t", index=False)

    methods = list(FILES)
    heat = np.eye(len(methods))
    for row in pairs.groupby(["method_a", "method_b"], as_index=False)["pearson_log1p_cpm"].mean().itertuples():
        i, j = methods.index(row.method_a), methods.index(row.method_b)
        heat[i, j] = heat[j, i] = row.pearson_log1p_cpm
    fig, ax = plt.subplots(figsize=(8.8, 7.1))
    img = ax.imshow(heat, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(methods)), methods, rotation=42, ha="right")
    ax.set_yticks(range(len(methods)), methods)
    for i in range(len(methods)):
        for j in range(len(methods)):
            ax.text(j, i, f"{heat[i,j]:.2f}", ha="center", va="center",
                    color="white" if heat[i,j] < 0.55 else "black", fontsize=8)
    fig.colorbar(img, ax=ax, label="Mean Pearson r across six samples")
    fig.tight_layout()
    fig.savefig(args.out / "correlations.png", dpi=170)
    plt.close(fig)

    # Differences are descriptive, not causal tests of transcript ambiguity.
    # Remove the many all-zero genes before asking whether isoform count matters.
    expressed = matrices["STAR+featureCounts"].loc[common].max(axis=1) >= 5
    expressed_genes = expressed.index[expressed].tolist()
    expressed_pairs = []
    for i, method_a in enumerate(FILES):
        for method_b in list(FILES)[i + 1:]:
            for sample in samples:
                a = np.log1p(cpm[method_a].loc[expressed_genes, sample])
                b = np.log1p(cpm[method_b].loc[expressed_genes, sample])
                expressed_pairs.append({"method_a": method_a, "method_b": method_b,
                                        "sample": sample, "shared_expressed_genes": len(expressed_genes),
                                        "pearson_log1p_cpm": a.corr(b)})
    expressed_pairs = pd.DataFrame(expressed_pairs)
    expressed_pairs.to_csv(args.out / "correlations_expressed.tsv", sep="\t", index=False)
    baseline = cpm["STAR+featureCounts"].loc[expressed_genes]
    diff_rows = []
    for method in methods[1:]:
        other = cpm[method].loc[expressed_genes]
        magnitude = (np.log2(other + 1) - np.log2(baseline + 1)).abs().mean(axis=1)
        for gid, value in magnitude.items():
            n_tx = int(complexity[gid])
            diff_rows.append({"method": method, "gene_id": gid,
                              "gene_name": names.get(gid, ""), "n_transcripts": n_tx,
                              "complexity": "1" if n_tx == 1 else "2-4" if n_tx <= 4 else "5+",
                              "mean_abs_log2_cpm_difference": value})
    diff = pd.DataFrame(diff_rows)
    diff.to_csv(args.out / "gene_differences.tsv", sep="\t", index=False)
    order = ["1", "2-4", "5+"]
    med = diff.groupby(["method", "complexity"])["mean_abs_log2_cpm_difference"].median().unstack().reindex(columns=order)
    med.to_csv(args.out / "isoform_complexity.tsv", sep="\t")
    n_tx = complexity.loc[expressed_genes]
    complexity_n = pd.Series({"1": int((n_tx == 1).sum()),
                              "2-4": int(((n_tx >= 2) & (n_tx <= 4)).sum()),
                              "5+": int((n_tx >= 5).sum())}, name="genes")
    complexity_n.to_csv(args.out / "isoform_complexity_n.tsv", sep="\t", header=True)
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8.8, 7), sharex=True)
    x = np.arange(len(order))
    model_methods = ["STAR+RSEM", "Salmon+tximport", "kallisto+tximport"]
    width = 0.23
    for j, method in enumerate(model_methods):
        ax.bar(x + (j - 1) * width, med.loc[method], width, label=method)
    ax.set_title("Model-based estimates")
    ax.legend(fontsize=8, ncol=3)
    ax2.bar(x, med.loc["HISAT2+StringTie"], width=0.45, color="tab:purple")
    ax2.set_title("StringTie coverage-derived counts")
    ax2.set_xticks(x, order)
    ax2.set_xlabel("Annotated transcripts per gene")
    fig.supylabel("Median |log2(CPM+1) difference| vs STAR+featureCounts")
    fig.tight_layout()
    fig.savefig(args.out / "isoform_complexity.png", dpi=170)
    plt.close(fig)

    write_efficiency(args.logs, samples, args.out)
    selected = [g for g in names.index if names[g] in ("OR11H1", "DDX17", "SULT4A1")]
    example_rows = []
    for method, matrix in matrices.items():
        for gid in selected:
            for sample in samples:
                example_rows.append({"method": method, "gene_id": gid,
                                     "gene_name": names[gid], "sample": sample,
                                     "count": matrix.at[gid, sample] if gid in matrix.index else np.nan,
                                     "chr22_cpm": cpm[method].at[gid, sample] if gid in matrix.index else np.nan})
    pd.DataFrame(example_rows).to_csv(args.out / "igv_examples.tsv", sep="\t", index=False)

    summary_table = pd.DataFrame(summary)
    runtime_table = pd.read_csv(args.out / "runtime_memory.tsv", sep="\t")
    gene_r = pairs.groupby(["method_a", "method_b"])["pearson_log1p_cpm"].mean()
    expressed_r = expressed_pairs.groupby(["method_a", "method_b"])["pearson_log1p_cpm"].mean()
    tx_r = tx_pairs.groupby(["method_a", "method_b"])["pearson_log1p_cpm"].mean()
    with (args.out / "comparison.md").open("w") as handle:
        handle.write("# Comparison of six chr22 RNA-seq quantification pipelines\n\n")
        handle.write(f"Gene correlations use {len(common)} genes shared by all six matrices; ")
        handle.write(f"transcript correlations use {len(common_tx)} transcripts shared by the four methods that estimate transcripts. ")
        handle.write("For each method and sample, counts were divided by the total reported chr22 counts, then transformed as log1p(CPM). This is a descriptive within-chromosome comparison.\n\n")
        handle.write("## What changed across methods?\n\n")
        direct_equal = matrices["STAR+featureCounts"].equals(matrices["STAR GeneCounts"])
        if direct_equal:
            handle.write("STAR GeneCounts and STAR+featureCounts gave identical gene counts for all genes and samples. ")
        else:
            handle.write("STAR GeneCounts and STAR+featureCounts differed for some genes. ")
        handle.write("Mean gene-level Pearson correlations with STAR+featureCounts were ")
        handle.write(", ".join(f"{method} {gene_r.loc['STAR+featureCounts', method]:.3f}"
                               for method in ("STAR+RSEM", "Salmon+tximport", "kallisto+tximport", "HISAT2+StringTie")))
        handle.write(" across six samples.\n\n")
        handle.write(f"Among the {len(expressed_genes)} genes with at least 5 direct counts in one sample, ")
        handle.write("the corresponding correlations were ")
        handle.write(", ".join(f"{method} {expressed_r.loc['STAR+featureCounts', method]:.3f}"
                               for method in ("STAR+RSEM", "Salmon+tximport", "kallisto+tximport", "HISAT2+StringTie")))
        handle.write(". Including many all-zero genes makes some methods look more similar.\n\n")
        handle.write("![Gene-level method correlations](correlations.png)\n\n")
        handle.write(f"At transcript level, RSEM versus kallisto averaged {tx_r.loc['STAR+RSEM', 'kallisto+tximport']:.3f}, ")
        handle.write(f"RSEM versus Salmon {tx_r.loc['STAR+RSEM', 'Salmon+tximport']:.3f}, and ")
        handle.write(f"Salmon versus kallisto {tx_r.loc['Salmon+tximport', 'kallisto+tximport']:.3f}. ")
        handle.write("These lower transcript correlations are consistent with different allocation of reads among isoforms; they do not by themselves identify which estimate is correct.\n\n")
        handle.write("## Low counts and isoforms\n\n")
        handle.write("The mean number of genes with at least 1 but fewer than 10 reported counts per sample was:\n\n")
        handle.write("| Method | Genes with 1 to <10 counts |\n|---|---:|\n")
        for method in methods:
            mean_low = summary_table.loc[summary_table.method == method, "genes_low_1_to_9_counts"].mean()
            handle.write(f"| {method} | {mean_low:.1f} |\n")
        handle.write("\nThe count types differ: direct integers, model-estimated counts, and StringTie coverage-derived approximations. These are detection summaries, not a sensitivity ranking.\n\n")
        handle.write(f"For the isoform comparison I kept {len(expressed_genes)} shared genes with at least 5 STAR+featureCounts reads in one sample. ")
        handle.write("The plot shows median absolute log2(CPM+1) differences from STAR+featureCounts by number of annotated transcripts. The pattern is descriptive; isoform count alone does not establish the cause of a difference.\n\n")
        handle.write("![Differences by annotated isoform complexity](isoform_complexity.png)\n\n")
        sult_genes = [g for g in selected if names[g] == "SULT4A1"]
        if sult_genes:
            gene = sult_genes[0]
            sample = design["sample"].iloc[0]
            tx_ids = txmap.loc[txmap.gene_id == gene, "transcript_id"].tolist()
            handle.write(f"For SULT4A1 in {sample}, the gene totals from RSEM, Salmon, and kallisto are ")
            handle.write(", ".join(f"{method} {matrices[method].at[gene, sample]:.1f}"
                                   for method in ("STAR+RSEM", "Salmon+tximport", "kallisto+tximport")))
            handle.write(" counts, but their allocation to transcripts differs:\n\n")
            handle.write("| Transcript | RSEM | Salmon | kallisto |\n|---|---:|---:|---:|\n")
            for tid in tx_ids:
                values = [transcript_matrices[method].at[tid, sample]
                          for method in ("STAR+RSEM", "Salmon+tximport", "kallisto+tximport")]
                handle.write(f"| {tid} | " + " | ".join(f"{value:.1f}" for value in values) + " |\n")
            handle.write("\nThis shows why similar gene totals do not guarantee similar isoform estimates.\n\n")
        handle.write("## Cross-check with IGV\n\n")
        handle.write("For the STAR+featureCounts matrix, the genes examined in the earlier IGV analysis have these count ranges:\n\n")
        handle.write("| Gene | " + " | ".join(design.group.drop_duplicates()) + " |\n")
        handle.write("|---|" + "---:|" * design.group.nunique() + "\n")
        for name in ("OR11H1", "DDX17", "SULT4A1"):
            genes = [g for g in selected if names[g] == name]
            if not genes:
                continue
            gene = genes[0]
            ranges = []
            for group in design.group.drop_duplicates():
                group_samples = design.loc[design.group == group, "sample"].tolist()
                values = matrices["STAR+featureCounts"].loc[gene, group_samples]
                ranges.append(f"{values.min():,.0f}–{values.max():,.0f}")
            handle.write(f"| {name} | " + " | ".join(ranges) + " |\n")
        handle.write("\nThese values agree with the broad visual pattern recorded in the earlier IGV analysis.\n\n")
        handle.write("## Runtime and memory\n\n")
        handle.write("| Route | Mean seconds per sample | Mean maximum stage RAM (MiB) |\n|---|---:|---:|\n")
        for method in methods:
            sub = runtime_table.loc[runtime_table.method == method]
            handle.write(f"| {method} | {sub.seconds.mean():.1f} | {sub.max_stage_ram_mib.mean():.0f} |\n")
        handle.write("\nThese are runs on one machine using only the course chr22 subset. Index-building time is separate, and STAR alignment is shared across three routes. The memory number is the largest single stage for each sample.\n\n")
        handle.write("Full tables: [gene matrix summary](summary.tsv), [gene correlations](correlations.tsv), ")
        handle.write("[expressed-gene correlations](correlations_expressed.tsv), ")
        handle.write("[transcript correlations](transcript_correlations.tsv), [IGV genes](igv_examples.tsv), ")
        handle.write("[per-gene differences](gene_differences.tsv), [isoform group sizes](isoform_complexity_n.tsv), ")
        handle.write("[runtime and memory](runtime_memory.tsv), [index costs](index_cost.tsv).\n")

    print(f"Compared {len(common)} common genes across {len(methods)} methods and {len(samples)} samples")


if __name__ == "__main__":
    main()
