# IL6-associated pSTAT3 binding and chromatin accessibility

Siersbæk et al. reported that IL6/STAT3 signaling remodels a subset of
regulatory regions in T47D breast cancer cells. This analysis tests a limited
part of that result: whether IL6 increases pSTAT3 chromatin occupancy, and
whether the gained pSTAT3 peaks occur in regions that were already accessible
or became accessible following IL6 treatment.

## Data and design

The analysis uses hg38 and two biological replicates per assay and condition.
ATAC-seq samples come from [GSE150131](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150131),
and pSTAT3 ChIP-seq and Input DNA samples come from
[GSE126004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126004).

| Sample | GEO sample | Assay | Condition | SRA runs |
|---|---|---|---|---|
| `atac_control_r1` | GSM4524017 | ATAC | Control | SRR11735697 |
| `atac_control_r2` | GSM4524021 | ATAC | Control | SRR11735701 |
| `atac_il6_r1` | GSM4524018 | ATAC | IL6 | SRR11735698 |
| `atac_il6_r2` | GSM4524022 | ATAC | IL6 | SRR11735702 |
| `pstat3_control_r1` | GSM3587801 | pSTAT3 | Control | SRR8520836–SRR8520838 |
| `pstat3_control_r2` | GSM3587805 | pSTAT3 | Control | SRR8520848–SRR8520850 |
| `pstat3_il6_r1` | GSM3587802 | pSTAT3 | IL6 | SRR8520839–SRR8520841 |
| `pstat3_il6_r2` | GSM3587806 | pSTAT3 | IL6 | SRR8520851–SRR8520853 |
| `input_control` | GSM3587833 | Input | Control | SRR8520932–SRR8520937 |
| `input_il6` | GSM3587835 | Input | IL6 | SRR8520944–SRR8520949 |

ATAC-seq is paired-end; ChIP-seq and Input DNA are single-end. Multiple SRA
runs from the same GSM are merged before alignment and remain one biological
sample. The machine-readable sample mapping and ChIP/Input pairing are in
[`design.csv`](design.csv).

## Reproducing the analysis

The isolated Pixi environment records all command-line dependencies in
`pixi.toml` and `pixi.lock`. The completed analysis used Bowtie2 2.5.5,
SAMtools 1.25, MACS3 3.1.0, BEDTools 2.31.1, deepTools 3.5.6, FastQC 0.12.1,
and SRA Toolkit 3.4.1.

From this directory, run:

```bash
pixi install
pixi run make all JOBS=2 THREADS=8
```

`JOBS` is the number of samples processed concurrently, and `THREADS` is the
number of threads available to each sample. Use `JOBS=1 THREADS=4` on a
smaller computer. The full run needs substantial temporary disk space;
approximately 100 GB free is recommended.

The default `MAX_SPOTS=0` processes every spot. A short technical test can be
run without colliding with full-run outputs:

```bash
pixi run make all MAX_SPOTS=10000 JOBS=2 THREADS=4
```

The workflow performs:

```text
SRA runs -> merged sample FASTQ -> FastQC -> Bowtie2/hg38
         -> MAPQ >= 30 sorted BAM -> CPM bigWig -> MACS3 peaks
         -> replicate overlap and accessibility classification
```

ATAC peaks are called from paired fragments with MACS3 `BAMPE`. pSTAT3 peaks
are called against the condition-matched Input with a fixed 147-bp extension.
The bigWigs show CPM-normalized coverage; Input correction is applied to the
ChIP peak calls, not subtracted from the displayed bigWig signal.

## Peak calls and replicate agreement

| Sample | Peaks |
|---|---:|
| ATAC control replicate 1 | 90,529 |
| ATAC control replicate 2 | 78,077 |
| ATAC IL6 replicate 1 | 85,950 |
| ATAC IL6 replicate 2 | 82,132 |
| pSTAT3 control replicate 1 | 36 |
| pSTAT3 control replicate 2 | 18 |
| pSTAT3 IL6 replicate 1 | 62,278 |
| pSTAT3 IL6 replicate 2 | 67,034 |

Replicate overlap was measured as the base-pair Jaccard index between the two
peak files in each group.

| Assay | Condition | Jaccard | Overlapping intervals |
|---|---|---:|---:|
| ATAC | Control | 0.614 | 65,867 |
| ATAC | IL6 | 0.614 | 66,205 |
| pSTAT3 | Control | 0.372 | 15 |
| pSTAT3 | IL6 | 0.706 | 56,756 |

The IL6 pSTAT3 replicates therefore show strong agreement. The control
pSTAT3 Jaccard is based on very few peaks and is less informative.

## Accessibility at IL6-gained pSTAT3 peaks

A consensus peak is supported by both biological replicates. An IL6-gained
pSTAT3 peak is an IL6 consensus peak that does not overlap a peak in either
control pSTAT3 replicate. Gained peaks overlapping the IL6 ATAC consensus were
classified as follows:

- **Pre-existing accessible:** overlaps the control ATAC consensus.
- **Candidate newly accessible:** overlaps neither control ATAC replicate.
- **Control ambiguous:** detected in only one control ATAC replicate.

| Category | pSTAT3 peaks | Percentage of IL6 ATAC-open peaks |
|---|---:|---:|
| IL6-gained pSTAT3 | 54,973 | — |
| IL6-gained and ATAC-open | 33,702 | 100% |
| Pre-existing accessible | 30,987 | 91.9% |
| Candidate newly accessible | 935 | 2.8% |
| Control ambiguous | 1,780 | 5.3% |

Thus, 61.3% of reproducible IL6-gained pSTAT3 peaks overlapped accessible
chromatin under IL6 treatment. Most were already accessible before treatment,
while a smaller subset represented candidate IL6-associated newly accessible
regions.

## Representative regions in IGV

Load the built-in `hg38` genome in IGV. Then load the eight CPM bigWigs from
`coverage/` and the eight experimental narrowPeak files from `peaks/` in this
order: ATAC control replicates, ATAC IL6 replicates, pSTAT3 control replicates,
and pSTAT3 IL6 replicates. Apply the same data range to all four ATAC bigWigs
and another shared range to all four pSTAT3 bigWigs so that track autoscaling
does not exaggerate condition differences.

### Candidate newly accessible region upstream of BCL3

Go to `chr19:44,740,000-44,755,000`. The selected peak is
`chr19:44,744,822-44,745,915`, approximately 2–4 kb upstream of the BCL3
transcription start site. Both IL6 replicates show increased ATAC coverage and
strong pSTAT3 occupancy, whereas the control ATAC samples have low coverage
and no called peak at this interval.

Save the screenshot as `figures/bcl3_candidate_new_accessibility.png`.

![Candidate IL6-associated newly accessible region upstream of BCL3](figures/bcl3_candidate_new_accessibility.png)

### Pre-existing accessible region at the MYC promoter

Go to `chr8:127,731,000-127,746,000`. The selected peak is
`chr8:127,734,758-127,736,623` at the MYC promoter. ATAC coverage is present in
both control replicates and remains accessible after IL6 treatment, while
pSTAT3 occupancy increases in both IL6 replicates.

Save the screenshot as `figures/myc_preexisting_accessibility.png`.

![Pre-existing accessible region at the MYC promoter](figures/myc_preexisting_accessibility.png)

### Region where ATAC replicates disagree

Go to `chr19:18,269,000-18,281,000`. The interval
`chr19:18,274,459-18,275,754` has substantially higher ATAC coverage in
replicate 1 than replicate 2 in both conditions. It is called as an IL6 ATAC
peak only in replicate 1 and has no corresponding IL6 pSTAT3 enrichment. This
example shows why replicate-specific peaks were not treated as reproducible
chromatin changes.

Save the screenshot as `figures/atac_replicate_disagreement.png`.

![Region where ATAC replicates disagree](figures/atac_replicate_disagreement.png)

The mean CPM values for all three examples are recorded in
[`results/representative_regions.tsv`](results/representative_regions.tsv).

## Conclusion

IL6 treatment produced a strong and reproducible increase in pSTAT3 chromatin
occupancy. Of 54,973 reproducible IL6-gained pSTAT3 peaks, 33,702 (61.3%)
overlapped accessible chromatin in IL6-treated cells. Most of these regions
were already accessible in control cells, whereas 935 (2.8%) were candidate
newly accessible regions. The MYC promoter and an upstream region near BCL3
illustrate these two patterns, respectively. These observations are
consistent with the localized IL6-associated chromatin remodeling reported by
Siersbæk et al.; peak presence or absence here is not a formal differential
accessibility test.

## Output files

- [`results/peak_counts.tsv`](results/peak_counts.tsv): peaks per sample.
- [`results/replicate_overlap.tsv`](results/replicate_overlap.tsv): replicate
  Jaccard statistics.
- [`results/accessibility_classes.tsv`](results/accessibility_classes.tsv):
  pSTAT3/ATAC classification counts.
- [`results/representative_regions.tsv`](results/representative_regions.tsv):
  signal values for the two IGV examples.

Large FASTQ, BAM, bigWig, peak, reference, and intermediate files are excluded
from GitHub by `.gitignore` and are regenerated by the Makefile.

## Reference

Siersbæk R. et al. (2020). [IL6/STAT3 Signaling Hijacks Estrogen Receptor α
Enhancers to Drive Breast Cancer Metastasis](https://pmc.ncbi.nlm.nih.gov/articles/PMC7116707/).
*Cancer Cell* 38:412–423.e9.
