# RNA-seq quantification with six methods

I generated transcript and gene count matrices for the same six human RNA-seq FASTQs used in [my IGV analysis](../week11/README.md), then compared six quantification workflows. The [course dataset](https://www.biostarhandbook.com/courses/2026/appbio/rnaseq/rnaseq-tutorial-data/) is a chromosome 22 teaching subset adapted from the tutorial described by [Griffith et al. (2015)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004393). The [design file](design.csv) groups three Human Brain Reference (HBR) and three Universal Human Reference (UHR) files, representing two commercially available reference RNA materials.

## Inputs

| Course archive file | What it contains and how I used it |
|---|---|
| `fastq/*_R1.fq` | The measured reads: six 100-bp, single-end FASTQs with 118,571–227,392 reads each. HBR is brain-derived total RNA from 23 donors; UHR combines RNA from ten cancer cell lines. |
| `fasta/chr22.genome.fa` | One 50,818,468-bp chromosome 22 reference. STAR and HISAT2 align to it; RSEM builds transcript sequences from it, and Salmon uses it as a decoy. |
| `fasta/chr22.transcripts.fa` | A fixed reference transcriptome: 4,506 annotated, spliced transcript sequences from 1,371 genes, including protein-coding isoforms, noncoding RNAs, and pseudogenes. Salmon and kallisto index these sequences directly. |
| `gtf/chr22.gtf` | Gene IDs, names, and exon coordinates linking those 4,506 transcripts to 1,371 genes. It guides counting and supplies the transcript-to-gene mapping. |

The transcript FASTA sequences match the exons in the genome and GTF; they are a reference, not transcripts measured from HBR or UHR. Some routes count reads directly at genes; others estimate transcript counts and sum them by gene. [Reference checksums](report/reference_sha256.txt) identify the exact files used here.

## Transcript and gene count matrices

For the two requested matrices, I used Salmon to estimate counts for each annotated transcript, then [tximport](src/import_tximport.R) to sum transcripts by gene using the GTF-derived mapping.

| Output in `results/matrices/` | Rows | Meaning |
|---|---:|---|
| `salmon.transcript.expected.tsv` | 4,506 | Estimated count for each transcript in each of the six samples. |
| `salmon.gene.expected.tsv` | 1,371 | Sum of its transcript estimates for each gene and sample. |

Both matrices include gene names and the six samples in `design.csv` order. These are estimated counts, which may be fractional; the workflow writes TPM separately.

Do the gene counts support the [IGV coverage comparison](../week11/images/01_wide_common_scale.png)? The Salmon gene matrix does:

| Gene | Earlier IGV observation | Salmon gene counts |
|---|---|---|
| [`OR11H1`](../week11/images/03_or11h1_unexpressed.png) | No visible coverage | 0 in all six samples |
| [`DDX17`](../week11/images/04_ddx17_gene.png) | Strong coverage in both groups | About 771–1,412 in every sample |
| [`SULT4A1`](../week11/images/06_sult4a1_gene.png) | Stronger in HBR | 284–327 in HBR; 3–5 in UHR |

## Six-method comparison

| Workflow | How counts are obtained |
|---|---|
| STAR → featureCounts | Count reads assigned to genes after genome alignment. |
| STAR GeneCounts | Count genes during STAR alignment. |
| STAR → RSEM | Estimate transcript and gene counts from transcriptome alignments. |
| Salmon → tximport | Estimate transcript counts by selective alignment; sum them by gene. |
| kallisto → tximport | Estimate transcript counts by pseudoalignment; sum them by gene. |
| HISAT2 → StringTie → prepDE | Estimate reference transcript coverage and derive count-like values. |

The workflow produced six gene matrices and four transcript matrices in `results/matrices/`. Direct, estimated, and coverage-derived counts have different meanings, so their filenames keep those types separate.

The two direct gene-count methods agreed exactly. Across 1,336 genes, RSEM, Salmon, and kallisto correlated strongly with featureCounts (mean Pearson `r` = 0.93–0.97 on log1p chr22 CPM); StringTie was lower (`r` = 0.79). The [comparison report](report/comparison.md) contains the transcript comparisons, figures, and runtime.

## Reproduce

From `week12/` on Linux with Pixi installed:

```bash
pixi install
pixi run make salmon THREADS=4
```

This produces the two Salmon matrices above. To run all six workflows:

```bash
mkdir -p src/vendor
curl -fL https://raw.githubusercontent.com/gpertea/stringtie/v2.2.3/prepDE.py3 -o src/vendor/prepDE.py3
sha256sum -c report/prepde_sha256.txt
pixi run make -j2 all THREADS=4
```

The Makefile downloads the course RNA-seq archive into `data/` if needed. To use an already extracted copy, set `DATA_DIR=/path/to/data` when running `make`. `-j2` runs up to two jobs at once; `THREADS=4` sets threads per tool. The library is reverse-stranded. For single-end quantification, I assumed a fragment length of 200 ± 20 bp; change `FRAG_MEAN` and `FRAG_SD` if measured values are available. Software versions are pinned in [pixi.lock](pixi.lock).

I used AI to draft the scripts, then checked the sample order, strand setting, and transcript-to-gene mapping. During validation I found that Salmon dropped five duplicate-sequence transcripts by default and that prepDE's gene table misgrouped some transcripts. The workflow keeps the duplicates and builds the StringTie gene matrix by summing prepDE transcript counts with the checked GTF mapping.
