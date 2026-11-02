# Parallel bacterial variant calling

This workflow takes paired-end Illumina **whole-genome sequencing (WGS)** runs from SRA,
maps them to a reference genome, calls variants for each sample, and combines
the calls into one multisample VCF.

The workflow is designed for **small haploid bacterial genomes**. It assumes that
the runs are **paired-end** and that all samples are analysed against the same
reference genome.

## Quick start

From the repository root, run:

```bash
cd week09
pixi run -m "$HOME/edu/bioinfo" make all
```

The default reference downloads automatically when it is missing. 

## Analysis

Each row in `design.csv` is processed independently:

```text
SRA run
  -> paired FASTQ
  -> fastp
  -> bwa mem alignment
  -> sorted and indexed BAM
  -> bcftools VCF
```

The sample VCFs are merged after all sample jobs finish. This is a merge of
independent sample calls, not joint genotyping. The workflow also accepts a
single row for a small test, although multiple samples are needed for a
meaningful multisample comparison.

The workflow writes generated files to:

```text
reads/       downloaded and trimmed FASTQ files
bam/         sorted BAM files and BAM indexes
vcf/         sample VCFs and the merged VCF
logs/        download logs and parallel.joblog
reports/     fastp reports and VCF summary statistics
```

These directories are ignored by Git. The workflow definition consists of
`Makefile`, `design.csv`, and this README.

## Input files

The workflow root contains the design file:

```text
week09/
├── design.csv
├── Makefile
├── refs/        generated default reference and indexes
└── README.md
```

The first two columns of `design.csv` are required:

```text
Run,Sample
SRR8369177,H_pylori_SRR8369177
SRR8369178,H_pylori_SRR8369178
SRR8369179,H_pylori_SRR8369179
```

`Run` is an SRA accession. `Sample` is used in output file names and as the
sample ID in the VCF. Sample IDs must be unique. Additional metadata columns
can follow these two columns.

By default, the workflow uses the exact RefSeq assembly
`GCF_000008525.1` for *H. pylori* 26695:

```text
refs/Helicobacter_pylori_26695.fa
```

If this default FASTA is absent, `make` downloads the compressed genomic FASTA
from NCBI RefSeq, decompresses it safely, and saves it at that path. Failed or
incomplete downloads are removed. `bwa index` and `samtools faidx` then create
indexes beside the FASTA.

The reference remains configurable. When `REF=/path/to/custom.fa` is supplied,
that file must already exist; the workflow never downloads or overwrites a
custom reference. All samples in a run must be compatible with the selected
reference. FASTQ files do not need to be downloaded manually because the
workflow obtains them from the SRA accessions.

## Environment

This workflow uses only tools included in the standard course `bioinfo`
environment:

```text
make parallel sra-tools fastp bwa samtools bcftools tabix wget gzip
```

No additional package is required for this workflow. The Makefile does not
install or update software. The commands below explicitly use the course
environment at `$HOME/edu/bioinfo`, so it does not need to be activated first.

## Running the workflow

Run from the workflow root:

```bash
cd week09
pixi run -m "$HOME/edu/bioinfo" make all
```

The default settings are:

```text
N=50000       paired-end spots downloaded per SRA run
JOBS=3        samples processed simultaneously
THREADS=1     threads used inside each sample job
DESIGN=design.csv
REF=refs/Helicobacter_pylori_26695.fa
PROJECT=H_pylori
```

For a small test:

```bash
pixi run -m "$HOME/edu/bioinfo" make all N=10000 JOBS=2 THREADS=1
```

`N` controls the amount of sequence processed from each run. It is a runtime
setting, not a biological replicate count. `JOBS` controls the number of
samples running at once. `THREADS` controls the CPU threads used within each
sample job. A practical starting point is `JOBS * THREADS` no greater than the
available CPU cores, with enough memory for all concurrent jobs.

For another compatible project, keep its design and reference outside this
directory:

```text
repository-root/
├── other_project/
│   ├── design.csv
│   └── reference.fa
└── week09/
```

Run from `week09/` with:

```bash
pixi run -m "$HOME/edu/bioinfo" make all \
  DESIGN=../other_project/design.csv \
  REF=../other_project/reference.fa \
  PROJECT=other_project \
  JOBS=4 THREADS=2
```

## Parallel execution

The workflow parallelizes independent samples rather than individual commands
inside one sample. The Makefile defines the complete single-sample workflow;
GNU Parallel launches that target once for each row in `design.csv`:

```bash
cat design.csv | parallel --colsep , --header : --jobs 3 \
  'make sample RUN={Run} SAMPLE={Sample} N=50000 THREADS=1'
```

`--colsep ,` reads the CSV format. `--header :` makes the column names
available as placeholders, so `{Run}` and `{Sample}` are replaced for each
row. The final merge runs only after all sample VCFs and indexes exist.

Preview the commands without running them:

```bash
pixi run -m "$HOME/edu/bioinfo" parallel --dry-run \
  --colsep , --header : --jobs 3 \
  'make sample RUN={Run} SAMPLE={Sample} N=50000 THREADS=1' \
  :::: design.csv
```

`logs/parallel.joblog` records the command, runtime, and exit status of each
parallel job.

## SRA access

SRA download requires network access. Test one small run before starting a full
analysis:

```bash
pixi run -m "$HOME/edu/bioinfo" fastq-dump --split-files \
  --skip-technical -X 1000 -O /tmp/sra-test SRR8369177
```

If the direct download fails, cache the run first and rerun the workflow:

```bash
pixi run -m "$HOME/edu/bioinfo" prefetch SRR8369177
pixi run -m "$HOME/edu/bioinfo" make all
```

The accession remains in `design.csv`; `prefetch` only makes the data available
in the local SRA cache.

## Current example

### Genome selected

I selected the complete RefSeq genome for **Helicobacter pylori 26695**.

- Organism: `Helicobacter pylori 26695`
- RefSeq assembly: `GCF_000008525.1_ASM852v1`
- Assembly name: `ASM852v1`
- Data source: NCBI RefSeq FTP
- FASTA: `GCF_000008525.1_ASM852v1_genomic.fna.gz`

### Samples selected

The example contains three independent _Helicobacter pylori_ isolates from
BioProject `PRJNA362473`:

| Run            | Sample                  |
| -------------- | ----------------------- |
| `SRR8369177` | `H_pylori_SRR8369177` |
| `SRR8369178` | `H_pylori_SRR8369178` |
| `SRR8369179` | `H_pylori_SRR8369179` |

The reference is _H. pylori_ 26695. With the default settings, the main output
is:

```text
vcf/H_pylori_multisample_N50000.vcf.gz
vcf/H_pylori_multisample_N50000.vcf.gz.tbi
```

Other useful outputs are:

```text
reports/H_pylori_sample_summary_N50000.tsv
reports/H_pylori_merged_N50000.stats
logs/parallel.joblog
```

The `PROJECT` and `N` values are included in report names so that results from
different projects or test sizes do not overwrite one another.

To remove generated files and start again:

```bash
pixi run -m "$HOME/edu/bioinfo" make clean
```
