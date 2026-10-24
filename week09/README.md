# Parallel bacterial variant calling

This workflow takes paired-end Illumina **whole-genome sequencing (WGS)** runs from SRA,
maps them to a reference genome, calls variants for each sample, and combines
the calls into one multisample VCF.

The workflow is designed for **small haploid bacterial genomes**. It assumes that
the runs are **paired-end** and that all samples are analysed against the same
reference genome.

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
independent sample calls, not joint genotyping.

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

The reference FASTA is supplied with `REF`. The current example uses:

```text
../week02/fasta/Helicobacter_pylori_26695.fa
```

Reference indexes are created beside the FASTA. FASTQ files do not need to be
downloaded manually; the workflow obtains them from the SRA accessions.

## Environment

The commands use the course Pixi environment at `$HOME/edu/bioinfo`. It
contains the tools used here:

```text
make parallel sra-tools fastp bwa samtools bcftools tabix
```

To create that environment:

```bash
mkdir -p "$HOME/edu/bioinfo"
cd "$HOME/edu/bioinfo"
pixi init
pixi workspace channel add bioconda
pixi project platform add osx-arm64 linux-64
pixi add make bwa samtools bcftools trimmomatic sra-tools fastqc seqkit fastp parallel \
  csvtk wget ncbi-datasets-cli ucsc-bedgraphtobigwig bedtools \
  jq python=3.12 'biopython<1.86' bio
```

The Makefile uses the existing environment; it does not install packages.

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
REF=../week02/fasta/Helicobacter_pylori_26695.fa
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

Other useful outputs are `reports/sample_summary.tsv`,
`reports/merged.stats`, and `logs/parallel.joblog`.

To remove generated files and start again:

```bash
pixi run -m "$HOME/edu/bioinfo" make clean
```
