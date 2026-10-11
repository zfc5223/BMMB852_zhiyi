# Week 7: Generate a VCF File

This assignment calls small variants from an alignment to the _Helicobacter pylori_ 26695 reference genome.

## Variant Caller Choice

I considered `bcftools`, FreeBayes, and GATK. I used `bcftools mpileup` and `bcftools call --ploidy 1` because it is already available in the course pixi environment, works directly with `samtools`/`htslib`, and is appropriate for a small bacterial short-read example. The `--ploidy 1` setting is used because _H. pylori_ is haploid. FreeBayes is also reasonable for microbial short-read data, but it is not installed in the current course environment. GATK is powerful but heavier than needed for this small assignment.

The commands are automated in the `Makefile`.

```bash
pixi run -m "$HOME/edu/bioinfo" make all
```

## Input Data

I did not reuse the previous Week 5 BAM as the main input because that file was made from only a tiny subset of transcriptomic reads from `SRR031127`. It had very low genome-wide coverage and was not a good basis for variant calling.

For this VCF workflow I used:

- Reference: `../week02/fasta/Helicobacter_pylori_26695.fa`
- Reference accession: `NC_000915.1`
- Reads: `SRR12829274`
- Organism metadata: _Helicobacter pylori_ 26695
- Library strategy: `WGS`
- Layout: paired-end
- Platform: Illumina HiSeq 2500
- Subset size: `N=100000` spots

The full run contains 8,750,301 read spots and 2,625,090,300 bases. Using 100,000 paired-end spots gives about 30,000,000 raw bases before trimming. For a 1,667,867 bp genome, the expected raw coverage is about:

```text
30,000,000 / 1,667,867 = 18.0x
```

This is a practical compromise: it is larger than the earlier tiny BAM and should produce a more informative VCF, but it avoids downloading the full SRA run.

## Workflow

The `Makefile` performs these steps:

1. Download a subset of paired-end reads with `fastq-dump`.
2. Trim reads with `fastp`.
3. Align reads to the reference with `bwa mem`.
4. Sort and index the BAM with `samtools`.
5. Call variants with `bcftools mpileup | bcftools call --ploidy 1`.
6. Index the VCF with `tabix`.
7. Summarize coverage and variant counts.

Main output files:

```text
bam/H_pylori_26695_SRR12829274_N100000.sorted.bam
bam/H_pylori_26695_SRR12829274_N100000.sorted.bam.bai
vcf/H_pylori_26695_SRR12829274_N100000.vcf.gz
vcf/H_pylori_26695_SRR12829274_N100000.vcf.gz.tbi
reports/coverage_summary.tsv
reports/variant_summary.tsv
reports/high_confidence_summary.tsv
reports/H_pylori_26695_SRR12829274_N100000.bcftools_stats.txt
```

Large generated files are ignored by Git.

## Results

The generated VCF contains small variants relative to the _H. pylori_ 26695 reference. The summary files report the final coverage and variant counts:

```bash
cat reports/coverage_summary.tsv
cat reports/variant_summary.tsv
cat reports/high_confidence_summary.tsv
```

For this run, the mean depth was `16.3862x`, but coverage was not uniform: `46.6917%` of the genome had at least one read and `27.3779%` had at least `10x` depth. The raw VCF contained `199` calls: `141` SNVs and `58` indels. A simple higher-confidence subset using `QUAL>=30` and `DP>=10` contained `70` calls: `39` SNVs and `31` indels.

The high-confidence subset is useful for discussion, but it is not a formal truth set. Low-depth calls such as variants with `DP=1` are much more likely to be errors or unsupported calls. Calls with high depth, high quality, and alternate alleles supported by many reads are more plausible true variants. Indels in homopolymer runs should be interpreted cautiously because alignment artifacts are common there.

Variant calls are computational predictions. To decide which calls look real, I would inspect the VCF and BAM together in IGV and check whether the alternate allele is supported by multiple reads on both strands, has reasonable depth and quality, and is not concentrated at read ends or in low-quality/alignment-problem regions.

## IGV

Load these files in IGV:

1. Genome: `../week02/fasta/Helicobacter_pylori_26695.fa`
2. Alignment: `bam/H_pylori_26695_SRR12829274_N100000.sorted.bam`
3. Variants: `vcf/H_pylori_26695_SRR12829274_N100000.vcf.gz`
4. Annotation: `../week02/gff/Helicobacter_pylori_26695.gff.gz`

Screenshot to add after visual inspection:

```text
images/vcf_igv_view.png
```
