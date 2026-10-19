# Week 8: Predict Variant Effects

Make a multisample VCF from three _Helicobacter pylori_ isolates and explain
the effect of coding variants in that VCF.

The point of the multisample VCF is comparison: at the same reference
coordinate, I can see which isolates carry the alternate allele.

## Samples

I used three independent isolates from the same _H. pylori_ WGS study:

| Sample                         | SRA run        |
| ------------------------------ | -------------- |
| `H_pylori_SRR8369177_N50000` | `SRR8369177` |
| `H_pylori_SRR8369178_N50000` | `SRR8369178` |
| `H_pylori_SRR8369179_N50000` | `SRR8369179` |

The three runs are listed in the same _H. pylori_ isolate study:
[https://pmc.ncbi.nlm.nih.gov/articles/PMC7545276/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7545276/).

All reads are aligned to the same _H. pylori_ 26695 reference, so the coordinates are comparable.

Here I used 50,000 paired-end spots per isolate. This produced about 72,500--73,600 called sites per sample;
these are reference-relative calls, not measurements of biological function.

## Software choice

The tools have separate jobs:

- `bcftools mpileup` and `bcftools call --ploidy 1` call variants from each
  BAM. _H. pylori_ is haploid, so the haploid calling model is appropriate.
- `bcftools merge` combines the three single-sample VCFs while preserving the
  sample genotypes. This is the multisample step taught in the course.
- `SnpEff` predicts coding effects. The course material specifically uses
  `SnpEff`, and its `ANN` output reports the gene, consequence, coding change,
  and protein change in one record.

The Makefile builds a small custom `HP26695` SnpEff database from the same 26695 FASTA and NCBI GFF used in the
alignment and in IGV. This keeps the reference, annotation, and effect
prediction consistent for _H. pylori_.

## Reproduce

```bash
pixi run -m "$HOME/edu/bioinfo" make -C week08 all
```

The Makefile downloads the read subsets, trims and aligns them, creates one
BAM and one VCF per isolate, merges the VCFs, builds the local SnpEff database,
and annotates the merged VCF.

Important outputs:

```text
vcf/H_pylori_multisample_N50000.vcf.gz
vcf/H_pylori_multisample_N50000.annotated.vcf.gz
reports/sample_summary.tsv
reports/effect_summary.tsv
```

The annotated VCF has these sample columns:

```text
H_pylori_SRR8369177_N50000
H_pylori_SRR8369178_N50000
H_pylori_SRR8369179_N50000
```

## Variant examples

SnpEff writes effects in the VCF `ANN` field. The two examples below show the
comparison use of the multisample VCF:

| Coordinate            | Gene     | Coding codon   | Protein change        | Isolates with ALT              | Effect               |
| --------------------- | -------- | -------------- | --------------------- | ------------------------------ | -------------------- |
| `NC_000915.1:21861` | `eptA` | `AAC -> AGC` | Asn286Ser (`N286S`) | `SRR8369177`, `SRR8369178` | `missense_variant` |
| `NC_000915.1:78906` | `glmM` | `AAG -> CAG` | Lys401Gln (`K401Q`) | `SRR8369178`, `SRR8369179` | `missense_variant` |

The first record is genomic `T>C`. Because `eptA` is on the reverse strand,
the coding-strand change is `AAC -> AGC`; SnpEff reports `c.857A>G` and
`p.Asn286Ser`. The genotypes are `1/1, 1/1, .` in the sample order above, so
the alternate allele is called in the first two isolates. SnpEff also flags the
`eptA` feature as `WARNING_TRANSCRIPT_NO_START_CODON`; the example therefore
focuses on the internal amino-acid change.

## IGV evidence

I visualized the `eptA` variant in IGV using the 26695 reference FASTA, the
`SRR8369177` BAM, the final annotated multisample VCF, and the matching GFF:

- `../week02/fasta/Helicobacter_pylori_26695.fa`
- `bam/H_pylori_SRR8369177_N50000.sorted.bam`
- `vcf/H_pylori_multisample_N50000.annotated.vcf.gz`
- `../week02/gff/Helicobacter_pylori_26695.gff.gz`

The view covers `NC_000915.1:21830-21890` and shows the variant marker at
position 21,861, the BAM coverage and reads, the reference sequence, and the
translated coding feature. The `eptA` feature is on the reverse strand, so the
coding-strand change is `AAC -> AGC` even though the genomic change is `T>C`.
The variant popup identifies the same site as a SNP with reference `T` and
alternate `C`, and displays the SnpEff `ANN` annotation.

![IGV view of the eptA variant](images/effect_igv_view.png)

The screenshot is the alignment evidence. The matching software result is:

```text
NC_000915.1  21861  T  C  DP=30
ANN: missense_variant | MODERATE | eptA | c.857A>G | p.Asn286Ser
GT:  1/1  1/1  .
```
