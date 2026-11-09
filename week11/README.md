# RNA-seq coverage in HBR and UHR samples

I compared RNA-seq coverage on chromosome 22 for three Human Brain Reference
samples (`HBR_1`-`HBR_3`) and three Universal Human Reference samples
(`UHR_1`-`UHR_3`). HBR is brain RNA, while UHR is a mixture of RNA from ten
cancer cell lines.

## IGV setup

I loaded `fasta/chr22.genome.fa`, the six coverage tracks in `bw/`, and
`gtf/chr22.gtf`. I colored HBR pink and UHR blue. All six tracks have the same
linear data range within each image, so their heights can be compared.

## Genome-wide patterns and replicate consistency

Region: `chr22:19,668,645-24,668,645`  
Data range: `0-150`

![Five-megabase overview with a common data range](images/01_wide_common_scale.png)

At this scale, coverage appears in separate transcribed regions rather than
continuously across the chromosome. The three replicates in each group have
similar peak locations. HBR and UHR are not identical; one clear example is
the strong signal near 22.9 Mb in all three UHR samples, which is much weaker
in HBR.

The `0-150` view keeps the strongest signals in perspective, but it hides many
smaller peaks. I also show the same region capped at 15. Values above 15 are
clipped in this image.

Region: `chr22:19,668,645-24,668,645`  
Data range: `0-15` (clipped)

![Five-megabase overview emphasizing lower coverage](images/02_wide_low_range.png)

The lower range makes the weaker peaks easier to see. Their positions are
fairly consistent among the three replicates in each group, although the peak
heights are not exactly the same. Overall, the replicates look reasonably
consistent.

## A gene not expressed in these samples

Gene: `OR11H1`  
Region: `chr22:15,527,500-15,530,000`  
Data range: `0-10`

![OR11H1 has no visible coverage](images/03_or11h1_unexpressed.png)

I chose `OR11H1` as an unexpressed example. It is present in the GTF, but none
of the six tracks has visible coverage over the gene.

## A gene strongly expressed in all samples

Gene: `DDX17`  
Region: `chr22:38,480,000-38,510,000`  
Data range: `0-30` (clipped)

![DDX17 gene and annotated transcripts](images/04_ddx17_gene.png)

I chose `DDX17` because all six samples have clear coverage across many of its
exons. The GTF contains several overlapping `DDX17` transcripts. Most of the
covered exons are shared, so I cannot choose one complete transcript from
these short-read tracks. I could only exclude a transcript if an exon unique
to that transcript had no coverage. I capped this overview at 30 so the exon
pattern is easier to see.

Representative exons: `chr22:38,494,200-38,495,100`  
Data range: `0-110`

![Coverage over representative DDX17 exons](images/05_ddx17_exon.png)

Across these exons, typical coverage is about 16-18 in HBR and 19-32 in UHR,
with local peaks of about 60-110. This supports strong `DDX17` expression in
both sample types.

## A gene differentially expressed between HBR and UHR

Gene: `SULT4A1`  
Region: `chr22:43,820,000-43,870,000`  
Data range: `0-30`

![SULT4A1 coverage across the complete gene](images/06_sult4a1_gene.png)

`SULT4A1` was the clearest difference I found between the two groups. All
three HBR samples have coverage across several exons, while the three UHR
samples have little or no coverage in the same places. I interpret this as
`SULT4A1` being upregulated in HBR relative to UHR.

Representative exon: `chr22:43,824,300-43,826,300`  
Data range: `0-30`

![Coverage over a representative SULT4A1 exon](images/07_sult4a1_exon.png)

Coverage over this exon is about 13-14 in the HBR samples and close to zero in
UHR, giving a visual estimate of about a 100-fold difference. I treated this
as an approximate fold change because this assignment did not include a
statistical differential-expression test.

I am fairly confident that this is not a random fluctuation. The same pattern
appears in all three replicates and across several exons: every HBR sample has
clear coverage, while every UHR sample stays low.

