# RNA-seq quantification with six methods

I quantified the same six human RNA-seq samples used in [my IGV analysis](../week11/README.md): three Human Brain Reference (HBR) and three Universal Human Reference (UHR) samples. All runs used the same single-end FASTQs, chromosome 22 genome, transcript FASTA, and GTF. The [sample sheet](design.csv) fixes the order in every matrix.

| Workflow | How counts are obtained |
|---|---|
| STAR → featureCounts | Count reads assigned to genes after genome alignment. |
| STAR GeneCounts | Count genes during STAR alignment. |
| STAR → RSEM | Estimate transcript and gene counts from transcriptome alignments. |
| Salmon → tximport | Estimate transcript counts by selective alignment; sum them by gene. |
| kallisto → tximport | Estimate transcript counts by pseudoalignment; sum them by gene. |
| HISAT2 → StringTie → prepDE | Estimate reference transcript coverage and derive count-like values. |

The workflow produced six gene matrices and four transcript matrices in `results/matrices/`. Direct, estimated, and coverage-derived counts have different meanings, so their filenames keep those types separate.

## Results

The two direct gene-count methods agreed exactly. Across 1,336 genes, RSEM, Salmon, and kallisto correlated strongly with featureCounts (mean Pearson `r` = 0.93–0.97 on log1p chr22 CPM); StringTie was lower (`r` = 0.79). For `SULT4A1`, featureCounts found 286–335 reads in HBR and 3–6 in UHR, matching the [IGV view](../week11/images/06_sult4a1_gene.png). RSEM, Salmon, and kallisto gave similar gene totals but different transcript estimates. The [comparison report](report/comparison.md) has the figures and full results.

## Reproduce

From `week12/` on Linux with Pixi installed:

```bash
pixi install
mkdir -p src/vendor
curl -fL https://raw.githubusercontent.com/gpertea/stringtie/v2.2.3/prepDE.py3 -o src/vendor/prepDE.py3
sha256sum -c report/prepde_sha256.txt
pixi run make -j2 all THREADS=4
```

The last command downloads the course RNA-seq archive into `data/` if needed. To reuse the existing files instead, add `DATA_DIR=../week11` to that command. `-j2` runs up to two jobs at once; `THREADS=4` sets threads per tool. The library is reverse-stranded. For single-end quantification, I assumed a fragment length of 200 ± 20 bp; change `FRAG_MEAN` and `FRAG_SD` if measured values are available. Software versions are pinned in [pixi.lock](pixi.lock), and the reference checksums are in [report/reference_sha256.txt](report/reference_sha256.txt).

I used AI to draft the scripts, then checked the sample order, strand setting, and transcript-to-gene mapping. During validation I found that Salmon dropped five duplicate-sequence transcripts by default and that prepDE's gene table misgrouped some transcripts. The workflow keeps the duplicates and builds the StringTie gene matrix by summing prepDE transcript counts with the checked GTF mapping.
