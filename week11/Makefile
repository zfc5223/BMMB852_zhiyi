# RNA-Seq analysis with HISAT2 (course tutorial workflow)

DATA_URL = http://data.biostarhandbook.com/data/rnaseq-tutorial-2026.tar.gz
FASTA = fasta/chr22.genome.fa
FASTA_IDX = $(FASTA).1.ht2
GTF = gtf/chr22.gtf
SAMPLE ?= HBR_1
R1 = fastq/$(SAMPLE)_R1.fq
BAM = bam/$(SAMPLE).bam
BW = bw/$(SAMPLE).bw
COUNT = counts/$(SAMPLE).txt
DESIGN = design.csv
PARALLEL = parallel --eta --colsep , --header : --jobs 2 --lb

$(FASTA) $(GTF) $(R1) &:
	curl -sSL $(DATA_URL) | tar xzvf -

data: $(FASTA) $(GTF) $(R1)
	ls -lh $(FASTA) $(GTF) $(R1)

$(FASTA_IDX): $(FASTA)
	hisat2-build $(FASTA) $(FASTA)

index: $(FASTA_IDX)
	ls -lh $(FASTA_IDX)

$(BAM): $(FASTA_IDX) $(R1)
	mkdir -p $(dir $(BAM))
	hisat2 --threads 4 --sensitive -x $(FASTA) -U $(R1) | \
		samtools sort --write-index -o $(BAM)

$(BW): $(BAM)
	mkdir -p $(dir $(BW))
	bamCoverage --binSize 1 -b $(BAM) -o $(BW)

align: $(BAM) $(BW)
	ls -lh $(BAM) $(BW)

$(COUNT): $(BW)
	mkdir -p $(dir $(COUNT))
	featureCounts -a $(GTF) -s 2 -o $(COUNT) $(BAM)

counts: $(COUNT)
	ls -lh $(COUNT)

all: counts

batch: index
	cat $(DESIGN) | $(PARALLEL) $(MAKE) SAMPLE={sample} all

.PHONY: all data index align counts batch

SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory
