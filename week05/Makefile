SHELL := /bin/bash
.DELETE_ON_ERROR:

SAMPLE ?= H_pylori_26695_SRR031127
REF ?= ../week02/fasta/Helicobacter_pylori_26695.fa
FASTQ ?= ../week04/reads/trimmed/$(SAMPLE).trimmed.fastq

BAM = bam/$(SAMPLE).sorted.bam

.PHONY: all check clean

all: check $(BAM).bai reports/$(SAMPLE).flagstat.txt reports/$(SAMPLE).idxstats.txt reports/$(SAMPLE).coverage.tsv reports/$(SAMPLE).coverage_summary.tsv

check:
	@test -s $(REF) || { echo "Missing reference FASTA: $(REF)"; exit 1; }
	@test -s $(FASTQ) || { echo "Missing FASTQ: $(FASTQ). Run week04 first."; exit 1; }

bam reports:
	mkdir -p $@

$(REF).bwt: $(REF)
	bwa index $(REF)

$(BAM): $(REF).bwt $(FASTQ) | bam
	bwa mem $(REF) $(FASTQ) | samtools sort -o $@

$(BAM).bai: $(BAM)
	samtools index $<

reports/$(SAMPLE).flagstat.txt: $(BAM) | reports
	samtools flagstat $< > $@

reports/$(SAMPLE).idxstats.txt: $(BAM) | reports
	samtools idxstats $< > $@

reports/$(SAMPLE).coverage.tsv: $(BAM) | reports
	samtools depth -a $< > $@

reports/$(SAMPLE).coverage_summary.tsv: reports/$(SAMPLE).coverage.tsv $(REF).fai | reports
	@awk 'BEGIN { OFS="\t"; print "metric","value" } FNR==NR { genome += $$2; next } { covered += ($$3 > 0); depth += $$3 } END { print "genome_bp", genome; print "covered_bp", covered; print "covered_fraction", covered / genome; print "mean_depth", depth / genome }' $(REF).fai $< > $@

clean:
	rm -rf bam reports
