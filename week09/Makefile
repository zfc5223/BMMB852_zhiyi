SHELL := /bin/bash
.DELETE_ON_ERROR:

N ?= 50000
JOBS ?= 3
THREADS ?= 1
DESIGN ?= design.csv
REF ?= ../week02/fasta/Helicobacter_pylori_26695.fa
PROJECT ?= H_pylori

SAMPLES := $(shell awk -F',' 'NR > 1 && $$2 != "" {print $$2}' $(DESIGN))
VCFS := $(foreach sample,$(SAMPLES),vcf/$(sample)_N$(N).vcf.gz)
MERGED := vcf/$(PROJECT)_multisample_N$(N).vcf.gz
SUMMARY := reports/$(PROJECT)_sample_summary_N$(N).tsv
STATS := reports/$(PROJECT)_merged_N$(N).stats

RUN ?=
SAMPLE ?=
PREFIX := $(SAMPLE)_N$(N)
RAW_R1 := reads/raw/$(PREFIX)_R1.fastq
RAW_R2 := reads/raw/$(PREFIX)_R2.fastq
RAW_DONE := reads/raw/$(PREFIX).done
TRIM_R1 := reads/trimmed/$(PREFIX)_R1.trimmed.fastq
TRIM_R2 := reads/trimmed/$(PREFIX)_R2.trimmed.fastq
BAM := bam/$(PREFIX).sorted.bam
VCF := vcf/$(PREFIX).vcf.gz

.PHONY: all check check-sample batch merge sample clean

all: check batch merge $(SUMMARY) $(STATS)

check:
	@test -s $(DESIGN) || { echo "Missing design file: $(DESIGN)"; exit 1; }
	@test "$$(head -n 1 $(DESIGN) | cut -d, -f1-2)" = "Run,Sample" || { echo "The first two design columns must be Run,Sample"; exit 1; }
	@test "$$(awk -F',' 'NR > 1 && (NF < 2 || $$1 == "" || $$2 == "") {bad=1} END {print bad+0}' $(DESIGN))" = "0" || { echo "Each design row needs a Run and Sample value"; exit 1; }
	@test "$$(awk -F',' 'NR > 1 {seen[$$2]++; if (seen[$$2] > 1) bad=1} END {print bad+0}' $(DESIGN))" = "0" || { echo "Sample names must be unique"; exit 1; }
	@test -s $(REF) || { echo "Missing reference FASTA: $(REF)"; exit 1; }
	@test -n "$(SAMPLES)" || { echo "Design file has no samples"; exit 1; }

check-sample:
	@test -n "$(RUN)" || { echo "RUN is required"; exit 1; }
	@test -n "$(SAMPLE)" || { echo "SAMPLE is required"; exit 1; }

reads/raw reads/trimmed bam vcf logs reports:
	mkdir -p $@

$(REF).bwt: $(REF)
	bwa index $(REF)

$(REF).fai: $(REF)
	samtools faidx $(REF)

$(RAW_DONE): | check-sample reads/raw logs
	fastq-dump --split-files --skip-technical -X $(N) -O reads/raw $(RUN) > logs/$(PREFIX).fastq-dump.log 2>&1
	@test -s reads/raw/$(RUN)_1.fastq || { echo "Missing R1 FASTQ from $(RUN)"; exit 1; }
	@test -s reads/raw/$(RUN)_2.fastq || { echo "Missing R2 FASTQ from $(RUN)"; exit 1; }
	mv reads/raw/$(RUN)_1.fastq $(RAW_R1)
	mv reads/raw/$(RUN)_2.fastq $(RAW_R2)
	touch $@

$(TRIM_R1) $(TRIM_R2): $(RAW_DONE) | reads/trimmed reports
	fastp -w $(THREADS) -i $(RAW_R1) -I $(RAW_R2) -o $(TRIM_R1) -O $(TRIM_R2) \
		-h reports/$(PREFIX).fastp.html -j reports/$(PREFIX).fastp.json

$(BAM): $(REF).bwt $(TRIM_R1) $(TRIM_R2) | bam
	bwa mem -t $(THREADS) $(REF) $(TRIM_R1) $(TRIM_R2) | samtools sort -@ $(THREADS) -o $@

$(BAM).bai: $(BAM)
	samtools index $<

$(VCF): $(BAM).bai $(REF).fai | vcf
	bcftools mpileup -Ou -f $(REF) $(BAM) | bcftools call --ploidy 1 -mv -Oz -o $@
	printf '%s\n' '$(SAMPLE)' > logs/$(PREFIX).sample
	bcftools reheader -s logs/$(PREFIX).sample -o $@.tmp $@
	mv $@.tmp $@

$(VCF).tbi: $(VCF)
	tabix -p vcf $<

sample: check-sample $(VCF).tbi

batch: check | logs
	@cat $(DESIGN) | parallel --colsep , --header : --jobs $(JOBS) --line-buffer \
		--joblog logs/parallel.joblog \
		'$(MAKE) --no-print-directory sample RUN={Run} SAMPLE={Sample} N=$(N) THREADS=$(THREADS)'

$(MERGED).tbi: batch | vcf
	@for vcf in $(VCFS); do test -s "$$vcf" && test -s "$$vcf.tbi" || { echo "Missing sample VCF: $$vcf"; exit 1; }; done
	@if test -s $(MERGED) && test -s $@; then \
		stale=0; \
		for vcf in $(VCFS); do \
			test "$$vcf" -nt $(MERGED) && stale=1; \
			test "$$vcf.tbi" -nt $(MERGED) && stale=1; \
		done; \
		if test "$$stale" -eq 0; then exit 0; fi; \
	fi; \
	bcftools merge --force-single -m all -Oz -o $(MERGED) $(VCFS) && \
	tabix -p vcf $(MERGED)

merge: $(MERGED).tbi

$(SUMMARY): $(MERGED).tbi | reports
	@{ echo -e "sample\tcalled_sites"; for sample in $(SAMPLES); do count=$$(bcftools query -s $$sample -i 'GT!="."' -f '%CHROM\n' $(MERGED) | wc -l); echo -e "$$sample\t$$count"; done; } > $@

$(STATS): $(MERGED).tbi | reports
	bcftools stats $(MERGED) > $@

clean:
	rm -rf reads bam vcf logs reports
