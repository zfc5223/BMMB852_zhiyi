SHELL := /bin/bash
.DELETE_ON_ERROR:

DESIGN ?= design.csv
THREADS ?= 8
JOBS ?= 2
MAX_SPOTS ?= 0
CHIP_EXTSIZE ?= 147
SUFFIX := $(if $(filter 0,$(MAX_SPOTS)),,_N$(MAX_SPOTS))

REF_URL ?= https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
REF ?= reference/hg38.fa
BT2_PREFIX ?= reference/hg38
BT2_INDEX := $(BT2_PREFIX).1.bt2

SAMPLES := $(shell awk -F, 'NR > 1 {print $$1}' $(DESIGN))
PEAK_SAMPLES := $(shell awk -F, 'NR > 1 && $$4 != "Input" {print $$1}' $(DESIGN))

SAMPLE ?=
RUNS ?=
ASSAY ?=
GROUP ?=
REPLICATE ?=
LAYOUT ?=
CONTROL ?=

READ1 := reads/$(SAMPLE)$(SUFFIX)_R1.fastq.gz
READ2 := reads/$(SAMPLE)$(SUFFIX)_R2.fastq.gz
READS := reads/$(SAMPLE)$(SUFFIX).fastq.gz
DOWNLOAD_DONE := reads/$(SAMPLE)$(SUFFIX).download.done
QC_DONE := qc/$(SAMPLE)$(SUFFIX).done
BAM := bam/$(SAMPLE)$(SUFFIX).sorted.bam
BAI := $(BAM).bai
CONTROL_BAM := bam/$(CONTROL)$(SUFFIX).sorted.bam
BIGWIG := coverage/$(SAMPLE)$(SUFFIX).CPM.bw
PEAK := peaks/$(SAMPLE)$(SUFFIX)_peaks.narrowPeak
PEAK_DEPS := $(BAM) $(if $(filter pSTAT3,$(ASSAY)),$(CONTROL_BAM))
PEAK_COUNTS := results/peak_counts$(SUFFIX).tsv
REPLICATE_OVERLAP := results/replicate_overlap$(SUFFIX).tsv
ACCESSIBILITY_CLASSES := results/accessibility_classes$(SUFFIX).tsv

.PHONY: all check check-sample genome samples sample call-peaks peak summaries clean

all: check genome samples call-peaks summaries

check:
	@test -s $(DESIGN) || { echo "Missing design file: $(DESIGN)"; exit 1; }
	@test "$$(head -n 1 $(DESIGN))" = "Sample,GSM,Runs,Assay,Group,Replicate,Layout,Control" || { echo "Unexpected design header"; exit 1; }
	@test "$$(awk -F, 'NR > 1 {seen[$$1]++; if (seen[$$1] > 1) bad=1} END {print bad+0}' $(DESIGN))" = 0 || { echo "Sample names must be unique"; exit 1; }

check-sample:
	@test -n "$(SAMPLE)" || { echo "SAMPLE is required"; exit 1; }
	@test -n "$(RUNS)" || { echo "RUNS is required"; exit 1; }
	@test "$(LAYOUT)" = "PAIRED" -o "$(LAYOUT)" = "SINGLE" || { echo "LAYOUT must be PAIRED or SINGLE"; exit 1; }

reference/ reads/ reads/runs/ qc/ bam/ coverage/ peaks/ logs/ results/ figures/:
	mkdir -p $@

$(REF): | reference/
	wget -q --show-progress -O $@.gz $(REF_URL)
	gzip -dc $@.gz > $@
	rm $@.gz

$(BT2_INDEX): $(REF)
	bowtie2-build --threads $(THREADS) $(REF) $(BT2_PREFIX)

genome: $(BT2_INDEX)

$(DOWNLOAD_DONE): $(DESIGN) | check-sample reads/ reads/runs/ logs/
	@set -euo pipefail; \
	run_list="$${RUNS//;/ }"; \
	for run in $$run_list; do \
		rm -f reads/runs/$${run}.fastq reads/runs/$${run}_1.fastq reads/runs/$${run}_2.fastq; \
		if test "$(MAX_SPOTS)" = 0; then \
			fasterq-dump --threads $(THREADS) --outdir reads/runs $$run \
				> logs/$$run.fasterq.log 2>&1; \
		else \
			fastq-dump --split-files --skip-technical -X $(MAX_SPOTS) \
				-O reads/runs $$run > logs/$$run.fastq-dump.log 2>&1; \
		fi; \
		if test "$(LAYOUT)" = SINGLE && test ! -s reads/runs/$${run}.fastq \
			&& test -s reads/runs/$${run}_1.fastq; then \
			mv reads/runs/$${run}_1.fastq reads/runs/$${run}.fastq; \
		fi; \
	done; \
	if test "$(LAYOUT)" = PAIRED; then \
		for run in $$run_list; do test -s reads/runs/$${run}_1.fastq; done; \
		for run in $$run_list; do cat reads/runs/$${run}_1.fastq; done | pigz -p $(THREADS) > $(READ1); \
		for run in $$run_list; do cat reads/runs/$${run}_2.fastq; done | pigz -p $(THREADS) > $(READ2); \
	else \
		for run in $$run_list; do test -s reads/runs/$${run}.fastq; done; \
		for run in $$run_list; do cat reads/runs/$${run}.fastq; done | pigz -p $(THREADS) > $(READS); \
	fi; \
	for run in $$run_list; do rm -f reads/runs/$${run}*.fastq; done; \
	touch $@

$(QC_DONE): $(DOWNLOAD_DONE) | qc/
	@if test "$(LAYOUT)" = PAIRED; then \
		fastqc --threads $(THREADS) --outdir qc $(READ1) $(READ2); \
	else \
		fastqc --threads $(THREADS) --outdir qc $(READS); \
	fi
	touch $@

$(BAM): $(DOWNLOAD_DONE) $(BT2_INDEX) | bam/
	@if test "$(LAYOUT)" = PAIRED; then \
		bowtie2 --very-sensitive --no-mixed --no-discordant -X 2000 \
			-p $(THREADS) -x $(BT2_PREFIX) -1 $(READ1) -2 $(READ2) | \
			samtools view -@ $(THREADS) -b -q 30 -f 2 -F 3844 | \
			samtools sort -@ $(THREADS) -o $@; \
	else \
		bowtie2 --very-sensitive -p $(THREADS) -x $(BT2_PREFIX) -U $(READS) | \
			samtools view -@ $(THREADS) -b -q 30 -F 3844 | \
			samtools sort -@ $(THREADS) -o $@; \
	fi

$(BAI): $(BAM)
	samtools index -@ $(THREADS) $<

$(BIGWIG): $(BAI) | coverage/
	bamCoverage --bam $(BAM) --outFileName $@ --outFileFormat bigwig \
		--normalizeUsing CPM --binSize 10 --skipNonCoveredRegions \
		--numberOfProcessors $(THREADS)

sample: check-sample $(QC_DONE) $(BAI) $(BIGWIG)

samples: check genome | logs/
	@awk -F, 'NR > 1 && $$1 != ""' $(DESIGN) | parallel --colsep , --jobs $(JOBS) --line-buffer \
		--joblog logs/sample_jobs$(SUFFIX).tsv \
		'$(MAKE) --no-print-directory sample SAMPLE={1} RUNS={3} ASSAY={4} GROUP={5} REPLICATE={6} LAYOUT={7} CONTROL={8} THREADS=$(THREADS) MAX_SPOTS=$(MAX_SPOTS)'

$(PEAK): $(PEAK_DEPS) $(BAI) | peaks/
	@if test "$(ASSAY)" = ATAC; then \
		macs3 callpeak -t $(BAM) -f BAMPE -g hs -n $(SAMPLE)$(SUFFIX) \
			--outdir peaks --keep-dup auto -q 0.01; \
	elif test "$(ASSAY)" = pSTAT3; then \
		test -s $(CONTROL_BAM) || { echo "Missing control BAM: $(CONTROL_BAM)"; exit 1; }; \
		macs3 callpeak -t $(BAM) -c $(CONTROL_BAM) -f BAM -g hs \
			-n $(SAMPLE)$(SUFFIX) --outdir peaks --keep-dup 1 \
			--nomodel --extsize $(CHIP_EXTSIZE) -q 0.01; \
	else \
		echo "Peak calling is not defined for assay $(ASSAY)"; exit 1; \
	fi

peak: check-sample $(PEAK)

call-peaks: samples | logs/
	@awk -F, 'NR > 1 && $$4 != "Input"' $(DESIGN) | \
	parallel --colsep , --jobs $(JOBS) --line-buffer --joblog logs/peak_jobs$(SUFFIX).tsv \
		'$(MAKE) --no-print-directory peak SAMPLE={1} RUNS={3} ASSAY={4} GROUP={5} REPLICATE={6} LAYOUT={7} CONTROL={8} THREADS=$(THREADS) MAX_SPOTS=$(MAX_SPOTS)'

$(PEAK_COUNTS) $(REPLICATE_OVERLAP) $(ACCESSIBILITY_CLASSES) &: call-peaks scripts/summarize.py $(DESIGN) | results/
	python scripts/summarize.py $(DESIGN) peaks results '$(SUFFIX)'

summaries: $(PEAK_COUNTS) $(REPLICATE_OVERLAP) $(ACCESSIBILITY_CLASSES)

clean:
	rm -rf reads qc bam coverage peaks logs results reference

