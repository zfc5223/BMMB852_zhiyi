SHELL := /bin/bash
.DELETE_ON_ERROR:

ACC ?= SRR031127
N ?= 1000
SAMPLE ?= H_pylori_26695_$(ACC)

ENA_URL = https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq(85962)&fields=run_accession,study_accession,experiment_accession,scientific_name,instrument_platform,instrument_model,library_strategy,library_layout,read_count,base_count&format=tsv&limit=0

META = metadata/ena_runs.tsv
SUMMARY = reports/ena_summary.tsv

.PHONY: all metadata fastq qc trim clean

all: metadata fastq qc trim

metadata: $(SUMMARY)

fastq: reads/raw/$(SAMPLE).done

qc: qc/raw/$(SAMPLE).done

trim: reads/trimmed/$(SAMPLE).done qc/trimmed/$(SAMPLE).done

metadata/.dir reports/.dir reads/raw/.dir reads/trimmed/.dir qc/raw/.dir qc/trimmed/.dir logs/.dir:
	mkdir -p $(@D)
	touch $@

$(META): | metadata/.dir
	curl -L --fail "$(ENA_URL)" -o $@

$(SUMMARY): $(META) | reports/.dir
	@echo -e "category\tvalue\truns" > $@
	@awk -F '\t' 'NR > 1 { platform[$$5]++; strategy[$$7]++; layout[$$8]++; runs++ } END { print "all\truns\t" runs; for (k in platform) print "platform\t" k "\t" platform[k]; for (k in strategy) print "strategy\t" k "\t" strategy[k]; for (k in layout) print "layout\t" k "\t" layout[k] }' $< | sort >> $@
	rm -f $@.tmp

reads/raw/$(SAMPLE).done: | reads/raw/.dir logs/.dir
	fastq-dump --split-3 --skip-technical -X $(N) -O reads/raw $(ACC) > logs/$(SAMPLE).fastq-dump.log 2>&1
	@if [ -s reads/raw/$(ACC)_1.fastq ]; then mv reads/raw/$(ACC)_1.fastq reads/raw/$(SAMPLE)_R1.fastq; fi
	@if [ -s reads/raw/$(ACC)_2.fastq ]; then mv reads/raw/$(ACC)_2.fastq reads/raw/$(SAMPLE)_R2.fastq; fi
	@if [ -s reads/raw/$(ACC).fastq ]; then mv reads/raw/$(ACC).fastq reads/raw/$(SAMPLE).fastq; fi
	touch $@

qc/raw/$(SAMPLE).done: reads/raw/$(SAMPLE).done | qc/raw/.dir
	fastqc -q -o qc/raw reads/raw/$(SAMPLE)*.fastq
	touch $@

reads/trimmed/$(SAMPLE).done: reads/raw/$(SAMPLE).done | reads/trimmed/.dir reports/.dir
	@if [ -s reads/raw/$(SAMPLE)_R1.fastq ] && [ -s reads/raw/$(SAMPLE)_R2.fastq ]; then \
		fastp -i reads/raw/$(SAMPLE)_R1.fastq -I reads/raw/$(SAMPLE)_R2.fastq \
			-o reads/trimmed/$(SAMPLE)_R1.trimmed.fastq -O reads/trimmed/$(SAMPLE)_R2.trimmed.fastq \
			-h reports/$(SAMPLE).fastp.html -j reports/$(SAMPLE).fastp.json; \
	else \
		fastp -i reads/raw/$(SAMPLE).fastq -o reads/trimmed/$(SAMPLE).trimmed.fastq \
			-h reports/$(SAMPLE).fastp.html -j reports/$(SAMPLE).fastp.json; \
	fi
	touch $@

qc/trimmed/$(SAMPLE).done: reads/trimmed/$(SAMPLE).done | qc/trimmed/.dir
	fastqc -q -o qc/trimmed reads/trimmed/$(SAMPLE)*.fastq
	touch $@

clean:
	rm -rf metadata reports reads qc logs
