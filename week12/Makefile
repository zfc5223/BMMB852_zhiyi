# Six quantification methods on the same six single-end chr22 RNA-seq samples.
DATA_DIR ?= data
DATA_URL ?= https://data.biostarhandbook.com/data/rnaseq-tutorial-2026.tar.gz
THREADS ?= 4
FRAG_MEAN ?= 200
FRAG_SD ?= 20
READ_LENGTH ?= 100
PREPDE_FILE ?= src/vendor/prepDE.py3
PREPDE ?= python $(PREPDE_FILE)

GENOME = $(DATA_DIR)/fasta/chr22.genome.fa
TRANSCRIPTS = $(DATA_DIR)/fasta/chr22.transcripts.fa
GTF = $(DATA_DIR)/gtf/chr22.gtf
SAMPLES := $(shell awk -F, 'NR>1 {print $$1}' design.csv)
FASTQS = $(addprefix $(DATA_DIR)/fastq/,$(addsuffix _R1.fq,$(SAMPLES)))
TXMAP = refs/tx2gene.tsv
GENTROME = refs/salmon_gentrome.fa
DECOYS = refs/decoys.txt
TIME = env time -v -o

STAR_BAMS = $(addprefix results/star/,$(addsuffix /Aligned.sortedByCoord.out.bam,$(SAMPLES)))
STAR_TXBAMS = $(addprefix results/star/,$(addsuffix /Aligned.toTranscriptome.out.bam,$(SAMPLES)))
STAR_COUNTS = $(addprefix results/star/,$(addsuffix /ReadsPerGene.out.tab,$(SAMPLES)))
FC_COUNTS = $(addprefix results/featurecounts/,$(addsuffix .txt,$(SAMPLES)))
RSEM_COUNTS = $(addprefix results/rsem/,$(addsuffix .genes.results,$(SAMPLES)))
SALMON_COUNTS = $(addprefix results/salmon/,$(addsuffix /quant.sf,$(SAMPLES)))
KALLISTO_COUNTS = $(addprefix results/kallisto/,$(addsuffix /abundance.tsv,$(SAMPLES)))
HISAT_BAMS = $(addprefix results/hisat2/,$(addsuffix .bam,$(SAMPLES)))
STRINGTIE_GTFS = $(addprefix results/stringtie/,$(addsuffix .gtf,$(SAMPLES)))

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.ONESHELL:
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c

# Set DATA_DIR to an already extracted copy of the course input data if desired.
$(GENOME) $(TRANSCRIPTS) $(GTF) $(FASTQS) &:
	if [ "$(DATA_DIR)" != "data" ]; then \
		echo "Missing input in DATA_DIR=$(DATA_DIR); provide all six FASTQs, FASTA and GTF" >&2; exit 1; \
	fi
	if [ -e $(DATA_DIR)/fasta ] || [ -e $(DATA_DIR)/gtf ] || [ -e $(DATA_DIR)/fastq ]; then \
		echo "Partial input directory; inspect $(DATA_DIR) before retrying" >&2; exit 1; \
	fi
	mkdir -p $(DATA_DIR)
	curl --fail --location --retry 3 --output $(DATA_DIR)/tutorial.tar.gz.part $(DATA_URL)
	tar -tzf $(DATA_DIR)/tutorial.tar.gz.part >/dev/null
	tar -xzf $(DATA_DIR)/tutorial.tar.gz.part -C $(DATA_DIR)
	mv $(DATA_DIR)/tutorial.tar.gz.part $(DATA_DIR)/tutorial.tar.gz
	for path in $(GENOME) $(TRANSCRIPTS) $(GTF) $(FASTQS); do test -s "$$path"; done

data: $(GENOME) $(TRANSCRIPTS) $(GTF) $(FASTQS)

$(TXMAP) $(GENTROME) $(DECOYS) &: $(GENOME) $(TRANSCRIPTS) $(GTF) src/prepare_reference.py
	python src/prepare_reference.py --gtf $(GTF) --transcripts $(TRANSCRIPTS) --genome $(GENOME) --out refs

refs/star/SA: $(GENOME) $(GTF)
	mkdir -p refs/star logs/index
	$(TIME) logs/index/star.time STAR --runMode genomeGenerate --runThreadN $(THREADS) \
		--genomeDir refs/star --genomeFastaFiles $(GENOME) --sjdbGTFfile $(GTF) \
		--sjdbOverhang 99 --genomeSAindexNbases 11 >logs/index/star.log 2>&1

refs/rsem/ref.grp: $(GENOME) $(GTF)
	mkdir -p refs/rsem logs/index
	$(TIME) logs/index/rsem.time rsem-prepare-reference --gtf $(GTF) \
		$(GENOME) refs/rsem/ref >logs/index/rsem.log 2>&1

refs/salmon/.complete: $(GENTROME) $(DECOYS)
	mkdir -p refs/salmon logs/index
	$(TIME) logs/index/salmon.time salmon index -t $(GENTROME) -d $(DECOYS) --keepDuplicates \
		-i refs/salmon -p $(THREADS) >logs/index/salmon.log 2>&1
	touch $@

refs/kallisto/transcripts.idx: $(TRANSCRIPTS) $(TXMAP)
	mkdir -p refs/kallisto logs/index
	$(TIME) logs/index/kallisto.time kallisto index -i $@ $(TRANSCRIPTS) >logs/index/kallisto.log 2>&1

refs/hisat2/chr22.1.ht2: $(GENOME)
	mkdir -p refs/hisat2 logs/index
	$(TIME) logs/index/hisat2.time hisat2-build -p $(THREADS) $(GENOME) refs/hisat2/chr22 \
		>logs/index/hisat2.log 2>&1

refs/splice_sites.txt: $(GTF)
	mkdir -p refs
	hisat2_extract_splice_sites.py $(GTF) >$@

results/star/%/ReadsPerGene.out.tab results/star/%/Aligned.sortedByCoord.out.bam results/star/%/Aligned.toTranscriptome.out.bam &: $(DATA_DIR)/fastq/%_R1.fq refs/star/SA
	mkdir -p results/star/$* logs/star
	$(TIME) logs/star/$*.time STAR --runThreadN $(THREADS) --genomeDir refs/star \
		--readFilesIn $(DATA_DIR)/fastq/$*_R1.fq --outFileNamePrefix results/star/$*/ \
		--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts TranscriptomeSAM \
		>logs/star/$*.log 2>&1

results/featurecounts/%.txt: results/star/%/Aligned.sortedByCoord.out.bam $(GTF)
	mkdir -p results/featurecounts logs/featurecounts
	$(TIME) logs/featurecounts/$*.time featureCounts -T $(THREADS) -s 2 -t exon -g gene_id \
		-a $(GTF) -o $@ $< >logs/featurecounts/$*.log 2>&1

results/rsem/%.genes.results: results/star/%/Aligned.toTranscriptome.out.bam refs/rsem/ref.grp
	mkdir -p results/rsem logs/rsem
	$(TIME) logs/rsem/$*.time rsem-calculate-expression --alignments \
		--strandedness reverse --fragment-length-mean $(FRAG_MEAN) \
		--fragment-length-sd $(FRAG_SD) --num-threads $(THREADS) --no-bam-output \
		$< refs/rsem/ref results/rsem/$* >logs/rsem/$*.log 2>&1

results/salmon/%/quant.sf: $(DATA_DIR)/fastq/%_R1.fq refs/salmon/.complete
	mkdir -p results/salmon/$* logs/salmon
	$(TIME) logs/salmon/$*.time salmon quant -i refs/salmon -l SR -r $< \
		--validateMappings --incompatPrior 0 --fldMean $(FRAG_MEAN) --fldSD $(FRAG_SD) \
		-p $(THREADS) -o results/salmon/$* >logs/salmon/$*.log 2>&1

results/kallisto/%/abundance.tsv: $(DATA_DIR)/fastq/%_R1.fq refs/kallisto/transcripts.idx
	mkdir -p results/kallisto/$* logs/kallisto
	$(TIME) logs/kallisto/$*.time kallisto quant -i refs/kallisto/transcripts.idx \
		-o results/kallisto/$* --single --rf-stranded -l $(FRAG_MEAN) -s $(FRAG_SD) \
		-t $(THREADS) $< >logs/kallisto/$*.log 2>&1

results/hisat2/%.bam: $(DATA_DIR)/fastq/%_R1.fq refs/hisat2/chr22.1.ht2 refs/splice_sites.txt
	mkdir -p results/hisat2 logs/hisat2
	$(TIME) logs/hisat2/$*.time hisat2 -p $(THREADS) --dta --rna-strandness R \
		--known-splicesite-infile refs/splice_sites.txt -x refs/hisat2/chr22 -U $< \
		-S results/hisat2/$*.sam 2>logs/hisat2/$*.log
	$(TIME) logs/hisat2/$*.sort.time samtools sort -@ $(THREADS) -o $@ results/hisat2/$*.sam
	samtools index $@
	rm results/hisat2/$*.sam

results/stringtie/%.gtf: results/hisat2/%.bam $(GTF)
	mkdir -p results/stringtie logs/stringtie
	$(TIME) logs/stringtie/$*.time stringtie --rf -e -G $(GTF) -p $(THREADS) \
		-o $@ $< >logs/stringtie/$*.log 2>&1

results/matrices/star_featurecounts.gene.direct.tsv: $(FC_COUNTS) $(TXMAP) design.csv src/merge_counts.py
	python src/merge_counts.py featurecounts --design design.csv --mapping $(TXMAP) \
		--input results/featurecounts --output $@

results/matrices/star_genecounts.gene.direct.tsv: $(STAR_COUNTS) $(TXMAP) design.csv src/merge_counts.py
	python src/merge_counts.py star --design design.csv --mapping $(TXMAP) \
		--input results/star --star-column 4 --output $@

results/matrices/rsem.gene.expected.tsv results/matrices/rsem.transcript.expected.tsv &: $(RSEM_COUNTS) $(TXMAP) design.csv src/merge_counts.py
	python src/merge_counts.py rsem --design design.csv --mapping $(TXMAP) \
		--input results/rsem --level gene --output results/matrices/rsem.gene.expected.tsv
	python src/merge_counts.py rsem --design design.csv --mapping $(TXMAP) \
		--input results/rsem --level transcript --output results/matrices/rsem.transcript.expected.tsv

results/matrices/salmon.gene.expected.tsv results/matrices/salmon.transcript.expected.tsv &: $(SALMON_COUNTS) $(TXMAP) design.csv src/import_tximport.R
	mkdir -p results/matrices
	Rscript src/import_tximport.R salmon design.csv $(TXMAP) results/salmon results/matrices

results/matrices/kallisto.gene.expected.tsv results/matrices/kallisto.transcript.expected.tsv &: $(KALLISTO_COUNTS) $(TXMAP) design.csv src/import_tximport.R
	mkdir -p results/matrices
	Rscript src/import_tximport.R kallisto design.csv $(TXMAP) results/kallisto results/matrices

results/stringtie/prepde.gene.csv results/stringtie/prepde.transcript.csv &: $(STRINGTIE_GTFS) design.csv
	mkdir -p logs/stringtie
	test -s $(PREPDE_FILE) || { echo "Missing $(PREPDE_FILE); see README" >&2; exit 1; }
	awk -F, 'NR>1 {print $$1 " results/stringtie/" $$1 ".gtf"}' design.csv >results/stringtie/sample_list.txt
	$(TIME) logs/stringtie/prepde.time $(PREPDE) -i results/stringtie/sample_list.txt \
		-g results/stringtie/prepde.gene.csv -t results/stringtie/prepde.transcript.csv \
		-l $(READ_LENGTH) >logs/stringtie/prepde.log 2>&1

results/matrices/stringtie.gene.derived.tsv results/matrices/stringtie.transcript.derived.tsv &: results/stringtie/prepde.gene.csv results/stringtie/prepde.transcript.csv $(TXMAP) src/merge_counts.py
	python src/merge_counts.py prepde --design design.csv --mapping $(TXMAP) \
		--input results/stringtie/prepde.transcript.csv --level gene --output results/matrices/stringtie.gene.derived.tsv
	python src/merge_counts.py prepde --design design.csv --mapping $(TXMAP) \
		--input results/stringtie/prepde.transcript.csv --level transcript --output results/matrices/stringtie.transcript.derived.tsv

GENE_MATRICES = results/matrices/star_featurecounts.gene.direct.tsv \
	results/matrices/star_genecounts.gene.direct.tsv \
	results/matrices/rsem.gene.expected.tsv \
	results/matrices/salmon.gene.expected.tsv \
	results/matrices/kallisto.gene.expected.tsv \
	results/matrices/stringtie.gene.derived.tsv
TX_MATRICES = results/matrices/rsem.transcript.expected.tsv \
	results/matrices/salmon.transcript.expected.tsv \
	results/matrices/kallisto.transcript.expected.tsv \
	results/matrices/stringtie.transcript.derived.tsv

matrices: $(GENE_MATRICES) $(TX_MATRICES)

report/summary.tsv: $(GENE_MATRICES) $(TXMAP) src/compare.py
	python src/compare.py --matrices results/matrices --mapping $(TXMAP) \
		--design design.csv --logs logs --out report

report: report/summary.tsv
all: matrices report versions

versions: data
	mkdir -p report
	sha256sum $(GENOME) $(TRANSCRIPTS) $(GTF) >report/reference_sha256.txt
	sha256sum $(PREPDE_FILE) >report/prepde_sha256.txt
	{ \
		STAR --version; featureCounts -v; rsem-calculate-expression --version; \
		salmon --version; kallisto version; hisat2 --version | sed -n '1p'; \
		stringtie --version; samtools --version | sed -n '1p'; env time --version | sed -n '1p'; Rscript --version; \
		Rscript -e 'cat("tximport", as.character(packageVersion("tximport")), "\n")'; \
	} >report/software_versions.txt 2>&1

star_fc: results/matrices/star_featurecounts.gene.direct.tsv
star_gc: results/matrices/star_genecounts.gene.direct.tsv
star_rsem: results/matrices/rsem.gene.expected.tsv results/matrices/rsem.transcript.expected.tsv
salmon: results/matrices/salmon.gene.expected.tsv results/matrices/salmon.transcript.expected.tsv
kallisto: results/matrices/kallisto.gene.expected.tsv results/matrices/kallisto.transcript.expected.tsv
stringtie: results/matrices/stringtie.gene.derived.tsv results/matrices/stringtie.transcript.derived.tsv

.PHONY: all data matrices report versions star_fc star_gc star_rsem salmon kallisto stringtie
