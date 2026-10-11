SHELL := /bin/bash
.DELETE_ON_ERROR:

ACC ?= SRR12829274
N ?= 100000
SAMPLE ?= H_pylori_26695_$(ACC)_N$(N)
REF ?= ../week02/fasta/Helicobacter_pylori_26695.fa

RAW_R1 = reads/raw/$(SAMPLE)_R1.fastq
RAW_R2 = reads/raw/$(SAMPLE)_R2.fastq
TRIM_R1 = reads/trimmed/$(SAMPLE)_R1.trimmed.fastq
TRIM_R2 = reads/trimmed/$(SAMPLE)_R2.trimmed.fastq
BAM = bam/$(SAMPLE).sorted.bam
VCF = vcf/$(SAMPLE).vcf.gz

.PHONY: all check clean

all: check $(BAM).bai $(VCF).tbi reports/$(SAMPLE).flagstat.txt reports/coverage_summary.tsv reports/variant_summary.tsv reports/high_confidence_summary.tsv reports/$(SAMPLE).bcftools_stats.txt

check:
	@test -s $(REF) || { echo "Missing reference FASTA: $(REF). Run week02 first."; exit 1; }

reads/raw reads/trimmed bam vcf reports logs:
	mkdir -p $@

$(REF).bwt: $(REF)
	bwa index $(REF)

$(REF).fai: $(REF)
	samtools faidx $(REF)

reads/raw/$(SAMPLE).done: | reads/raw logs
	fastq-dump --split-files --skip-technical -X $(N) -O reads/raw $(ACC) > logs/$(SAMPLE).fastq-dump.log 2>&1
	@test -s reads/raw/$(ACC)_1.fastq || { echo "Missing R1 FASTQ from $(ACC)"; exit 1; }
	@test -s reads/raw/$(ACC)_2.fastq || { echo "Missing R2 FASTQ from $(ACC)"; exit 1; }
	mv reads/raw/$(ACC)_1.fastq $(RAW_R1)
	mv reads/raw/$(ACC)_2.fastq $(RAW_R2)
	touch $@

$(TRIM_R1) $(TRIM_R2): reads/raw/$(SAMPLE).done | reads/trimmed reports
	fastp -i $(RAW_R1) -I $(RAW_R2) -o $(TRIM_R1) -O $(TRIM_R2) \
		-h reports/$(SAMPLE).fastp.html -j reports/$(SAMPLE).fastp.json

$(BAM): $(REF).bwt $(TRIM_R1) $(TRIM_R2) | bam
	bwa mem $(REF) $(TRIM_R1) $(TRIM_R2) | samtools sort -o $@

$(BAM).bai: $(BAM)
	samtools index $<

$(VCF): $(BAM).bai $(REF).fai | vcf
	bcftools mpileup -Ou -f $(REF) $(BAM) | bcftools call --ploidy 1 -mv -Oz -o $@

$(VCF).tbi: $(VCF)
	tabix -p vcf $<

reports/coverage.tsv: $(BAM).bai | reports
	samtools depth -a $(BAM) > $@

reports/$(SAMPLE).flagstat.txt: $(BAM).bai | reports
	samtools flagstat $(BAM) > $@

reports/coverage_summary.tsv: reports/coverage.tsv $(REF).fai | reports
	@awk 'BEGIN { OFS="\t"; print "metric","value" } FNR==NR { genome += $$2; next } { covered += ($$3 > 0); depth += $$3; if ($$3 >= 10) covered10++ } END { print "genome_bp", genome; print "covered_bp", covered; print "covered_fraction", covered / genome; print "covered_10x_bp", covered10; print "covered_10x_fraction", covered10 / genome; print "mean_depth", depth / genome }' $(REF).fai $< > $@

reports/variant_summary.tsv: $(VCF).tbi | reports
	@echo -e "metric\tvalue" > $@
	@bcftools view -H $(VCF) | awk 'BEGIN { snv=0; indel=0; other=0 } { if (length($$4)==1 && length($$5)==1) snv++; else if (length($$4)!=length($$5)) indel++; else other++ } END { print "total_variants\t" snv + indel + other; print "snv\t" snv; print "indel\t" indel; print "other\t" other }' >> $@

reports/high_confidence_summary.tsv: $(VCF).tbi | reports
	@echo -e "metric\tvalue" > $@
	@bcftools view -H -i 'QUAL>=30 && INFO/DP>=10' $(VCF) | awk 'BEGIN { snv=0; indel=0; other=0 } { if (length($$4)==1 && length($$5)==1) snv++; else if (length($$4)!=length($$5)) indel++; else other++ } END { print "filter\tQUAL>=30 and DP>=10"; print "total_variants\t" snv + indel + other; print "snv\t" snv; print "indel\t" indel; print "other\t" other }' >> $@

reports/$(SAMPLE).bcftools_stats.txt: $(VCF).tbi | reports
	bcftools stats $(VCF) > $@

clean:
	rm -rf reads bam vcf logs reports
