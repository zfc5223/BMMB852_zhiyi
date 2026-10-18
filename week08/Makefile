SHELL := /bin/bash
.DELETE_ON_ERROR:

N ?= 50000
REF ?= ../week02/fasta/Helicobacter_pylori_26695.fa
GFF ?= ../week02/gff/Helicobacter_pylori_26695.gff.gz
CSQ_GFF := annotation/Helicobacter_pylori_26695.csq.gff3

ACCESSIONS := SRR12829274 SRR8369179 SRR27606652
SAMPLES := $(foreach acc,$(ACCESSIONS),H_pylori_$(acc)_N$(N))
BAMS := $(foreach sample,$(SAMPLES),bam/$(sample).sorted.bam)
VCFS := $(foreach sample,$(SAMPLES),vcf/$(sample).vcf.gz)
MERGED := vcf/H_pylori_multisample_N$(N).vcf.gz
ANNOTATED := vcf/H_pylori_multisample_N$(N).annotated.vcf.gz

.PHONY: all check clean

all: check $(ANNOTATED).tbi reports/sample_summary.tsv reports/effect_summary.tsv

check:
	@test -s $(REF) || { echo "Missing reference FASTA: $(REF)"; exit 1; }
	@test -s $(GFF) || { echo "Missing GFF annotation: $(GFF)"; exit 1; }

reads/raw reads/trimmed bam vcf reports logs:
	mkdir -p $@

$(REF).bwt: $(REF)
	bwa index $(REF)

$(REF).fai: $(REF)
	samtools faidx $(REF)

annotation:
	mkdir -p $@

$(CSQ_GFF): $(GFF) | annotation
	@{ echo '##gff-version 3'; echo '##sequence-region NC_000915.1 1 1667867'; \
	 zcat $(GFF) | awk -F '\t' 'BEGIN{OFS="\t"} $$3=="CDS" {id="cds:" NR; tx="transcript:" NR; gene="gene:" NR; name="unknown"; n=split($$9,a,";"); for(i=1;i<=n;i++){if(a[i] ~ /^gene=/){sub(/^gene=/,"",a[i]); name=a[i]}}; print $$1,"bcftools","gene",$$4,$$5,".",$$7,".","ID=" gene ";Name=" name ";biotype=protein_coding"; print $$1,"bcftools","transcript",$$4,$$5,".",$$7,".","ID=" tx ";Parent=" gene ";Name=" name ";biotype=protein_coding"; print $$1,$$2,"CDS",$$4,$$5,$$6,$$7,$$8,"ID=" id ";Parent=" tx ";gene=" name}' ; } > $@

define SAMPLE_RULES
RAW_R1_$(1) = reads/raw/$(1)_R1.fastq
RAW_R2_$(1) = reads/raw/$(1)_R2.fastq
TRIM_R1_$(1) = reads/trimmed/$(1)_R1.trimmed.fastq
TRIM_R2_$(1) = reads/trimmed/$(1)_R2.trimmed.fastq
BAM_$(1) = bam/$(1).sorted.bam
VCF_$(1) = vcf/$(1).vcf.gz

reads/raw/$(1).done: | reads/raw logs
	fastq-dump --split-files --skip-technical -X $(N) -O reads/raw $(2) > logs/$(1).fastq-dump.log 2>&1
	@test -s reads/raw/$(2)_1.fastq || { echo "Missing R1 FASTQ from $(2)"; exit 1; }
	@test -s reads/raw/$(2)_2.fastq || { echo "Missing R2 FASTQ from $(2)"; exit 1; }
	mv reads/raw/$(2)_1.fastq $$(RAW_R1_$(1))
	mv reads/raw/$(2)_2.fastq $$(RAW_R2_$(1))
	touch $$@

$$(TRIM_R1_$(1)) $$(TRIM_R2_$(1)): reads/raw/$(1).done | reads/trimmed reports
	fastp -i $$(RAW_R1_$(1)) -I $$(RAW_R2_$(1)) -o $$(TRIM_R1_$(1)) -O $$(TRIM_R2_$(1)) \
		-h reports/$(1).fastp.html -j reports/$(1).fastp.json

$$(BAM_$(1)): $(REF).bwt $$(TRIM_R1_$(1)) $$(TRIM_R2_$(1)) | bam
	bwa mem $(REF) $$(TRIM_R1_$(1)) $$(TRIM_R2_$(1)) | samtools sort -o $$@

$$(BAM_$(1)).bai: $$(BAM_$(1))
	samtools index $$<

$$(VCF_$(1)): $$(BAM_$(1)).bai $(REF).fai | vcf
	bcftools mpileup -Ou -f $(REF) $$(BAM_$(1)) | bcftools call --ploidy 1 -mv -Oz -o $$@
	printf '%s\n' '$(1)' > logs/$(1).sample
	bcftools reheader -s logs/$(1).sample -o $$@.tmp $$@
	mv $$@.tmp $$@

$$(VCF_$(1)).tbi: $$(VCF_$(1))
	tabix -p vcf $$<
endef

$(foreach pair,$(join $(SAMPLES),$(addprefix |,$(ACCESSIONS))),$(eval $(call SAMPLE_RULES,$(word 1,$(subst |, ,$(pair))),$(word 2,$(subst |, ,$(pair))))))

$(MERGED).tbi: $(VCFS:%=%.tbi) | vcf
	bcftools merge -m all -Oz -o $(MERGED) $(VCFS)
	tabix -p vcf $(MERGED)

$(ANNOTATED): $(MERGED).tbi $(REF).fai $(CSQ_GFF) | vcf
	bcftools csq --no-version -f $(REF) -g $(CSQ_GFF) -Oz -o $@ $(MERGED)

$(ANNOTATED).tbi: $(ANNOTATED)
	tabix -p vcf $<

reports/sample_summary.tsv: $(ANNOTATED).tbi | reports
	@{ echo -e "sample\tcalled_sites"; for sample in $(SAMPLES); do count=$$(bcftools query -s $$sample -i 'GT!="."' -f '%CHROM\n' $(ANNOTATED) | wc -l); echo -e "$$sample\t$$count"; done; } > $@

reports/effect_summary.tsv: $(ANNOTATED).tbi | reports
	@{ echo -e "effect\tcount"; bcftools query -f '%INFO/BCSQ\n' $(ANNOTATED) | tr ',' '\n' | awk -F'|' '$$1 ~ /^(missense|synonymous|stop_gained|stop_lost|start_lost|frameshift|inframe_insertion|inframe_deletion|inframe_altering)$$/ {count[$$1]++} END {for (effect in count) print effect"\t"count[effect]}' | sort; } > $@

clean:
	rm -rf reads bam vcf logs reports
