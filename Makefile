OUTPUT_DIR  ?= /mnt/treasure/asc25/jiazhaopeng/rna/output
INPUT_DIR   ?= /mnt/treasure/asc25/rna/input

HISAT_BUILD ?= ./hisat-3n/hisat-3n-build
SAMTOOLS    ?= ./samtools/samtools

# Stage 0: prepare index，这一部分不计时

${OUTPUT_DIR}/Homo_sapiens.GRCh38.%.fa: ${INPUT_DIR}/Homo_sapiens.GRCh38.%.fa.gz
	gunzip -c $< > $@

# extract index fasta
${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.index_phony: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	${HISAT_BUILD} -p 32 --base-change C,T $^ $^
	touch $@

${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	${SAMTOOLS} faidx $^

# 这里文档里写的是 OFS="\\t", 但是 bash 里应该用 \t，fish 里用 \\t
${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
	awk 'BEGIN{{OFS="\t"}}{{print $$1,$$1,0,$$2,"+"}}' $^ > $@

${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa.index_phony: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa
	${HISAT_BUILD} -p 16 --base-change C,T $^ $^
	touch $@

${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa.fai: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa
	${SAMTOOLS} faidx $^

all: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.index_phony ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa.index_phony ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa.fai ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf
