OUTPUT_DIR  ?= /mnt/treasure/zhangyang/asc-rna/output
INPUT_DIR   ?= /mnt/treasure/zhangyang/asc-rna/input

CUTSEQ      ?= ./cutseq/bin/cutseq
HISAT       ?= ./hisat-3n/hisat-3n
HISAT_BUILD ?= ./hisat-3n/hisat-3n-build
HISAT_TABLE ?= ./hisat-3n/hisat-3n-table
UMICOLLAPSE ?= java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=/mnt/treasure/zhangyang/asc-rna/tmp -jar ./umicollapse.jar
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

# Stage 1:

${OUTPUT_DIR}/SRR23538290_genome.arrow:
	make -f Makefile_case CASE_ID=SRR23538290 $@
${OUTPUT_DIR}/SRR23538291_genome.arrow:
	make -f Makefile_case CASE_ID=SRR23538291 $@
${OUTPUT_DIR}/SRR23538292_genome.arrow:
	make -f Makefile_case CASE_ID=SRR23538292 $@

# Stage 2:

${OUTPUT_DIR}/WT.arrow: ${OUTPUT_DIR}/SRR23538290_genome.arrow ${OUTPUT_DIR}/SRR23538291_genome.arrow ${OUTPUT_DIR}/SRR23538292_genome.arrow
	./m5C-UBSseq/bin/python m5C-UBSseq/bin/group_pileup.py -i $^ -o $@

${OUTPUT_DIR}/WT.prefilter.tsv: ${OUTPUT_DIR}/WT.arrow
	./m5C-UBSseq/bin/python m5C-UBSseq/bin/select_sites.py -i $^ -o $@

${OUTPUT_DIR}/%.res_phony: ${OUTPUT_DIR}/%_genome.arrow ${OUTPUT_DIR}/WT.prefilter.tsv
	./m5C-UBSseq/bin/python ./m5C-UBSseq/bin/filter_sites.py -i $< -m ${OUTPUT_DIR}/WT.prefilter.tsv -b $(basename $@).bg.tsv -o $(basename $@).filtered.tsv
	touch $@ # touch phony target

all: ${OUTPUT_DIR}/SRR23538290.res_phony ${OUTPUT_DIR}/SRR23538291.res_phony ${OUTPUT_DIR}/SRR23538292.res_phony
