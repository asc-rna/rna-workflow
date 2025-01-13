OUTPUT_DIR  := /mnt/treasure/zhangyang/asc-rna/output
INPUT_DIR   := /mnt/treasure/zhangyang/asc-rna/input
CASE_ID     := SRR23538290

CUTSEQ      := ./cutseq/bin/cutseq
HISAT       := ./hisat-3n/hisat-3n
HISAT_BUILD := ./hisat-3n/hisat-3n-build
UMICOLLAPSE := java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=./tmp/ -jar ./umicollapse.jar
SAMTOOLS    := ./samtools/samtools


# Stage 0: prepare index，这一部分不计时

# extract index fasta
${OUTPUT_DIR}/%.fa: ${INPUT_DIR}/%.fa.gz
	gunzip -k $^

${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.3n.%.ht2: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	${HISAT_BUILD} -p 32 --base-change C,T $^ $^

${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	${SAMTOOLS} faidx $^

# 这里文档里写的是 OFS="\\t", 但是 bash 里应该用 \t，fish 里用 \\t
${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
	awk 'BEGIN{{OFS="\t"}}{{print $$1,$$1,0,$$2,"+"}}' $^ > $@

Homo_sapiens.GRCh38.ncrna.fa.3n.%.ht2: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa
	${HISAT_BUILD} -p 16 --base-change C,T $^ $^

${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa.fai: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa
	${SAMTOOLS} faidx $^

# Stage 1: 

# extract fastq

${OUTPUT_DIR}/${CASE_ID}.fastq: ${INPUT_DIR}/${CASE_ID}/${CASE_ID}.sra # 这一步也不计时
	./sra-tools/bin/fasterq-dump $^ -o $@ 
