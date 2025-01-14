OUTPUT_DIR  ?= /mnt/treasure/zhangyang/asc-rna/output
INPUT_DIR   ?= /mnt/treasure/zhangyang/asc-rna/input
CASE_ID     ?= SRR23538314

CUTSEQ      ?= ./cutseq/bin/cutseq
HISAT       ?= ./hisat-3n/hisat-3n
HISAT_BUILD ?= ./hisat-3n/hisat-3n-build
HISAT_TABLE ?= ./hisat-3n/hisat-3n-table
UMICOLLAPSE ?= java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=/mnt/treasure/zhangyang/asc-rna/tmp -jar ./umicollapse.jar
SAMTOOLS    ?= ./samtools/samtools


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
	./sra-tools/bin/fasterq-dump $^ --concatenate-reads -o $@ 


# FIXME Makefile 没法处理多个目标共用一条生成命令，-j2 的时候会对两个目标同时跑两遍，需要引入一个 dummy output 解决
CUTSEQ_TARGET := ${OUTPUT_DIR}/${CASE_ID}.fastq_cut ${OUTPUT_DIR}/${CASE_ID}.fastq_tooshort ${OUTPUT_DIR}/${CASE_ID}.fastq_untrimmed
${CUTSEQ_TARGET}: ${OUTPUT_DIR}/${CASE_ID}.fastq
	${CUTSEQ} $^ -t 20 -A INLINE -m 20 --trim-polyA --ensure-inline-barcode -o ${OUTPUT_DIR}/${CASE_ID}.fastq_cut -s ${OUTPUT_DIR}/${CASE_ID}.fastq_tooshort -u ${OUTPUT_DIR}/${CASE_ID}.fastq_untrimmed

NCRNA_BAM_TARGET := ${OUTPUT_DIR}/${CASE_ID}.ncrna.unmapped.bam ${OUTPUT_DIR}/${CASE_ID}.ncrna.mapped.bam
${NCRNA_BAM_TARGET}: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa ${OUTPUT_DIR}/${CASE_ID}.fastq_cut
	${HISAT} --index ${OUTPUT_DIR}/Homo_sapiens.GRCh38.ncrna.fa --summary-file ${OUTPUT_DIR}/${CASE_ID}_map2ncrna.output.summary --new-summary -q -U ${OUTPUT_DIR}/${CASE_ID}.fastq_cut -p 16 --base-change C,T --mp 8,2 --no-spliced-alignment --directional-mapping | ${SAMTOOLS} view -@ 16 -e '!flag.unmap' -O BAM -U ${OUTPUT_DIR}/${CASE_ID}.ncrna.unmapped.bam -o ${OUTPUT_DIR}/${CASE_ID}.ncrna.mapped.bam
	# 中间会有一个 Warning: Unsupported file format

${OUTPUT_DIR}/${CASE_ID}.mRNA.fastq: ${OUTPUT_DIR}/${CASE_ID}.ncrna.unmapped.bam
	${SAMTOOLS} fastq -@ 16 -O $^ > $@

# FIXME 同上
MRNA_BAM_TARGET := ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.unmapped.bam ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.bam
${MRNA_BAM_TARGET}: ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${OUTPUT_DIR}/${CASE_ID}.mRNA.fastq
	${HISAT} --index ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa -p 16 --summary-file ${OUTPUT_DIR}/${CASE_ID}_map2genome.output.summary --new-summary -q -U ${OUTPUT_DIR}/${CASE_ID}.mRNA.fastq --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 | ${SAMTOOLS} view -@ 16 -e '!flag.unmap' -O BAM -U ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.unmapped.bam -o ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.bam
	# 可以适当改大 hisat 的并行数

${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.bam: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.bam
	${SAMTOOLS} sort -@ 16 --write-index -O BAM -o $@ $^

${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.bam.tsv: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.bam
	${SAMTOOLS} view -@ 20 -F 3980 -c $^ > $@

${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.bam
	${UMICOLLAPSE} bam -t 2 -T 16 --data naive --merge avgqual --two-pass -i $^ -o $@ > ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.log

${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam.bai: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam
	${SAMTOOLS} index -@ 8 $^ $@

${OUTPUT_DIR}/${CASE_ID}_unfiltered_uniq.tsv.gz: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam
	${SAMTOOLS} view -e "rlen<100000" -h ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam | ${HISAT_TABLE} -p 16 -u --alignments - --ref ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > ${OUTPUT_DIR}/${CASE_ID}_unfiltered_uniq.tsv.gz

${OUTPUT_DIR}/${CASE_ID}_unfiltered_multi.tsv.gz: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam
	${SAMTOOLS} view -e "rlen<100000" -h ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam | ${HISAT_TABLE} -p 16 -m --alignments - --ref ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > ${OUTPUT_DIR}/${CASE_ID}_unfiltered_multi.tsv.gz
	# hisat-table 部分很慢，每个 thread cpu 占用率不高，代码不长，可能要重写

${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam
	${SAMTOOLS} view -@ 8 -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" $^ -O BAM -o $@

${OUTPUT_DIR}/${CASE_ID}_filtered_uniq.tsv.gz: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam
	${SAMTOOLS} view -e "rlen<100000" -h ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam | ${HISAT_TABLE} -p 16 -u --alignments - --ref ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > ${OUTPUT_DIR}/${CASE_ID}_filtered_uniq.tsv.gz

${OUTPUT_DIR}/${CASE_ID}_filtered_multi.tsv.gz: ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam
	${SAMTOOLS} view -e "rlen<100000" -h ${OUTPUT_DIR}/${CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam | ${HISAT_TABLE} -p 16 -m --alignments - --ref ${OUTPUT_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > ${OUTPUT_DIR}/${CASE_ID}_filtered_multi.tsv.gz

${OUTPUT_DIR}/${CASE_ID}_genome.arrow: ${OUTPUT_DIR}/${CASE_ID}_unfiltered_uniq.tsv.gz ${OUTPUT_DIR}/${CASE_ID}_unfiltered_multi.tsv.gz ${OUTPUT_DIR}/${CASE_ID}_filtered_uniq.tsv.gz ${OUTPUT_DIR}/${CASE_ID}_filtered_multi.tsv.gz
	./m5C-UBSseq/bin/python m5C-UBSseq/bin/join_pileup.py -i $^ -o $@

# Stage 2:

${OUTPUT_DIR}/WT.arrow: ${OUTPUT_DIR}/SRR23538312_genome.arrow ${OUTPUT_DIR}/SRR23538313_genome.arrow ${OUTPUT_DIR}/SRR23538314_genome.arrow
	./m5C-UBSseq/bin/python m5C-UBSseq/bin/group_pileup.py -i $^ -o $@

${OUTPUT_DIR}/WT.prefilter.tsv: ${OUTPUT_DIR}/WT.arrow
	./m5C-UBSseq/bin/python m5C-UBSseq/bin/select_sites.py -i $^ -o $@

${OUTPUT_DIR}/%.res_phony: ${OUTPUT_DIR}/%_genome.arrow ${OUTPUT_DIR}/WT.prefilter.tsv
	./m5C-UBSseq/bin/python ./m5C-UBSseq/bin/filter_sites.py -i $< -m ${OUTPUT_DIR}/WT.prefilter.tsv -b $(basename $@).bg.tsv -o $(basename $@).filtered.tsv
	touch $@ # touch phony target

all: ${OUTPUT_DIR}/SRR23538312.res_phony ${OUTPUT_DIR}/SRR23538313.res_phony ${OUTPUT_DIR}/SRR23538314.res_phony
