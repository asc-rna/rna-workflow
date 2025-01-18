# Snakefile
import os

# 配置路径
OUTPUT_DIR = "/mnt/treasure/scy24/output"
INPUT_DIR = "/mnt/treasure/asc25/rna/input"
REF_DIR = "/home/scy24/asc-rna/ref"

CUTSEQ = "cutseq"
HISAT = "./hisat-3n/hisat-3n"
HISAT_TABLE = "./hisat-3n/hisat-3n-table"
UMICOLLAPSE = "java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=/mnt/treasure/scy24/tmp -jar ./umicollapse.jar"
SAMTOOLS = "./samtools/bin/samtools"

CASE_ID = config['case_id']


# 最终目标
rule all:
    input:
        f"{OUTPUT_DIR}/{CASE_ID}_genome.arrow"

# 1. SRA 转 FASTQ
rule sra_to_fastq:
    input:
        sra_file = f"{INPUT_DIR}/{CASE_ID}/{CASE_ID}.sra"
    output:
        fastq = f"{OUTPUT_DIR}/{CASE_ID}.fastq"
    shell:
        "./sra-tools/bin/fasterq-dump {input.sra_file} --concatenate-reads -o {output.fastq}"

# 2. 修剪 FASTQ
rule cut_fastq:
    input:
        fastq = f"{OUTPUT_DIR}/{CASE_ID}.fastq"
    output:
        cut = f"{OUTPUT_DIR}/{CASE_ID}.fastq_cut",
        tooshort = f"{OUTPUT_DIR}/{CASE_ID}.fastq_tooshort",
        untrimmed = f"{OUTPUT_DIR}/{CASE_ID}.fastq_untrimmed"
    shell:
        """
        {CUTSEQ} {input.fastq} -t 20 -A INLINE -m 20 --trim-polyA --ensure-inline-barcode \
        -o {output.cut} -s {output.tooshort} -u {output.untrimmed}
        """

# 3. 比对到 ncRNA
rule align_to_ncrna:
    input:
        fastq_cut = f"{OUTPUT_DIR}/{CASE_ID}.fastq_cut",
        ref_fa_ncrna = f"{REF_DIR}/Homo_sapiens.GRCh38.ncrna.fa",
    output:
        unmapped_bam = f"{OUTPUT_DIR}/{CASE_ID}.ncrna.unmapped.bam",
        mapped_bam = f"{OUTPUT_DIR}/{CASE_ID}.ncrna.mapped.bam",
        summary = f"{OUTPUT_DIR}/{CASE_ID}_map2ncrna.output.summary"
    shell:
        """
        {HISAT} --index {input.ref_fa_ncrna} --summary-file {output.summary} \
        --new-summary -q -U {input.fastq_cut} -p 16 --base-change C,T --mp 8,2 --no-spliced-alignment \
        --directional-mapping | {SAMTOOLS} view -@ 16 -e '!flag.unmap' -O BAM \
        -U {output.unmapped_bam} -o {output.mapped_bam}
        """

# 4. 从 BAM 转 FASTQ
rule bam_to_fastq:
    input:
        unmapped_bam = f"{OUTPUT_DIR}/{CASE_ID}.ncrna.unmapped.bam"
    output:
        mRNA_fastq = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.fastq"
    shell:
        "{SAMTOOLS} fastq -@ 16 -O {input.unmapped_bam} > {output.mRNA_fastq}"

# 5. 比对到基因组
rule align_to_genome:
    input:
        fastq = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.fastq"
    output:
        unmapped_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.unmapped.bam",
        mapped_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.bam",
        summary = f"{OUTPUT_DIR}/{CASE_ID}_map2genome.output.summary"
    shell:
        """
        {HISAT} --index {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa -p 16 --summary-file {output.summary} \
        --new-summary -q -U {input.fastq} --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 | \
        {SAMTOOLS} view -@ 16 -e '!flag.unmap' -O BAM -U {output.unmapped_bam} -o {output.mapped_bam}
        """

# 6. 排序 BAM
rule sort_bam:
    input:
        mapped_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.bam"
    output:
        sorted_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.bam"
    shell:
        "{SAMTOOLS} sort -@ 16 --write-index -O BAM -o {output.sorted_bam} {input.mapped_bam}"

# 7. UMI 去重
rule umi_dedup:
    input:
        sorted_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.bam"
    output:
        dedup_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam",
        log = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.log"
    shell:
        """
        {UMICOLLAPSE} bam -t 2 -T 16 --data naive --merge avgqual --two-pass \
        -i {input.sorted_bam} -o {output.dedup_bam} > {output.log}
        """

# # 8. 转换为表格
# rule bam_to_tsv:
#     input:
#         dedup_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam"
#     output:
#         tsv = f"{OUTPUT_DIR}/{CASE_ID}_unfiltered_uniq.tsv.gz"
#     shell:
#         """
#         {SAMTOOLS} view -e "rlen<100000" -h {input.dedup_bam} | \
#         {HISAT_TABLE} -p 16 -u --alignments - --ref {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#         --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > {output.tsv}
#         """

# 8.1：生成 unfiltered_uniq.tsv.gz ==========
rule unfiltered_uniq:
    input:
        dedup_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam",
        ref_fa = f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        unfiltered_uniq = f"{OUTPUT_DIR}/{CASE_ID}_unfiltered_uniq.tsv.gz"
    shell:
        """
        {SAMTOOLS} view -e "rlen<100000" -h {input.dedup_bam} |
        {HISAT_TABLE} -p 16 -u --alignments - \
            --ref {input.ref_fa} \
            --output-name /dev/stdout --base-change C,T |
        cut -f 1,2,3,5,7 |
        gzip -c > {output.unfiltered_uniq}
        """

# 8.2：生成 unfiltered_multi.tsv.gz ==========
rule unfiltered_multi:
    input:
        dedup_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam",
        ref_fa = f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        unfiltered_multi = f"{OUTPUT_DIR}/{CASE_ID}_unfiltered_multi.tsv.gz"
    shell:
        """
        {SAMTOOLS} view -e "rlen<100000" -h {input.dedup_bam} |
        {HISAT_TABLE} -p 16 -m --alignments - \
            --ref {input.ref_fa} \
            --output-name /dev/stdout --base-change C,T |
        cut -f 1,2,3,5,7 |
        gzip -c > {output.unfiltered_multi}
        """

# 8.3：生成 filtered BAM ==========
rule filter_bam:
    """
    生成 {CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam
    主要根据 [XM], [Zf], [Yf], qlen, sclen 等 TAG 进行过滤
    """
    input:
        dedup_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.bam",
    output:
        filtered_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam"
    shell:
        """
        {SAMTOOLS} view -@ 8 \
            -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" \
            {input.dedup_bam} -O BAM -o {output.filtered_bam}
        """

# 8.4：生成 filtered_uniq.tsv.gz ==========
rule filtered_uniq:
    input:
        filtered_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam",
        ref_fa = f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        filtered_uniq = f"{OUTPUT_DIR}/{CASE_ID}_filtered_uniq.tsv.gz"
    shell:
        """
        {SAMTOOLS} view -e "rlen<100000" -h {input.filtered_bam} |
        {HISAT_TABLE} -p 16 -u --alignments - \
            --ref {input.ref_fa} \
            --output-name /dev/stdout --base-change C,T |
        cut -f 1,2,3,5,7 |
        gzip -c > {output.filtered_uniq}
        """

# 8.5：生成 filtered_multi.tsv.gz ==========
rule filtered_multi:
    input:
        filtered_bam = f"{OUTPUT_DIR}/{CASE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam",
        ref_fa = f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        filtered_multi = f"{OUTPUT_DIR}/{CASE_ID}_filtered_multi.tsv.gz"
    shell:
        """
        {SAMTOOLS} view -e "rlen<100000" -h {input.filtered_bam} |
        {HISAT_TABLE} -p 16 -m --alignments - \
            --ref {input.ref_fa} \
            --output-name /dev/stdout --base-change C,T |
        cut -f 1,2,3,5,7 |
        gzip -c > {output.filtered_multi}
        """


# 9. 合并 TSV 并生成最终输出
rule join_pileup:
    input:
        unfiltered_uniq = f"{OUTPUT_DIR}/{CASE_ID}_unfiltered_uniq.tsv.gz",
        unfiltered_multi = f"{OUTPUT_DIR}/{CASE_ID}_unfiltered_multi.tsv.gz",
        filtered_uniq = f"{OUTPUT_DIR}/{CASE_ID}_filtered_uniq.tsv.gz",
        filtered_multi = f"{OUTPUT_DIR}/{CASE_ID}_filtered_multi.tsv.gz"
    output:
        arrow = f"{OUTPUT_DIR}/{CASE_ID}_genome.arrow"
    shell:
        """
        python m5C-UBSseq/bin/join_pileup.py -i {input.unfiltered_uniq} {input.unfiltered_multi} \
        {input.filtered_uniq} {input.filtered_multi} -o {output.arrow}
        """