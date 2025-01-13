<!-- （以下均为对Problem_Description，下称文档，的提取）

分为5个阶段：构建索引、清洗数据、基因组比对、排序去重、检测和过滤，
每个阶段用不同软件去处理


初始输入文件在 文档中P8给了链接，建议用SRA Toolkit下载和提取
TODO1:下载一下输入文件，用自己的话描述一下这是什么。（是不是文本文件，能不能直接读？）

P8-P12详细介绍了每个环节的输入输出，用到的软件和常用的命令。
P14-P17 供参考的软件使用的脚本代码
疑似除了最后一步的软件在文档中给出链接的m5C-UBSseq仓库里有外，别的都软件都只给出了名字，得自己去找。

TODO2:下载一下工具，跑一下，用自己的话讲一个每个阶段在做什么，输入输出是啥，记录脚本 -->

## 依赖软件

+ 下载 sra tools: <https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz>

+ hisat3n: <http://daehwankimlab.github.io/hisat2/hisat-3n/>
```
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout -b hisat-3n origin/hisat-3n
make
```

+ samtools: <https://github.com/samtools/samtools>
```
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
# make install
```

+ umicollapse: <https://raw.githubusercontent.com/Daniel-Liu-c0deb0t/UMICollapse/refs/heads/master/umicollapse.jar>

## 下载

+ 下载输入文件

```
prefetch <ID>
fasterq-dump --fasta <ID>
```

输入文件：ascii 字符串，都是 ACGTU 状物，fastq 格式 120+G，这题共有 3 个。

小样例：SRR000001，但不知道是否能得到预期的结果。

## 各步骤

根据题目末尾的 Workflow Example

TODO: 替换里面的文件名，写成脚本/Makefile跑起来。

+ Stage 0

```
hisat-3n/hisat-3n-build -p 32 --base-change C,T /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools-1.21/samtools faidx /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai >Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf
hisat-3n/hisat-3n-build -p 16 --base-change C,T /asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa /asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa
samtools-1.21/samtools faidx /asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa
```

根据参考基因和非编码rna建立索引，这一步应该和数据集没什么关系。

TODO: 找到这里面需要的 GRCh38 dna he ncrna 文件。

+ Stage 1

```
cutseq /asc25/SRR23538290/SRR23538290.fastq -t 20 -A INLINE -m 20 --trim-polyA --ensure-inline -barcode -o /asc25/SRR23538290/SRR23538290.fastq_cut -s /asc25/SRR23538290/SRR23538290.fastq_tooshort -u /asc25/SRR23538290/SRR23538290.fastq_untrimmed

hisat-3n/hisat-3n --index /asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa --summary-file /asc25/SRR23538290/map2ncrna.output.summary --new-summary -q -U /asc25/SRR23538290/SRR23538290.fastq_cut -p 16 --base-change C,T --mp 8,2 --no-spliced-alignment --directional -mapping | /asc25/samtools-1.21/samtools view -@ 16 -e '!flag.unmap' -O BAM -U /asc25/SRR23538290/SRR23538290.ncrna.unmapped.bam -o /asc25/SRR23538290/SRR23538290.ncrna.mapped.bam

samtools-1.21/samtools fastq -@ 16 -O /asc25/SRR23538290/SRR23538290.ncrna.unmapped.bam >/asc25/SRR23538290/SRR23538290.mRNA.fastq

hisat-3n/hisat-3n --index /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -p 16 --summary-file /asc25/SRR23538290/map2genome.output.summary --new-summary -q -U /asc25/SRR23538290/SRR23538290.mRNA.fastq --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 | samtools-1.21/samtools view -@ 16 -e '!flag.unmap' -O BAM -U /asc25/SRR23538290/SRR23538290.mRNA.genome.unmapped.bam -o /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.bam

samtools-1.21/samtools sort -@ 16 --write-index -O BAM -o /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.bam

samtools-1.21/samtools view -@ 20 -F 3980 -c /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam >/asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam.tsv

java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=/asc25/SRR23538290 -jar /asc25/UMICollapse-1.0.0/umicollapse.jar bam -t 2 -T 16 --data naive --merge avgqual --two-pass -i /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam -o /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam > /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.log

samtools-1.21/samtools index -@ 8 /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam.bai

samtools-1.21/samtools view -e "rlen<100000" -h /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam | hisat-3n/hisat-3n-table -p 16 -u --alignments - --ref /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > /asc25/SRR23538290/SRR23538290_unfiltered_uniq.tsv.gz

samtools-1.21/samtools view -e "rlen<100000" -h /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam | hisat-3n/hisat-3n-table -p 16 -m --alignments - --ref /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut - f 1,2,3,5,7 | gzip -c > /asc25/SRR23538290/SRR23538290_unfiltered_multi.tsv.gz

samtools-1.21/samtools view -@ 8 -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam -O BAM -o /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam

samtools-1.21/samtools view -e "rlen<100000" -h /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam | /asc25/hisat-3n/hisat-3n-table -p 16 -u --alignments - --ref /mnt/nvme2n1/asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > /asc25/SRR23538290/SRR23538290_filtered_uniq.tsv.gz

samtools-1.21/samtools view -e "rlen<100000" -h /asc25/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam | hisat-3n/hisat-3n-table -p 16 -m --alignments - --ref /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > /asc25/SRR23538290/SRR23538290_filtered_multi.tsv.gz

python m5C-UBSseq-main/bin/join_pileup.py -i /asc25/SRR23538290/SRR23538290_unfiltered_uniq.tsv.gz /asc25/SRR23538290/SRR23538290_unfiltered_multi.tsv.gz /asc25/SRR23538290/SRR23538290_filtered_uniq.tsv.gz /asc25/SRR23538290/SRR23538290_filtered_multi.tsv.gz -o /asc25/SRR23538290/SRR23538290_genome.arrow
```

对每个数据分别运行

+ Stage 2

```
python /asc25/m5C-UBSseq-main/bin/group_pileup.py -i ./SRR23538290/SRR23538290_genome.arrow ./SRR23538291/SRR23538291_genome.arrow ./SRR23538292/SRR23538292_genome.arrow -o WT.arrow
python m5C-UBSseq-main/bin/select_sites.py -i ./WT.arrow -o ./WT.prefilter.tsv 
python ./m5C-UBSseq-main/bin/filter_sites.py -i ./SRR23538290/SRR23538290_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538290/SRR23538290.bg.tsv -o ./SRR23538290/SRR23538290.filtere
d.tsv
python ./m5C-UBSseq-main/bin/filter_sites.py -i ./SRR23538291/SRR23538291_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538291/SRR23538291.bg.tsv -o ./SRR23538291/SRR23538291.filtere
d.tsv
python ./m5C-UBSseq-main/bin/filter_sites.py -i ./SRR23538292/SRR23538292_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538292/SRR23538292.bg.tsv -o ./SRR23538292/SRR23538292.filtere
d.tsv
```

得到最终结果。
