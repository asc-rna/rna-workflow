CUTSEQ := ./cutseq/bin/cutseq
HISAT  := ./hisat-3n/hisat-3n
HISAT_BUILD := ./hisat-3n/hisat-3n-build
UMICOLLAPSE := java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=./tmp/ -jar ./umicollapse.jar
SAMTOOLS := ./samtools/samtools

# Stage 0: prepare index，这一部分不计时

# extract index fasta
data/%.fa: data/%.fa.gz
	gunzip -k $^

data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.3n.%.ht2: data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	${HISAT_BUILD} -p 32 --base-change C,T $^ $^

data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai: data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	${SAMTOOLS} faidx $^

# 这里文档里写的是 OFS="\\t", 但是 bash 里应该用 \t，fish 里用 \\t
data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf: data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
	awk 'BEGIN{{OFS="\t"}}{{print $$1,$$1,0,$$2,"+"}}' $^ > $@

Homo_sapiens.GRCh38.ncrna.fa.3n.%.ht2: data/Homo_sapiens.GRCh38.ncrna.fa
	${HISAT_BUILD} -p 16 --base-change C,T $^ $^

data/Homo_sapiens.GRCh38.ncrna.fa.fai: data/Homo_sapiens.GRCh38.ncrna.fa
	${SAMTOOLS} faidx $^

# Stage 1: 
