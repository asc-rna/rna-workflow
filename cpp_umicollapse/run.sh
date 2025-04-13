OUTPUT_DIR="/mnt/ramdisk/rna/yyyyxh/output"
TEMP_DIR="/mnt/ramdisk/rna/yyyyxh/tmp"
CASE_ID="SRR23538290"
export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH
./umicollapse \
    "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.bam" \
    "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.dedup.bam" \
    > "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.dedup.log"