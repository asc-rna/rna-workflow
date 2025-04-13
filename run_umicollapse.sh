OUTPUT_DIR="/mnt/ramdisk/rna/yyyyxh/output"
TEMP_DIR="/mnt/ramdisk/rna/yyyyxh/tmp"
CASE_ID="SRR23538292"

cd Umicollapse_upd

rm -f umicollapse.jar
rm -rf bin

mkdir bin

javac --release 11 -cp "lib/htsjdk-4.1.3-imp.jar:lib/snappy-java-1.1.7.3.jar" -d bin Main.java
if [ $? -ne 0 ]; then
    echo "Java compilation failed. Exiting."
    exit 1
else
    echo "Java compilation succeeded."
fi

cd bin

jar -c -m ../Manifest.txt -f ../umicollapse.jar *.class
if [ $? -ne 0 ]; then
    echo "Jar creation failed. Exiting."
    exit 1
else
    echo "Jar creation succeeded."
fi

cd ..

cd ..


rm -f "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.dedup.bam"
rm -f "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.dedup.log"

# exit 0

java -server -Xms10G -Xmx40G -Xss256K -Djava.io.tmpdir="$TEMP_DIR" \
    -XX:+UseZGC \
    -Dthread.pool.size=16 \
    -Dsamjdk.sort_col_threads=2 \
    -Dsamjdk.use_async_io_read_samtools=true \
    -Dsamjdk.compression_level=0 \
    -jar ./Umicollapse_upd/umicollapse.jar \
    "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.bam" \
    "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.dedup.bam" \
    > "$OUTPUT_DIR/$CASE_ID.mRNA.genome.mapped.sorted.dedup.log"
