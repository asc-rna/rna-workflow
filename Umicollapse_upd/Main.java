import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFileHeader;

import java.util.Set;
import java.util.Map;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.stream.Stream;
import java.util.stream.Collectors;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import java.io.File;

class BitSet implements Comparable{
    private static final int CHUNK_SIZE = 64;

    private long[] bits;
    private long[] nBits;
    private boolean recalcHash;
    private int hash;

    public BitSet(int length){
        this.bits = new long[length / CHUNK_SIZE + (length % CHUNK_SIZE == 0 ? 0 : 1)];
        this.recalcHash = true;
    }

    private BitSet(long[] bits){
        this.bits = bits;
        this.recalcHash = true;
    }

    private BitSet(long[] bits, int hash){
        this.bits = bits;
        this.recalcHash = false;
        this.hash = hash;
    }

    public boolean get(int idx){
        return (bits[idx / CHUNK_SIZE] & (1L << (idx % CHUNK_SIZE))) != 0L;
    }

    // does not set the nBits array, so distance calculations could be wrong if not careful!
    public void set(int idx, boolean bit){
        recalcHash = true;
        int i = idx / CHUNK_SIZE;
        int j = idx % CHUNK_SIZE;
        bits[i] = bit ? (bits[i] | (1L << j)) : (bits[i] & ~(1L << j));
    }

    public void setNBit(int idx, boolean bit){
        if(nBits == null)
            nBits = new long[bits.length];

        int i = idx / CHUNK_SIZE;
        int j = idx % CHUNK_SIZE;
        nBits[i] = bit ? (nBits[i] | (1L << j)) : (nBits[i] & ~(1L << j));
    }

    public int bitCountXOR(BitSet o){
        int res = 0;

        for(int i = 0; i < bits.length; i++){
            long xor = (nBits == null ? 0L : nBits[i]) ^ (o.nBits == null ? 0L : o.nBits[i]);
            res += Long.bitCount(xor | (bits[i] ^ o.bits[i])) - Long.bitCount(xor) / Read.ENCODING_LENGTH;
        }

        return res;
    }

    @Override
    public boolean equals(Object obj){
        if(!(obj instanceof BitSet))
            return false;

        BitSet o = (BitSet)obj;

        if(this == o)
            return true;

        if(bits.length != o.bits.length)
            return false;

        for(int i = 0; i < bits.length; i++){
            if(bits[i] != o.bits[i])
                return false;
        }

        return true;
    }

    @Override
    public int compareTo(Object o){
        BitSet other = (BitSet)o;

        if(bits.length != other.bits.length)
            return bits.length - other.bits.length;

        for(int i = 0; i < bits.length; i++){
            if(bits[i] != other.bits[i])
                return Long.compare(bits[i], other.bits[i]);
        }

        return 0;
    }

    public BitSet clone(){
        if(recalcHash)
            return new BitSet(Arrays.copyOf(bits, bits.length));
        else
            return new BitSet(Arrays.copyOf(bits, bits.length), hash);
    }

    @Override
    public int hashCode(){
        if(recalcHash){
            long h = 1234L; // same as Java's built-in BitSet hash function

            for(int i = bits.length; --i >= 0;)
                h ^= bits[i] * (i + 1L);

            hash = (int)((h >> 32) ^ h);
            recalcHash = false;
        }

        return hash;
    }

    @Override
    public String toString(){
        StringBuilder res = new StringBuilder();

        for(int i = 0; i < bits.length; i++){
            String s = Long.toBinaryString(bits[i]);
            res.append(reverse(s));
            res.append(make('0', CHUNK_SIZE - s.length()));
        }

        return res.toString();
    }

    private String make(char c, int n){
        char[] res = new char[n];

        for(int i = 0; i < n; i++)
            res[i] = c;

        return new String(res);
    }

    private String reverse(String s){
        char[] res = new char[s.length()];

        for(int i = 0; i < s.length(); i++)
            res[i] = s.charAt(s.length() - 1 - i);

        return new String(res);
    }
}

class Read{
    public static final int ENCODING_DIST = 2;
    public static final int ENCODING_LENGTH = 3;
    public static final Map<Character, Integer> ENCODING_MAP = new HashMap<>();
    public static final Map<Integer, Integer> ENCODING_IDX = new HashMap<>();
    public static final char[] ALPHABET = {'A', 'T', 'C', 'G', 'N'};
    public static final int UNDETERMINED = 0b100;
    public static final char UNDETERMINED_CHAR = 'N';
    public static final int ANY = 0b111;

    static{
        ENCODING_MAP.put('A', 0b000);
        ENCODING_MAP.put('T', 0b101);
        ENCODING_MAP.put('C', 0b110);
        ENCODING_MAP.put('G', 0b011);
        ENCODING_MAP.put(UNDETERMINED_CHAR, UNDETERMINED);

        ENCODING_IDX.put(0b000, 0);
        ENCODING_IDX.put(0b101, 1);
        ENCODING_IDX.put(0b110, 2);
        ENCODING_IDX.put(0b011, 3);
        ENCODING_IDX.put(UNDETERMINED, 4);
    }

    private static Pattern defaultUMIPattern;
    private SAMRecord record;
    private int avgQual;
    private String umi;

    private static BitSet toBitSet(String s){
        BitSet res = new BitSet(s.length() * Read.ENCODING_LENGTH);

        for(int i = 0; i < s.length(); i++){
            charSet(res, i, Read.ENCODING_MAP.get(s.charAt(i)));

            if(s.charAt(i) == Read.UNDETERMINED_CHAR)
                charSetNBit(res, i);
        }

        return res;
    }

    private static BitSet charSet(BitSet a, int idx, int b){
        for(int i = 0; i < Read.ENCODING_LENGTH; i++)
            a.set(idx * Read.ENCODING_LENGTH + i, ((b & (1 << i)) != 0));

        return a;
    }

    private static BitSet charSetNBit(BitSet a, int idx){
        for(int i = 0; i < Read.ENCODING_LENGTH; i++)
            a.setNBit(idx * Read.ENCODING_LENGTH + i, true);

        return a;
    }

    public Read(SAMRecord record){
        this.record = record;
        int idx = record.getReadName().lastIndexOf('_');
        this.umi = (idx == -1 || idx == record.getReadName().length() - 1) ? "" : record.getReadName().substring(idx + 1);

        float avg = 0.0f;

        for(byte b : record.getBaseQualities())
            avg += b;

        this.avgQual = (int)(avg / record.getReadLength());
    }

    public BitSet getUMI(int maxLength){
        String umi_cut = this.umi;
        if (maxLength >= 0 && umi.length() > maxLength) {
            umi_cut = this.umi.substring(0, maxLength);
        }

        return toBitSet(umi_cut);
    }

    public int getUMILength(){
        return umi.length();
    }

    public int getAvgQual(){
        return avgQual;
    }

    @Override
    public boolean equals(Object o){
        Read r = (Read)o;
        return record.equals(r.record);
    }

    public int getMapQual(){
        return record.getMappingQuality();
    }

    public SAMRecord toSAMRecord(){
        return record;
    }
}

class ReadFreq{
    public Read read;
    public int freq;

    public ReadFreq(Read read, int freq){
        this.read = read;
        this.freq = freq;
    }
}

class UmiFreq{
    public BitSet umi;
    public ReadFreq readFreq;

    public UmiFreq(BitSet umi, ReadFreq readFreq){
        this.umi = umi;
        this.readFreq = readFreq;
    }
}

class ClusterTracker{
    private boolean track;
    private int offset;

    private List<BitSet> temp;
    private int tempFreq;

    private Map<BitSet, Integer> toUniqueIdx;
    private List<ClusterStats> clusters;
    private int idx;

    public ClusterTracker(boolean track){
        this.track = track;
        this.offset = 0;

        this.temp = new ArrayList<BitSet>();
        this.tempFreq = 0;

        this.toUniqueIdx = new HashMap<BitSet, Integer>();
        this.clusters = new ArrayList<ClusterStats>();
        this.idx = 0;
    }

    public boolean shouldTrack(){
        return this.track;
    }

    public void setOffset(int offset){
        this.offset = offset;
    }

    public int getOffset(){
        return this.offset;
    }

    public void addAll(Set<BitSet> s, Map<BitSet, ReadFreq> reads){
        if(this.track){
            this.temp.addAll(s);

            for(BitSet umi : s)
                this.tempFreq += reads.get(umi).freq;
        }
    }

    public void addAllImp(List<BitSet> s, Map<BitSet, ReadFreq> reads){
        if (this.track) {
            // 遍历数组，把所有 BitSet 加入 temp
            for (BitSet umi : s) {
                this.temp.add(umi);
    
                // 防止 reads.get(umi) 为空导致 NullPointerException
                ReadFreq rf = reads.get(umi);
                this.tempFreq += rf.freq;
                
            }
        }
    }

    public void track(BitSet unique, Read read){
        if(this.track){
            for(BitSet s : this.temp)
                this.toUniqueIdx.put(s, idx);

            clusters.add(new ClusterStats(unique, tempFreq, read));

            this.temp.clear();
            this.tempFreq = 0;
            this.idx++;
        }
    }

    public int getId(BitSet umi){
        return toUniqueIdx.get(umi);
    }

    public ClusterStats getStats(int id){
        return clusters.get(id);
    }

    public static class ClusterStats{
        private BitSet umi;
        private int freq;
        private Read read;

        public ClusterStats(BitSet umi, int freq, Read read){
            this.umi = umi;
            this.freq = freq;
            this.read = read;
        }

        public BitSet getUMI(){
            return this.umi;
        }

        public int getFreq(){
            return this.freq;
        }

        public Read getRead(){
            return this.read;
        }
    }
}

class Data{
    private Map<BitSet, Integer> umiFreq;

    public static int umiDist(BitSet a, BitSet b){
        // divide by the pairwise Hamming distance in the encoding
        return a.bitCountXOR(b) / Read.ENCODING_DIST;
    }

    public void init(Map<BitSet, Integer> umiFreq, int umiLength, int maxEdits){
        this.umiFreq = umiFreq;
    }

    public Set<BitSet> removeNear(BitSet umi, int k, int maxFreq){
        Set<BitSet> res = new HashSet<>();

        for(Iterator<Map.Entry<BitSet, Integer>> it = umiFreq.entrySet().iterator(); it.hasNext();){
            Map.Entry<BitSet, Integer> e = it.next();
            BitSet o = e.getKey();
            int f = e.getValue();
            int dist = umiDist(umi, o);

            if(dist <= k && (dist == 0 || f <= maxFreq)){
                res.add(o);
                it.remove();
            }
        }
        return res;
        
    }

    public List<BitSet> removeNearImp(BitSet umi, int maxFreq) {
        List<BitSet> res = new ArrayList<>(1);  // 用 ArrayList 替代 HashSet

        Iterator<Map.Entry<BitSet, Integer>> it = umiFreq.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<BitSet, Integer> e = it.next();
            BitSet o = e.getKey();
            int f = e.getValue();
            int dist = umiDist(umi, o);

            if (dist <= 1 && (dist == 0 || f <= maxFreq)) {
                res.add(o);
                it.remove();
            }
        }

        return res;
    }   

    public boolean contains(BitSet umi){
        return umiFreq.containsKey(umi);
    }

    public Map<String, Float> stats(){
        Map<String, Float> res = new HashMap<>();
        return res;
    }
}

class Merge{
    public Read merge(Read a, Read b){
        if(a.getAvgQual() >= b.getAvgQual())
            return a;
        else
            return b;
    }
}

class Algo{
    public List<Read> apply(Map<BitSet, ReadFreq> reads, Data data, ClusterTracker tracker, int umiLength){
        UmiFreq[] freq = new UmiFreq[reads.size()];
        List<Read> res = new ArrayList<>();
        Map<BitSet, Integer> m = new HashMap<>();

        int idx = 0;

        for(Map.Entry<BitSet, ReadFreq> e : reads.entrySet()){
            freq[idx] = new UmiFreq(e.getKey(), e.getValue());
            m.put(e.getKey(), e.getValue().freq);
            idx++;
        }


        Arrays.sort(freq, (a, b) -> b.readFreq.freq - a.readFreq.freq);
        data.init(m, umiLength, 1);

        for(int i = 0; i < freq.length; i++){
            if(data.contains(freq[i].umi)){
                visitAndRemove(freq[i].umi, reads, data, tracker);
                tracker.track(freq[i].umi, freq[i].readFreq.read);
                res.add(freq[i].readFreq.read);
            }
        }

        return res;
    }

    private void visitAndRemove(BitSet u, Map<BitSet, ReadFreq> reads, Data data, ClusterTracker tracker){
        List<BitSet> ci = data.removeNearImp(u, (int)(0.5f * (reads.get(u).freq + 1)));
        tracker.addAllImp(ci, reads);

        for (BitSet v : ci) {
            if (u.equals(v))
                continue;

            visitAndRemove(v, reads, data, tracker);
        }
    }
}

class DeduplicateSAM{
    private int avgUMICount;
    private int maxUMICount;
    private int dedupedCount;
    private int umiLength;
    public static final int HASH_CONST = 31;

    // trade off speed for lower memory usage
    // input should be sorted based on alignment for best results
    public void deduplicateAndMergeTwoPass(File in, File out) {
        Algo algo = new Algo();
        Merge merge = new Merge();
        System.out.println("[" + LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss")) + "] Start the first pass!");


        SamReader firstPass = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(in);
        Writer writer = new Writer(in, out, firstPass);
        Map<Alignment, AlignReads> align = new HashMap<>(1 << 29);
        int totalReadCount = 0;
        int unmapped = 0;
        int unpaired = 0;
        int chimeric = 0;
        int readCount = 0;

        // first pass to figure out where each alignment position ends
        for(SAMRecord record : firstPass){
            // always skip the reversed read

            totalReadCount++;

            if(record.getReadUnmappedFlag()){ // discard unmapped reads
                unmapped++;
                continue;
            }


            Alignment alignment = new Alignment(
                record.getReadNegativeStrandFlag(),
                record.getReadNegativeStrandFlag() ? record.getUnclippedEnd() : record.getUnclippedStart(),
                record.getReferenceName()
            );

            if(!align.containsKey(alignment))
                align.put(alignment, new AlignReads());

            align.get(alignment).latest = readCount;
            readCount++;
        }

        try{
            firstPass.close();
        }catch(Exception e){
            e.printStackTrace();
        }

        firstPass = null;

        System.gc(); // attempt to clear up memory before second pass

        System.out.println("[" + LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss")) + "] Done with the first pass!");

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(in);

        int idx = 0;
        int alignPosCount = align.size();
        avgUMICount = 0;
        maxUMICount = 0;
        dedupedCount = 0;

        for(SAMRecord record : reader){
            if(record.getReadUnmappedFlag()) // discard unmapped reads
                continue;

            Alignment alignment = null;

            alignment = new Alignment(
                    record.getReadNegativeStrandFlag(),
                    record.getReadNegativeStrandFlag() ? record.getUnclippedEnd() : record.getUnclippedStart(),
                    record.getReferenceName()
            );

            AlignReads alignReads = align.get(alignment);

            if(alignReads.umiRead == null)
                alignReads.umiRead = new HashMap<BitSet, ReadFreq>(4);

            Read read = new Read(record);
            BitSet umi = read.getUMI(-1);
            umiLength = read.getUMILength();

            if(alignReads.umiRead.containsKey(umi)){
                ReadFreq prev = alignReads.umiRead.get(umi);
                prev.read = merge.merge(read, prev.read);
                prev.freq++;
            }else{
                alignReads.umiRead.put(umi, new ReadFreq(read, 1));
            }

            if(idx >= alignReads.latest){
                List<Read> deduped;
                Data data = new Data();

                deduped = algo.apply(alignReads.umiRead, data, new ClusterTracker(false), umiLength);

                avgUMICount += alignReads.umiRead.size();
                maxUMICount = Math.max(maxUMICount, alignReads.umiRead.size());
                dedupedCount += deduped.size();

                for(Read r : deduped)
                    writer.write(r.toSAMRecord());

                // done with the current alignment position, so free up memory
                align.remove(alignment);
            }

            idx++;
        }


        System.out.println("[" + LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss")) + "] Done with the second pass!");
        try{
            reader.close();
        }catch(Exception e){
            e.printStackTrace();
        }

        writer.close();

        System.out.println("Number of input reads\t" + totalReadCount);
        System.out.println("Number of removed unmapped reads\t" + unmapped);

        System.out.println("Number of unremoved reads\t" + readCount);
        System.out.println("Number of unique alignment positions\t" + alignPosCount);
        System.out.println("Average number of UMIs per alignment position\t" + ((double)avgUMICount / alignPosCount));
        System.out.println("Max number of UMIs over all alignment positions\t" + maxUMICount);
        System.out.println("Number of reads after deduplicating\t" + dedupedCount);
    }

    // heavily inspired by TwoPassPairWriter from UMI-tools
    private static class Writer{
        private SAMFileWriter writer;
        private File in;
        private String ref = null;

        public Writer(File in, File out, SamReader r){
            SAMFileHeader header = r.getFileHeader();
            this.writer = new SAMFileWriterFactory()
            .makeSAMOrBAMWriter(header, false, out);
        } 

        public void write(SAMRecord record){
            writer.addAlignment(record);
        }

        public void close(){
            writer.close();
        }

    }

    private static class AlignReads{
        public int latest;
        public Map<BitSet, ReadFreq> umiRead;

        public AlignReads(){
            this.latest = 0;
            this.umiRead = null;
        }
    }

    private static class Alignment implements Comparable{
        private boolean strand;
        private int coord;
        private String ref;

        public Alignment(boolean strand, int coord, String ref){
            this.strand = strand;
            this.coord = coord;
            this.ref = ref.intern();
        }

        public String getRef(){
            return ref;
        }

        @Override
        public boolean equals(Object o){
            if(!(o instanceof Alignment))
                return false;

            Alignment a = (Alignment)o;

            if(this == a)
                return true;

            if(strand != a.strand)
                return false;

            if(coord != a.coord)
                return false;

            if(ref != a.ref) // can directly compare interned strings
                return false;

            return true;
        }

        @Override
        public int hashCode(){
            int hash = strand ? 1231 : 1237;
            hash = hash * HASH_CONST + coord;
            hash = hash * HASH_CONST + ref.hashCode();
            return hash;
        }

        @Override
        public int compareTo(Object o){
            Alignment other = (Alignment)o;

            if(strand != other.strand)
                return Boolean.compare(strand, other.strand);

            if(coord != other.coord)
                return coord - other.coord;

            return ref.compareTo(other.ref);
        }
    }
}

public class Main{
    public static void main(String[] args){
        System.out.println("Arguments\t" + Arrays.toString(args));

        long startTime = System.currentTimeMillis();

        if(args.length != 2)
            throw new IllegalArgumentException("Please provide input and output files!");

        File in = new File(args[0]);
        File out = new File(args[1]);

        System.out.println("input: " + args[0]);
        System.out.println("output: " + args[1]);

        DeduplicateSAM dedup = new DeduplicateSAM();
        dedup.deduplicateAndMergeTwoPass(in, out);

        System.out.println("UMI collapsing finished in " + ((System.currentTimeMillis() - startTime) / 1000.0) + " seconds!");
    }
}