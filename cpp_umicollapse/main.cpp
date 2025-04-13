#include <ctime>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

const size_t HASH_BASE = 329;
const int CODE_LENGTH = 5;
const int CODE_DIST = 2;
const float PERCENTAGE = 0.5f;
const int HAMMING_LIM = 1;

uint16_t ChartoBit(char c) {
    switch (c) {
        case 'A' : return 0;
        case 'C' : return 3;
        case 'T' : return 5;
        case 'G' : return 9;
        default : return 17;
    }
    __builtin_unreachable();
}

struct BitSet {
    uint64_t* bits;
    size_t UMI_length, length;

    BitSet() : bits(NULL), UMI_length(0), length(0) {
        bits = NULL;
    }

    BitSet(size_t UMI_len, char *UMI_str) : bits(NULL) {
        UMI_length = UMI_len;
        length = (UMI_len * CODE_LENGTH + 63) >> 6;
        bits = (uint64_t*)calloc(length, sizeof(uint64_t));
        if (bits == NULL) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            throw std::runtime_error("Memory allocation failed");
        }
        for (size_t i = 0, pos = 0; i < UMI_len; ++i) {
            uint16_t x = ChartoBit(UMI_str[i]);
            for (int j = 0; j < CODE_LENGTH; ++j) {
                if (x >> j & 1)
                    bits[pos >> 6] |= 1ULL << (pos & 63);
                ++pos;
            }
        }
    }

    BitSet(const BitSet& other) : bits(NULL), UMI_length(other.UMI_length), length(other.length) {
        if (other.bits) {
            bits = (uint64_t*)calloc(length, sizeof(uint64_t));
            if (bits == NULL) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                throw std::runtime_error("Memory allocation failed");
            }
            memcpy(bits, other.bits, length * sizeof(uint64_t));
        }
    }

    void operator=(const BitSet& other) {
        if (this != &other) {
            if (bits != NULL) free(bits);
            bits = NULL;
            UMI_length = other.UMI_length;
            length = other.length;
            if (other.bits) {
                bits = (uint64_t*)calloc(length, sizeof(uint64_t));
                if (bits == NULL) {
                    fprintf(stderr, "Error: Memory allocation failed\n");
                    throw std::runtime_error("Memory allocation failed");
                }
                memcpy(bits, other.bits, length * sizeof(uint64_t));
            }
        }
    }

    ~BitSet() {
        if (bits != NULL)
            free(bits);
    }
    
    size_t HammingDist(const BitSet &other) {
        if (UMI_length != other.UMI_length) {
            fprintf(stderr, "Error: UMI lengths do not match\n");
            throw std::runtime_error("UMI lengths do not match");
        }
        size_t dist = 0;
        for (size_t i = 0; i < length; ++i)
            dist += __builtin_popcountll(bits[i] ^ other.bits[i]);
        return dist / CODE_DIST;
    }

    bool operator!=(const BitSet &other) const {
        if (UMI_length != other.UMI_length) {
            fprintf(stderr, "Error: UMI lengths do not match\n");
            throw std::runtime_error("UMI lengths do not match");
        }
        for (size_t i = 0; i < length; ++i)
            if (bits[i] != other.bits[i])
                return true;
        return false;
    }

    bool operator<(const BitSet &other) const {
        if (UMI_length != other.UMI_length) {
            fprintf(stderr, "Error: UMI lengths do not match\n");
            throw std::runtime_error("UMI lengths do not match");
        }
        for (size_t i = 0; i < length; ++i)
            if (bits[i] != other.bits[i])
                return bits[i] < other.bits[i];
        return false;
    }
};


bool IsUnmapped(const bam1_t *record) {
    return (record->core.flag & BAM_FUNMAP) != 0;
}

bool IsNegativeStrand(const bam1_t *record) {
    return (record->core.flag & BAM_FREVERSE) != 0;
}

// AnalyseUMI function to extract UMI from the record

BitSet AnalyseUMI(const bam1_t *record) {
    char *qname = bam_get_qname(record);
    char *UMI_str = strrchr(qname, '_');

    if (UMI_str == NULL) {
        fprintf(stderr, "Error: UMI not found in read name %s\n", qname);
        throw std::runtime_error("UMI not found");
    }

    size_t UMI_len = strlen(++UMI_str);

    if (UMI_len == 0) {
        fprintf(stderr, "Error: UMI length is zero in read name %s\n", qname);
        throw std::runtime_error("UMI length is zero");
    }

    return BitSet(UMI_len, UMI_str);
}

int32_t getAlignPos(const bam1_t *record) {
    int32_t align_pos = record->core.pos;
    uint32_t *cigar = bam_get_cigar(record);
    int32_t n_cigar = record->core.n_cigar;
    if (IsNegativeStrand(record)) {
        for (int i = 0; i < n_cigar; ++i) {
            int32_t op = bam_cigar_op(cigar[i]);
            int32_t len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
                op == BAM_CEQUAL || op == BAM_CDIFF) {
                align_pos += len;
            }
        }
        for (int i = n_cigar - 1; i >= 0; --i) {
            int32_t op = bam_cigar_op(cigar[i]);
            int32_t len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
                align_pos += len;
            else
                break;
        }
        return align_pos;
    }
    else {
        for (int i = 0; i < n_cigar; ++i) {
            int32_t op = bam_cigar_op(cigar[i]);
            int32_t len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
                align_pos -= len;
            else
                break;
        }
        return align_pos + 1;
    }
}

struct Alignment {
    bool neg_strand;
    int32_t align_pos, ref_tid;

    Alignment() : neg_strand(false), align_pos(0), ref_tid(0) {
    }

    Alignment(const bam1_t *record) :
        neg_strand(IsNegativeStrand(record)),
        align_pos(getAlignPos(record)),
        ref_tid(record->core.tid) {
        if (ref_tid < 0) {
            fprintf(stderr, "Error: Invalid reference ID\n");
            throw std::runtime_error("Invalid reference ID");
        }
    }

    bool operator==(const Alignment &other) const {
        return neg_strand == other.neg_strand &&
                align_pos == other.align_pos &&
                  ref_tid == other.ref_tid;
    }
};

struct AlignmentHasher {
    size_t operator()(const Alignment &align) const {
        // hash involve neg_strand align_pos ref_tid
        return std::hash<bool>()(align.neg_strand) ^ std::hash<int32_t>()(align.align_pos) ^ std::hash<int32_t>()(align.ref_tid);
    }
};

struct Read {
    BitSet UMI;
    int origin_pos, avgQual;
    Read(int id, const bam1_t *record) : UMI(AnalyseUMI(record)), origin_pos(id), avgQual(0) {
        if (record == NULL) {
            fprintf(stderr, "Error: Record is NULL\n");
            throw std::runtime_error("Record is NULL");
        }
        // Caculate Average Quality
        uint8_t *qual = bam_get_qual(record);
        int32_t len = record->core.l_qseq;
        float avg = 0.0f;

        for (int32_t i = 0; i < len; ++i)
            avg += qual[i];

        avgQual = (int)(avg / len);
    }
};

struct FreqRead {
    int freq, origin_pos, jumper;
    bool exist;
    BitSet UMI;
    FreqRead(int freq_val, int id, const BitSet &umi) :
        freq(freq_val), origin_pos(id), jumper(), exist(true), UMI(umi) {
    }
};

struct AlignRead {
    int latest;
    std::vector<Read> reads;
};

std::unordered_map<Alignment, AlignRead, AlignmentHasher> align_map;
std::vector<bool> exist_pos;
std::vector<FreqRead> freq_vec;
int freq_length;
std::vector<int> jump_stack;

int JumpToExist(int x) {
    if (x < 0)
        return -1;
    if (freq_vec[x].exist)
        return x;
    return freq_vec[x].jumper = JumpToExist(freq_vec[x].jumper);
}

void RemoveNear(int x) {
    freq_vec[x].exist = false;
    freq_vec[x].jumper = x - 1;
    int limit = ceil(freq_vec[x].freq * PERCENTAGE);
    int y = JumpToExist(freq_length - 1);

    while (y >= 0) {
        if (freq_vec[y].freq > limit)
            break;
        if (freq_vec[x].UMI.HammingDist(freq_vec[y].UMI) <= HAMMING_LIM)
            RemoveNear(y);
        y = JumpToExist(y - 1);
    }
}

void DedupFreqVec() {
    std::sort(freq_vec.begin(), freq_vec.end(), [](const FreqRead &a, const FreqRead &b) {
        return a.freq > b.freq;
    });

    freq_length = freq_vec.size();

    for (int i = 0; i < freq_length; ++i)
        freq_vec[i].jumper = i;

    for (int i = 0; i < freq_length; ++i)
        if (freq_vec[i].exist) {
            exist_pos[freq_vec[i].origin_pos] = true;
            RemoveNear(i);
        }
}


void MarkDedupThreePass(htsFile *fp, htsFile *out_fp) {
    align_map.clear();
    exist_pos.clear();

    bam_hdr_t *header = sam_hdr_read(fp);
    if (header == NULL) {
        fprintf(stderr, "Error: Unable to read header\n");
        throw std::runtime_error("Unable to read header");
    }

    bam1_t *record = bam_init1();
    if (record == NULL) {
        fprintf(stderr, "Error: Unable to initialize BAM record\n");
        sam_hdr_destroy(header);
        hts_close(fp);
        hts_close(out_fp);
        throw std::runtime_error("Unable to initialize BAM record");
    }

    BGZF *bgzfp = fp->fp.bgzf;

    if (bgzfp == NULL) {
        fprintf(stderr, "Error: Unable to get BGZF pointer\n");
        bam_destroy1(record);
        sam_hdr_destroy(header);
        hts_close(fp);
        hts_close(out_fp);
        throw std::runtime_error("Unable to get BGZF pointer");
        return;
    }
    

    int64_t data_start = bgzf_tell(bgzfp);

    // Begin First Pass

    int count_unmapped = 0;
    int record_length = 0;

    fprintf(stderr, "Start First Pass\n");

    while (sam_read1(fp, header, record) >= 0) {
        if (IsUnmapped(record)) {
            ++count_unmapped;
            continue;
        }
        align_map[Alignment(record)].latest = record_length++;
    }

    exist_pos.resize(record_length, false);

    // End First Pass

    // Reset File Pointer

    if (bgzf_seek(fp->fp.bgzf, data_start, SEEK_SET) < 0) {
        fprintf(stderr, "Error: Unable to reset file pointer\n");
        bam_destroy1(record);
        sam_hdr_destroy(header);
        hts_close(fp);
        hts_close(out_fp);
        throw std::runtime_error("Unable to reset file pointer");
        return;
    }

    // Begin Second Pass

    int id = 0, UMIDedup = 0;

    fprintf(stderr, "Start Second Pass\n");

    while (sam_read1(fp, header, record) >= 0) {
        if (IsUnmapped(record))
            continue;
        std::unordered_map<Alignment, AlignRead, AlignmentHasher>::iterator
            it = align_map.find(Alignment(record));

        if (it == align_map.end()) {
            fprintf(stderr, "Error: Alignment not found in map\n");
            bam_destroy1(record);
            sam_hdr_destroy(header);
            hts_close(fp);
            hts_close(out_fp);
            throw std::runtime_error("Alignment not found in map");
        }

        #define r_vec it->second.reads

        r_vec.emplace_back(id, record);

        if (id >= it->second.latest) {
            std::sort(r_vec.begin(), r_vec.end(), [](const Read &a, const Read &b) {
                if (a.UMI != b.UMI)
                    return a.UMI < b.UMI;
                else
                    return a.avgQual < b.avgQual;
            });
            size_t reads_length = r_vec.size();
            freq_vec.clear();
            int freq = 0;
            for (size_t i = 0; i < reads_length; ++i) {
                ++freq;
                if (i == reads_length - 1 || r_vec[i].UMI != r_vec[i + 1].UMI) {
                    ++UMIDedup;
                    freq_vec.emplace_back(freq, r_vec[i].origin_pos, r_vec[i].UMI);
                    freq = 0;
                }
            }
            DedupFreqVec();
            align_map.erase(it);
        }

        #undef r_vec
        
        ++id;
    }

    // End Second Pass

    // Reset File Pointer
    
    fprintf(stderr, "Start Third Pass\n");

    if (bgzf_seek(fp->fp.bgzf, data_start, SEEK_SET) < 0) {
        fprintf(stderr, "Error: Unable to reset file pointer\n");
        bam_destroy1(record);
        sam_hdr_destroy(header);
        hts_close(fp);
        hts_close(out_fp);
        throw std::runtime_error("Unable to reset file pointer");
        return;
    }

    // Begin Third Pass - Write Pass

    id = 0;

    // Write Header

    if (sam_hdr_write(out_fp, header) < 0) {
        fprintf(stderr, "Error: Unable to write header\n");
        bam_destroy1(record);
        sam_hdr_destroy(header);
        hts_close(fp);
        hts_close(out_fp);
        throw std::runtime_error("Unable to write header");
    }
    
    int count = 0;

    while (sam_read1(fp, header, record) >= 0) {
        if (IsUnmapped(record))
            continue;
        if (exist_pos[id]) {
            ++count;
            // Write This Record
            if (sam_write1(out_fp, header, record) < 0) {
                fprintf(stderr, "Error: Unable to write record\n");
                bam_destroy1(record);
                sam_hdr_destroy(header);
                hts_close(fp);
                hts_close(out_fp);
                throw std::runtime_error("Unable to write record");
            }
        }
        ++id;
    }

    // End Third Pass

    printf("Number of removed unmapped reads %d\n", count_unmapped);
    printf("Number of reads after deduplicating %d\n", count);
    
    bam_destroy1(record);
    sam_hdr_destroy(header);
    hts_close(fp);
    hts_close(out_fp);
}

int main(int argc, char *argv[]) {
    // Record Start Time
    time_t start_time = time(NULL);

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.bam> <output.bam>\n", argv[0]);
        return 1;
    }
    
    htsFile *fp = sam_open(argv[1], "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error: Unable to open %s\n", argv[1]);
        return 1;
    }
    
    htsFile *out_fp = sam_open(argv[2], "w");
    if (out_fp == NULL) {
        fprintf(stderr, "Error: Unable to open %s\n", argv[2]);
        hts_close(fp);
        return 1;
    }

    MarkDedupThreePass(fp, out_fp);

    time_t end_time = time(NULL);

    printf("UMICollapse finished in %.2f seconds\n", difftime(end_time, start_time));
    return 0;
}