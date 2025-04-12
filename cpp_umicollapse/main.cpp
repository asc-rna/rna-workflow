#include <ctime>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <functional>
#include <unordered_map>
#include "htslib/sam.h"

const size_t HASH_BASE = 329;
const int CODE_LENGTH = 5;
const int CODE_DIST = 2;

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
    size_t UMI_length, length, hash;

    BitSet() : bits(NULL), UMI_length(0), length(0), hash(0) {
        bits = NULL;
    }

    BitSet(size_t UMI_len, char *UMI_str) : bits(NULL), hash(0) {
        UMI_length = UMI_len;
        length = (UMI_len * CODE_LENGTH + 63) >> 6;
        bits = (uint64_t*)calloc(length, sizeof(uint64_t));
        if (bits == NULL) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            throw std::runtime_error("Memory allocation failed");
        }
        hash = 0;
        for (int i = 0, pos = 0; i < UMI_len; ++i) {
            uint16_t x = ChartoBit(UMI_str[i]);
            hash = hash * HASH_BASE + UMI_str[i];
            for (int j = 0; j < CODE_LENGTH; ++j) {
                if (x >> j & 1) {
                    int pos = i * CODE_LENGTH + j;
                    bits[pos >> 6] |= 1ULL << (pos & 63);
                }
                ++pos;
            }
        }
    }

    BitSet(const BitSet& other) : bits(NULL), UMI_length(other.UMI_length), length(other.length), hash(other.hash) {
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
            hash = other.hash;
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

    bool operator==(const BitSet &other) const {
        if (UMI_length != other.UMI_length) {
            fprintf(stderr, "Error: UMI lengths do not match\n");
            throw std::runtime_error("UMI lengths do not match");
        }
        return hash == other.hash && memcmp(bits, other.bits, length * sizeof(uint64_t)) == 0;
    }
};

struct BitSetHasher {
    size_t operator()(const BitSet &bitset) const {
        return bitset.hash;
    }
};

bool IsUnmapped(const bam1_t *record) {
    return (record->core.flag & BAM_FUNMAP) != 0;
}

bool IsNegativeStrand(const bam1_t *record) {
    return (record->core.flag & BAM_FREVERSE) != 0;
}

void AnalyseUMI(const bam1_t *record, BitSet &UMI) {
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

    UMI = BitSet(UMI_len, UMI_str);
}

struct Alignment {
    bool neg_strand;
    int32_t align_pos, ref_tid;
    size_t origin_pos;

    Alignment() : neg_strand(false), align_pos(0), origin_pos(0), ref_tid(0) {
    }

    bool operator==(const Alignment &other) const {
        return neg_strand == other.neg_strand && align_pos == other.align_pos && ref_tid == other.ref_tid;
    }
};

struct AlignmentHasher {
    size_t operator()(const Alignment &align) const {
        // hash involve neg_strand align_pos ref_tid
        return std::hash<bool>()(align.neg_strand) ^ std::hash<int32_t>()(align.align_pos) ^ std::hash<int32_t>()(align.ref_tid);
    }
};

unordered_map<Aligment, BitSet, AlignmentHasher> UMI_map;

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

std::vector<bam1_t*> orgin_record;

void MarkDedup(htsFile *fp, htsFile *out_fp) {
    fflush(stderr);
    bam_hdr_t *header = sam_hdr_read(fp);
    if (header == NULL) {
        fprintf(stderr, "Error: Unable to read header\n");
        throw std::runtime_error("Unable to read header");
    }

    bam1_t *record = bam_init1();
    if (record == NULL) {
        fprintf(stderr, "Error: Unable to initialize BAM record\n");
        throw std::runtime_error("Unable to initialize BAM record");
    }

    int count_unmapped = 0;

    while (sam_read1(fp, header, record) >= 0) {
        if (IsUnmapped(record)) {
            ++count_unmapped;
            continue;
        }
        orgin_record.push_back(bam_dup1(record));
    }

    for (size_t i = 0; i < orgin_record.size(); ++i) {
        Alignment align;
        align.neg_strand = IsNegativeStrand(orgin_record[i]);
        align.align_pos = getAlignPos(orgin_record[i]);
        align.origin_pos = i;
        align.ref_tid = orgin_record[i]->core.tid;
    }

    for (size_t i = 0; i < orgin_record.size(); ++i)
        bam_destroy1(orgin_record[i]);
    bam_destroy1(record);
    sam_hdr_destroy(header);
    orgin_record.clear();
    orgin_record.shrink_to_fit();
}

int main(int argc, char *argv[]) {
    // Record Start Time
    time_t start_time = time(NULL);

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.bam> <output.bam>\n", argv[0]);
        return 1;
    }
    
    htsFile *fp = hts_open(argv[1], "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: Unable to open %s\n", argv[1]);
        return 1;
    }
    
    htsFile *out_fp = hts_open(argv[2], "wb");
    if (out_fp == NULL) {
        fprintf(stderr, "Error: Unable to open %s\n", argv[2]);
        hts_close(fp);
        return 1;
    }

    MarkDedup(fp, out_fp);

    hts_close(fp);
    hts_close(out_fp);

    // Record End Time
    time_t end_time = time(NULL);

    printf("UMICollapse finished in %.2f seconds\n", difftime(end_time, start_time));
    return 0;
}