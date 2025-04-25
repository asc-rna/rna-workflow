/*

merge_split_bam:

do merge sort among input files(BAM)
write to output files according to chromosome(SAM)

*/

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
#include <ctime>
#include <fstream>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#define PII std::pair<int, int>
#define MP std::make_pair

#define MAX_FILE_NUM 23
#define INF_ 1234567890

struct InputFile {
    htsFile *fp;
    bam_hdr_t *header;
    bam1_t *record;
    PII chr_pos;
    inline void next_record() {
        int ret = sam_read1(fp, header, record);
        chr_pos = MP(INF_, INF_);
        if (ret != 0) chr_pos = MP(record->core.tid, record->core.pos);
    }
    inline void init(const char *name_prefix, int id) {
        char filename[100];
        sprintf(filename, "%s.%d.bam", name_prefix, id);  // input file example: xxx.0.bam
        fp = sam_open(filename, "rb");
        if (fp == NULL) {
            fprintf(stderr, "Error: Unable to open %s\n", filename);
            exit(-1);
        }
        header = sam_hdr_read(fp);
        record = bam_init1();
        next_record();
    }
    inline bam_hdr_t* get_header() {
        return header;
    }
    ~InputFile() {
        hts_close(fp);
        bam_destroy1(record);
        sam_hdr_destroy(header);
    }
} input_files[MAX_FILE_NUM];

int input_num, output_num;

struct OutputFile {
    htsFile *fp;
    bam_hdr_t *header;
    int num_rec;
    inline void init(const char *name_prefix, int id, bam_hdr_t* hdr) {
        header = hdr;
        char filename[100];
        sprintf(filename, "%s.split.%d.bam", name_prefix, id);  // output file example: xxx.split.0.bam
        fp = sam_open(filename, "w");
        if (fp == NULL) {
            fprintf(stderr, "Error: Unable to open %s\n", filename);
            exit(-1);
        }
        if (sam_hdr_write(fp, hdr) < 0) {
            fprintf(stderr, "Error: Unable to write header\n");
            sam_hdr_destroy(hdr);
            throw std::runtime_error("Unable to write header");
        }
        num_rec = 0;
    }
    inline void wirte_record(bam1_t *record) {
        num_rec++;
        if (sam_write1(fp, header, record) < 0) {
            fprintf(stderr, "Error: Unable to write record\n");
            exit(-1);
        }
    }
    ~OutputFile() {
        hts_close(fp);
        sam_hdr_destroy(header);
    }
} output_files[MAX_FILE_NUM];

inline int get_min_input_id() {
    PII pr = input_files[0].chr_pos;
    int min_id = 0;
    for (int i = 1; i < input_num; ++i) if (input_files[i].chr_pos < pr)
        pr = input_files[i].chr_pos, min_id = i;
    return min_id;
}

inline int get_min_output_id() {
    int min_id = 0, mn = output_files[0].num_rec;
    for (int i = 1; i < output_num; ++i) if (output_files[i].num_rec < mn)
        min_id = i, mn = output_files[i].num_rec;
    return min_id;
}

int main(int argc, char *argv[]) {
    // Record Start Time
    time_t start_time = time(NULL);

    if (argc != 4 || strcmp(argv[1], "-h") == 0) {
        fprintf(stderr, "Usage: %s <input prefix> <input num> <output prefix> <output num>\n", argv[0]);
        return 1;
    }

    input_num = atoi(argv[2]);
    output_num = atoi(argv[3]);
    
    htsFile *fp[MAX_FILE_NUM];
    for (int i = 0; i < input_num; ++i) {
        input_files[i].init(argv[1], i);
    }
    bam_hdr_t *header = input_files[0].get_header();
    
    htsFile *out_fp[MAX_FILE_NUM];
    for (int i = 0; i < output_num; ++i) {
        output_files[i].init(argv[3], i, header);
    }

    int cur_in = get_min_input_id();
    int cur_chr = input_files[cur_in].chr_pos.first, lst_chr = cur_chr;
    int cur_out = get_min_output_id();
    while (cur_chr != INF_) {
        if (lst_chr != cur_chr) cur_out = get_min_output_id();
        output_files[cur_out].wirte_record(input_files[cur_in].record);
        input_files[cur_in].next_record();
        
        cur_in = get_min_input_id();
        cur_chr = input_files[cur_in].chr_pos.first;
    }

    time_t end_time = time(NULL);

    printf("BAM splitting finished in %.2f seconds\n", difftime(end_time, start_time));
    return 0;
}
