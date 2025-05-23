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
#include <omp.h>

#define MAX_FILE_NUM 37
#define INF_ 1234567890

#define CHR_NUM 194
#define INPUT_PARALLEL_SIZE 4

struct data_t {
    bam1_t *record;
    int64_t pos;
    // int32_t tid;
    // data_t(bam1_t *rec) : record(bam_dup1(rec)), pos(rec->core.pos), tid(rec->core.tid) {}
    data_t(bam1_t *rec) : record(bam_dup1(rec)), pos(rec->core.pos) {}
    bool operator<(const data_t &other) const {
        return pos < other.pos;
        // return tid < other.tid || (tid == other.tid && pos < other.pos);
    }
};
std::vector<data_t> data[CHR_NUM+3];
namespace InputNameSpace {
    htsFile *fp;
    bam_hdr_t *header;
    bam1_t *record;
    void scan_record() {
        while (sam_read1(fp, header, record) >= 0) {
            data[record->core.tid].emplace_back(record);
        }
        // printf("[JZPDEBUG] jzpcnt = %d\n", jzpcnt);
    }
    void closeInputFile() {
        if (fp) hts_close(fp);
        if (record) bam_destroy1(record);
    }
    void init(const char *name_prefix, int id) {
        // printf("--------------------------\n");
        char filename[200];
        sprintf(filename, "%s.%d.bam", name_prefix, id);  // input file example: xxx.0.bam
        fp = sam_open(filename, "rb");
        if (fp == NULL) {
            fprintf(stderr, "Error: Unable to open %s\n", filename);
            exit(-1);
        }
        if (header) sam_hdr_destroy(header);
        header = sam_hdr_read(fp);
        record = bam_init1();
        scan_record();
        closeInputFile();
    }
}

int input_num, output_num;

struct OutputFile {
    htsFile *fp;
    bam_hdr_t *header;
    int num_rec;
    void init(const char *name_prefix, int id, bam_hdr_t* hdr) {
        header = hdr;
        char filename[200];
        sprintf(filename, "%s.sorted.split.%d.sam", name_prefix, id);  // output file example: xxx.split.0.sam
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
    void wirte_record(bam1_t *record) {
        if (sam_write1(fp, header, record) < 0) {
            fprintf(stderr, "Error: Unable to write record\n");
            exit(-1);
        }
    }
} output_files[MAX_FILE_NUM];

int main(int argc, char *argv[]) {
    // Record Start Time
    time_t start_time = clock();

    if (argc != 4 || strcmp(argv[1], "-h") == 0) {
        fprintf(stderr, "Usage: %s <input prefix> <input num> <output num>\n", argv[0]);
        return 1;
    }

    input_num = atoi(argv[2]);
    output_num = atoi(argv[3]);
    
    for (int i = 0; i < input_num; ++i) {
        InputNameSpace::init(argv[1], i);
        // printf("[JZPDEBUG] finished init %d\n", i);
        // fflush(stdout);
    }

    // printf("[JZPDEBUG] finished init input files, time: %.3f sec\n", 
    //        1.0*(clock()-start_time)/CLOCKS_PER_SEC);
    // fflush(stdout);


    htsFile *out_fp[MAX_FILE_NUM];
    for (int i = 0; i < output_num; ++i) {
        output_files[i].init(argv[1], i, InputNameSpace::header);
    }

    fprintf(stderr, "[JZPDEBUG] finished init output files, time: %.3f sec\n", 
        1.0*(clock()-start_time)/CLOCKS_PER_SEC);

    std::sort(data, data+CHR_NUM, [](auto a, auto b) { return a.size() > b.size(); });
    #pragma omp parallel for schedule(dynamic, 1) 
    for (int i = 0; i < CHR_NUM; ++i) {
	    int id = omp_get_thread_num();
        // fprintf(stderr, "[JZPDEBUG] data[%d].size() = %d, chromosome scheduled for process %d\n", i, data[i].size(), id); fflush(stdout);
        std::sort(data[i].begin(), data[i].end());
	    for (const auto j: data[i]) {
            output_files[id].wirte_record(j.record);
        }
    }
    for (int i = 0; i < output_num; ++i) hts_close(output_files[i].fp);

    fprintf(stderr, "[JZPDEBUG] finished writing output files, time: %.3f sec\n", 
        1.0*(clock()-start_time)/CLOCKS_PER_SEC);

    sam_hdr_destroy(InputNameSpace::header);
    // for (int i = 0; i < CHR_NUM; ++i)
    //     for (const auto j: data[i])
    //         bam_destroy1(j.record);

    fprintf(stderr, "BAM splitting finished in %.3f seconds\n", 1.0*(clock()-start_time)/CLOCKS_PER_SEC);
    return 0;
}
