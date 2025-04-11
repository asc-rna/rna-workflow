#include <cstdio>
#include <ctime>
#include "htslib/sam.h"

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

    

    hts_close(fp);
    hts_close(out_fp);


    // Record End Time
    time_t end_time = time(NULL);

    printf("UMICollapse finished in %.2f seconds\n", difftime(end_time, start_time));

    return 0;
}