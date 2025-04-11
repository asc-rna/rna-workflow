#include <cstdio>
#include <ctime>
#include "htslib/sam.h"


bool is_read_unmapped(const bam1_t *record) {
    return (record->core.flag & BAM_FUNMAP) != 0;
}

bool is_read_negative_strand(const bam1_t *record) {
    return (record->core.flag & BAM_FREVERSE) != 0;
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

    bam_hdr_t *header = sam_hdr_read(fp);

    if (header == NULL) {
        fprintf(stderr, "Error: Unable to read header from %s\n", argv[1]);
        hts_close(fp);
        hts_close(out_fp);
        return 1;
    }

    if (sam_hdr_write(out_fp, header) < 0) {
        fprintf(stderr, "Error: Unable to write header to %s\n", argv[2]);
        bam_hdr_destroy(header);
        hts_close(fp);
        hts_close(out_fp);
        return 1;
    }

    bam1_t *record = bam_init1();

    if (record == NULL) {
        fprintf(stderr, "Error: Unable to initialize BAM record\n");
        bam_hdr_destroy(header);
        hts_close(fp);
        hts_close(out_fp);
        return 1;
    }
    int cnt = 0;
    const int lim = 1000;

    while (sam_read1(fp, header, record) >= 0) {
        if (is_read_unmapped(record)) {
            fprintf(stderr, "Read %s is unmapped\n", bam_get_qname(record));
            continue;
        }

        // Analyze CIGAR
        int cigar_len = record->core.n_cigar;
        uint32_t *cigar = bam_get_cigar(record);

        if (cnt < lim)
            fprintf(stderr, "POS %d NEG %d NAME %s : ", record->core.pos, is_read_negative_strand(record), bam_get_qname(record));

        for (int i = 0; i < cigar_len; i++) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            if (cnt < lim)
                fprintf(stderr, "%c%d ", "MIDNSHP=XB"[op], len);
        }
        if (cnt < lim)
            fprintf(stderr, "\n");
        ++cnt;
    }

    bam_destroy1(record);
    bam_hdr_destroy(header);

    hts_close(fp);
    hts_close(out_fp);
    

    // Record End Time
    time_t end_time = time(NULL);

    printf("UMICollapse finished in %.2f seconds\n", difftime(end_time, start_time));

    return 0;
}