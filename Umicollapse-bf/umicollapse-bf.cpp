# include <iostream>
# include <ctime>
# include <string>
# include <unordered_map>
# include "htslib/sam.h"

using namespace std;

inline int getUnclippedStart(bam1_t *b) {
    int alignmentStart = b->core.pos;
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    for (int i = 0; i < b->core.n_cigar; ++i) {                                   //no. of cigar data entries
        int len = bam_cigar_oplen(cigar[i]);
        char ch = bam_cigar_op(cigar[i]);       //the macros gives the count of operation and the symbol of operation for given cigar entry
        if (ch != BAM_CSOFT_CLIP && ch != BAM_CHARD_CLIP) {
            break;
        }
        alignmentStart -= len;
    }
    return alignmentStart;
}

inline int getReferenceLength(bam1_t *b) {
    // TODO: 和 htsjdk.samtools/Cigar.class 里面的 getReferenceLength 对齐。但并没有搞懂为啥里面是对 1,2,3,4,5 进行加法，感觉应该是 M, D, N, =, X 才对
    return -1;
}

inline int getUnclippedEnd(bam1_t *b) {
    int alignmentEnd = b->core.pos + getReferenceLength(b) - 1;
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    for (int i = b->core.n_cigar - 1; i >= 0; --i) {                                   //no. of cigar data entries
        int len = bam_cigar_oplen(cigar[i]);
        char ch = bam_cigar_op(cigar[i]);       //the macros gives the count of operation and the symbol of operation for given cigar entry
        if (ch != BAM_CSOFT_CLIP && ch != BAM_CHARD_CLIP) {
            break;
        }
        alignmentEnd += len;
    }
    return alignmentEnd;
}

struct Alignment {
    bool strand;
    int coord;
    int ref;
    Alignment (bam1_t *bamdata) {
        strand = (bamdata->core.flag & BAM_FREVERSE) ? true : false;
        coord = strand ? getUnclippedEnd(bamdata) : getUnclippedStart(bamdata);
        ref = bamdata->core.tid;
    }
    inline long long get_hash() {
        return 1ll * coord << 32 | ref << 1 | strand;
    }
};

const char *input_path;
const char *output_path;
void parseOptions(int argc, const char **argv) {
    /*
        很简陋的 parser
        只支持两个参数，且第一个是输入文件，第二个是输出文件
    */
    if (argc != 3) {
        cerr << "Usage: ./umicollapse_bf /path/to/input.bam /path/to/output.bam" << endl;
        throw (1);
    }
    input_path = argv[1];
    output_path = argv[2];
}

void deduplicateAndMerge(sam_hdr_t *in_samhdr) {
    bam1_t *bamdata = NULL;
    if (!(bamdata = bam_init1())) {
        printf("Failed to allocate data memory!\n");
        throw(1);
    }

    printf("Start the first pass!\n");
    clock_t time = clock();

    unordered_map<long long, int> latest;
    int totalReadCount = 0, unmapped = 0, unpaired = 0, chimeric = 0, readCount = 0;
    samFile *firstPass = NULL;
    //open input file
    if (!(firstPass = sam_open(input_path, "r"))) {
        printf("Could not open %s\n", input_path);
        throw(2);
    }
    while (sam_read1(firstPass, in_samhdr, bamdata) >= 0) {
        totalReadCount++;
        if (bamdata->core.flag & BAM_FUNMAP) {  // discard unmapped reads
            unmapped++;
            continue;
        }
        Alignment alignment = Alignment(bamdata);
        latest[alignment.get_hash()] = readCount;
        readCount++;
    }
}

void umicollapse_bf() {
    const char *flags = NULL, *tidname = NULL;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    uint8_t *data = NULL;
    uint32_t *cigar = NULL;

    //open input file
    if (!(infile = sam_open(input_path, "r"))) {
        printf("Could not open %s\n", input_path);
        throw(2);
    }
    // read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        throw(3);
    }

    deduplicateAndMerge(in_samhdr);

    // int ret_r = 0;
    // while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
    //     // QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAG:TYPE:VALUE]…
    //     printf("NAME: %s\n", bam_get_qname(bamdata));                                   //get the query name using the macro
    //     flags = bam_flag2str(bamdata->core.flag);                                       //flags as string
    //     printf("FLG: %d - %s\n", bamdata->core.flag, flags);                            //flag is available in core structure
    //     free((void*)flags);
    //     tidname = sam_hdr_tid2name(in_samhdr, bamdata->core.tid);
    //     printf("RNAME/TID: %d - %s\n", bamdata->core.tid, tidname? tidname: "" );       //retrieves the target name using the value in bam and by referring the header
    //     printf("POS: %d\n", bamdata->core.pos + 1);                          //internally position is 0 based and on text output / SAM it is 1 based
    //     printf("MQUAL: %d\n", bamdata->core.qual);                                      //map quality value

    //     cigar = bam_get_cigar(bamdata);                                                 //retrieves the cigar data
    //     printf("CGR: ");
    //     for (int i = 0; i < bamdata->core.n_cigar; ++i) {                                   //no. of cigar data entries
    //         printf("%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));       //the macros gives the count of operation and the symbol of operation for given cigar entry
    //     }
    //     printf("\nTLEN/ISIZE: %d\n", bamdata->core.isize);

    //     data = bam_get_seq(bamdata);                                                    //get the sequence data
    //     if (bamdata->core.l_qseq != bam_cigar2qlen(bamdata->core.n_cigar, cigar)) {     //checks the length with CIGAR and query
    //         printf("\nLength doesnt matches to cigar data\n");
    //         throw(4);
    //     }

    //     printf("SEQ: ");
    //     for (int i = 0; i < bamdata->core.l_qseq ; ++i) {                                   //sequence length
    //         printf("%c", seq_nt16_str[bam_seqi(data, i)]);                              //retrieves the base from (internal compressed) sequence data
    //     }
    //     printf("\nQUAL: ");
    //     for (int i = 0; i < bamdata->core.l_qseq ; ++i) {
    //         printf("%c", bam_get_qual(bamdata)[i]+33);                                  //retrives the quality value
    //     }
    //     printf("\n\n");
    // }
}

int main(int argc, const char** argv) {
    try {
        parseOptions(argc, argv);
        umicollapse_bf();
    } catch(std::exception& e) {
        cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
        cerr << "Command: ";
        for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
        cerr << endl;
        return 1;
    } catch(int e) {
        if (e != 0) {
            cerr << "Error: Encountered internal Umicollapse-bf exception (#" << e << ")" << endl;
            cerr << "Command: ";
            for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
            cerr << endl;
        }
        return e;
    }
    return 0;
}