/*
example of the logs:

Number of input reads	8532151
Number of removed unmapped reads	0
Number of unremoved reads	8532151
Number of unique alignment positions	1206022
Average number of UMIs per alignment position	5.6625069858
Max number of UMIs over all alignment positions	2146
Number of reads after deduplicating	6485295
*/

#include <cstdio>
#include <cstring>
#include <cctype>
#include <cstdlib>
using namespace std;
inline int read() {
    int x = 0; char c = getchar();
    while (!isdigit(c)) c = getchar();
    while (isdigit(c)) {
        x = x * 10 + c - '0';
        c = getchar();
    }
    return x;
}
inline double read_double() {
    double x = read(), b = 0.1;
    char c = getchar();
    while (!isdigit(c)) c = getchar();
    while (isdigit(c)) {
        x += (c - '0') * b;
        b *= 0.1;
        c = getchar();
    }
    return x;
}
int main(int argc, char *argv[]) {
    if (argc != 3 || strcmp(argv[1], "-h") == 0) {
        fprintf(stderr, "Usage: %s <log prefix> <log num>\n", argv[0]);
        return 1;
    }
    int log_num = atoi(argv[2]);
    char filename[200];
    
    int tot = 0, removed = 0, unremoved = 0, unique = 0, max_umi = 0, remained_num = 0;
    double avg_umi = 0;
    for (int i = 0; i < log_num; i++) {
        sprintf(filename, "%s.%d.log", argv[1], i);  // input file example: xxx.0.log
        freopen(filename, "r", stdin);
        tot += read();
        removed += read();
        unremoved += read();
        int unique_alignment = read();
        unique += unique_alignment;
        avg_umi += read_double() * unique_alignment;
        int tmp = read();
        if (tmp > max_umi) max_umi = tmp;
        remained_num += read();
    }
    printf("Number of input reads: %d\n", tot);
    printf("Number of removed unmapped reads: %d\n", removed);
    printf("Number of unremoved reads: %d\n", unremoved);
    printf("Number of unique alignment positions: %d\n", unique);
    printf("Average number of UMIs per alignment position: %.10f\n", avg_umi / unique);
    printf("Max number of UMIs over all alignment positions: %d\n", max_umi);
    printf("Number of reads after deduplicating: %d\n", remained_num);
    return 0;
}