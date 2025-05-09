/*
example of the logs:

HISAT2 summary stats:
	Total reads: 76051745
		Aligned 0 time: 6445095 (8.47%)
		Aligned 1 time: 60736402 (79.86%)
		Aligned >1 times: 8870248 (11.66%)
	Overall alignment rate: 91.53%
*/

#include <cstdio>
#include <cstring>
#include <cctype>
#include <cstdlib>
using namespace std;
inline int read() {
    int x = 0; char c = getchar();
    while (!isdigit(c) && c != '.') c = getchar();
    while (isdigit(c) || c == '.') {
        if (c != '.') x = x * 10 + c - '0';
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
    int tot = 0, aligned_0 = 0, aligned_1 = 0, aligned_2 = 0;
    for (int i = 0; i < log_num; i++) {
        sprintf(filename, "%s.%d.summary", argv[1], i);  // input file example: xxx.0.summary
        freopen(filename, "r", stdin);
        read();
        tot += read();
        read();
        aligned_0 += read();
        read(); read();
        aligned_1 += read();
        read(); read();
        aligned_2 += read();
    }
    printf("HISAT2 summary stats:\n");
    printf("\tTotal reads: %d\n", tot);
    printf("\t\tAligned 0 time: %d (%.2f%%)\n", aligned_0, 100.0 * aligned_0 / tot);
    printf("\t\tAligned 1 time: %d (%.2f%%)\n", aligned_1, 100.0 * aligned_1 / tot);
    printf("\t\tAligned >1 times: %d (%.2f%%)\n", aligned_2, 100.0 * aligned_2 / tot);
    printf("\tOverall alignment rate: %.2f%%\n", 100.0 * (aligned_1 + aligned_2) / tot);
    return 0;
}