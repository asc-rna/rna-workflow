#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 3 || strcmp(argv[1], "-h") == 0) {
        puts("Usage: split <output prefix> <split num>");
        exit(1);
    }

    int file_num = stoi(argv[2]);
    FILE *fp[23];
    for (int i = 0; i < file_num; ++i) {
        char filename[100];
        sprintf(filename, "%s.%d.fastq_cut", argv[1], i);
        fp[i] = fopen(filename, "w");
        if (fp[i] == nullptr) {
            perror("fopen error");
            return 1;
        }
    }

    static char buff[10000];
    int cnt = 0, cur = 0;
    while (true) {
        if (fgets(buff, sizeof(buff), stdin) == NULL) break;
        string line(buff);
        if (line.empty()) continue;
        fprintf(fp[cur], "%s", line.c_str());
        ++cnt;
        if (cnt == 1000) {
            cnt = 0;
            cur = (cur + 1) % file_num;
        }
    }
    return 0;
}