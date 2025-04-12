#include "scanner.hpp"
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>

template <class _BufferReader>
class FastqReader : protected _BufferReader {
  FastqReader() = delete;
  using _BufferReader::_BufferReader;
  using _BufferReader::end;
  using _BufferReader::get;
  using _BufferReader::next;

  bool _is_cr(char c) {
    return c == '\n' || c == '\r';
  }

public:
  char *next_read(char *s) {
    for (int cr_cnt = 0; cr_cnt < 4 && !end(); ) {
      cr_cnt += _is_cr(get());
      *s++ = next();
    }
    if (end()) {
      return nullptr;
    } else {
      *s = '\0';
      return s;
    }
  }
};

int main(int argc, char *argv[]) {
  if (argc < 4 || strcmp(argv[1], "-h") == 0) {
    puts("Usage: par-run <process num> <output prefix> <command-file>");
    exit(1);
  }

  int np = 1;
  sscanf(argv[1], "%d", &np);

  FILE *cmd_f = fopen(argv[3], "r");
  if (cmd_f == nullptr) {
    puts("unable to read command-file");
    return 1;
  }
  cm::scanner<cm::optimal_reader> cmd_sc(cmd_f);

  int np_in[np];
  for (int i = 0; i < np; i++) {
    int pipefd[2];
    if (pipe(pipefd) == -1) {
      perror("pipe2 error");
      return 1;
    }
    np_in[i] = pipefd[1];

    char cmd[4096];
    cmd_sc.next_line(cmd);

    if (fork() == 0) { // child process
      std::cout << "subprocess " << i << std::endl;
      close(pipefd[1]);
      // set stdin as pipe
      dup2(pipefd[0], STDIN_FILENO);
      // set stdout as prefix.i
      char out_filename[4096];
      sprintf(out_filename, "%s.%d", argv[2], i);
      std::cout << "output filename: " << out_filename << std::endl;
      freopen(out_filename, "w", stdout);

      if (execl("/bin/sh", "sh", "-c", cmd, NULL) == -1) {
        perror("unable to start subprocess");
        exit(1);
      }
    }
    else { // parent process
      close(pipefd[0]);
    }
  }

  FastqReader<cm::optimal_reader> sc(stdin);
  // BUFF_SIZE depends on kernel pipe size, and is determined by experiment
  const int BUFF_SIZE = 64 * 1024;
  static char buf[BUFF_SIZE];
  char *p = buf;
  for (int i = 0, c = 0; ; i++) {
    p = sc.next_read(p);
    if (p == nullptr) {
      write(np_in[c], buf, p - buf);
      break;
    }
    
    if ((p - buf) > BUFF_SIZE - 2000) {
      write(np_in[c], buf, p - buf);
      c = (c + 1) % np;
      p = buf;
    }
  }

  for (int i = 0; i < np; i++)
    close(np_in[i]);
  while (wait(NULL) > 0);

  return 0;
}
