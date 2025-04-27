CXX := g++
CPPFLAGS := -O3 -g

SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp, %, ${SRC})

all: $(OBJ)

split_fastq: split_fastq.cpp
	$(CXX) $(CPPFLAGS) -o $@ $^

merge_split_bam: merge_split_bam.cpp
	$(CXX) $^ -o $@ $(CPPFLAGS) -I${HOME}/local/include -L${HOME}/local/lib -lhts -lpthread -fopenmp
