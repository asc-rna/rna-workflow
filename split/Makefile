CXX := g++
CPPFLAGS := -O3 -g

SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp, %, ${SRC})

all: $(OBJ)

%: %.cpp
	$(CXX) $(CPPFLAGS) -o $@ $^
