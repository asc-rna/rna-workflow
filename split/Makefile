CXX := g++
CPPFLAGS := -O2 -g

SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp, %, ${SRC})

all: $(OBJ)

%: %.cpp
	$(CXX) $(CPPFLAGS) -o $@ $^
