
CXX = g++

CXXFLAGS = -O3 -I${HOME}/local/include

LDFLAGS = -L${HOME}/local/lib -lhts

all: umicollapse

umicollapse: main.cpp
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDFLAGS) -lpthread

example: example.cpp
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDFLAGS) -lpthread

clean:
	rm -f example umicollapse
