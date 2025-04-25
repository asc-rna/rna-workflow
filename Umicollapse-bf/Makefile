
CXX = g++

CXXFLAGS = -O3 -I${HOME}/local/include

LDFLAGS = -L${HOME}/local/lib -lhts

SOURCES = umicollapse-bf.cpp

TARGET = umicollapse-bf

all: $(TARGET) umicollapse

$(TARGET):$(SOURCES)
	$(CXX) $(SOURCES) -o $(TARGET) $(CXXFLAGS) $(LDFLAGS)

umicollapse: main.cpp
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDFLAGS) -lpthread

clean:
	rm -f $(TARGET) umicollapse
