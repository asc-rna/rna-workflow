
CXX = g++

CXXFLAGS = -O3 -I${HOME}/local/include

LDFLAGS = -L${HOME}/local/lib -lhts

SOURCES = umicollapse-bf.cpp

TARGET = umicollapse-bf

all: $(TARGET)

$(TARGET):$(SOURCES)
	$(CXX) $(SOURCES) -o $(TARGET) $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f $(TARGET)
