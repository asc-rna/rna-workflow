CXX := g++
CPPFLAGS += -std=c++17
CPPFLAGS += -g -Wall -Wextra -Wconversion -Wno-sign-conversion 
CPPFLAGS += -fsanitize=address,undefined -fno-sanitize-recover 
CPPFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC 
CPPFLAGS += -DCM_DEBUG
REMFLAGS := -O2 -g

SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp, %, ${SRC})
RUN := $(patsubst %.cpp, %_run, ${SRC})
REM := $(patsubst %.cpp, %_rem, ${SRC})

PCHS := $(patsubst %, -include-pch %, $(wildcard ~/cm/*.hpp.gch))


all: $(OBJ)

%: %.cpp
	$(CXX) $(CPPFLAGS) -o $@ $^


in:
	@touch in

.PHONY:run
run: in $(RUN)

%_run: %
	@echo -e "\033[1;34m--- $^ ---\033[0m"
	@time ./$^ < in > $^.out


.PHONY:remote
remote: $(REM)


%_rem: %.cpp
	$(CXX) $(REMFLAGS) -o $(basename $^) $^
