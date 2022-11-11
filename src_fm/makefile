EXECUTABLE=gMap_fm

CXX = g++

CFLAGS = -lpthread -msse4 -std=c++11

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
endif
ifeq ($(UNAME_S),Darwin)
	CFLAGS += -I/opt/local/include
endif

$(EXECUTABLE): Main.cpp basic.o B+tree.o build.o ExactMatch.o print.o read.o 
	$(CXX) $^ -o $@ $(CFLAGS) 

basic.o: basic.cpp
	$(CXX) -c $< -o $@

B+tree.o: B+tree.cpp 
	$(CXX) -c $< -o $@

build.o:build.cpp 
	$(CXX) -c $< -o $@

ExactMatch.o: ExactMatch.cpp 
	$(CXX) -c $< -o $@

print.o: print.cpp 
	$(CXX) -c $< -o $@

read.o: read.cpp
	$(CXX) -c $< -o $@
 
.PHONY : clean

clean:
	rm *.o $(EXECUTABLE)
