EXECUTABLE=gMap

CXX = g++

CFLAGS = -lpthread -msse4 -std=c++11

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
endif
ifeq ($(UNAME_S),Darwin)
	CFLAGS += -I/opt/local/include
endif

$(EXECUTABLE): Main.cpp basic.o bitVec.o BplusTreeBit.o countMap.o DMindex.o exactMatchFMindex.o FreqSeeds.o Hash.o hashT.o inputFMindex.o inputRef.o Locate.o posList.o DMindex_Oper.o 
	$(CXX) $^ -o $@ $(CFLAGS) 

basic.o: basic.cpp
	$(CXX) -c $< -o $@

bitVec.o: bitVec.cpp #inputFMindex.h DMindex_Oper.h
	$(CXX) -c $< -o $@

BplusTreeBit.o: BplusTreeBit.cpp
	$(CXX) -c $< -o $@

countMap.o: countMap.cpp #seeding.h
	$(CXX) -c $< -o $@

DMindex_Oper.o: DMindex_Oper.cpp
	$(CXX) -c $< -o $@

DMindex.o: DMindex.cpp 
	$(CXX) -c $< -o $@

exactMatchFMindex.o:exactMatchFMindex.cpp 
	$(CXX) -c $< -o $@

FreqSeeds.o:FreqSeeds.cpp 
	$(CXX) -c $< -o $@

Hash.o:Hash.cpp 
	$(CXX) -c $< -o $@

hashT.o: hashT.cpp 
	$(CXX) -c $< -o $@

inputFMindex.o: inputFMindex.cpp 
	$(CXX) -c $< -o $@

inputRef.o:inputRef.cpp 
	$(CXX) -c $< -o $@

Locate.o: Locate.cpp 
	$(CXX) -c $< -o $@

posList.o: posList.cpp 
	$(CXX) -c $< -o $@

.PHONY : clean

clean:
	rm *.o $(EXECUTABLE)
