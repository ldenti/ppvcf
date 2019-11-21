CXX=g++
CXXFLAGS=-DNDEBUG -Wall -O3 -std=c++11 -fopenmp
CPPFLAGS=-I. -I./htslib/htslib
LDFLAGS=-L./htslib
LDLIBS=-lhts

.PHONY: all

all: main

main: main.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $<

clean:
	rm -rf main.o main
