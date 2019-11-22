CXX=g++
CXXFLAGS=-DNDEBUG -Wall -O3 -std=c++11 -fopenmp
CXXDEBUGFLAGS=-Wall -Wpedantic -O0 -g -std=c++11 -fopenmp
CPPFLAGS=-I. -I./htslib/htslib
LDFLAGS=-L./htslib
LDLIBS=-lhts

.PHONY: all

all: main

main: main.o variant.o vcf_file.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $<

debug: main.debug

main.debug: main.o.debug variant.o.debug vcf_file.o.debug
	@echo "* Linking $@"
	$(CXX) $(CXXDEBUGFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o.debug: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXDEBUGFLAGS) $(CPPFLAGS) -o $@ -c $<

clean:
	rm -rf *.o *.o.debug main main.debug
