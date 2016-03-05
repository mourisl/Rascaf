CXX = g++
CXXFLAGS= -Wformat -O3
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS =

all: rascaf join

rascaf: main.o 
	if [ ! -f ./samtools-0.1.19/libbam.a ] ; \
	then \
		cd samtools-0.1.19 ; make ;\
	fi ; 
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

join: join.o
	$(CXX) -o rascaf-join $(LINKPATH) $(CXXFLAGS) $(OBJECTS) join.o $(LINKFLAGS)
	
main.o: main.cpp alignments.hpp blocks.hpp scaffold.hpp support.hpp genome.hpp KmerCode.hpp defs.h ContigGraph.hpp
join.o: join.cpp alignments.hpp blocks.hpp support.hpp genome.hpp KmerCode.hpp defs.h ContigGraph.hpp

clean:
	rm -f *.o *.gch rascaf rascaf-join
