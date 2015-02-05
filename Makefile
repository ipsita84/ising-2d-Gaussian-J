# Creates 4 binaries by compiling 4 different source files:
# normal:E-vs-beta-sim-anl.cc
# replicaA: EmA-vs-beta.cc
# replicaB: EmB-vs-beta.cc
# mutualinfo: Mutual-info-vs-beta.cc

# Example:
# make normal will compile E-vs-beta-Normal.cc and generate the
# executable "normal", which can then be run as ./normal and so on.

# Special commands:
# "make clean" will remove all generated binaries and .o object
# files;
# "make all" will compile all 4 source files and generate all
# 4 binaries in one step.

all: disorder normal replicaA replicaB mutualinfo

disorder: disorder-gauss.cc
	g++ -Wall -O3 disorder-gauss.cc -o disorder

normal: E-vs-beta-Normal.cc
	g++ -Wall -O3 E-vs-beta-sim-anl.cc -o normal
	
replicaA: EmA-vs-beta.cc
	g++ -Wall -O3 EmB-vs-beta.cc -o replicaA

replicaB: EmB-vs-beta.cc
	g++ -Wall -O3 EmA-vs-beta.cc -o replicaB

mutualinfo: Mutual-info-vs-beta.cc
	g++ -Wall -O3 \
	`pkg-config --cflags --libs gsl tabdatrw-0.4 interp2dpp` \
#	g++ -Wall -O3 \
#	`pkg-config --cflags --libs gsl tabdatrw interp2dpp` \
#	Mutual-info-vs-beta.cc -o mutualinfo
	
.PHONY: clean

clean:
	rm -f normal replicaA replicaB mutualinfo *.o 
