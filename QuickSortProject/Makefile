######################################################################
# Makefile for runSorts
######################################################################
CFLAGS= -Wall -O3 -std=c++11 -fopenmp
CC = g++

######################################################################
runSorts: runSorts.cpp standardSortAlgorithms.o listGenerator.o util.o betterQuick.o betterMerge.o PQuick.o FJQuick.o
	$(CC) $(CFLAGS) -o runSorts runSorts.cpp standardSortAlgorithms.o listGenerator.o util.o betterQuick.o betterMerge.o PQuick.o FJQuick.o

parCode: parallelCode.cpp util.o
	$(CC) $(CFLAGS) -o parCode parallelCode.cpp util.o

testPart: testPart.cpp util.o yourSorts.o
	$(CC) $(CFLAGS) -o testPart testPart.cpp yourSorts.o util.o
	
listGenerator.o: listGenerator.h listGenerator.cpp
	$(CC) $(CFLAGS) -c listGenerator.cpp
util.o: util.h util.cpp
	$(CC) $(CFLAGS) -c util.cpp
 
standardSortAlgorithms.o: standardSortAlgorithms.h standardSortAlgorithms.cpp
	$(CC) $(CFLAGS) -c standardSortAlgorithms.cpp
 
betterMerge.o: betterMerge.h betterMerge.cpp
	$(CC) $(CFLAGS) -c betterMerge.cpp
betterQuick.o: betterQuick.h betterQuick.cpp
	$(CC) $(CFLAGS) -c betterQuick.cpp
PQuick.o: PQuick.h PQuick.cpp
	$(CC) $(CFLAGS) -c PQuick.cpp
FJQuick.o: FJQuick.h FJQuick.cpp
	$(CC) $(CFLAGS) -c FJQuick.cpp
 
#yourSorts.o: yourSorts.h yourSorts.cpp
#	$(CC) $(CFLAGS) -c yourSorts.cpp

######################################################################
