OBJS = solveTask6.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -O3 -std=c++11 -msse3
LFLAGS = -Wall

task678 : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o task678

memAlign.o: memAlignOS.h
	$(CC) $(CFLAGS) memAlignOS.cpp memAlignOS.h

solveTask6.o : memAlign.o memAlignOS.h
	$(CC) $(CFLAGS) solveTask6.cpp memAlignOS.h

main.o : main.cpp 
	$(CC) $(CFLAGS) main.cpp

clean:
	\rm *.o *~ task5.exe *.dat
