OBJS = jacobiIterationSolver.o solveTask5.o main.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -O3 -std=c++11
LFLAGS = -Wall

task5 : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o task5

jacobiIterationSolver.o : jacobiIterationSolver.cpp jacobiIterationSolver.h
	$(CC) $(CFLAGS) jacobiIterationSolver.cpp 

solveTask5.o : solveTask5.cpp executionData.h jacobiIterationSolver.h
	$(CC) $(CFLAGS) solveTask5.cpp

main.o : main.cpp 
	$(CC) $(CFLAGS) main.cpp

clean:
	\rm *.o *~ task5.exe *.dat
