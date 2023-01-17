CC:=g++
FLAG:=-Wall -Wextra
OPT:=-O3

debug: main.cpp
	$(CC) $(FLAG) $(OPT) -D DEBUGing -g -c main.cpp
record: main.cpp
	$(CC) $(FLAG) $(OPT) -D RECORD -g -c main.cpp
dep: main.cpp
	$(CC) $(FLAG) $(OPT) -c main.cpp
all: main.o
	$(CC) -o main main.o
plot: plot.txt
	gnuplot plot.gp
clean:main.o
	-rm -rf main.o
