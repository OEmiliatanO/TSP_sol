CC:=g++
FLAG:=-Wall -Wextra

debug:main.cpp
	$(CC) $(FLAG) -D DEBUGing -g -c main.cpp
dep:main.cpp
	$(CC) -c main.cpp
all:main.o
	$(CC) -o main main.o
clean:main.o
	-rm -rf main.o
