.PHONY: all clean

all: main

main: main.cpp
	g++ -O3 -Wall -fopenmp -std=c++11 -o $@ $<

clean:
	rm -rf main
