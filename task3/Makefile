CXXFLAGS = -O2 -std=c++11

.PHONY: all clean

all: main

main: main.cpp
	mpicxx $(CXXFLAGS) -o $@ $<
	
clean:
	rm -rf main