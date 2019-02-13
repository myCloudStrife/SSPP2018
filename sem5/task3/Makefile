CXXFLAGS = -O2 -std=c++98

.PHONY: all clean

all: main pthreads

main: main.cpp
	mpicxx $(CXXFLAGS) -o $@ $<
	
pthreads: pthreads.cpp
	g++ -O3 -std=c++11 -o $@ $< -lpthread
	
clean:
	rm -rf main pthreads