.PHONY: all clean

N = 20
K = 7

all: main

main: main.cpp
	g++ -O3 -Wall -fopenmp -std=c++11 -o $@ $<

report: main
	bsub -o time_n$(N)_k$(K)_t1 ./main $(N) $(K) 1
	bsub -o time_n$(N)_k$(K)_t2 ./main $(N) $(K) 2
	bsub -o time_n$(N)_k$(K)_t4 ./main $(N) $(K) 4
	bsub -o time_n$(N)_k$(K)_t8 ./main $(N) $(K) 8
	bsub -o time_n$(N)_k$(K)_t16 ./main $(N) $(K) 16
	bsub -o time_n$(N)_k$(K)_t32 ./main $(N) $(K) 32


clean:
	rm -rf main
