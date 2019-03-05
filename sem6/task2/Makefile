.PHONY: all clean report

N = 20
K = 7

all: main

main: main.cpp
	mpicxx -O3 -Wall -std=c++98 -o $@ $<

report:
	mpisubmit.pl -n 1 --stdout time$(N)n_$(K)k_1p ./main $(N) $(K)
	mpisubmit.pl -n 2 --stdout time$(N)n_$(K)k_2p ./main $(N) $(K)
	mpisubmit.pl -n 4 --stdout time$(N)n_$(K)k_4p ./main $(N) $(K)
	mpisubmit.pl -n 8 --stdout time$(N)n_$(K)k_8p ./main $(N) $(K)

clean:
	rm -rf main