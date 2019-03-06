.PHONY: all clean report

N = 20
K = 7

all: main

main: main.cpp
	mpicxx -O3 -Wall -std=c++98 -o $@ $<

report4:
	make report N=20
	make report N=24
	make report N=28
	make report N=32

report: main
	mpisubmit.pl -p 1 --stdout time$(N)n_$(K)k_1p --stderr errors ./main $(N) $(K)
	mpisubmit.pl -p 2 --stdout time$(N)n_$(K)k_2p --stderr errors ./main $(N) $(K)
	mpisubmit.pl -p 4 --stdout time$(N)n_$(K)k_4p --stderr errors ./main $(N) $(K)
	mpisubmit.pl -p 8 --stdout time$(N)n_$(K)k_8p --stderr errors ./main $(N) $(K)

clean:
	rm -rf main