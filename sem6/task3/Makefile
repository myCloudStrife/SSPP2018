.PHONY: all clean report

N = 28
T = 8

all: main

main: main.cpp
	mpicxx -O3 -fopenmp -Wall -std=c++11 -o $@ $<

fullreport:
	make report T=8
	sleep 300
	make report T=4
	sleep 300
	make report T=2
	sleep 300
	make report T=1

report: main
	mpisubmit.pl -p 1 -t $(T) --stdout time$(N)n_$(T)t_1p --stderr errors ./main $(N) 0.01 $(T)
	sleep 300
	mpisubmit.pl -p 2 -t $(T) --stdout time$(N)n_$(T)t_2p --stderr errors ./main $(N) 0.01 $(T)
	sleep 300
	mpisubmit.pl -p 4 -t $(T) --stdout time$(N)n_$(T)t_4p --stderr errors ./main $(N) 0.01 $(T)

clean:
	rm -rf main
