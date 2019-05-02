.PHONY: all clean fullreport report test lint

FLAGS = -O3 -Wall -Werror

N=28
T=8

all: main print random

lint: *.cpp
	cpplint --extensions=c,cc,h,hpp,c++,h++,hh,cu,cpp,hxx,cxx,cuh --filter=-legal/copyright $^

test:

random: ../task4/random.cpp
	g++ $(FLAGS) -o $@ $^

print: ../task4/print.c
	gcc $(FLAGS) -o $@ $^

main: main.cpp
	mpicxx $(FLAGS) -fopenmp -o $@ $^

fullreport:
	make report T=8
	sleep 300
	make report T=4
	sleep 300
	make report T=2
	sleep 300
	make report T=1

report: main
	mpisubmit.pl -p 1 -t $(T) --stdout time$(N)n_$(T)t_1p --stderr errors ./main $(N) $(T)
	sleep 300
	mpisubmit.pl -p 2 -t $(T) --stdout time$(N)n_$(T)t_2p --stderr errors ./main $(N) $(T)
	sleep 300
	mpisubmit.pl -p 4 -t $(T) --stdout time$(N)n_$(T)t_4p --stderr errors ./main $(N) $(T)


clean:
	rm -rf print random main
