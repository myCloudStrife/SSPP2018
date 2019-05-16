.PHONY: all clean fullreport report test lint

FLAGS = -O3 -Wall -Werror

N=28
T=8

all: main print random

lint: *.cpp
	cpplint --extensions=c,cc,h,hpp,c++,h++,hh,cu,cpp,hxx,cxx,cuh --filter=-legal/copyright $^

test: main print
	mpirun --oversubscribe -n 1 ./main tests/test1in.data 1 tmp1.data
	mpirun --oversubscribe -n 1 ./main tests/test1in.data 2 tmp2.data
	mpirun --oversubscribe -n 2 ./main tests/test1in.data 1 tmp3.data
	mpirun --oversubscribe -n 4 ./main tests/test1in.data 4 tmp4.data
	./print tmp1.data > tmp1.txt
	./print tmp2.data > tmp2.txt
	./print tmp3.data > tmp3.txt
	./print tmp4.data > tmp4.txt
	diff tmp1.txt tests/test1out.txt
	diff tmp2.txt tests/test1out.txt
	diff tmp3.txt tests/test1out.txt
	diff tmp4.txt tests/test1out.txt
	rm tmp*

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
