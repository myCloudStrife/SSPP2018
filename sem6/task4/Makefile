.PHONY: all clean test test_not test_cnot test_adamar lint

FLAGS = -O3 -Wall -Werror

TEST_Q = 7
TEST_K = 1
TEST_L = 2

all: test_my_not

lint: *.cpp *.c *.hpp tests/*.cpp
	cpplint --extensions=c,cc,h,hpp,c++,h++,hh,cu,cpp,hxx,cxx,cuh --filter=-legal/copyright $^

test: random
	make test_not TEST_Q=7 TEST_K=3
	make test_not TEST_Q=8 TEST_K=2
	make test_cnot TEST_Q=5 TEST_K=1 TEST_L=2
	make test_cnot TEST_Q=5 TEST_K=4 TEST_L=3
	make test_adamar TEST_Q=6 TEST_K=1
	make test_adamar TEST_Q=6 TEST_K=6

test_not: test_correct_not test_my_not random
	./random $(TEST_Q) src
	./test_correct_not src $(TEST_K) res
	mpirun --oversubscribe -n 1 ./test_my_not src $(TEST_K) myres 1
	diff res myres
	mpirun --oversubscribe -n 1 ./test_my_not src $(TEST_K) myres 2
	diff res myres
	mpirun --oversubscribe -n 1 ./test_my_not src $(TEST_K) myres 4
	diff res myres
	mpirun --oversubscribe -n 2 ./test_my_not src $(TEST_K) myres 1
	diff res myres
	mpirun --oversubscribe -n 4 ./test_my_not src $(TEST_K) myres 1
	diff res myres
	mpirun --oversubscribe -n 4 ./test_my_not src $(TEST_K) myres 4
	diff res myres
	rm src res myres

test_cnot: test_correct_cnot test_my_cnot random
	./random $(TEST_Q) src
	./test_correct_cnot src $(TEST_K) $(TEST_L) res
	mpirun --oversubscribe -n 1 ./test_my_cnot src $(TEST_K) $(TEST_L) myres 1
	diff res myres
	mpirun --oversubscribe -n 1 ./test_my_cnot src $(TEST_K) $(TEST_L) myres 2
	diff res myres
	mpirun --oversubscribe -n 1 ./test_my_cnot src $(TEST_K) $(TEST_L) myres 4
	diff res myres
	mpirun --oversubscribe -n 2 ./test_my_cnot src $(TEST_K) $(TEST_L) myres 1
	diff res myres
	mpirun --oversubscribe -n 4 ./test_my_cnot src $(TEST_K) $(TEST_L) myres 1
	diff res myres
	mpirun --oversubscribe -n 4 ./test_my_cnot src $(TEST_K) $(TEST_L) myres 4
	diff res myres
	rm src res myres

test_adamar: test_my_adamar random print
	./random $(TEST_Q) src
	./print src > src.txt
	mpirun --oversubscribe -n 1 ./test_my_adamar src $(TEST_K) tmp 1
	mpirun --oversubscribe -n 1 ./test_my_adamar tmp $(TEST_K) myres 1
	./print myres > myres.txt
	diff src.txt myres.txt
	mpirun --oversubscribe -n 2 ./test_my_adamar src $(TEST_K) tmp 2
	mpirun --oversubscribe -n 2 ./test_my_adamar tmp $(TEST_K) myres 2
	./print myres > myres.txt
	diff src.txt myres.txt
	mpirun --oversubscribe -n 4 ./test_my_adamar src $(TEST_K) tmp 1
	mpirun --oversubscribe -n 4 ./test_my_adamar tmp $(TEST_K) myres 1
	./print myres > myres.txt
	diff src.txt myres.txt
	rm src src.txt tmp myres myres.txt


random: random.cpp
	g++ $(FLAGS) -o random random.cpp

print: print.c
	gcc $(FLAGS) -o print print.c

test_correct_not: tests/test_not.cpp
	g++ $(FLAGS) -o $@ $^

test_my_not: tests/test_mynot.cpp
	mpicxx $(FLAGS) -fopenmp -o $@ $^

test_correct_cnot: tests/test_cnot.cpp
	g++ $(FLAGS) -o $@ $^

test_my_cnot: tests/test_mycnot.cpp
	mpicxx $(FLAGS) -fopenmp -o $@ $^

test_my_adamar: tests/test_myadamar.cpp
	mpicxx $(FLAGS) -fopenmp -o $@ $^

clean:
	rm -rf print random test_*
