.PHONY: all clean test

TEST_MATRIX_A = ../task1/tests/A.dat
TEST_MATRIX_B = ../task1/tests/B.dat
TEST_MATRIX_C = ../task1/tests/C.dat

SIZE = 1024

all: main gen compare print mapping

main: main.cpp
	mpicxx -O3 -Wall -std=c++98 -o $@ $<

gen: ../task1/gen.cpp
	g++ -O3 -Wall -std=c++98 -o $@ $<
	
mapping: mapping.cpp
	g++ -O3 -Wall -std=c++98 -o $@ $<

compare: ../task1/compare.cpp
	g++ -O3 -Wall -std=c++11 -o $@ $<

print: ../task1/print.cpp
	g++ -O3 -Wall -std=c++11 -o $@ $<

test: $(TEST_MATRIX_A) $(TEST_MATRIX_B) $(TEST_MATRIX_C) main compare
	mpirun -n 8 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C1.dat
	./compare $(TEST_MATRIX_C) C1.dat
	mpirun -n 27 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C2.dat
	./compare $(TEST_MATRIX_C) C2.dat
	mpirun -n 125 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C3.dat
	./compare $(TEST_MATRIX_C) C3.dat
	rm C?.dat

report: main gen mapping
	./gen d $(SIZE) $(SIZE) A.dat
	./gen d $(SIZE) $(SIZE) B.dat
	mpisubmit.bg -n 1 --stdout time$(SIZE)_1n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 8 --stdout time$(SIZE)_8n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 27 --stdout time$(SIZE)_27n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 125 --stdout time$(SIZE)_125n ./main A.dat B.dat tmp.dat
	./mapping > my.map
	./mympisubmit.bg -n 125 --stdout time$(SIZE)_125n_mapped ./main A.dat B.dat tmp.dat

clean:
	rm -rf main gen compare print mapping