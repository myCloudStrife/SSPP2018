CXXFLAGS = -O3 -Wall -std=c++11

TEST_MATRIX_A = tests/A.dat
TEST_MATRIX_B = tests/B.dat
TEST_MATRIX_C = tests/C.dat

SIZE_M = 512
SIZE_N = 512

.PHONY: all clean test report full_report

all: main gen compare print mapping

main: main.cpp
	mpicxx -O3 -Wall -std=c++98 -o $@ $<

gen: ../task1/gen.cpp
	g++ -O3 -Wall -std=c++98 -o $@ $<
	
mapping: mapping.cpp
	g++ -O3 -Wall -std=c++98 -o $@ $<

compare: ../task1/compare.cpp
	g++ $(CXXFLAGS) -o $@ $<

print: ../task1/print.cpp
	g++ $(CXXFLAGS) -o $@ $<

test: $(TEST_MATRIX_A) $(TEST_MATRIX_B) $(TEST_MATRIX_C) main compare
	mpirun -n 1 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C1.dat
	mpirun -n 2 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C2.dat
	mpirun -n 3 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C3.dat
	mpirun -n 4 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C4.dat
	mpirun -n 50 ./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C5.dat
	./compare $(TEST_MATRIX_C) C1.dat
	./compare $(TEST_MATRIX_C) C2.dat
	./compare $(TEST_MATRIX_C) C3.dat
	./compare $(TEST_MATRIX_C) C4.dat
	./compare $(TEST_MATRIX_C) C5.dat
	rm C?.dat

full_report:
	echo "Sorry, it doesn't work on Bluegene, try separately"
	false
	make report SIZE_M=512 SIZE_N=512
	make report SIZE_M=1024 SIZE_N=1024
	make report SIZE_M=2048 SIZE_N=2048
	make report SIZE_M=4096 SIZE_N=4096
	make report SIZE_M=4096 SIZE_N=1024
	make report SIZE_M=1024 SIZE_N=4096

report: main gen mapping
	./gen d $(SIZE_M) $(SIZE_N) A.dat
	./gen d $(SIZE_N) 1 B.dat
	mpisubmit.bg -n 1 --stdout time$(SIZE_M)x$(SIZE_N)_1n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 32 --stdout time$(SIZE_M)x$(SIZE_N)_32n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 64 --stdout time$(SIZE_M)x$(SIZE_N)_64n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 128 --stdout time$(SIZE_M)x$(SIZE_N)_128n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 256 -w 00:00:45 --stdout time$(SIZE_M)x$(SIZE_N)_256n ./main A.dat B.dat tmp.dat
	mpisubmit.bg -n 512 -w 00:00:45 --stdout time$(SIZE_M)x$(SIZE_N)_512n ./main A.dat B.dat tmp.dat
	./mapping > my.map
	./mympisubmit.bg -n 512 -w 00:00:45 --stdout time$(SIZE_M)x$(SIZE_N)_512n_mapped ./main A.dat B.dat tmp.dat
	


clean:
	rm -rf main gen compare print mapping
