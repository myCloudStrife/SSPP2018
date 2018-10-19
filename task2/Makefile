CXXFLAGS = -O2 -std=c++11
OBJS = main gen

TEST_MATRIX_A = ../task1/tests/A.dat
TEST_MATRIX_B = ../task1/tests/B.dat
TEST_MATRIX_C = ../task1/tests/C.dat

SIZE = 500

.PHONY: all clean test report

all: $(OBJS)

main: main.cpp
	g++ $(CXXFLAGS) -o $@ $< /usr/local/lib/libpapi.a

gen: ../task1/gen.cpp
	g++ $(CXXFLAGS) -o $@ $<

compare: ../task1/compare.cpp
	g++ $(CXXFLAGS) -o $@ $<

test: $(TEST_MATRIX_A) $(TEST_MATRIX_B) $(TEST_MATRIX_C) main compare
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C0.dat 0
	./compare $(TEST_MATRIX_C) C0.dat
	rm C0.dat compare

report: main gen
	./gen f ${SIZE} ${SIZE} A.dat
	./gen f ${SIZE} ${SIZE} B.dat
	./main A.dat B.dat /dev/null 0 10 > plot_${SIZE}_time.dat
	./main A.dat B.dat /dev/null 1 10 > plot_${SIZE}_L1.dat
	./main A.dat B.dat /dev/null 2 10 > plot_${SIZE}_L2.dat
	./main A.dat B.dat /dev/null 3 10 > plot_${SIZE}_totcyc.dat
	gnuplot plot.gnu -e "set output 'Time.svg';\
						set ylabel 'Time';\
						plot 'plot_${SIZE}_time.dat' with linespoints title '${SIZE}x${SIZE} matrix'"
	gnuplot plot.gnu -e "set output 'L1_misses.svg';\
						set ylabel 'L1 misses';\
						plot 'plot_${SIZE}_L1.dat' with linespoints title '${SIZE}x${SIZE} matrix'"
	gnuplot plot.gnu -e "set output 'L2_misses.svg';\
						set ylabel 'L2 misses';\
						plot 'plot_${SIZE}_L2.dat' with linespoints title '${SIZE}x${SIZE} matrix'"
	gnuplot plot.gnu -e "set output 'Total_cycles.svg';\
						set ylabel 'Total cycles';\
						plot 'plot_${SIZE}_totcyc.dat' with linespoints title '${SIZE}x${SIZE} matrix'"
	rm A.dat B.dat plot*.dat

clean:
	rm -rf $(OBJS)
	