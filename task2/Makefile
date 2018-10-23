CXXFLAGS = -O2 -std=c++11
OBJS = main gen

TEST_MATRIX_A = tests/A.dat
TEST_MATRIX_B = tests/B.dat
TEST_MATRIX_C = tests/C.dat

.PHONY: all clean test report

all: $(OBJS)

main: main.cpp
	g++ $(CXXFLAGS) -o $@ $< -lpapi

gen: ../task1/gen.cpp
	g++ $(CXXFLAGS) -o $@ $<

compare: ../task1/compare.cpp
	g++ $(CXXFLAGS) -o $@ $<

test: $(TEST_MATRIX_A) $(TEST_MATRIX_B) $(TEST_MATRIX_C) main compare
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C0.dat 0
	./compare $(TEST_MATRIX_C) C0.dat
	rm C0.dat compare

report: main gen
	./gen f 1000 1000 A.dat
	./gen f 1000 1000 B.dat
	./main A.dat B.dat /dev/null 0 10 > plot_ijk.dat
	./main A.dat B.dat /dev/null 1 10 > plot_ikj.dat
	./main A.dat B.dat /dev/null 2 10 > plot_ikj+.dat
	./gen f 2000 2000 A.dat
	./gen f 2000 2000 B.dat
	./main A.dat B.dat /dev/null 0 5 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 5 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 5 >> plot_ikj+.dat
	./gen f 3000 3000 A.dat
	./gen f 3000 3000 B.dat
	./main A.dat B.dat /dev/null 0 2 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 2 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 2 >> plot_ikj+.dat
	./gen f 4000 4000 A.dat
	./gen f 4000 4000 B.dat
	./main A.dat B.dat /dev/null 0 2 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 2 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 2 >> plot_ikj+.dat
	./gen f 5000 5000 A.dat
	./gen f 5000 5000 B.dat
	./main A.dat B.dat /dev/null 0 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 >> plot_ikj+.dat
	gnuplot plot.gnu -e "set output 'Time.svg';\
						set ylabel 'Time, microseconds';\
						plot 'plot_ijk.dat' using 1:2 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:2 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:2 with linespoints title 'ikj, 72x72'"
	gnuplot plot.gnu -e "set output 'L1_misses.svg';\
						set ylabel 'L1 misses';\
						plot 'plot_ijk.dat' using 1:3 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:3 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:3 with linespoints title 'ikj, 72x72'"
	gnuplot plot.gnu -e "set output 'L2_misses.svg';\
						set ylabel 'L2 misses';\
						plot 'plot_ijk.dat' using 1:4 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:4 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:4 with linespoints title 'ikj, 72x72'"
	gnuplot plot.gnu -e "set output 'Total_cycles.svg';\
						set ylabel 'Total cycles';\
						plot 'plot_ijk.dat' using 1:5 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:5 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:5 with linespoints title 'ikj, 72x72'"
	gnuplot plot.gnu -e "set output 'Flop.svg';\
						set ylabel 'Floating point operation';\
						plot 'plot_ijk.dat' using 1:6 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:6 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:6 with linespoints title 'ikj, 72x72'"
	gnuplot plot.gnu -e "set output 'TLB_misses.svg';\
						set ylabel 'TLB_misses';\
						plot 'plot_ijk.dat' using 1:7 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:7 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:7 with linespoints title 'ikj, 72x72'"
	rm A.dat B.dat plot*.dat

clean:
	rm -rf $(OBJS)
	