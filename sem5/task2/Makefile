CXXFLAGS = -O2 -std=c++11
OBJS = main gen compare

TEST_MATRIX_A = tests/A.dat
TEST_MATRIX_B = tests/B.dat
TEST_MATRIX_C = tests/C.dat

.PHONY: all clean test report plot

all: $(OBJS)

main: main.cpp
	g++ $(CXXFLAGS) -o $@ $< -lpapi

gen: ../task1/gen.cpp
	g++ $(CXXFLAGS) -o $@ $<

compare: ../task1/compare.cpp
	g++ $(CXXFLAGS) -o $@ $<

test: $(TEST_MATRIX_A) $(TEST_MATRIX_B) $(TEST_MATRIX_C) main compare
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C0.dat 0
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C1.dat 1
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C2.dat 2
	./compare $(TEST_MATRIX_C) C0.dat
	./compare $(TEST_MATRIX_C) C1.dat
	./compare $(TEST_MATRIX_C) C2.dat
	rm C?.dat compare

report: main gen
	./gen f 1000 1000 A.dat
	./gen f 1000 1000 B.dat
	./main A.dat B.dat /dev/null 0 > plot_ijk.dat
	./main A.dat B.dat /dev/null 1 > plot_ikj.dat
	./main A.dat B.dat /dev/null 2 > plot_ikj+.dat
	./gen f 2000 2000 A.dat
	./gen f 2000 2000 B.dat
	./main A.dat B.dat /dev/null 0 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 >> plot_ikj+.dat
	./gen f 3000 3000 A.dat
	./gen f 3000 3000 B.dat
	./main A.dat B.dat /dev/null 0 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 >> plot_ikj+.dat
	./gen f 4000 4000 A.dat
	./gen f 4000 4000 B.dat
	./main A.dat B.dat /dev/null 0 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 >> plot_ikj+.dat
	./gen f 5000 5000 A.dat
	./gen f 5000 5000 B.dat
	./main A.dat B.dat /dev/null 0 >> plot_ijk.dat
	./main A.dat B.dat /dev/null 1 >> plot_ikj.dat
	./main A.dat B.dat /dev/null 2 >> plot_ikj+.dat
	make plot
	rm A.dat B.dat plot*.dat
	
plot: 
	gnuplot plot.gnu -e "set output 'Time.png';\
						set ylabel 'Time, microseconds';\
						plot 'plot_ijk.dat' using 1:2 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:2 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:2 with linespoints title 'ikj, 52x52'"
	gnuplot plot.gnu -e "set output 'L1_misses.png';\
						set ylabel 'L1 misses';\
						plot 'plot_ijk.dat' using 1:3 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:3 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:3 with linespoints title 'ikj, 52x52'"
	gnuplot plot.gnu -e "set output 'L2_misses.png';\
						set ylabel 'L2 misses';\
						plot 'plot_ijk.dat' using 1:4 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:4 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:4 with linespoints title 'ikj, 52x52'"
	gnuplot plot.gnu -e "set output 'Total_cycles.png';\
						set ylabel 'Total cycles';\
						plot 'plot_ijk.dat' using 1:5 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:5 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:5 with linespoints title 'ikj, 52x52'"
	gnuplot plot.gnu -e "set output 'Flop.png';\
						set ylabel 'Floating point operation';\
						plot 'plot_ijk.dat' using 1:6 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:6 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:6 with linespoints title 'ikj, 52x52'"
	gnuplot plot.gnu -e "set output 'TLB_misses.png';\
						set ylabel 'TLB misses';\
						plot 'plot_ijk.dat' using 1:7 with linespoints title 'ijk, 32x32',\
							 'plot_ikj.dat' using 1:7 with linespoints title 'ikj, 32x32',\
							 'plot_ikj+.dat' using 1:7 with linespoints title 'ikj, 52x52'"

clean:
	rm -rf $(OBJS)
	
