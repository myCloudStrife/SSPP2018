CXXFLAGS = -O2 -std=c++11
OBJS = main gen print compare

TEST_MATRIX_A = tests/A.dat
TEST_MATRIX_B = tests/B.dat
TEST_MATRIX_C = tests/C.dat

SIZE = 500
REPORT_MATRIX_A = tmpA.dat
REPORT_MATRIX_B = tmpB.dat
REPORT_SCRIPT = report.gnu
REPORT_PLOT_DATA = plot.data
REPORT_OUTPUT = plot.svg

.PHONY: all test report clean

all: $(OBJS)

main: main.cpp
	g++ $(CXXFLAGS) -o $@ $<

gen: gen.cpp
	g++ $(CXXFLAGS) -o $@ $<
	
print: print.cpp
	g++ $(CXXFLAGS) -o $@ $<
	
compare: compare.cpp
	g++ $(CXXFLAGS) -o $@ $<

test: $(TEST_MATRIX_A) $(TEST_MATRIX_B) $(TEST_MATRIX_C) main compare
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C0.dat 0
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C1.dat 1
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C2.dat 2
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C3.dat 3
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C4.dat 4
	./main $(TEST_MATRIX_A) $(TEST_MATRIX_B) C5.dat 5
	./compare $(TEST_MATRIX_C) C0.dat
	./compare $(TEST_MATRIX_C) C1.dat
	./compare $(TEST_MATRIX_C) C2.dat
	./compare $(TEST_MATRIX_C) C3.dat
	./compare $(TEST_MATRIX_C) C4.dat
	./compare $(TEST_MATRIX_C) C5.dat
	rm C?.dat 

report: main gen
	./gen d ${SIZE} ${SIZE} $(REPORT_MATRIX_A)
	./gen d ${SIZE} ${SIZE} $(REPORT_MATRIX_B)
	echo "# mode time" > $(REPORT_PLOT_DATA)
	./main $(REPORT_MATRIX_A) $(REPORT_MATRIX_B) tmp.dat 0 10 >> $(REPORT_PLOT_DATA)
	./main $(REPORT_MATRIX_A) $(REPORT_MATRIX_B) tmp.dat 1 10 >> $(REPORT_PLOT_DATA)
	./main $(REPORT_MATRIX_A) $(REPORT_MATRIX_B) tmp.dat 2 10 >> $(REPORT_PLOT_DATA)
	./main $(REPORT_MATRIX_A) $(REPORT_MATRIX_B) tmp.dat 3 10 >> $(REPORT_PLOT_DATA)
	./main $(REPORT_MATRIX_A) $(REPORT_MATRIX_B) tmp.dat 4 10 >> $(REPORT_PLOT_DATA)
	./main $(REPORT_MATRIX_A) $(REPORT_MATRIX_B) tmp.dat 5 10 >> $(REPORT_PLOT_DATA)
	gnuplot $(REPORT_SCRIPT) -e "set output '$(REPORT_OUTPUT)'; plot '$(REPORT_PLOT_DATA)' with linespoints"
	rm tmp.dat $(REPORT_MATRIX_A) $(REPORT_MATRIX_B)

clean:
	rm -rf $(OBJS) $(REPORT_PLOT_DATA) $(REPORT_OUTPUT)
