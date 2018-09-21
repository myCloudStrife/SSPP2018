set terminal svg size 600,500
set output 'plots/plot.svg'
set title "Task1"
set xlabel 'Mode'
set ylabel 'Time, sec'
set xtics ("ijk" 0, "ikj" 1, "kij" 2, "jik" 3, "jki" 4, "kji" 5)
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "white" behind
plot 'plots/plot.dat' with linespoints