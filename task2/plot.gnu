set terminal svg size 600,500
set title "Task2"
set xlabel 'Mode'
set xtics ("ijk,32x32" 0, "ikj,32x32" 1, "ikj,45x45" 2)
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "white" behind
