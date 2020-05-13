#!/bin/bash
let l=$1-1

echo "
set title 'Uncertainty until 04.05: Germany, Richard ODE, #deaths' font ',11'
set xlabel 'days [1 < April < 30]'
set ylabel 'cases'
set ytic 2000
set xtic 4
set grid
set xrange[0:$2]
set key left top
set key font ',10'
plot '../datasets/deceased/germany.txt' with lines lc rgb 'blue' lw 2 title 'Real data: interpolated from day 0 to $l', '$1-$2-best.txt' with lines lc rgb 'brown' lw 2 title 'Best scenario', '$1-$2-worst.txt' with lines lc rgb 'red' lw 2 title 'Worst scenario.', '$1-$2-expected.txt' with lines lw 2 title 'Second best scenario'" | gnuplot -p
