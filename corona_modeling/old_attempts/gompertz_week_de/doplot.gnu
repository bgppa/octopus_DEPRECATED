#!/bin/bash

echo "
set title 'Weekly prediction: Germany, Gompertz, deceased' font ',11'
set xlabel 'days [1 < April < 30]'
set ylabel 'cases'
set ytic 2000
set xtic 6
set grid
set xrange[0:$2]
set key left top
set key font ',11'
plot '../datasets/deceased/complete_germany.txt' with lines lc rgb 'blue' lw 2 title 'Real data: interpolated from day 0 to $1', '$1-$2-best.txt' with lines lc rgb 'brown' lw 2 title 'Best', '$1-$2-worst.txt' with lines lc rgb 'red' lw 2 title 'Worst'" | gnuplot -p
