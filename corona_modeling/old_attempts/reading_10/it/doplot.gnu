#!/bin/bash
let l=$1-1

echo "
set title 'Uncertainty until 05.05: Italy, Richard ODE, #deaths' font ',11'
set xlabel 'days [7 < April < 30]'
set ylabel 'cases'
set ytic 2000
set xtic 4
set grid
set key left top
set key font ',10'
plot '../../datasets/deceased/italy.txt' with points lc rgb 'blue' lw 2 title 'Real data: interpolated from day 7 to 16', 'best.txt' with lines lc rgb 'brown' lw 2 title 'Best scenario', 'worst.txt' with lines lc rgb 'red' lw 2 title 'Worst scenario: 99% confidence interval.', 'worst80.txt' with lines lc rgb 'purple' lw 2 title 'Worst scenario: 80% confidence interval'" | gnuplot -p
