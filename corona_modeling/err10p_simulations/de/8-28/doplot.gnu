#!/bin/bash

echo "
set title 'Germany: prediction until 01.06 using 3 weeks of data.' font ',11'
set xlabel 'days [8 < April < 30]'
set ylabel 'deceased people'
set ytic 3000
set xtic 10
set grid
set xrange[8:62]
set key right bottom
set key font ',12'
plot '../../../datasets/deceased/germany.txt' with points lc rgb 'blue' lw 2 title 'Real data: interpolated from 8.04 to 28.04', 'best.txt' with lines lc rgb 'brown' lw 2 title 'Best', 'worst.txt' with lines lc rgb 'red' lw 2 title 'Worst', 'exp.txt' with lines lc rgb 'purple' lw 2 title 'Expected" | gnuplot -p
