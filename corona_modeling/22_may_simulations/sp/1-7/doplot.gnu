#!/bin/bash

echo "
set title 'Spain: prediction until 18.05 using 1 week of data.' font ',11'
set xlabel 'days [1 < April < 30]'
set ylabel 'deceased people'
set ytic 5000
set xtic 4
set grid
set xrange[1:48]
set key right bottom
set key font ',13'
plot '../../../datasets/deceased/spain.txt' with points lc rgb 'blue' lw 2 title 'Real data: interpolated from 1.04 to 7.04', 'best.txt' with lines lc rgb 'brown' lw 2 title 'Best', 'worst.txt' with lines lc rgb 'red' lw 2 title 'Worst', 'exp.txt' with lines lc rgb 'purple' lw 2 title 'Expected" | gnuplot -p
