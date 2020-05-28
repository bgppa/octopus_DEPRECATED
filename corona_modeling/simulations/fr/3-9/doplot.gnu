#!/bin/bash

echo "
set title 'France: prediction until 18.05 using 1 week of data.' font ',11'
set xlabel 'days [3 < April < 30]'
set ylabel 'deceased people'
set ytic 20000
set xtic 4
set grid
set xrange[3:48]
set key left top
set key font ',13'
plot '../../../datasets/deceased/france.txt' with points lc rgb 'blue' lw 2 title 'Real data: interpolated from 3.04 to 9.04', 'best.txt' with lines lc rgb 'brown' lw 2 title 'Best', 'worst.txt' with lines lc rgb 'red' lw 2 title 'Worst', 'exp.txt' with lines lc rgb 'purple' lw 2 title 'Expected" | gnuplot -p
