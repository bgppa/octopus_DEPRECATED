#!/bin/bash

echo "
set title 'France: prediction until 01.06 using 3 weeks of data.' font ',11'
set xlabel 'days [3 < April < 30]'
set ylabel 'deceased people'
set ytic 5000
set xtic 4
set grid
set xrange[3:62]
set key right bottom
set key font ',13'
plot '../../../datasets/deceased/france.txt' with points lc rgb 'blue' lw 2 title 'Real data: interpolated from 3.04 to 23.04', 'best.txt' with lines lc rgb 'brown' lw 2 title 'Best', 'worst.txt' with lines lc rgb 'red' lw 2 title 'Worst', 'exp.txt' with lines lc rgb 'purple' lw 2 title 'Expected" | gnuplot -p
