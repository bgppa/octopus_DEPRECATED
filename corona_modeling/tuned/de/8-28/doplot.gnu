#!/bin/bash

echo "
set title 'Germany: prediction until 18.05 using 3 weeks of data.' font ',11'
set xlabel 'days [8 < April < 30]'
set ylabel 'deceased people'
set ytic 2000
set xtic 4
set grid
set xrange[8:48]
set key left top
set key font ',13'
plot '../../../datasets/deceased/germany.txt' with points lc rgb 'blue' lw 2 title 'Real data: interpolated from 8.04 to 28.04', 'best.txt' with lines lc rgb 'brown' lw 2 title 'Best scenario. Prob: $1%', 'worst.txt' with lines lc rgb 'red' lw 2 title 'Worst scenario. Prob: $2%', 'map.txt' with lines lc rgb 'grey' lw 2 title 'Most probably, prob: $3%" | gnuplot -p
