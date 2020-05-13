#!/bin/bash
let l=$1-1

echo "
set title '10 days prediction: Italy, Richard, deceased' font ',11'
set xlabel 'days [1 < April < 30]'
set ylabel 'cases'
set ytic 2000
set xtic 4
set grid
set xrange[0:$2]
set key left top
set key font ',11'
plot '../datasets/deceased/italy.txt' with lines lc rgb 'blue' lw 2 title 'Real data: interpolated from day 0 to $l', '$1-$2-best.txt' with points lc rgb 'brown' lw 1 title 'Best scenario. Fitting error: $3%, Prob: $4%', '$1-$2-worst.txt' with points lc rgb 'red' lw 1 title 'Worst scenario. Fitting error: $5%, Prob: $6%', '$1-$2-expected.txt' with lines title 'Expected' lc rgb 'violet' lw 2" | gnuplot -p

