#!/bin/bash

for f in ./*.txt; do
	echo "set terminal png
	set output '$f.png'
	set title 'Bayesian posterior distribution'
	set ytic 1000
	set xtic 0.05
	set xrange [0:0.2]
	set yrange [4000:10000]
	set grid
	set xlabel 'q (growth speed)'
	set ylabel 'Q (max infected)'
	plot '$f' with dots title 'DE Gompertz posterior $f-32'" | gnuplot
done
