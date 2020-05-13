#!/bin/bash

for f in ./*.txt; do
	echo "set terminal png
	set output '$f.png'
	set title 'Bayesian posterior distribution'
	set ytic 5000
	set xtic 0.1
	set yrange[18000:40000]
	set xrange[0:0.5]
	set grid
	set xlabel 'q (growth speed)'
	set ylabel 'Q (max infected)'
	plot '$f' with dots title 'DE gompertz days 10-32'" | gnuplot
done
