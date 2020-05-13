#!/bin/bash

for f in ./*.txt; do
	echo "set terminal png
	set output '$f.png'
	set title 'Bayesian posterior distribution'
	set ytic 2500
	set xtic 0.1
	set xrange [0:0.3]
	set yrange [1:22000]
	set grid
	set xlabel 'r (growth speed)'
	set ylabel 'Q (max infected)'
	plot '$f' with dots title 'gompertz toy model'" | gnuplot
done
