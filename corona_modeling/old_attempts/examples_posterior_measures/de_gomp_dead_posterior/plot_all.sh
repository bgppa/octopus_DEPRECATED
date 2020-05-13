#!/bin/bash

for f in ./*.txt; do
	echo "set terminal png
	set output '$f.png'
	set title 'Bayesian posterior distribution'
	set ytic 1000
	set xtic 0.05
	set yrange [1000:15000]
	set xrange [0:0.2]
	set grid
	set xlabel 'r (growth speed)'
	set ylabel 'Q (max infected)'
	plot '$f' with dots" | gnuplot
done
