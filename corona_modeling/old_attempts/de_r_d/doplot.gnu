set title "Prediction until 15.05: Germany, RichardODE, deaths" font ",11"
set xlabel "days [1 < April < 30]"
set ylabel "cases"
set ytic 2000
set xtic 5
set grid
set xrange[0:45]
set key left top
set key font ",11"
plot "../datasets/deceased/germany.txt" with lines lc rgb "blue" lw 2 title "Real data", "best.txt" with lines lc rgb "brown" lw 2 title "Best", "worst.txt" with lines lc rgb "red" lw 2 title "Worst", "../datasets/deceased/new_germany.txt" with points lc rgb "blue" lw 2 title "New data (not interpolated)"
