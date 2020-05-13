set title "Prediction until 31.05: Germany, Logistic, NEW DAILY cases of death" font ",11"
set xlabel "days [1 < April < 30]"
set ylabel "cases"
set ytic 100
set xtic 10
set grid
set xrange[0:61]
set key font ",11"
plot "peak_germany.txt" with lines lc rgb "blue" title "Real data", "best.txt" with lines lc rgb "brown" lw 2 title "Best", "worst.txt" with lines lc rgb "red" lw 2 title "Worst", "new_data.txt" with points lc rgb "blue" lw 2 title "New data (not interpolated)"
