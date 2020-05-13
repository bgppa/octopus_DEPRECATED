set title "Prediction until 31.05: Germany, RichardsODE, total infected" font ",11"
set xlabel "days [1 < April < 30]"
set ylabel "cases"
set ytic 15000
set xtic 10
set grid
set yrange[70000:170000]
set xrange[0:61]
set key bottom right
set key font ",10"
plot "germany.txt" with lines lc rgb "blue" lw 2 title "Real data", "worst.txt" with yerrorbars lc rgb "red" title "Uncertainty around the worst possible case", "new_data.txt" with points lc rgb "blue" lw 2 title "New data (not interpolated)"
