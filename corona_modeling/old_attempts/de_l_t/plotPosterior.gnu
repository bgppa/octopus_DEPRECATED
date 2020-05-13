set title "Posterior distribution: Germany, Logistic, total infected" font ",11"
set xlabel "q"
set ylabel "Q"
set ytic 10000
set xtic 0.1
set grid
set key left top
set key font ",11"
plot "posterior.txt" with dots title "Apprx. interpolation error for most of the points: 2%"
