set title "Posterior measure: Italy, Richards ODE, total infected" font ",13"
set xlabel "q"
set ylabel "Q"
set zlabel "v"
set ytic 15000
set xtic 0.1
set ztic 0.1
set grid
set key right top
set key font ",10"
splot "posterior.txt" with dots title "Apprx. interpolation error for most of the points: 3%"
