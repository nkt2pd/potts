set xrange [0:1.6]
set yrange [0:1.6]

set title "Specific Heat vs. T"
set xlabel "T"
set ylabel "Specific Heat"
set key outside

plot "Potts4_heat.dat" using 1:2 with lp title "L = 4", \
     "Potts12_heat.dat" using 1:2 with lp title "L = 12", \
     "Potts32_heat.dat" using 1:2 with lp title "L = 32"

set term png
set output "PottsHeat.png"
replot
unset output
set term qt