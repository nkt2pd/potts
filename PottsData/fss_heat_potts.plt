set xrange [0:1.5]
set yrange [0:2]

set key outside

plot "Potts4_heat.dat" using 1:2 with lp pt 1 lc rgb "black" title "L = 4", \
     "Potts12_heat.dat" using 1:2 with lp pt 1 lc rgb "black" title "L = 12", \
     "Potts32_heat.dat" using 1:2 with lp pt 1 lc rgb "black" title "L = 32", \
     "Potts48_heat.dat" using 1:2 with lp pt 1 lc rgb "black" title "L = 48", \
     "Potts72_heat.dat" using 1:2 with lp pt 1 lc rgb "black" title "L = 72", \

set term png
set output "Tvs.C_FSS.png"
replot
unset output
set term qt