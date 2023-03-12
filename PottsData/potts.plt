set xrange [0:1.5]
set yrange [0:2]

plot "Potts4_heat.dat" using 1:2, \
     "Potts8_heat.dat" using 1:2, \
     "Potts12_heat.dat" using 1:2, \
     "Potts16_heat.dat" using 1:2, \
     "Potts20_heat.dat" using 1:2, \
     "Potts32_heat.dat" using 1:2