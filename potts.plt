set xrange [0:1.5]
set yrange [0:.15]

plot "Potts8_m.dat" using 1:2 with p, \
     "Potts12_m.dat" using 1:2 with p, \
     "Potts16_m.dat" using 1:2 with p, \
     "Potts20_m.dat" using 1:2 with p