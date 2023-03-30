set xrange [3.5:100]
set yrange [.4:1]
set logscale x
set logscale y
set key outside


plot "t=1.05_m.vs.L_fss.dat" using 1:2 with lp pt 1 lc rgb "black" title "t = 1.05", \
     "t=1.0_m.vs.L_fss.dat" using 1:2 with lp pt 2 lc rgb "black" title "t = 1.0", \
     "t=0.95_m.vs.L_fss.dat" using 1:2 with lp pt 3 lc rgb "black" title "t = .95", \
     "t=0.9_m.vs.L_fss.dat" using 1:2 with lp pt 4 lc rgb "black" title "t = 0.90", \
     "t=0.8_m.vs.L_fss.dat" using 1:2 with lp pt 5 lc rgb "black" title "t = 0.80", \
     "t=0.7_m.vs.L_fss.dat" using 1:2 with lp pt 6 lc rgb "black" title "t = 0.70", \
     "t=0.65_m.vs.L_fss.dat" using 1:2 with lp pt 7 lc rgb "black" title "t = 0.65", \
     "t=0.55_m.vs.L_fss.dat" using 1:2 with lp pt 8 lc rgb "black" title "t = 0.55"

set term png
set output "Lvs.T_FSS.png"
replot
unset output
set term qt