set xrange [0:100]
set yrange [0:1]
set logscale x
set logscale y


plot for [i=6:9] "t=0.".i."_m.vs.L_fss.dat" using 1:2 with lp pt i lc rgb "black", \
    "t=1.0_m.vs.L_fss.dat" using 1:2 with lp pt 10 lc rgb "black"