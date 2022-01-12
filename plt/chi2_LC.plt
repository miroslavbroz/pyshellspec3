#!/usr/bin/gnuplot

set colors classic
set term x11

shift = 0.2

set xl "JD - 2400000"
set yl "magnitude [mag]"

set yr [:] reverse
#set ytics shift
set grid ytics
set key left

p \
  "<awk '{ print $1,$4,ARGIND; print $1,$6,ARGIND; print null; }' *.lc.syn.dat" u ($1-2400000):($2+$3*shift)   t "residua"   w l lt 1 lw 3,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.lc.syn.dat"             u ($1-2400000):($4+$9*shift)   not           w l lt 3,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.lc.syn.dat"             u ($1-2400000):($4+$9*shift):5 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.lc.syn.dat"             u ($1-2400000):($6+$9*shift)   t "synthetic" w lp pt 1 lt 7
pa -1

set term png small
set out "chi2_LC.png"
rep

q

  "<awk '(FNR==2){ print $1,$2,$6,ARGIND,FILENAME; }' *.cp.syn.dat" u (2.3e8):($4*360):5 not w labels left font "Helvetica,8"
  "<awk '($NF+0>100){ print $1,$2,$6,$10,ARGIND; print $1,$2,$6,$12,ARGIND; print null; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*360) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\


