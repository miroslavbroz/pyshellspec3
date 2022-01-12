#!/usr/bin/gnuplot

set colors classic
set term x11

rad = 180./pi

set xl "B/lambda [cycles]"
set yl "triple product amplitude |T_3| [] (shifted by dataset number)"

shift = 2
set yr [0:60]
set ytics shift
set grid ytics ytics
set zeroaxis

p \
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$6):($9+$15*shift)   t "synthetic" w p lt 7 pt 1,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$6):($7+$15*shift):8 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '{ print $1,$2,$6,$7,ARGIND; print $1,$2,$6,$9,ARGIND; print null; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*shift) t "residua"  w l lt 1 lw 3,\
  "<awk '($(NF-1)+0>100){ print $1,$2,$6,$7,ARGIND; print $1,$2,$6,$9,ARGIND; print null; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*shift) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '(FNR==2){ print $1,$2,$7,ARGIND,FILENAME; }' *.cp.syn.dat" u (2.3e8):($4*shift):5 not w labels left font "Helvetica,10"
pa -1

set term png small
set out "chi2_T3.png"
rep



