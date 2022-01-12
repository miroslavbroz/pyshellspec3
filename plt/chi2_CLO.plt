#!/usr/bin/gnuplot

set colors classic
set term x11

rad = 180./pi

set xl "B/lambda [cycles]"
set yl "closure phase arg T_3 [deg] (shifted by dataset number)"

set yr [0:30*360]
set ytics 360
set mytics 2
set grid ytics mytics
set zeroaxis

p \
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$6):($12+$15*360) t "synthetic" w p lt 7 pt 1,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$6):($10+$15*360):11 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '{ print $1,$2,$6,$10,ARGIND; print $1,$2,$6,$12,ARGIND; print null; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*360) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$6,$10,ARGIND; print $1,$2,$6,$12,ARGIND; print null; }' *.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*360) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '(FNR==2){ print $1,$2,$6,ARGIND,FILENAME; }' *.cp.syn.dat" u (2.3e8):($4*360):5 not w labels left font "Helvetica,8"

pa -1

set term png small size 1920,1080
set out "chi2_CLO.png"
rep



