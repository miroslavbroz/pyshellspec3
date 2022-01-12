#!/usr/bin/gnuplot

set colors classic
set term x11

rad = 180./pi

set xl "B/lambda [cycles]"
set yl "visibility phase arg V [] (shifted by dataset number)"

#set xr [:1.65e8]
#set yr [0:40]
set ytics 360
set grid ytics
set zeroaxis

p \
  "<awk '{ print $1,$2,$4,$9,ARGIND; print $1,$2,$4,$11,ARGIND; print null; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*360-360) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$4,$9,ARGIND; print $1,$2,$4,$11,ARGIND; print null; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*360-360) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$4):($11+$17*360-360)   t "synthetic" w p   lt 7 pt 1,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$4):($9 +$17*360-360):10 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '(FNR==2){ print $1,$2,$4,$9,ARGIND,FILENAME; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*360-360):(sprintf("   %s",stringcolumn(6))) not w labels left font "Helvetica,10"

pa -1

set term png small
set out "chi2_VPHI.png"
rep



