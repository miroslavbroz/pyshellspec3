#!/usr/bin/gnuplot

set colors classic
set term x11

rad = 180./pi

set xl "B/lambda [cycles]"
set yl "visibility amplitude |V| [] (shifted by dataset number)"

#set xr [:1.15e8]
#set yr [0:40]
set ytics 1
set grid ytics
set zeroaxis

p \
  "<awk '{ print $1,$2,$4,$6,ARGIND; print $1,$2,$4,$8,ARGIND; print null; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5) t "residua"  w l lt 1 lw 3,\
  "<awk '($(NF-1)+0>100){ print $1,$2,$4,$6,ARGIND; print $1,$2,$4,$8,ARGIND; print null; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '($(NF-1)+0>1e4){ print $1,$2,$4,$6,ARGIND; print $1,$2,$4,$8,ARGIND; print null; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5) t "chi^2 > 10^4" w p lc 'red' pt 2 ps 3.5 lw 1,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$4):($8+$17)   t "synthetic" w p   lt 7 pt 1,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$4):($6+$17):7 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '(FNR==2){ print $1,$2,$4,ARGIND,FILENAME; }' *.vis.syn.dat" u (sqrt($1**2+$2**2)/$3):($3+$4):(sprintf("   %s",stringcolumn(5))) not w labels left font "Helvetica,10"

pa -1

set term png small
set out "chi2_VAMP.png"
rep



