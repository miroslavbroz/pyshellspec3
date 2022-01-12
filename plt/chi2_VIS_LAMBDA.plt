#!/usr/bin/gnuplot

set colors classic
set term x11

rad = 180./pi
nm = 1.e-9
ang = 1.e-10/nm

set xl "lambda [nm]"
set yl "squared visbility V^2 [] (shifted by dataset number)"

set yr [0:40]
set ytics 1
set grid ytics
set zeroaxis

call "line.plt" "Halpha" 6563
call "line.plt" "Hbeta"  4861
call "line.plt" "Hgamma" 4341
call "line.plt" "Hdelta" 4102
call "line.plt" "HeI"    4009
call "line.plt" "HeI"    4026
call "line.plt" "HeI"    4120
call "line.plt" "HeI"    4143
call "line.plt" "HeI"    4387
call "line.plt" "HeI"    4471
call "line.plt" "HeI"    4713
call "line.plt" "HeI"    4922
call "line.plt" "HeI"    5016
call "line.plt" "HeI"    5047
call "line.plt" "HeI"    5876
call "line.plt" "HeI"    6678
call "line.plt" "CII"    4267
call "line.plt" "MgII"   4481
call "line.plt" "SiII"   4128
call "line.plt" "SiII"   4130
call "line.plt" "SiII"   6347
call "line.plt" "SiII"   6371
call "line.plt" "NeI"    6402

p \
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.vis2.syn.dat" u ($4/nm):($7+$9)   t "synthetic" w p   lt 7 pt 1,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.vis2.syn.dat" u ($4/nm):($5+$9):6 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '{ print $1,$2,$4,$5,ARGIND; print $1,$2,$4,$7,ARGIND; print null; }' *.vis2.syn.dat" u ($3/nm):($4+$5) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$4,$5,ARGIND; print $1,$2,$4,$7,ARGIND; print null; }' *.vis2.syn.dat" u ($3/nm):($4+$5) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '(FNR==2){ print $1,$2,$4,ARGIND,FILENAME; }' *.vis2.syn.dat" u (2000.):($4):5 not w labels left font "Helvetica,8"

pa -1

set term png small
set out "chi2_VIS_LAMBDA.png"
rep



