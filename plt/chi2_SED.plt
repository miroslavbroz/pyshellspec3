#!/usr/bin/gnuplot

nm = 1.e-9

set colors classic
set term x11

set xl "lambda [A]"
set yl "F_{lambda} [J s^-1 m^-2 m^-1] (shifted by dataset number)"

set yr [0:]
set zeroaxis

fac = 0.002
fac = 0.0
ang = 1

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
call "line.plt" "CII"    6578
call "line.plt" "CII"    6583
call "line.plt" "MgII"   4481
call "line.plt" "SiII"   4128
call "line.plt" "SiII"   4130
call "line.plt" "SiII"   6347
call "line.plt" "SiII"   6371
call "line.plt" "NeI"    6402

ang = 1.e-10

p \
  "<awk '{ print $2,$4,$7; print $2,$6,$7; print null; }' *.sed.syn.dat" u ($1/ang):($2+$3*fac) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$4,$6,$7; print $1,$2,$4,$6,$7; print null; }' *.sed.syn.dat" u ($2/ang):($4+$5*fac) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.sed.syn.dat" u ($2/ang):($6+$7*fac) t "synthetic" w lp lt 7 pt 1,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.sed.syn.dat" u ($2/ang):($4+$7*fac):5 t "observed" w err lt 3 pt 1 ps 0.5,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.sed.syn.dat" u ($2/ang):($4+$7*fac):5 not w lp  lt 3 pt 1 ps 0.5,\
  "<awk '(FNR==2){ print $1,$2,$4,FILENAME; }' *.sed.syn.dat" u (8000.):($4*360):5 not w labels left font "Helvetica,8"

pa -1

set term png small
set out "chi2_SED.png"
rep



