#!/usr/bin/gnuplot

nm = 1.e-9
ang = 1.e-10/nm

fac = 2.0

set colors classic
set term x11

set xl "lambda [nm]"
set yl "F_{lambda} [] (shifted by dataset number)"

#set xr [654:658]
set yr [0:]
set xtics 1.0
set ytics 0.5
set grid ytics
set zeroaxis
set errorbars small

call "line2.plt" "Halpha" 6562.8
call "line.plt" "+187 km" (6562.8*(1+187e3/3e8))
call "line.plt" "-187"    (6562.8*(1-187e3/3e8))
call "line.plt" "+41"     (6562.8*(1+41e3/3e8))
call "line.plt" "-41"     (6562.8*(1-41e3/3e8))
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
call "line.plt" "-187 km" (6347*(1-187e3/3e8))
call "line.plt" "+187 km" (6347*(1+187e3/3e8))
call "line.plt" "SiII"   6371
call "line.plt" "NeI"    6402

p \
  "<awk '{ print $2,$4,$7; print $2,$6,$7; print null; }' *.spe.syn.dat" u ($1/nm):($2+($3-1)*fac) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$4,$6,$7; print $1,$2,$4,$6,$7; print null; }' *.spe.syn.dat" u ($2/nm):($4+($5-1)*fac) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.spe.syn.dat" u ($2/nm):($6+($7-1)*fac) t "synthetic" w lp lt 7 pt 1,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.spe.syn.dat" u ($2/nm):($4+($7-1)*fac):5 t "observed" w err lt 3 pt 1 ps 0.5,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.spe.syn.dat" u ($2/nm):($4+($7-1)*fac):5 not w lp  lt 3 pt 1 ps 0.5,\
  1.0 lt 0

pa -1

set term png small size 1920,1080
set out "chi2_SPE.png"
rep

q

  "<awk '(FNR==2){ print $1,$2,$4,FILENAME; }' *.spe.syn.dat" u (660.):(($4-1)*fac):5 not w labels left font "Helvetica,8",\


