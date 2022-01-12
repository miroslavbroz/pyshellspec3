#!/usr/bin/gnuplot

set colors classic
set term x11

nm = 1.e-9
ang = 1.e-10/nm
rad = 180./pi

set xl "lambda [nm]"
set yl "visibility amplitude |V| [] (normalized to 0)"

#set yr [0:40]
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

#load "wavelengths.lab"

p \
  "<./null.awk 4 6 7 *.vis.syn.dat" u ($1/nm):2:3 not w err lc 'gray' ps 0,\
  "<./null.awk 4 8 7 *.vis.syn.dat" u ($1/nm):2   t "synthetic" w l lt 7,\
  "<./null.awk 4 6 7 *.vis.syn.dat" u ($1/nm):2   t "observed"  w l lt 3,\

pa -1

set term png small
set out "chi2_VAMP_null.png"
rep


q

