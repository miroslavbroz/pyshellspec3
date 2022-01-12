#!/usr/bin/gnuplot

set colors classic

set term post eps enh color solid "Helvetica" 18
set out "chi2_VAMP_null.eps"
set size 0.8,0.6

nm = 1.e-9
ang = 1.e-10/nm
rad = 180./pi

set xl "{/Symbol-Oblique l} [nm]"
set yl "|d{/Helvetica-Oblique V}| [] (normalized to 1)"

set yr [0:2]
set xtics 1
set mxtics 5
set ytics 1
set mytics 2
set grid ytics
set zeroaxis
set key at graph 0.975,graph 0.95 samplen 1.5 font "Helvetica,14"

set lmargin 5.8
set rmargin 1.6
set tmargin 1.0

tmp=1.05
call "line.plt" "H_{/Symbol a}" 6563
call "line.plt" ""              6563
call "line.plt" "H_{/Symbol b}" 4861
call "line.plt" "H_{/Symbol g}" 4341
call "line.plt" "H_{/Symbol d}" 4102
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
  "<./null.awk 4 6 7 ../*.vis.syn.dat" u ($1/nm):($2+1):3 not w err lc 'gray' ps 0,\
  "<./null.awk 4 8 7 ../*.vis.syn.dat" u ($1/nm):($2+1)   t "synthetic" w l lt 6,\
  "<./null.awk 4 6 7 ../*.vis.syn.dat" u ($1/nm):($2+1)   t "observed"  w l lt 3,\
  1.0 w l not lt 0

system("patch chi2_VAMP_null.eps patch.diff")

q

