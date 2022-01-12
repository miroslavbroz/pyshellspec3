#!/usr/bin/gnuplot

cm = 1.e-2
R_S = 6.955e8
deg = pi/180.
arcsec = deg/3600.
mas = 1.e-3*arcsec
AU = 1.49597870700e11
pc = 648000.0/pi*AU

sma = 58.349
qq = 0.223
omega_an = 253.7409264632917
dd = 316.2810415258063

omega_an = omega_an*deg
dd = dd*pc
xb = -qq*sma/(1.0+qq)

e(x) = x*cm/R_S
f(x,y) = (e(x)-xb)*cos(-omega_an-90*deg) - e(y)*sin(-omega_an-90*deg)
g(x,y) = -((e(x)-xb)*sin(-omega_an-90*deg) + e(y)*cos(-omega_an-90*deg))

########################################################################

set term post eps enh color solid "Helvetica" 18
set out "templc04_fort_96.eps"

set xl "{/Symbol Da} [mas]"
set yl "{/Symbol Dd} [mas]"
set cbl "{/Helvetica-Oblique I}_{/Symbol n} [10^{-4} erg cm^{-2} s^{-1} Hz^{-1} sr^{-1}]" offset +2,0

set xr [-1.0:0.8]
set yr [-0.6:0.6]
set cbr [0:0.0006]

set cbtics (\
  "0" 0.e-4 0,\
  "1" 1.e-4 0,\
  "2" 2.e-4 0,\
  "3" 3.e-4 0,\
  "4" 4.e-4 0,\
  "5" 5.e-4 0,\
  "6" 6.e-4 0,\
  "7" 7.e-4 0,\
  "8" 8.e-4 0,\
  )

set pm3d map
set size ratio -1
set nokey
set border lc rgb "white"
set tics tc rgb "black"
set cbtics tc rgb "black"

set label "{/Symbol l} = 155 nm (FUV)" at graph 0.65,0.9 tc rgb "white" front

splot \
  0,\
  "<cat ../templc04*/fort.96" u (f($1,$2)*R_S/dd/mas):(g($1,$2)*R_S/dd/mas):3


