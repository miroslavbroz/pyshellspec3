#!/usr/bin/gnuplot

set colors classic
set term x11

set xl "u [cycles]"
set yl "v [cycles]"
set cbl "V^2 []"

tmp=5.e8
set xr [-tmp:tmp]
set yr [-tmp:tmp]

set size ratio -1
set zeroaxis
set palette rgbformulae 31,13,10

p \
  "<awk '{ print; }' *.vis2.syn.dat" u ($1/$4):($2/$4):($7) t "all" w p lt 7 pt 1 lc palette z,\

pa -1

set term png small
set out "uv.png"
rep

q

  "<awk '($NF+0>100){ print; }' *.vis2.syn.dat" u ($1/$4):($2/$4) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\


