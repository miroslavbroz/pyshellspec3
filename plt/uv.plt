#!/usr/bin/gnuplot

set colors classic
set term x11

set xl "u [cycles]"
set yl "v [cycles]"

tmp=2.e8
#set xr [-tmp:tmp]
#set yr [-tmp:tmp]

set size ratio -1
set zeroaxis

p \
  "<awk '{ print; }' *.vis.syn.dat" u ($1/$4):($2/$4) t "all" w p lt 7 pt 1,\
  "<awk '($NF+0>100){ print; }' *.vis.syn.dat" u ($1/$4):($2/$4) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\

pa -1

set term png small
set out "uv.png"
rep

q


