#!/usr/bin/gnuplot

set term x11

set xl "iunt"
set yl "dens"

p "<awk '/akonc/{ i=0; }(i==1) && ($5<1.e15){ print; }/vel./{ i=1; }' shellspec.out" u 1:5 w lp
pa -1



