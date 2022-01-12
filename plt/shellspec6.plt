#!/usr/bin/gnuplot

set term x11

set xl "iunt"
set yl "opdep"

set yr [0:]

p \
  "<awk '/opac2/{ i=0; }(i==1){ print; }/opdep/{ i=1; }' shellspec.out" u 1:2 w lp
pa -1



