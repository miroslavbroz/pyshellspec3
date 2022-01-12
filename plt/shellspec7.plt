#!/usr/bin/gnuplot

set term x11

set xl "iunt"
set yl "contrib.f"

#set yr [0:]
set logscale y

p \
  "<awk '/opac2/{ i=0; }(i==1){ print; }/opdep/{ i=1; }' shellspec.out" u 1:3 w lp
pa -1



