#!/usr/bin/gnuplot

set term x11

set xl "iunt"
set yl "dusttt, temp"

p \
  "<awk '/akonc/{ i=0; }(i==1) && ($5<3.e4){ print; }/vel./{ i=1; }' shellspec.out" u 1:7 w lp,\
  "<awk '/AKONC/{ i=0; }(i==1){ print; }/akonc/{ i=1; }' shellspec.out" u 1:4 w lp
pa -1



