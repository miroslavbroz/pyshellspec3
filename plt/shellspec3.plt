#!/usr/bin/gnuplot

set term x11

set xl "iunt"
set yl "akonc, ekonc"

set logscale y

p \
  "<awk '/AKONC/{ i=0; }(i==1){ print; }/akonc/{ i=1; }' shellspec.out" u 1:2 w lp,\
  "<awk '/AKONC/{ i=0; }(i==1){ print; }/akonc/{ i=1; }' shellspec.out" u 1:3 w lp
pa -1



