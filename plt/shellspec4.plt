#!/usr/bin/gnuplot

set term x11

set xl "iunt"
set yl "AKONC, HI, H-, H2, HEI, HIPF, HIIPF"

set yr [1.e-20:]
set logscale y

p \
  "<awk '/ilin/{ i=0; }(i==1){ print; }/AKONC/{ i=1; }' shellspec.out" u 1:2 t "AKONC" w lp,\
  "<awk '/ilin/{ i=0; }(i==1){ print; }/AKONC/{ i=1; }' shellspec.out" u 1:3 t "HI"    w lp,\
  "<awk '/ilin/{ i=0; }(i==1){ print; }/AKONC/{ i=1; }' shellspec.out" u 1:4 t "H-"    w lp,\
  "<awk '/ilin/{ i=0; }(i==1){ print; }/AKONC/{ i=1; }' shellspec.out" u 1:5 t "H2"    w lp,\
  "<awk '/ilin/{ i=0; }(i==1){ print; }/AKONC/{ i=1; }' shellspec.out" u 1:6 t "HEI"   w lp,\
  "<awk '/ilin/{ i=0; }(i==1){ print; }/AKONC/{ i=1; }' shellspec.out" u 1:7 t "HIPF"  w lp,\
  "<awk '/ilin/{ i=0; }(i==1){ print; }/AKONC/{ i=1; }' shellspec.out" u 1:8 t "HIIPF" w lp
pa -1



