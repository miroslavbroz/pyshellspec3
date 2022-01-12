#!/usr/bin/gnuplot

set term x11

set xl "iunt"
#set yl "log10 ophbf, ophff, opth, ophrs, ophn, opmieb, opmiej, oplin, opmol"
set yl "log10 opacity [cgs]"

p \
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($2))  t "ophbf"  w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($3))  t "ophff"  w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($4))  t "opth"   w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($5))  t "ophrs"  w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($6))  t "ophn"   w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($7))  t "opmieb" w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($8))  t "opmiej" w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($9))  t "oplin"  w lp,\
  "<awk '/opdep/{ i=0; }(i==1){ print; }/ophbf/{ i=1; }' shellspec.out" u 1:(log10($10)) t "opmol"  w lp
pa -1

set term png small
set out "shellspec5.png"
rep


