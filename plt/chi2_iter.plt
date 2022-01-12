#!/usr/bin/gnuplot

chi2 = `awk 'BEGIN{ m=1e38; }($(NF-9)<m){ m=$(NF-9); }END{ print m; }' fit.log`

set colors classic

set xl "iter"
set yl "chi^2"

set term x11

set logscale y
#set yr [100000:200000]; set nologscale y
set yr [:1.e7]

p \
  "<awk '{ print $1,$(NF- 9); }' fit.log" u 1:2 t "chi^2" w lp lt 1 ps 0.5,\
  "<awk '{ print $1,$(NF-17); }' fit.log" u 1:2 t "LC"    w l,\
  "<awk '{ print $1,$(NF-16); }' fit.log" u 1:2 t "VIS"   w l,\
  "<awk '{ print $1,$(NF-15); }' fit.log" u 1:2 t "CLO"   w l,\
  "<awk '{ print $1,$(NF-14); }' fit.log" u 1:2 t "T3"    w l,\
  "<awk '{ print $1,$(NF-13); }' fit.log" u 1:2 t "VAMP"  w l,\
  "<awk '{ print $1,$(NF-12); }' fit.log" u 1:2 t "VPHI"  w l,\
  "<awk '{ print $1,$(NF-11); }' fit.log" u 1:2 t "SED"   w lp pt 2 ps 0.5,\
  "<awk '{ print $1,$(NF-10); }' fit.log" u 1:2 t "SPE"   w lp pt 3 ps 0.5,\
  chi2 w l t sprintf("%f",chi2) lt 0
pa -1 

set term png small
set out "chi2_iter.png"
rep


