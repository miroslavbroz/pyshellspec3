#!/usr/bin/gnuplot

set colors classic

set term post eps enh color dashed "Helvetica" 18
set out "chi2_iter.eps"
set size 0.725,0.8

chi2 = `awk 'BEGIN{ m=1e38; }($(NF-9)<m){ m=$(NF-9); }END{ print m; }' fit.log`

set xl "{/Helvetica-Oblique N}_{iter}"
set yl "{/Symbol-Oblique c}^2"

set xr [0:2000]
set logscale y
set yr [1e4:1e7]
set ytics (\
  "10^4" 1e4 0,\
  ""     2e4 1,\
  ""     3e4 1,\
  ""     4e4 1,\
  ""     5e4 1,\
  ""     6e4 1,\
  ""     7e4 1,\
  ""     8e4 1,\
  ""     9e4 1,\
  "10^5" 1e5 0,\
  ""     2e5 1,\
  ""     3e5 1,\
  ""     4e5 1,\
  ""     5e5 1,\
  ""     6e5 1,\
  ""     7e5 1,\
  ""     8e5 1,\
  ""     9e5 1,\
  "10^6" 1e6 0,\
  ""     2e6 1,\
  ""     3e6 1,\
  ""     4e6 1,\
  ""     5e6 1,\
  ""     6e6 1,\
  ""     7e6 1,\
  ""     8e6 1,\
  ""     9e6 1,\
  "10^7" 1e7 0,\
)
set key at graph 1.24,graph 1.0 spacing 1.3 samplen 1.5 font "Helvetica,14"
set grid noxtics noytics front

set lmargin 6.0
set rmargin 7.5
set bmargin 3.1
set tmargin 0.7

p \
  "<awk '{ print $1,$(NF- 9); }' fit.log" u 1:2 t "{/Symbol-Oblique c}^2" w lp lc 'red' lt 1 ps 0.5,\
  "<awk '{ print $1,$(NF-17); }' fit.log" u 1:2 t "LC"    w l  lc '#00cc66',\
  "<awk '{ print $1,$(NF-16); }' fit.log" u 1:2 t "VIS"   w l  lc '#3333ff',\
  "<awk '{ print $1,$(NF-15); }' fit.log" u 1:2 t "CLO"   w l  lc '#cc00cc',\
  "<awk '{ print $1,$(NF-14); }' fit.log" u 1:2 t "T3"    w l  lc 'cyan',\
  "<awk '{ print $1,$(NF-11); }' fit.log" u 1:2 t "SED"   w lp lc '#0000aa' pt 2 ps 0.5,\
  "<awk '{ print $1,$(NF-10); }' fit.log" u 1:2 t "SPE"   w lp lc 'orange' pt 3 ps 0.5,\
  "<awk '{ print $1,$(NF-13); }' fit.log" u 1:2 t "VAMP"  w l  lc '#cc0000',\
  "<awk '{ print $1,$(NF-12); }' fit.log" u 1:2 t "VPHI"  w l  lc '#cccc00',\
  chi2 w l t sprintf("{/=10 %.0f}",chi2) lt 0



