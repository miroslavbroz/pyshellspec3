#!/usr/bin/gnuplot

set colors classic

set term post eps enh color solid "Helvetica" 18
set out "uv.eps"
set size 0.65,0.9

set xl "{/Helvetica-Oblique u} [10^8 cycles]"
set yl "{/Helvetica-Oblique v} [10^8 cycles]"

set size ratio -1
set zeroaxis
set key left samplen 0.5
set bmargin 3.2
set tmargin 0.5
set rmargin 0.5

tmp=3.8
set xr [-tmp:tmp]
set yr [-tmp:tmp]
set mxtics 2
set mytics 2

f(x)=x/1.e8

set style line 1 lc '#00ff00' pt 7 ps 0.3
set style line 2 lc '#999999' pt 7 ps 0.4
set style line 3 lc '#00ffff' pt 7 ps 0.6
set style line 4 lc 'black'   pt 6 ps 0.8 lw 0.8

p \
  "<awk '{ print; }' ../npoi*.vis2.syn.dat"  u (f($1/$4)):(f($2/$4))   t "NPOI"  w p ls 1,\
  "<awk '{ print; }' ../npoi*.vis2.syn.dat"  u (-f($1/$4)):(-f($2/$4)) not       w p ls 1,\
  "<awk '{ print; }' ../npoi*.cp.syn.dat"    u (f($1/$6)):(f($2/$6))   not       w p ls 1,\
  "<awk '{ print; }' ../npoi*.cp.syn.dat"    u (f($3/$6)):(f($4/$6))   not       w p ls 1,\
  "<awk '{ print; }' ../MIRC*.vis2.syn.dat"  u (f($1/$4)):(f($2/$4))   t "MIRC"  w p ls 2,\
  "<awk '{ print; }' ../MIRC*.vis2.syn.dat"  u (-f($1/$4)):(-f($2/$4)) not       w p ls 2,\
  "<awk '{ print; }' ../MIRC*.cp.syn.dat"    u (f($1/$6)):(f($2/$6))   not       w p ls 2,\
  "<awk '{ print; }' ../MIRC*.cp.syn.dat"    u (f($3/$6)):(f($4/$6))   not       w p ls 2,\
  "<awk '{ print; }' ../vega*.vis2.syn.dat"  u (f($1/$4)):(f($2/$4))   t " VEGA" w p ls 3,\
  "<awk '{ print; }' ../vega*.vis2.syn.dat"  u (-f($1/$4)):(-f($2/$4)) not       w p ls 3,\
  "<awk '{ print; }' ../IC2013*.vis.syn.dat" u (f($1/$4)):(f($2/$4))   t " |dV|" w p ls 4,\
  "<awk '{ print; }' ../IC2013*.vis.syn.dat" u (-f($1/$4)):(-f($2/$4)) not       w p ls 4,\

system("patch uv.eps patch.diff")

q

  "<awk '{ print; }' ../vega*.cp.syn.dat"   u (f($1/$6)):(f($2/$6)) not      w p lt 4 pt 7,\
  "<awk '{ print; }' ../vega*.cp.syn.dat"   u (f($3/$6)):(f($4/$6)) not      w p lt 4 pt 7,\
