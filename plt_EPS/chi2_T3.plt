#!/usr/bin/gnuplot

set colors classic

set term post eps enh color solid "Helvetica" 18
set out "chi2_T3.eps"
set size 1.2,1.2

rad = 180./pi

set xl "{/Helvetica-Oblique B}/{/Symbol l} [cycles]"
set yl "triple product amplitude |{/Helvetica-Oblique T}_3| [] (shifted by dataset number)"

shift = 2
set xr [0.0:2.3e8]
set yr [0:12.25*shift]
set xtics (\
  "0"               0.0e8 0,\
  ""                0.1e8 1,\
  ""                0.2e8 1,\
  ""                0.3e8 1,\
  ""                0.4e8 1,\
  "0.5{{\264}}10^8" 0.5e8 0,\
  ""                0.6e8 1,\
  ""                0.7e8 1,\
  ""                0.8e8 1,\
  ""                0.9e8 1,\
  "1.0{{\264}}10^8" 1.0e8 0,\
  ""                1.1e8 1,\
  ""                1.2e8 1,\
  ""                1.3e8 1,\
  ""                1.4e8 1,\
  "1.5{{\264}}10^8" 1.5e8 0,\
  ""                1.6e8 1,\
  ""                1.7e8 1,\
  ""                1.8e8 1,\
  ""                1.9e8 1,\
  "2.0{{\264}}10^8" 2.0e8 0,\
  ""                2.1e8 1,\
  ""                2.2e8 1,\
  ""                2.3e8 1,\
  ""                2.4e8 1,\
  "2.5{{\264}}10^8" 2.5e8 0,\
  )
set ytics shift
set grid ytics ytics
set zeroaxis
set key samplen 1.0 width -3 font "Helvetica,14"

set lmargin 9.0
set rmargin 0.5
set bmargin 3.1
set tmargin 0.5

p \
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' ../*.cp.syn.dat" u (sqrt($1**2+$2**2)/$6):($9+$15*shift)   t "synthetic" w p lt 6 pt 1,\
  "<awk '(FNR==1){ print null; }($7<2.0) && ($8<1.0){ print $0,ARGIND; }' ../*.cp.syn.dat" u (sqrt($1**2+$2**2)/$6):($7+$15*shift):8 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '($7<2.0) && ($8<1.0){ print $1,$2,$6,$7,ARGIND; print $1,$2,$6,$9,ARGIND; print null; }' ../*.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*shift) t "residua"  w l lt 1 lw 3,\
  "<awk '($(NF-1)+0>100){ print $1,$2,$6,$7,ARGIND; print $1,$2,$6,$9,ARGIND; print null; }' ../*.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*shift) t "{/Symbol c}^2 > 100" w p lt 1 pt 6 ps 1.5,\

system("patch chi2_T3.eps patch.diff")

q

  "<awk '(FNR==2){ print $1,$2,$7,ARGIND,FILENAME; }' ../*.cp.syn.dat" u (2.3e8):($4*shift):5 not w labels left font "Helvetica,8"

  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' ../*.cp.syn.dat" u (sqrt($1**2+$2**2)/$6):($7+$15*shift):8 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '{ print $1,$2,$6,$7,ARGIND; print $1,$2,$6,$9,ARGIND; print null; }' ../*.cp.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5*shift) t "residua"  w l lt 1 lw 3,\

