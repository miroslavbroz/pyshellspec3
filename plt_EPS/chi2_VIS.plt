#!/usr/bin/gnuplot

set colors classic

set term post eps enh color solid "Helvetica" 18
set out "chi2_VIS.eps"
set size 1.2,1.2

rad = 180./pi

set xl "{/Helvetica-Oblique B}/{/Symbol l} [cycles]"
set yl "squared visbility {/Helvetica-Oblique V }^{2} [] (shifted by dataset number)"

set xr [0:4.0e8]
set yr [0:30]
set xtics (\
  "0"               0.0e8 0,\
  "0.5{{\264}}10^8" 0.5e8 0,\
  "1.0{{\264}}10^8" 1.0e8 0,\
  "1.5{{\264}}10^8" 1.5e8 0,\
  "2.0{{\264}}10^8" 2.0e8 0,\
  "2.5{{\264}}10^8" 2.5e8 0,\
  "3.0{{\264}}10^8" 3.0e8 0,\
  "3.5{{\264}}10^8" 3.5e8 0,\
  "4.0{{\264}}10^8" 4.0e8 0,\
  )
set ytics 1
set grid ytics
set zeroaxis
set key samplen 1.0 width -3 font "Helvetica,14"

set lmargin 7.0
set rmargin 2.7
set bmargin 3.1
set tmargin 0.5

p \
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' ../*.vis2.syn.dat" u (sqrt($1**2+$2**2)/$4):($7+$9)   t "synthetic" w p   lt 6 pt 1,\
  "<awk '(FNR==1){ print null; }($5<2.0) && ($6<1.0){ print $0,ARGIND; }' ../*.vis2.syn.dat" u (sqrt($1**2+$2**2)/$4):($5+$9):6 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '($5<2.0) && ($6<1.0){ print $1,$2,$4,$5,ARGIND; print $1,$2,$4,$7,ARGIND; print null; }' ../*.vis2.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$4,$5,ARGIND; print $1,$2,$4,$7,ARGIND; print null; }' ../*.vis2.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5) t "{/Symbol c}^2 > 100" w p lt 1 pt 6 ps 1.5,\

system("patch chi2_VIS.eps patch.diff")

q

  "<awk '(FNR==2){ print $1,$2,$4,ARGIND,FILENAME; }' ../*.vis2.syn.dat" u (4.0e8):($4):5 not w labels left font "Helvetica,10"

  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' ../*.vis2.syn.dat" u (sqrt($1**2+$2**2)/$4):($5+$9):6 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '{ print $1,$2,$4,$5,ARGIND; print $1,$2,$4,$7,ARGIND; print null; }' ../*.vis2.syn.dat" u (sqrt($1**2+$2**2)/$3):($4+$5) t "residua"  w l lt 1 lw 3,\


