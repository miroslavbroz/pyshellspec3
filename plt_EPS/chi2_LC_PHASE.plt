#!/usr/bin/gnuplot

# Eq. (1)
P0 = 12.913779
JD0 = 2408247.968
dotP = 3.87265e-6

# JD = JD0 + P0 E + dotP E^2
# dotP E^2 + P0 E + JD0-JD = 0
# a x^2    + b x  + c      = 0
# x = (-b +- sqrt(b^2 - 4ac)) / 2a

x1(a,b,c) = (-b + sqrt(b**2 - 4.*a*c)) / (2.*a)

E(JD) = x1(dotP, P0, JD0-JD)
zero_(E) = E > 0 ? E : E+1.0
phase(JD) = zero_(E(JD)-int(E(JD)))

########################################################################

set colors classic

set term post eps enh color solid "Helvetica" 18
set out "chi2_LC_PHASE.eps"
set size 0.8,2.0

shift = 2.0

set xl "JD {/Symbol -} 2400000"
set yl "magnitude [mag] (shifted by dataset number)" offset 0.5,0

set xr [-0.05:1.05]
set yr [46:4]
set grid ytics
set key at graph -0.01,graph 0.96 width -2 samplen 0.5 font "Helvetica,14"

set lmargin 9.0
set rmargin 5.5
set bmargin 3.1
set tmargin 0.5

p \
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' ../*.lc.syn.dat"             u (phase($1)):($6+$9*shift)   t "synthetic" w p pt 1 lt 6,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' ../*.lc.syn.dat"             u (phase($1)):($4+$9*shift):5 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '{ print $1,$4,ARGIND; print $1,$6,ARGIND; print null; }' ../*.lc.syn.dat" u (phase($1)):($2+$3*shift)   t "residua"   w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$4,ARGIND; print $1,$6,ARGIND; print null; }' ../*.lc.syn.dat" u (phase($1)):($2+$3*shift) t "{/Symbol c}^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '(FNR==2){ s=gensub(\"../\(.*\).phot.lc.syn.dat\",\"{\\220}\\\\1\",1,FILENAME); print $1,$4,ARGIND,s; }' ../*.lc.syn.dat" u (1.05):($2+$3*shift):4 not w labels left font "Helvetica,10"

system("patch chi2_LC_PHASE.eps patch.diff")

q



