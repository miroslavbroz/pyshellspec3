#!/usr/bin/gnuplot

# Eq. (1)
P0 = 12.913779
JD0 = 2408247.968  # i.e. as in Ak etal. (2007); NOT shifted by 1/2 P
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
set term x11

shift = 2.0

set xl "JD - 2400000"
set yl "magnitude [mag]"

set xr [-0.2:1.2]
set yr [:] reverse
set yr [15:5]
#set ytics shift
set grid ytics
set key left width -2 samplen 0.5

p \
  "<awk '{ print $1,$4,ARGIND; print $1,$6,ARGIND; print null; }' *.lc.syn.dat" u (phase($1)):($2+$3*shift)   t "residua"   w l lt 1 lw 3,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.lc.syn.dat"             u (phase($1)):($4+$9*shift):5 t "observed"  w err lt 3 pt 1 ps 0.5,\
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.lc.syn.dat"             u (phase($1)):($6+$9*shift)   t "synthetic" w p pt 1 lt 7,\
  "<awk '($NF+0>100){ print $1,$4,ARGIND; print $1,$6,ARGIND; print null; }' *.lc.syn.dat" u (phase($1)):($2+$3*shift) t "chi^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '(FNR==2){ print $1,$4,ARGIND,FILENAME; }' *.lc.syn.dat"                u (1.0):($2+$3*shift):4       not           w labels left font "Helvetica,10"
pa -1

set term png small size 1920,1080
set out "chi2_LC_PHASE.png"
rep

q



