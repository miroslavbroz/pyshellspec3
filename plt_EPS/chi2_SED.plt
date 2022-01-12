#!/usr/bin/gnuplot

set colors classic
set term post eps enh color solid "Helvetica" 18
set out "chi2_SED.eps"

nm = 1.e-9

set xl "{/Symbol-Oblique l} [nm]"
set yl "{/Helvetica-Oblique F}_{/Symbol-Oblique l} [J s^{-1} m^{-2} m^{-1}] (shifted by dataset number)" offset 1,0

set xr [300:750]
set yr [0:0.00575]
set zeroaxis
set key at graph 0.985,graph 0.975 spacing 1.2 samplen 1.5 font "Helvetica,14"

set lmargin 9.0
set rmargin 1.6

fac = 0.002
fac = 0.0
ang = 0.1

call "line2.plt" "H_{/Symbol {\245}}" 3646
call "line2.plt" "H_{/Symbol a}" 6563
call "line2.plt" "H_{/Symbol b}" 4861
call "line2.plt" "H_{/Symbol g}" 4341
call "line2.plt" "H_{/Symbol d}" 4102
#call "line3.plt" "HeI"    4009
#call "line3.plt" "HeI"    4026
#call "line3.plt" "HeI"    4120
#call "line3.plt" "HeI"    4143
#call "line3.plt" "HeI"    4387
#call "line3.plt" "HeI"    4471
#call "line3.plt" "HeI"    4713
#call "line3.plt" "HeI"    4922
#call "line3.plt" "HeI"    5016
#call "line3.plt" "HeI"    5047
call "line3.plt" "HeI"    5876
call "line3.plt" "  HeI"    6678
#call "line3.plt" "CII"    4267
#call "line3.plt" "CII"    6578
#call "line3.plt" "CII"    6583
#call "line3.plt" "MgII"   4481
#call "line3.plt" "SiII"   4128
#call "line3.plt" "SiII"   4130
call "line3.plt" "SiII"   6347
#call "line3.plt" "SiII"   6371
#call "line3.plt" "NeI"    6402

ang = 1.e-10

p \
  "<awk '{ print $2,$4,$7; print $2,$6,$7; print null; }' *.sed.syn.dat" u ($1/nm):($2+$3*fac) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$4,$6,$7; print $1,$2,$4,$6,$7; print null; }' *.sed.syn.dat" u ($2/nm):($4+$5*fac) t "{/Symbol-Oblique c}^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.sed.syn.dat" u ($2/nm):($6+$7*fac) t "synthetic" w lp lt 6 pt 1,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.sed.syn.dat" u ($2/nm):($4+$7*fac):5 t "observed" w err lt 3 pt 1 ps 0.5,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.sed.syn.dat" u ($2/nm):($4+$7*fac):5 not w lp  lt 3 pt 1 ps 0.5,\
  "<awk '(FNR==2){ print $1,$2,$4,FILENAME; }' *.sed.syn.dat" u (800.):($4*360):5 not w labels left font "Helvetica,8"

system("patch chi2_SED.eps patch.diff")



