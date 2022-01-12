#!/usr/bin/gnuplot

set colors classic

set term post eps enh color solid "Helvetica" 18
set out "chi2_SPE_.eps"
set size 0.9,0.7

nm = 1.e-9
ang = 1.e-10/nm
fac = 0.0

set xl "{/Symbol-Oblique l} [nm]"
set yl "{/Helvetica-Oblique I}_{/Symbol-Oblique l} []"

set xr [655:660]
set yr [0.5:3.5]
#set xtics 1.0
set mxtics 5
set ytics 1.0
set grid ytics
set zeroaxis
set errorbars 0.5
set key at graph 0.9,graph 0.750 samplen 1.5 font "Helvetica,14"

set rmargin 1.6

tmp=1.040
call "line.plt" "H_{/Symbol a}" 6562.8
#call "line.plt" "+187 km" (6562.8*(1+187e3/3e8))
#call "line.plt" "-187"    (6562.8*(1-187e3/3e8))
#call "line.plt" "+41"     (6562.8*(1+41e3/3e8))
#call "line.plt" "-41"     (6562.8*(1-41e3/3e8))
#call "line2.plt" "-80"   (6562.8*(1-80e3/3e8))
#call "line2.plt" "+374"   (6562.8*(1+374e3/3e8))
call "line.plt" "H_{/Symbol b}"  4861
call "line.plt" "H_{/Symbol g}" 4341
call "line.plt" "H_{/Symbol d}" 4102
call "line.plt" "HeI"    4009
call "line.plt" "HeI"    4026
call "line.plt" "HeI"    4120
call "line.plt" "HeI"    4143
call "line.plt" "HeI"    4387
call "line.plt" "HeI"    4471
call "line.plt" "HeI"    4713
call "line.plt" "HeI"    4922
call "line.plt" "HeI"    5016
call "line.plt" "HeI"    5047
call "line.plt" "HeI"    5876
call "line.plt" "HeI"    6678
call "line.plt" "CII"    4267
call "line.plt" "CII"    6578
call "line.plt" "CII"    6583
call "line.plt" "MgII"   4481
call "line.plt" "SiII"   4128
call "line.plt" "SiII"   4130
call "line.plt" "SiII"   6347
#call "line.plt" "-187 km" (6347*(1-187e3/3e8))
#call "line.plt" "+187 km" (6347*(1+187e3/3e8))
call "line.plt" "SiII"   6371
call "line.plt" "NeI"    6402

p \
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' ../../fitting_shell7_ABUNDANCES1__65/*.spe.syn.dat" u ($2/nm):($6+($7-1)*fac) t "solar" w l lt 6,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.spe.syn.dat" u ($2/nm):($6+($7-1)*fac) t "low C" w l lc '#333333',\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.spe.syn.dat" u ($2/nm):($4+($7-1)*fac):5 t "observed" w lp lc 'cyan' pt 1 ps 0.0,\
  "<awk '(FNR==2){ print $1,$2,$4,FILENAME; }' *.spe.syn.dat" u (660.):(($4-1)*fac):5 not w labels left font "Helvetica,8",\
  1.0 not lt 0

system("patch chi2_SPE_.eps patch.diff")

q
  "<awk '{ print $2,$4,$7; print $2,$6,$7; print null; }' *.spe.syn.dat" u ($1/nm):($2+($3-1)*fac) t "residua"  w l lt 1 lw 3,\
  "<awk '($NF+0>100){ print $1,$2,$4,$6,$7; print $1,$2,$4,$6,$7; print null; }' *.spe.syn.dat" u ($2/nm):($4+($5-1)*fac) t "{/Symbol-Oblique c}^2 > 100" w p lt 1 pt 6 ps 1.5,\
  "<awk '($7!=last){ print null; }{ print $0,ARGIND; last=$7; }' *.spe.syn.dat" u ($2/nm):($4+($7-1)*fac):5 t "observed" w err lt 3 pt 1 ps 0.0,\

  "<cat *.spe.syn.dat | awk '($7!=last){ i++; print $0,i; }{ last=$7; }'" u ($2/nm):($6+$10*fac-2*fac+0.30):(sprintf("{/=12 %.3f }",$8)) not w labels right,\
  "<cat *.spe.syn.dat | awk '($7!=last){ i++; print $0,i; }{ last=$7; }'" u ($2/nm):($6+$10*fac-2*fac-0.30):(sprintf("{/=12 %.0f }",$7)) not w labels right,\
