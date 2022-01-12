#!/usr/bin/gnuplot

# BEWARE! directories do N

set colors classic

set term post eps enh color solid "Helvetica" 18
set out "chi2_SPE_PHASE_.eps"
set size 1,1.8

nm = 1.e-9
ang = 1.e-10/nm

set xl "{/Symbol-Oblique l} [nm]"
set yl "{/Helvetica-Oblique I}_{/Symbol-Oblique l} [] (shifted by dataset number)" offset 1,0
set cbl "phase"

set yr [0:24.0]
set cbr [0:1.1]
set zeroaxis
set colorbox
set key at graph 0.9,graph 0.990 samplen 1.5 font "Helvetica,14"

set lmargin 6.0
set rmargin 1.6

tmp=1.01
call "line.plt" "H_{/Symbol a}" 6563
call "line.plt" "H_{/Symbol b}" 4861
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
call "line.plt" "MgII"   4481
call "line.plt" "SiII"   4128
call "line.plt" "SiII"   4130
call "line.plt" "SiII"   6347
call "line.plt" "SiII"   6371
call "line.plt" "NeI"    6402

fac = 2.0
f(x) = x-2
g(x) = x

p \
  "<./0.5.awk *.spe.syn.dat | sort -k8 -g | awk '($7!=last){ i++; print null; }{ print $0,i; last=$7; }'" u ($2/nm):($4+f($10)*fac):8 t "observed"  w l lc palette z lw 3,\
  "<./0.5.awk *.spe.syn.dat | sort -k8 -g | awk '($7!=last){ i++; print null; }{ print $0,i; last=$7; }'" u ($2/nm):($6+f($10)*fac):8 t "synthetic" w l lc palette z,\
  "<./0.5.awk *.spe.syn.dat | sort -k8 -g | awk '($7!=last){ i++; print $0,i; }{ last=$7; }'" u ($2/nm):($6+f($10)*fac+0.30):(sprintf("{/=12 %.3f }",g($8))) not w labels right,\
  "<./0.5.awk *.spe.syn.dat | sort -k8 -g | awk '($7!=last){ i++; print $0,i; }{ last=$7; }'" u ($2/nm):($6+f($10)*fac-0.30):(sprintf("{/=12 %.0f }",$7)) not w labels right,\

q


