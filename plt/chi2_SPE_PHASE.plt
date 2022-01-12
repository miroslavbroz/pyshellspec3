#!/usr/bin/gnuplot

# BEWARE! directories do N

nm = 1.e-9
ang = 1.e-10/nm

set colors classic
set term x11

set xl "lambda [nm]"
set yl "F_{lambda} [] (shifted by dataset number)"
set cbl "phase"

set yr [0:]
set cbr [0:1.1]
set zeroaxis
set colorbox

call "line.plt" "Halpha" 6563
call "line.plt" "Hbeta"  4861
call "line.plt" "Hgamma" 4341
call "line.plt" "Hdelta" 4102
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

p \
  "<sort -k8 -g *.spe.syn.dat | awk '($7!=last){ i++; print null; }{ print $0,i; last=$7; }'" u ($2/nm):($4+$10*fac-2*fac):8 t "observed"  w l lc palette z lw 3,\
  "<sort -k8 -g *.spe.syn.dat | awk '($7!=last){ i++; print null; }{ print $0,i; last=$7; }'" u ($2/nm):($6+$10*fac-2*fac):8 t "synthetic" w l lc palette z,\
  "<sort -k8 -g *.spe.syn.dat | awk '($7!=last){ i++; print $0,i; }{ last=$7; }'" u ($2/nm):($6+$10*fac-2*fac):(sprintf("%.3f ",$8)) not w labels right,\
  "<sort -k8 -g *.spe.syn.dat | awk '($7!=last){ i++; print $0,i; }{ last=$7; }'" u ($2/nm):($6+$10*fac-2*fac-0.30):(sprintf("%.0f ",$7)) not w labels right,\

pa -1

set term png small size 1920,1080
set out "chi2_SPE_PHASE.png"
rep

q


