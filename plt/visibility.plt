#!/usr/bin/gnuplot

set colors classic
set term x11

rad = 180./pi

set xl "B/lambda [cycles]"
set yl "squared visbility V^2 [] (shifted by dataset number)"

set yr [-0.5:2.5]
set ytics 1
set grid ytics
set zeroaxis

p \
  "<awk '(FNR==1){ print null; }{ print $0,ARGIND; }' *.vis2.syn.dat" u (sqrt($1**2+$2**2)/$4):($7+$9-1)   t "synthetic" w p   lt 7 pt 1,\
  "<awk '(FNR==2){ print $1,$2,$4,ARGIND,FILENAME; }' *.vis2.syn.dat" u (3.9e8):($4-1+1):5 not w labels left font "Helvetica,10"

pa -1

set term png small
set out "visibility.png"
rep

q



