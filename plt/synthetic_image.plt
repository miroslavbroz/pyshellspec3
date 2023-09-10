#!/usr/bin/gnuplot

cm = 1.e-2  # m
R_S = 6.9634e8  # m

set xl "x [R_S]"
set yl "y [R_S]"
set cbl "I_{nu} [erg cm^{-2} s^{-1} Hz^{-1} sr^{-1}]" offset 1,0

set cbr [:]
set logscale cb

set pm3d map
set size ratio -1
set nokey
set zeroaxis
#set palette gray

splot "tempspe00/2Dimage_001" u ($1*cm/R_S):($2*cm/R_S):3

pa -1

set term png small
set out "synthetic_image.png"
rep


