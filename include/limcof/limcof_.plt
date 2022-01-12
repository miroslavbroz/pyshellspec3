#!/usr/bin/gnuplot

set term x11

set xl "T_eff [K]"
set yl "log g [cm s^{-2}]"
set zl "x_lin []"

set ticslevel 0

p "limcof.dat" u 2:3 w p
pa -1


