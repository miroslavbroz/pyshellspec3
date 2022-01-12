#!/bin/sh

make

Teff=14600
logg=2.31357006
Z=0.0
lambda1=0
lambda2=-1e-10
dlambda=1e-10

./limcof_test $Teff $logg $Z $lambda1 $lambda2 $dlambda

#awk '/^ *#.* = / && !/Warning/{ print $2 " = " $4; }' < limcof_test.out > limcof_test.lab
#./limcof_test.plt

