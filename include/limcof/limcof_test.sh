#!/bin/sh

make

Teff=13300
logg=4.25
Z=0.0
lambda1=1250
lambda2=47500
dlambda=1000

./limcof_test $Teff $logg $Z $lambda1 $lambda2 $dlambda > limcof_test.out

awk '/^ *#.* = / && !/Warning/{ print $2 " = " $4; }' < limcof_test.out > limcof_test.lab
./limcof_test.plt

