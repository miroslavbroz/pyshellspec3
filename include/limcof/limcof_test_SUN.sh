#!/bin/sh

make

Teff=5780.
logg=4.25
Z=0.0
lambda1=3500
lambda2=-1
dlambda=100

./limcof_test $Teff $logg $Z $lambda1 $lambda2 $dlambda > limcof_test_SUN.out

