#!/usr/bin/gnuplot

f(x) = log10(abs(x))

array s[99]
`awk 'BEGIN{ for (i=1; i<=99; i++){ print "s[" i "] = \"\"; "; } }' /dev/null`
`awk '(NR==1){ for (i=1; i<=NF; i++){ print "s[" i "] = \"" $i "\"; "; }}' fit.log`

set xl "Niter"
set yl "params, f(x) = log10(abs(x))"

set xyplane 0.0
set key outside

tmp=135; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0

p \
  "fit.log" u 1:(f($2))  t s[ 2] w p ,\
  "fit.log" u 1:(f($3))  t s[ 3] w p ,\
  "fit.log" u 1:(f($4))  t s[ 4] w p ,\
  "fit.log" u 1:(f($5))  t s[ 5] w p ,\
  "fit.log" u 1:(f($6))  t s[ 6] w p ,\
  "fit.log" u 1:(f($7))  t s[ 7] w p ,\
  "fit.log" u 1:(f($8))  t s[ 8] w p ,\
  "fit.log" u 1:(f($9))  t s[ 9] w p ,\
  "fit.log" u 1:(f($10)) t s[10] w l ,\
  "fit.log" u 1:(f($11)) t s[11] w l ,\
  "fit.log" u 1:(f($12)) t s[12] w l ,\
  "fit.log" u 1:(f($13)) t s[13] w l ,\
  "fit.log" u 1:(f($14)) t s[14] w l ,\
  "fit.log" u 1:(f($15)) t s[15] w l ,\
  "fit.log" u 1:(f($16)) t s[16] w l ,\
  "fit.log" u 1:(f($17)) t s[17] w l ,\
  "fit.log" u 1:(f($18)) t s[18] w lp,\
  "fit.log" u 1:(f($19)) t s[19] w lp,\
  "fit.log" u 1:(f($20)) t s[20] w lp,\
  "fit.log" u 1:(f($21)) t s[21] w lp,\
  "fit.log" u 1:(f($22)) t s[22] w lp,\
  "fit.log" u 1:(f($23)) t s[23] w lp,\
  "fit.log" u 1:(f($24)) t s[24] w lp,\
  "fit.log" u 1:(f($25)) t s[25] w lp,\
  "fit.log" u 1:(f($26)) t s[26] w p ps 2,\
  "fit.log" u 1:(f($27)) t s[27] w p ps 2,\
  "fit.log" u 1:(f($28)) t s[28] w p ps 2,\
  "fit.log" u 1:(f($29)) t s[29] w p ps 2,\
  "fit.log" u 1:(f($30)) t s[30] w p ps 2,\
  "fit.log" u 1:(f($31)) t s[31] w p ps 2,\
  "fit.log" u 1:(f($32)) t s[32] w p ps 2,\
  "fit.log" u 1:(f($33)) t s[33] w p ps 2,\
  "fit.log" u 1:(f($34)) t s[34] w p ps 2,\
  "fit.log" u 1:(f($35)) t s[35] w lp,\
  "fit.log" u 1:(f($36)) t s[36] w lp,\
  "fit.log" u 1:(f($37)) t s[37] w lp,\
#  "fit.log" u 1:(f($38)) t s[38] w lp,\
#  "fit.log" u 1:(f($39)) t s[39] w lp,\

pa -1

set term png small
set out "chi2_params.png"
rep


q
