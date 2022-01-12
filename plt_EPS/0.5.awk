#!/usr/bin/gawk -f

{
  print $1,$2,$3,$4,$5,$6,$7,f($8),$9;
}

func f(x) {
  if (x>0.5) {
    return x-0.5;
  } else {
    return x+0.5;
  }
}

