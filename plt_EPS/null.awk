#!/usr/bin/gawk -f

BEGIN{
  X=ARGV[1];
  Y=ARGV[2];
  E=ARGV[3];
  ARGV[1]=ARGV[2]=ARGV[3]="";
  n=0;
}
(FNR==1) && (n>0){
  prn();
}
!/^#/{
  n++;
  x[n]=$X;
  y[n]=$Y;
  e[n]=$E;
  filename=FILENAME;
#  e[n]=0.0;
}
END{
  ARGIND++;
  prn();
}

func prn() {
  s=0.0;
  m=0;
  for (i=1; i<=n; i++){
#    if ((i<=10) || (i>30)) {
#    if (i<=10) {
    if ((x[i]<=655.0e-9) || (x[i]>658.0e-9)) {
      m++;
      s=s+y[i];
      print x[i],y[i],e[i] >"null.tmp";
    }
  }
  ym=s/m;
  print "# m = ", m;
  print "# ym = ", ym;
  print "" >"null.tmp";

  for (i=1; i<=n; i++){
    print x[i],y[i]-ym,e[i],ARGIND-4,filename;
  }
  print "";
  
  n=0;
}

func abs(x) {
  if (x > 0) {
    return x;
  } else {
    return -x;
  }
}


