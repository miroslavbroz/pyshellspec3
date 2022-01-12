#!/usr/bin/env python

import sys
import subprocess
import numpy as np

def f_G(x,mu,sigma):
    """Gauss function"""
    return 1.0/(np.sqrt(2.0*np.pi)*sigma) * np.exp(-0.5*((x-mu)/sigma)**2)

def convolve(x, y, sigma=1.0, decim=1):
    """Convolve array y(x) with function f_G"""
    wd = 3.0*sigma
    i = 1
    j = 0
    while x[i]-x[0] < wd and i < len(x):
        i += 1

    xc = []
    yc = []
    while x[-1]-x[i] > wd and i < len(x):
        s = 0.0
        norm = 0.0
        k = j
        print "--"
        print "i = ", i
        print "j = ", j
        while x[k]-x[i] < wd:
            print "k = ", k
            tmp = f_G(x[k]-x[i],0.0,sigma)
            s += y[k] * tmp
            norm += tmp
            k += 1
        print "--"
        xc.append(x[i])
        yc.append(s/norm)
        i += decim
        while x[i]-x[j] > wd:
            j += 1

    return xc, yc

def main():
    """Test of a spectrum convolution with an intrumental profile."""

    sigma = 1.0
    s = 0.0
    n = 10
    x1 = -5.0*sigma
    x2 = -x1
    dx = (x2-x1)/n
    for i in xrange(0,n+1):
        x = x1+(x2-x1)*i/n
        tmp = f_G(x,0.0,sigma)
        s += tmp*dx
        print "x = ", x, ", f_G = ", tmp
    print "s = ", s
    sys.exit(0)

    data = np.loadtxt("shellspectrum", usecols=[0,3]).transpose()
    x = data[0]; y = data[1]
    xc, yc = convolve(x, y, sigma=5.0, decim=1)
    np.savetxt("convolve.out", np.transpose([xc,yc]))

    subprocess.check_output("./convolve.plt", shell=True)

if __name__ == "__main__":
    main()


