#!/usr/bin/env python

import numpy as np

def hermite(xp,x,f):
    ier=1
    io=0
    iup=0
    n=len(x)
    if x[1]<x[0]:
        iup=1
    n1=n-2
    if (xp>=x[n-1] and iup==0) or (xp<=x[n-1] and iup==1):
        p=f[n-1]
        ier=2
        return p,ier
    elif (xp<=x[0] and iup==0) or (xp>=x[0] and iup==1):
        p=f[0]
        ier=2
        return p,ier
    i=io
    while i<n and xp>x[i] and iup==0:
        i+=1
    while i<n and xp<x[i] and iup==1:
        i+=1
    i-=1
    io=i+1
    lp1=1.0/(x[i]-x[i+1])
    lp2=1.0/(x[i+1]-x[i])
    if i==0:
        fp1=(f[1]-f[0])/(x[1]-x[0])
    else:
        fp1=(f[i+1]-f[i-1])/(x[i+1]-x[i-1])
    if i>=n1-1:
        fp2=(f[n-1]-f[n-2])/(x[n-1]-x[n-2])
    else:
        fp2=(f[i+2]-f[i])/(x[i+2]-x[i])
    xpi1=xp-x[i+1]
    xpi=xp-x[i]
    l1=xpi1*lp1
    l2=xpi*lp2
    p=f[i]*(1.0-2.0*lp1*xpi)*l1*l1 + f[i+1]*(1.0-2.0*lp2*xpi1)*l2*l2 + fp2*xpi1*l2*l2 + fp1*xpi*l1*l1
    return p,ier

def main():

    x = np.linspace(0.0,2.0*np.pi,10)
    y = np.sin(x)
    
    xnew = np.linspace(0.0,2.0*np.pi,100)
    ynew = np.sin(xnew)
    yher = np.zeros(len(xnew))
    for i in xrange(0,len(xnew)):
        yher[i],ier = hermite(xnew[i],x,y)
        print xnew[i],yher[i],ier
    
    import matplotlib.pyplot as plt

    ax = plt.subplot(111)
    ax.plot(x,y,'+-')
    ax.plot(xnew,ynew,'-')
    ax.plot(xnew,yher,'+-')
    plt.show()

if __name__ == "__main__":
    main()

