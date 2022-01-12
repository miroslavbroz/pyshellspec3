#!/usr/bin/env python

import math
import os
import re
import subprocess

__author__ = "Miroslav Broz (miroslav.broz@email.cz)"
__version__ = "Dec 7th 2017"

path = os.path.realpath(__file__)
shsp_limcof = os.path.join('/'.join(path.split('/')[:-1]), 'include/limcof/limcof_test')

class Limcof(object):
    """Interpolation of linear limb-darkening coefficients u(lambda) of Van Hamme (1993)."""

    def __init__(self, Teff=None, logg=None, Z=None, R=None, M=None):
        """Call Fortran program to get u(lambda) for given temperature Teff, gravity logg and metallicity Z [M/H units]"""

        print "Teff = ", Teff  # dbg
        print "logg = ", logg  # dbg
        print "Z = ", Z  # dbg
        print "R = ", R  # dbg
        print "M = ", M  # dbg

        if logg == None:
            G = 6.67408e-11  # SI
            M_Sun = 1.9891e30  # kg
            R_Sun = 6.957e8  # m
            logg = math.log10(G*M*M_Sun/math.pow(R*R_Sun,2) * 1.e2)
            print "logg = ", logg  # dbg
        if Z == None:
            Z = 0.0  # i.e. solar metallicity
            print "Warning: Assuming solar metallicity in limcof.__init__."

        lambda1 = 0.0
        lambda2 = -1.0
        dlambda = 1.0

        try:
            output = subprocess.check_output("%s %.8e %.8e %.8e %.8e %.8e %.8e" % (shsp_limcof, Teff, logg, Z, lambda1, lambda2, dlambda), shell=True)
        except subprocess.CalledProcessError as ex:
            print "Error calling limcof_test! Is limcof.dat symlinked in the current directory?"
            raise ex

        # parse output
        self.wvl = []
        self.u = []
        for line in output.split("\n"):
            if (not re.match("^ *#", line)) and len(line) > 1:
                l = line.split()
                try:
                    self.wvl.append(float(l[0]))
                    self.u.append(float(l[1]))
                except ValueError as ex:
                    print output
                    print ex

        print "wvl = ", self.wvl  # dbg
        print "u = ", self.u  # dbg

    def interp(self, wvl):
        """Interpolate limb darkening coef. for given wavelength wvl"""
        if len(self.wvl) == 0:
            print "Warning: No limb-darkening data for interpolation! Returning 0.0..."
            u = 0.0
            return u

        j = 1
        while self.wvl[j] < wvl and j < len(self.wvl)-1:
            j += 1
        if wvl < self.wvl[j-1]:
            print "Warning: Extrapolation of limb-darkening coef. for lambda = ", wvl, " m"
            u = self.u[j-1]
        elif wvl > self.wvl[j]:
            print "Warning: Extrapolation of limb-darkening coef. for lambda = ", wvl, " m"
            u = self.u[j]
        else:
            u = self.__interp(self.wvl[j-1], self.wvl[j], self.u[j-1], self.u[j], wvl)
        return u

    def __interp(self, x1, x2, y1, y2, x):
        return y1 + (y2-y1)*(x-x1)/(x2-x1)

def main():
    """Test the init and interpolation methods."""
    limcof = Limcof(Teff=13300.0, logg=4.25, Z=0.0)

    Ang = 1.0e-10
    print "# lambda [Ang] & limb-darkening coef. u []"
    for wvl in [1250., 1365, 1430, 1550, 1910, 2460, 2980, 3320, 3670, 4360, 5450 ]:
        u = limcof.interp(wvl*Ang)
        print wvl, u

if __name__ == "__main__":
    main()


