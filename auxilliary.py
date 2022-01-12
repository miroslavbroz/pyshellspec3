#!/usr/bin/env python3

import numpy as np
from astropy import units
from scipy.interpolate import splrep
from scipy.interpolate import splev
from scipy.optimize import minimize_scalar

from . import hermite

def add_wing_phase(ph, data):
    """
    Reflects phases
    :param ph:
    :param data:
    :return:
    """

    ind = np.where((ph >= 0.8) & (ph <= 1.0))
    data = np.insert(data, 0, data[ind])
    ph = np.insert(ph, 0, ph[ind] - 1.0)

    ind	= np.where((ph >= 0.0) & (ph <= 0.2))
    data = np.append(data, data[ind])
    ph = np.append(ph, ph[ind] + 1.0)

    return ph, data


def Bisection(func, a, b, args=()):
    """
    Finds root of an equation in an interval <a,b>.
    :param func: the equation
    :param a: lower boundary of the search interval
    :param b: upper boundary of the search interval
    :param args: arguments passed to the equation
    """
    while (abs(a-b) > 1e-10):
        c = (a + b)/2.
        if func(a, *args)*func(c, *args) < 0:
            b = c
        else:
            a = c

    return (a + b)/2

def DFT2d(img, x, y, fu, fv):
    """
    :param img:
    :param x:
    :param y:
    :param fu:
    :param fv:
    :return:
    """
    vis = np.exp(-1j * 2 * np.pi * y * fv).dot(img.dot(np.exp(-1j * 2 * np.pi * x * fu)))

    return vis

def filling_factor_to_rpoint(ffact, l1, component='primary'):
    """
    Computes rpoint for a given filling factor and
    position of the l1 with respect to primary.
    :param ffact:
    :param l1:
    :return:
    """
    if component.lower() == 'primary':
        return ffact * l1
    elif component.lower() == 'secondary':
        return ffact * (1. - l1)
    else:
        raise ValueError('Component can be either \'primary\' or \'secondary\'.')

def get_offset(a, b):
    """
    Computes offset between two arrays.
    :param a:
    :param b:
    :return:
    """

    def eq(offset):
        return np.sum(np.absolute(a - b - offset))

    return minimize_scalar(eq)

def get_mulfac(a, b):
    """
    Computes multiplicative factor between two arrays.
    :param a:
    :param b:
    :return:
    """

    def eq(mulfac):
        return np.sum(np.absolute(a - mulfac*b))

    return minimize_scalar(eq)


def has_element_float(e, l, precision=1e-12):
    """
    Check that element is present in a list.
    
    :param e: element
    :param l: list
    """
    for i in range(len(l)):
        if abs(e - l[i]) < precision:
            return True
    
    return False

def kepler_equation(M, ecc, tol=1e-6):
    """
    Solves the Kepler equation
    :param M:
    :param e:
    :return:
    """
    # nothing to compute
    if ecc < tol:
        return M

    # first iteration
    E0 = M.copy()
    E1 = M + ecc * np.sin(M)
    ind = np.where(np.abs(E0 - E1) > tol)[0]

    # iterate only ver those that do not satisfy condition
    while len(ind) > 0:
        E0[ind] = E1[ind]
        E1[ind] = M[ind] + ecc * np.sin(E0[ind])
        ind = np.where(np.abs(E0 - E1) > tol)[0]

        print(E1)

    return E1




def interpolate1d(xnew, x, y):
    """
    Interpolates in y = f(x) at xnew.
    :param xnew:
    :param x:
    :param y:
    :return:
    """
    tck = splrep(x, y)
    ynew = splev(xnew, tck)

    return ynew


def interpolate1d_hermite(xnew, x, y):
    ynew = np.zeros(len(xnew))
    for i in range(0,len(xnew)):
        ynew[i], ier = hermite.hermite(xnew[i], x, y)
    return ynew


def load_ascii(f, cols, delimiter=None, comment='#'):
    """
    Loads a matrix from a file.
    """

    # reads the file
    ifile = open(f, 'r')
    lines = ifile.readlines()
    ifile.close()
    
    # empty data block
    block = []
    
    # which columns are stored
    colind = [i for i in range(len(cols)) if i != None]
    
    # processes the file
    for l in lines:
        d = l.split(delimiter)
        if d[0].find(comment) > -1:
            continue
        else:
            #print l  # dbg
            #print d
            #print [d[i] for i in colind]
            d = list(map(float, [d[i] for i in colind]))
            block.append(d)
    
    return np.array(block)

def l1_position(q):
    """
    Position through iteration.
    :param q: mass ratio
    """

    def eq(x, q):
        """
        """
        y = - 1. / x ** 2 + q / (1 - x) ** 2 + (1 + q) * x - q
        return y

    L1 = Bisection(eq, 1e-6, 1.0, args=[q])

    return L1

def mean_anomaly(t, period, hjd0, dperiod):
    """
    Computes mean anomaly according to
    Hadrava, P.\ 2004, Publications of the Astronomical Institute, 92, 1
    http://www.asu.cas.cz/~had/fotel.pdf
    Eq. (12)
    :param t:
    :param period:
    :param hjd0:
    :param dperiod:
    :return:
    """
    tau = t - hjd0
    M = 4. * np.pi * tau * (1. / (1. + np.sqrt(1. + 2. * dperiod * tau / period))) / period

    return M


def polar_radius(pot, q, e=0.0, component='primary'):
    """
    Determines polar radiues
    :param pot:
    :param q:
    :param e:
    :return:
    """
    if component.lower() == 'secondary':
        q = 1. / q
        pot = pot * q + 0.5 * (1. / q - 1) * q
    elif component.lower() == 'primary':
        pass
    else:
        raise ValueError('Component can be either \'primary\' or \'secondary\'.')

    def eq(rpole, pot, q, e):
        return pot - potential_polar(rpole, q, e)

    rpole = Bisection(eq, 1e-6, 0.99999, args=(pot, q, e))

    return rpole

def potential_point(rho, q, ecc=0.0, F=1.0, component='primary'):
    """
    Computes potrential at point
    :param rho: relative radius
    :param q:
    :param ecc:
    :param F:
    :param component:
    :return:
    """
    delta = 1-ecc
    if component.lower() == 'primary':
        y = 1./rho + q*(1./((delta-rho)**2)**0.5 - rho/delta**2) + 0.5*(1.+q)*(F*rho)**2
    elif component.lower() == 'secondary':
        q = 1. / q
        y = 1./rho + q*(1./((delta-rho)**2)**0.5 - rho/delta**2) + 0.5*(1.+q)*(F*rho)**2
        y = y / q + 0.5 * (q - 1) / q
    else:
        raise ValueError('Component can be either \'primary\' or \'secondary\'.')
    return y

def potential_polar(rpole, q, e=0.0):
    """
    :param rpole:
    :param q:
    :param e:
    :return:
    """
    delta = 1.0 - e
    pot = 1. / rpole + q * (1. / np.sqrt(delta ** 2 + rpole ** 2))

    return pot

def quad_phase(time, T0, Per0, a=0.0):
    """
    Computes phase from a quadratic ephemeris.
    """
  
    phase = root_quad_eq(a, Per0, T0-time)[0];
    return phase

def semiamplitude(ecc, q, sma, incl, period, component):
    """
    :return:
    """
    if component.lower() == 'primary':
        K = 2 * np.pi * sma * np.sin(incl) * q / (period * np.sqrt(1 - ecc * ecc) * (1 + q))
    elif component.lower() == 'secondary':
        K = -2 * np.pi * sma * np.sin(incl) / (period * np.sqrt(1 - ecc * ecc) * (1 + q))
    else:
        raise ValueError('Component can be either \'primary\' or \'secondary\'.')
    K *= units.meter / units.second
    return K

def radial_velocity(t, ecc, q, sma, incl, omega, period, dperiod, t0, component):
    """
    :return:
    """
    # get mean anomaly
    M = mean_anomaly(t, period, t0, dperiod)

    # get eccentric anomaly
    E = kepler_equation(M, ecc)

    # get true anomaly
    V = true_anomaly(E, ecc)

    K = semiamplitude(ecc, q, sma, incl, period, component)

#    print "K = ", K * units.meter / units.second  # dbg
    return K * (np.cos(V + omega) + ecc * np.cos(omega)) * units.meter / units.second

def root_quad_eq(a, b, c):
    """
    Computes roots of a quadratic Eq.
    """

    if a > 1e-6:
        D = np.sqrt(b**2 - 4*a*c)
        x1 = -b/(2.*a) + D/(2.*a)
        x2 = -b/(2.*a) - D/(2.*a)
    else:
        x1 = -c/b
        x2 = 0.0
    
    return x1, x2

def rotateX(x, y, z, ang):
    """
    Rotates vector x, y, z around the z-axis
    :param x: 1st Cartesian coordinate
    :param y: 2nd Cartesian coordinate
    :param z: 3rd Cartesian coordinate
    :param ang: the rotation angle in radians
    :return:
    """
    X = x
    Y = y * np.cos(ang) + z * np.sin(ang)
    Z = -y * np.sin(ang) + z * np.cos(ang)

    return X, Y, Z

def rotateZ(x, y, z, ang):
    """
    Rotates vector x, y, z around the z-axis
    :param x: 1st Cartesian coordinate
    :param y: 2nd Cartesian coordinate
    :param z: 3rd Cartesian coordinate
    :param ang: the rotation angle in radians
    :return:
    """
    X = x * np.cos(ang) + y * np.sin(ang)
    Y = -x * np.sin(ang) + y * np.cos(ang)
    Z = z

    return X, Y, Z

def true_anomaly(E, ecc, tol=1e-6):
    """
    Computes true anomaly.
    :param E: eccentric anomaly
    :param ecc: eccentricity
    :return:
    """
    if ecc < tol:
        return E
    return 2 * np.arctan2(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E/2.), 1)

def write_file(f, s):
    """
    Writes list of strings to file.
    """

    ofile = open(f, 'w')
    ofile.writelines(s)
    ofile.close()
    
def append_file(f, s):
    """
    Appends list of strings to file.
    """

    ofile = open(f, 'a')
    ofile.writelines(s)
    ofile.close()
    
