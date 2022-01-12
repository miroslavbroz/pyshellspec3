#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np
import pyshellspec

def read_lc(lfile):
    """
    """
    data = np.loadtxt(lfile, dtype=str)

    names = data[:,0]
    pbands = data[:,1]
    errs = data[:,2].astype('float64')

    return names, pbands, errs

def read_if(ifile):
    """

    :param ifile:
    :return:
    """
    return np.loadtxt(ifile, dtype=str)


# constants
ra = (18 + 50 / 60. + 4.79525 / 3600.)/24.*360.
dec = 33 + 21 / 60. + 45.6100 / 3600.
#dir_data = '/home/mira/a/betalyr2/data_20200415_MOURARD'
dir_data = '/scratch/mira/a/betalyr2/data_20200415_MOURARD'

def main():
    """
    Main loop
    :return:
    """
    # read the data
    obs = []

    # first photometry
    dir = os.path.join(dir_data, 'photometry')
    infile = os.path.join(dir, 'lc.lis')
    names, pbands, errs = read_lc(infile)
    for i in range(0, len(names)):
        print names[i]
        obs.append(pyshellspec.LCData(filename=os.path.join(dir, names[i]),
                                      passband=pbands[i],
                                      global_error=errs[i],
                                      delimiter=','))

    # data from MIRC
    dir = os.path.join(dir_data, 'mirc2017')
    infile = os.path.join(dir, 'mirc.ascii.lis')
    names = read_if(infile)
    for i in range(0, len(names)):
        print names[i]
        obs.append(pyshellspec.IFData(filename=os.path.join(dir, names[i]),
                                      location='chara',
                                      ra=ra, dec=dec, 
                                      format='ascii', weight_vis2=0.5, weight_t3amp=0.5))

    # data from NPOI
    dir = os.path.join(dir_data, 'npoi')
    infile = os.path.join(dir, 'npoi.ascii.lis')
    names = read_if(infile)
    for i in range(0, len(names)):
        print names[i]
        obs.append(pyshellspec.IFData(filename=os.path.join(dir, names[i]),
                                      location='npoi',
                                      ra=ra, dec=dec,
                                      format='ascii', exclude_t3amp=True))

    # data from VEGA
    dir = os.path.join(dir_data, 'vega')
    infile = os.path.join(dir, 'vega.ascii.lis')
    d = read_if(infile)
    fnames = d[:, 0]
    names = d[:, 1]
    for i in range(0, len(names)):
        print names[i]
        obs.append(pyshellspec.IFData(filename=os.path.join(dir, fnames[i]),
                                      location='chara',
                                      ra=ra, dec=dec,
                                      format='ascii'))
        obs[-1].set_filename(names[i])

    # ICFITS data
    dir = os.path.join(dir_data, 'ICFITSAPRIL2020_6540')
    infile = os.path.join(dir, 'icfits_HA_unique.ascii.lis')
    names = read_if(infile)
    for i in range(0, len(names)):
        print names[i]
        #obs.append(pyshellspec.DFData(filename=os.path.join(dir, names[i]), location='chara', ra=ra, dec=dec, format='fits'))
        obs.append(pyshellspec.DFData(filename=os.path.join(dir, names[i]), location='chara', ra=ra, dec=dec, format='ascii'))

    # SED data
    dir = os.path.join(dir_data, 'burnasev')
    obs.append(pyshellspec.SEDData(filename=os.path.join(dir,'Sed.dat')))

    # spectral data
    dir = os.path.join(dir_data, 'ondrejov_6670')
    obs.append(pyshellspec.SPEData(filename=os.path.join(dir,'Spectra.dat')))

    # construct data class
    data = pyshellspec.Data(obs)

    # construct the model
    central = pyshellspec.CentralObject()
    companion = pyshellspec.Companion()
    nebula = pyshellspec.Nebula()
    jet = pyshellspec.Jet()
    spot = pyshellspec.Spot()
    envelope = pyshellspec.Envelope()
    shell = pyshellspec.Shell()
    orbit = pyshellspec.Orbit()
    objs = [central, companion, nebula, jet, spot, envelope, shell, orbit]
    model = pyshellspec.Model(objects=objs)

    # construct the Interface
    itf = pyshellspec.Interface(model=model, data=data, ncpu=8, image_size=180,
        if_phase_precision=3,
        df_phase_precision=3,
        lc_phase_precision=2,
        sed_phase_precision=2,
        spe_phase_precision=2, 
        if_ew_precision=7,
        df_ew_precision=10,
        lc_ew_precision=9,
        sed_ew_precision=10,
        spe_ew_precision=10,
        shellspec_template="template.in",
        shellspec_abundance="/scratch/mira/a/betalyr2/fitting_shell8_20200415/abundances",
        use_offset=True, use_differential=True,
        exclude_visphi=False, dry_run=False)

    # set grid resolution
    itf.set_parameter('rmdfx1', value=-81.0001)
    itf.set_parameter('rmdfx2', value=+81.0001)
    itf.set_parameter('rmdfy1', value=-81.0001)
    itf.set_parameter('rmdfy2', value=+81.0001)
    itf.set_parameter('rmdfz1', value=-81.0001)
    itf.set_parameter('rmdfz2', value=+81.0001)
    itf.set_parameter('stepfx', value=2.0)
    itf.set_parameter('stepfy', value=2.0)
    itf.set_parameter('stepfz', value=2.0)

    itf.set_parameter('rmdx1',  value=-81.0)
    itf.set_parameter('rmdx2',  value=+81.0)
    itf.set_parameter('rmdy1',  value=-81.0)
    itf.set_parameter('rmdy2',  value=+81.0)
    itf.set_parameter('rmdz1',  value=-81.0)
    itf.set_parameter('rmdz2',  value=+81.0)
    itf.set_parameter('stepx', value=2.0)
    itf.set_parameter('stepy', value=2.0)
    itf.set_parameter('stepz', value=2.0)

    # set fitted parameter
#    itf.set_parameter('istar', value=0)  # dbg
#    itf.set_parameter('icomp', value=0)  # dbg
    itf.set_parameter('inebl', value=1)  # nebula
#    itf.set_parameter('inebl', value=0)  # no nebula
    itf.set_parameter('itnb', value=3)  # power-law
    itf.set_parameter('ijet', value=2)
#    itf.set_parameter('ijet', value=0)  # no jet
    itf.set_parameter('ispot', value=0)  # no spot
#    itf.set_parameter('ienv', value=3) # evelope Roche contact
    itf.set_parameter('ienv', value=0) # no envelope
    itf.set_parameter('ishell', value=3)

    itf.set_parameter('ichemc', value=1)
    itf.set_parameter('ielnd', value=1)
    itf.set_parameter('ithom', value=1) # scattering
    itf.set_parameter('irayl', value=1)
    itf.set_parameter('imie', value=0)
    itf.set_parameter('imiepf', value=0)
    itf.set_parameter('iline', value=1)
    itf.set_parameter('iinvnb', value=1)  # inversion of T
    itf.set_parameter('ivelnb', value=1)  # radial wind
    itf.set_parameter('ishdnb', value=1)  # shadow

    itf.set_parameter('tempcp'  , value=14347.75381356163    , fitted=True , vmin=12000., vmax=14600.)
    itf.set_parameter('rinnb'   , value=8.650280108545585    , fitted=True , vmin=6.0, vmax=12.0)
    itf.set_parameter('routnb'  , value=30.84997931234177    , fitted=True , vmin=26.0, vmax=35.5)
    itf.set_parameter('hinvnb'  , value=3.243474108298643    , fitted=True , vmin=1.0, vmax=5.0)
    itf.set_parameter('tinvnb'  , value=1.590521297561623    , fitted=True , vmin=1.0, vmax=2.0)
    itf.set_parameter('hwindnb' , value=3.04272005472429     , fitted=True , vmin=3.0, vmax=5.0)
    itf.set_parameter('hcnb'    , value=4.168576968569478    , fitted=True , vmin=1.0, vmax=15.0)
    itf.set_parameter('hvelnb'  , value=1.528437310641364    , fitted=True , vmin=0.1, vmax=5.0)
    itf.set_parameter('vnb'     , value=114.2229824701649    , fitted=True , vmin=0.0, vmax=200.0)
    itf.set_parameter('evelnb'  , value=1.946709742926056    , fitted=True , vmin=0.0, vmax=2.0)
    itf.set_parameter('hshdnb'  , value=4.688917840654381    , fitted=True , vmin=1.0, vmax=5.0)
    itf.set_parameter('tempnb'  , value=30619.50770794726    , fitted=True , vmin=23000., vmax=34000.)
    itf.set_parameter('densnb'  , value=9.069408113108799e-10, fitted=True , vmin=1e-12, vmax=5e-9)
    itf.set_parameter('vtrbnb'  , value=12.24015734155999    , fitted=True , vmin=0.0, vmax=100.0)
    itf.set_parameter('edennb'  , value=-0.6760859260358993  , fitted=True , vmin=-3.0, vmax=-0.5)
    itf.set_parameter('etmpnb'  , value=-0.7087622446736188  , fitted=True , vmin=-1.1, vmax=-0.70)
    itf.set_parameter('ajet'    , value=28.70592426352512    , fitted=True , vmin=2.0, vmax=60.0)
    itf.set_parameter('rinjt'   , value=5.372699547592712    , fitted=True , vmin=3.0, vmax=10.0)
    itf.set_parameter('routjt'  , value=33.29831283973326    , fitted=True , vmin=15.0, vmax=45.0)
    itf.set_parameter('vjt'     , value=708.6812052715633    , fitted=True , vmin=0.0, vmax=1500.0)
    itf.set_parameter('eveljt'  , value=1.307882435170993    , fitted=True , vmin=0.0, vmax=2.0)
    itf.set_parameter('tempjt'  , value=15307.40334685057    , fitted=True , vmin=10000.0, vmax=35000.0)
    itf.set_parameter('densjt'  , value=4.99657349796472e-12 , fitted=True , vmin=1.0e-13, vmax=1.0e-10)
    itf.set_parameter('vtrbjt'  , value=65.1959012296722     , fitted=True , vmin=0.0, vmax=300.0)
    itf.set_parameter('rpoljt'  , value=32.4122061349182     , fitted=True , vmin=0.0, vmax=35.0)
    itf.set_parameter('vpoljt'  , value=8.506349520551186    , fitted=True , vmin=0.0, vmax=200.0)
    itf.set_parameter('pangjt'  , value=-43.79714298832457   , fitted=True , vmin=-360.0, vmax=360.0)
    itf.set_parameter('rinsh'   , value=7.297210039265366    , fitted=True , vmin=7.0, vmax=30.0)
    itf.set_parameter('routsh'  , value=73.64725970096735    , fitted=True , vmin=30.0, vmax=120.0)
    itf.set_parameter('vsh'     , value=76.74011548543764    , fitted=True , vmin=0.0, vmax=100.0)
    itf.set_parameter('evelsh'  , value=1.950760467077645    , fitted=True , vmin=0.0, vmax=2.0)
    itf.set_parameter('vysh'    , value=-25.65998874824725   , fitted=True , vmin=-100.0, vmax=100.0)
    itf.set_parameter('tempsh'  , value=5660.61091036055     , fitted=True , vmin=5000.0, vmax=35000.0)
    itf.set_parameter('denssh'  , value=2.773518035640456e-11, fitted=True , vmin=1.0e-13, vmax=5.0e-11)
    itf.set_parameter('vtrbsh'  , value=99.58202753289451    , fitted=True , vmin=0.0, vmax=200.0)
    itf.set_parameter('etmpsh'  , value=-0.005836973987931269, fitted=True , vmin=-1.0, vmax=0.0)
    itf.set_parameter('dinc'    , value=96.21590901592278    , fitted=True , vmin=91., vmax=97.)
    itf.set_parameter('omega_an', value=254.557924277449     , fitted=True , vmin=252., vmax=255.)
    itf.set_parameter('dd'      , value=330                  , fitted=True , vmin=305., vmax=330.)

    itf.set_parameter('rcjt'    , value=1.9                  , fitted=False, vmin=1.0, vmax=5.0)
    itf.set_parameter('asymjt'  , value=0.0                  , fitted=False, vmin=-1.0, vmax=1.0)

    itf.set_parameter('aneb'    , value=5.0                  , fitted=False, vmin=3.0, vmax=12.0)
    itf.set_parameter('asini'   , value=58.19                , fitted=False, vmin=53., vmax=63.)

    itf.set_parameter('emstar'  , value=13.048               , fitted=False, vmin=10.0, vmax=16.0)
    itf.set_parameter('q'       , value=0.2230               , fitted=False, vmin=0.10, vmax=0.30)

    itf.set_parameter('rcsh'    , value=6.0                  , fitted=False, vmin=6.0, vmax=6.0)

#    itf.set_parameter('emen'        , value=13.048            , fitted=False, vmin=0., vmax=10000.)
#    itf.set_parameter('qqen'        , value=0.2230            , fitted=False, vmin=0., vmax=1.)
#    itf.set_parameter('aen'         , value=58.349            , fitted=False, vmin=0., vmax=10000.)
#    itf.set_parameter('ffen'        , value=2.0               , fitted=False, vmin=1., vmax=2.)
#    itf.set_parameter('hen'         , value=100.0             , fitted=False, vmin=0., vmax=10000.)
#    itf.set_parameter('tempen'      , value=6000.0            , fitted=False, vmin=5500., vmax=13500.)
#    itf.set_parameter('densen'      , value=1.0e-11           , fitted=False, vmin=0., vmax=10000.)
#    itf.set_parameter('vtrben'      , value=0.0               , fitted=False, vmin=0., vmax=10000.)

#    itf.set_parameter('rsp'     , value=5.116504353388351    , fitted=False, vmin=4., vmax=2.0*6.36485904820971)
#    itf.set_parameter('tempsp'  , value=7131.412355683008    , fitted=False, vmin=6000., vmax=25000.)
#    itf.set_parameter('denssp'  , value=8.668325439973066e-09, fitted=False, vmin=1.e-11, vmax=1.e-7)
#    itf.set_parameter('rpolsp'  , value=29.82593138502884    , fitted=False, vmin=20., vmax=31.09436071132663)
#    itf.set_parameter('vpolsp'  , value=0.0                  , fitted=False, vmin=0.0, vmax=200.0)
#    itf.set_parameter('pangsp'  , value=4.683571066671123    , fitted=False, vmin=-180., vmax=180.)

    itf.set_parameter('vgamma'  , value=-18.0                , fitted=False, vmin=-100., vmax=100.)

    # compute one/fit
#    itf.compute_chi2(verbose=True)
#    itf.run_fit(fitter='sp_diff_evol', tol=1e-2, maxiter=10000)
    itf.run_fit(fitter='nlopt_nelder_mead', ftol=1e-6, maxiter=10000)
#    itf.run_fit(fitter='nlopt_sbplx', ftol=1e-3, maxiter=10000)

    itf.write_iterations()
    itf.set_model_to_shellspec()
    itf.write_template('final.in')
    itf.write_model()

    print "Note: fit.py ended successfully."
    sys.exit(0) 

    files = [] 
    for o in obs:
        files.append(o.get_filename().split('/')[-1])
    for f in list(set(files)):
        itf.plot_comparison(filename=f)

if __name__ == '__main__':
    main()


