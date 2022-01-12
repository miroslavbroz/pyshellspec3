#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
import pyshellspec3

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
dir_data = 'data'

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
        print(names[i])
        obs.append(pyshellspec3.LCData(filename=os.path.join(dir, names[i]),
                                      passband=pbands[i],
                                      global_error=errs[i],
                                      delimiter=','))

    # data from MIRC
    dir = os.path.join(dir_data, 'mirc2017')
    infile = os.path.join(dir, 'mirc.ascii.lis')
    names = read_if(infile)
    for i in range(0, len(names)):
        print(names[i])
        obs.append(pyshellspec3.IFData(filename=os.path.join(dir, names[i]),
                                      location='chara',
                                      ra=ra, dec=dec, 
                                      format='ascii', weight_vis2=0.5, weight_t3amp=0.5))

    # data from NPOI
    dir = os.path.join(dir_data, 'npoi')
    infile = os.path.join(dir, 'npoi.ascii.lis')
    names = read_if(infile)
    for i in range(0, len(names)):
        print(names[i])
        obs.append(pyshellspec3.IFData(filename=os.path.join(dir, names[i]),
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
        print(names[i])
        obs.append(pyshellspec3.IFData(filename=os.path.join(dir, fnames[i]),
                                      location='chara',
                                      ra=ra, dec=dec,
                                      format='ascii'))
        obs[-1].set_filename(names[i])

    # ICFITS data
    dir = os.path.join(dir_data, 'GOODPHI_6540')
    infile = os.path.join(dir, 'icfits_HA_unique.ascii.lis')
    names = read_if(infile)
    for i in range(0, len(names)):
        print(names[i])
        #obs.append(pyshellspec3.DFData(filename=os.path.join(dir, names[i]), location='chara', ra=ra, dec=dec, format='fits'))
        obs.append(pyshellspec3.DFData(filename=os.path.join(dir, names[i]), location='chara', ra=ra, dec=dec, format='ascii'))

    # SED data
    dir = os.path.join(dir_data, 'burnasev')
    obs.append(pyshellspec3.SEDData(filename=os.path.join(dir,'Sed.dat')))

    # spectral data
    dir = os.path.join(dir_data, 'ondrejov_6670')
    obs.append(pyshellspec3.SPEData(filename=os.path.join(dir,'Spectra.dat')))

    # construct data class
    data = pyshellspec3.Data(obs)

    # construct the model
    central = pyshellspec3.CentralObject()
    companion = pyshellspec3.Companion()
    nebula = pyshellspec3.Nebula()
    jet = pyshellspec3.Jet()
    spot = pyshellspec3.Spot()
    envelope = pyshellspec3.Envelope()
    shell = pyshellspec3.Shell()
    orbit = pyshellspec3.Orbit()
    objs = [central, companion, nebula, jet, spot, envelope, shell, orbit]
    model = pyshellspec3.Model(objects=objs)

    # construct the Interface
    itf = pyshellspec3.Interface(model=model, data=data, ncpu=8, image_size=180,
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
        shellspec_abundance=os.path.join(os.getcwd(), "abundances"),
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

    itf.set_parameter('tempcp'  , value=14334.08417484313    , fitted=True , vmin=12000., vmax=14600.)
    itf.set_parameter('rinnb'   , value=8.65881274820674     , fitted=True , vmin=6.0, vmax=12.0)
    itf.set_parameter('routnb'  , value=31.49478694096131    , fitted=True , vmin=26.0, vmax=35.5)
    itf.set_parameter('hinvnb'  , value=3.457352142937204    , fitted=True , vmin=1.0, vmax=5.0)
    itf.set_parameter('tinvnb'  , value=1.487603388182587    , fitted=True , vmin=1.0, vmax=2.0)
    itf.set_parameter('hwindnb' , value=3.045547309727574    , fitted=True , vmin=3.0, vmax=5.0)
    itf.set_parameter('hcnb'    , value=3.774743840155588    , fitted=True , vmin=1.0, vmax=15.0)
    itf.set_parameter('hvelnb'  , value=1.359190591983155    , fitted=True , vmin=0.1, vmax=5.0)
    itf.set_parameter('vnb'     , value=112.1285187618496    , fitted=True , vmin=0.0, vmax=200.0)
    itf.set_parameter('evelnb'  , value=1.905800947120465    , fitted=True , vmin=0.0, vmax=2.0)
    itf.set_parameter('hshdnb'  , value=4.963444632686908    , fitted=True , vmin=1.0, vmax=5.0)
    itf.set_parameter('tempnb'  , value=30345.1138311543     , fitted=True , vmin=23000., vmax=34000.)
    itf.set_parameter('densnb'  , value=1.205858367874966e-09, fitted=True , vmin=1e-12, vmax=5e-9)
    itf.set_parameter('vtrbnb'  , value=11.23662689479421    , fitted=True , vmin=0.0, vmax=100.0)
    itf.set_parameter('edennb'  , value=-0.5718663190976636  , fitted=True , vmin=-3.0, vmax=-0.5)
    itf.set_parameter('etmpnb'  , value=-0.7294933009989126  , fitted=True , vmin=-1.1, vmax=-0.70)
    itf.set_parameter('ajet'    , value=28.80312785138627    , fitted=True , vmin=2.0, vmax=60.0)
    itf.set_parameter('rinjt'   , value=5.620846116858949    , fitted=True , vmin=3.0, vmax=10.0)
    itf.set_parameter('routjt'  , value=35.85681868576612    , fitted=True , vmin=15.0, vmax=45.0)
    itf.set_parameter('vjt'     , value=675.6149964998236    , fitted=True , vmin=0.0, vmax=1500.0)
    itf.set_parameter('eveljt'  , value=1.265897054885929    , fitted=True , vmin=0.0, vmax=2.0)
    itf.set_parameter('tempjt'  , value=15089.07468847655    , fitted=True , vmin=10000.0, vmax=35000.0)
    itf.set_parameter('densjt'  , value=5.519452158846459e-12, fitted=True , vmin=1.0e-13, vmax=1.0e-10)
    itf.set_parameter('vtrbjt'  , value=66.14241437776984    , fitted=True , vmin=0.0, vmax=300.0)
    itf.set_parameter('rpoljt'  , value=32.96133162700006    , fitted=True , vmin=0.0, vmax=35.0)
    itf.set_parameter('vpoljt'  , value=10.0429239150491     , fitted=True , vmin=0.0, vmax=200.0)
    itf.set_parameter('pangjt'  , value=-70.15492190383611   , fitted=True , vmin=-360.0, vmax=360.0)
    itf.set_parameter('rinsh'   , value=7.41104754334615     , fitted=True , vmin=7.0, vmax=30.0)
    itf.set_parameter('routsh'  , value=72.92941657397402    , fitted=True , vmin=30.0, vmax=120.0)
    itf.set_parameter('vsh'     , value=78.96020166701454    , fitted=True , vmin=0.0, vmax=100.0)
    itf.set_parameter('evelsh'  , value=1.896922837722265    , fitted=True , vmin=0.0, vmax=2.0)
    itf.set_parameter('vysh'    , value=-5.240026676722604   , fitted=True , vmin=-100.0, vmax=100.0)
    itf.set_parameter('tempsh'  , value=5631.482711034088    , fitted=True , vmin=5000.0, vmax=35000.0)
    itf.set_parameter('denssh'  , value=2.863947911386624e-11, fitted=True , vmin=1.0e-13, vmax=5.0e-11)
    itf.set_parameter('vtrbsh'  , value=102.2237200068673    , fitted=True , vmin=0.0, vmax=200.0)
    itf.set_parameter('etmpsh'  , value=-0.00664352380023825 , fitted=True , vmin=-1.0, vmax=0.0)
    itf.set_parameter('dinc'    , value=96.32899543197529    , fitted=True , vmin=91., vmax=97.)
    itf.set_parameter('omega_an', value=254.5553839912507    , fitted=True , vmin=252., vmax=255.)
    itf.set_parameter('dd'      , value=328.3736904956343    , fitted=True , vmin=305., vmax=330.)

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
    itf.compute_chi2(verbose=True)
#    itf.run_fit(fitter='sp_diff_evol', tol=1e-2, maxiter=10000)
#    itf.run_fit(fitter='nlopt_nelder_mead', ftol=1e-6, maxiter=10000)
#    itf.run_fit(fitter='nlopt_sbplx', ftol=1e-3, maxiter=10000)

    itf.write_iterations()
    itf.set_model_to_shellspec()
    itf.write_template('final.in')
    itf.write_model()

    print("Note: fit.py ended successfully.")
    sys.exit(0) 

    files = [] 
    for o in obs:
        files.append(o.get_filename().split('/')[-1])
    for f in list(set(files)):
        itf.plot_comparison(filename=f)

if __name__ == '__main__':
    main()


