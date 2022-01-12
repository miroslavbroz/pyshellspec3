!    *******************************************************************
!
      subroutine partfn (idno,theta,u,eplow)
        IMPLICIT REAL*8 (A-H, O-Z)
!
!     calculates the partition functions for element idno
!
!     completely new version by B. Smalley October 1992
!
!     braket function implemented by B. Smalley April 1993
!
!      real g(92,5)
!      real dele(92,5)
 	dimension g(92,5),dele(92,5)
!-----------------------------------------------------------------------
!      include 'params.inc'
!       begin 'params.inc'
!
!i    melmln         - size of element abreviations
!
      integer    melmln
      parameter (melmln =     2)
!
!i    melmns         - number of elements
!
      integer    melmns
      parameter (melmns =    92)
!
!i    mlayer         - maximum number of model depths
!
      integer    mlayer
      parameter (mlayer =   200)
!
!i    mlines         - maximum number of lines in line data buffer
!
      integer    mlines
      parameter (mlines =  1000)
!
!i    msynth         - size of "synth" arrays
!
      integer    msynth
      parameter (msynth = 10000)
!
!i    mbroad         - size of "broad" arrays
!
      integer    mbroad
      parameter (mbroad = 10000)
!
!i    mtitle         - size of character strings for titles and filenames
!
      integer    mtitle
      parameter (mtitle =    80)
!
!i    mdipso         - size of "dipso" arrays
!
      integer    mdipso
      parameter (mdipso = 100000)
!
!i    mfpram         - maximum number of free parameters for chifit
!
      integer    mfpram
      parameter (mfpram =     6)
!
!i    mchmax         -
!
      integer    mchmax
      parameter (mchmax =    50)
!
!i    mdskin         - file input unit
!
      integer    mdskin
      parameter (mdskin =     2)
!
!i    mdskot         - file output unit
!
      integer    mdskot
      parameter (mdskot =     3)
!
!i    mchlin         - default input unit
!
      integer    mchlin
      parameter (mchlin =     5)
!
!i    mchlot         - default output unit
!
      integer    mchlot
      parameter (mchlot =     6)
!
!i    mchses         - unit for session log file
!
      integer    mchses
      parameter (mchses =    10)
!
!r    alge10         - value of log(10)
!
      real alge10
      parameter (alge10 = 2.302585093)
!
!r    around         - value by which to round-up errors
!
      real around
      parameter (around = 0.01)
!
!i    minstp         - maximum number of points in instrumental profile
!
      integer minstp
      parameter (minstp = 1000)
!	end 'params.inc' 
!-----------------------------------------------------------------------
!     include 'commons.inc'
!     begin 'commons.inc'
!
!     /data/           element data
!
!c    elmchr(melmns) - element symbols
!r    gmabds(melmns) - element number fractions  n(el)/n(h)
!r    atmass(melmns) - atomic masses (amu)
!r    aipot1(melmns) - 1st ionization potention (ev)
!r    aipot2(melmns) - 2nd ionization potention (ev)
!r    aipot3(melmns) - 3rd ionization potention (ev)
!r    aipot4(melmns) - 4th ionization potention (ev)
!r    aipot5(melmns) - 5th ionization potention (ev)
!
      character*(melmln) elmchr
      real gmabds,atmass,aipot1,aipot2,aipot3,aipot4,aipot5
      common /atdata/elmchr(melmns),gmabds(melmns),atmass(melmns),      &
     &               aipot1(melmns),aipot2(melmns),aipot3(melmns),      &
     &               aipot4(melmns),aipot5(melmns)
!
!     /model/          model deck parameters
!
!r    altau0(mlayer) - log(tau)
!r    aatau0(mlayer) - tau
!r    atheta(mlayer) - theta
!r    apgass(mlayer) - gas pressure
!r    apelec(mlayer) - electron pressure
!r    akapp0(mlayer) - kappa
!r    atelec(mlayer) - electron temperature
!r    anelec(mlayer) - electron density
!i    nlayer         - number of layers in model atmosphere
!c    cmodel         - name of model either filename or teff and log g
!
      character*(mtitle) cmodel
      real altau0,aatau0,atheta,apgass,apelec,akapp0,atelec,anelec
      integer nlayer
      common /model/ altau0(mlayer), aatau0(mlayer), atheta(mlayer),    &
     &               apgass(mlayer), apelec(mlayer), akapp0(mlayer),    &
     &               atelec(mlayer), anelec(mlayer),                    &
     &               nlayer,cmodel
!
!     /xlines/         line buffer
!
!d    xwaves(mlines) - wavelength (a)
!i    idents(mlines) - element number
!i    ionsts(mlines) - ionization stage (ev)
!r    expots(mlines) - lower level excitation potential (ev)
!r    gfvals(mlines) - log gf
!r    g1vals(mlines) - radiative damping constant
!r    g2vals(mlines) - van der vaals damping constant
!r    g3vals(mlines) - stark damping constant
!r    abunds(mlines) - log a
!r    eqwids(mlines) - equivalent width (a)
!c    catref(mlines) - atomic data reference
!i    nlines         - number of lines in buffer
!
      double precision xwaves
      real expots,gfvals,g1vals,g2vals,g3vals,abunds,eqwids
      integer idents,ionsts,nlines
      character*4        catref
      common /xlines/ xwaves(mlines),idents(mlines),ionsts(mlines),     &
     &                expots(mlines),gfvals(mlines),g1vals(mlines),     &
     &                g2vals(mlines),g3vals(mlines),abunds(mlines),     &
     &                eqwids(mlines),nlines,catref(mlines)
!
!     /radial/         cloud data
!
!r    ataucl(mlayer) - log tau values of cloud
!r    afclou(mlayer) - d(log a) values for each log tau value
!i    ncloud         - number of cloud points
!r    fcloud(mlayer) - abundance scaling factors for each altauo value
!i    identc         - identification of cloud element
!
      real ataucl,afclou,fcloud
      integer ncloud, identc
      common /radial/ ataucl(mlayer),afclou(mlayer),fcloud(mlayer),     &
     &                ncloud, identc
!
!     /dumb/           dumps
!
!i    idumps(10)     - dump flag array
!
      integer idumps
      common /dumb/ idumps(10)
!
!     /freq/
!
!r    xwavlh         -
!r    aehvkt(mlayer) -
!r    akappa(mlayer) -
!r    arhonu(mlayer) -
!r    ataunu(mlayer) -
!
      real xwavlh,aehvkt,akappa,arhonu,ataunu
      common /freq/ xwavlh,aehvkt(mlayer),akappa(mlayer),               &
     &                     arhonu(mlayer),ataunu(mlayer)
!
!     /line/           current line data
!
!i    identl         - element atomic number
!i    ionstl         - ionization stage
!r    xwavel         - wavelength (a)
!r    expotl         - lower level excitation potential (ev)
!r    gfvall         - log gf
!r    g1vall         - radiative damping constant
!r    g2vall         - van der waals damping constant
!r    g3vall         - stark damping constant
!r    eqwidl         - equivalent width (a)
!r    abundl         - log a
!i    ifewid         -
!
      real xwavel,expotl,gfvall,g1vall,g2vall,g3vall,eqwidl,abundl
      integer identl,ionstl,ifewid
      common /line/ identl,ionstl,xwavel,expotl,gfvall,g1vall,          &
     &              g2vall,g3vall,eqwidl,abundl,ifewid
!
!     /oldval/         last line used by sahan
!
!i    idold          - last element atomic number
!i    ionold         - last ionization stage
!r    epold          - last lower level excitation potential
!
      real epold
      integer idold,ionold
      common /oldval/idold,ionold,epold
!
!     /linop/          line opacity
!
!r    alinop(mlayer) -
!r    alfhln(mlayer) - hydrogen line opacity
!r    alphal(mlayer) -
!r    avogta(mlayer) -
!r    avogtv(mlayer) -
!r    agamma(mlayer) -
!r    axisqu(mlayer) -
!i    modesy         - synthesis mode (-1,0,1)          
!
      real alinop,alphal,avogta,avogtv,agamma,axisqu,alfhln
      integer modesy
      common /linop/ alinop(mlayer),alphal(mlayer),avogta(mlayer),      &
     &               avogtv(mlayer),agamma(mlayer),axisqu(mlayer),      &
     &               alfhln(mlayer),modesy
!
!     /sauce/
!
!r    abnufn(mlayer) -
!r    asnufn(mlayer) -
!r    ajmsnu(mlayer) -
!r    ajaynu(mlayer) -
!i    ifscat         -
!i    iflime         -
!
      real abnufn,asnufn,ajmsnu,ajaynu
      integer ifscat,iflime
      common /sauce/ abnufn(mlayer),asnufn(mlayer),ajmsnu(mlayer),      &
     &               ajaynu(mlayer),ifscat        ,iflime
!
!     /saha/
!
!r    anrs(mlayer)   -
!r    anrnel(mlayer) -
!r    potlow(mlayer) -
!
      real anrs,anrnel,potlow
      common /saha/ anrs(mlayer),anrnel(mlayer),potlow(mlayer)
!
!     /fract/          abundance fractions
!
!r    frah0(mlayer)  - neutral hydrogen (h i)
!r    frah2(mlayer)  - protons (h ii)
!r    framg0(mlayer) - mg i
!r    frasi0(mlayer) - si i
!r    fraal0(mlayer) - al i
!r    frac0(mlayer)  - c i
!r    framg1(mlayer) - mg ii
!r    frasi1(mlayer) - si ii
!r    fraal1(mlayer) - al ii
!r    frac1(mlayer)  - c ii
!r    frane(mlayer)  - electrons
!r    frahe0(mlayer) - neutral helium (he i)
!r    frahe2(mlayer) - he ii
!r    frahe3(mlayer) - he iii
!
      real frah0,frah2,framg0,frasi0,fraal0,frac0
      real framg1,frasi1,fraal1,frac1
      real frane,frahe0,frahe2,frahe3
      common /fract/ frah0(mlayer) ,frah2(mlayer) ,framg0(mlayer),      &
     &               frasi0(mlayer),fraal0(mlayer),frac0(mlayer) ,      &
     &               framg1(mlayer),frasi1(mlayer),fraal1(mlayer),      &
     &               frac1(mlayer) ,frane(mlayer) ,frahe0(mlayer),      &
     &               frahe2(mlayer),frahe3(mlayer)
!
!     /kappas/         continuous opacity
!
!r    alfh(mlayer)   -
!r    alfhm(mlayer)  -
!r    alfh2p(mlayer) -
!r    alfh2m(mlayer) -
!r    alfmg0(mlayer) -
!r    alfsi0(mlayer) -
!r    alfal0(mlayer) -
!r    alfc0(mlayer)  -
!r    alfmg1(mlayer) -
!r    alfsi1(mlayer) -
!r    alfal1(mlayer) -
!r    alfc1(mlayer)  -
!r    alfhe(mlayer)  -
!r    alfhep(mlayer) -
!r    alfhem(mlayer) -
!
      real alfh,alfhm,alfh2p,alfh2m
      real alfmg0,alfsi0,alfal0,alfc0,alfmg1,alfsi1,alfal1,alfc1
      real alfhe,alfhep,alfhem
      common /kappas/ alfh(mlayer)  ,alfhm(mlayer) ,alfh2p(mlayer),     &
     &                alfh2m(mlayer),alfmg0(mlayer),alfsi0(mlayer),     &
     &                alfal0(mlayer),alfc0(mlayer) ,alfmg1(mlayer),     &
     &                alfsi1(mlayer),alfal1(mlayer),alfc1(mlayer) ,     &
     &                alfhe(mlayer) ,alfhep(mlayer),alfhem(mlayer)
!
!     /finess/         quality values
!
!i    nemerg         -
!i    ntau           -
!i    nsaha          -
!
      integer nemerg, ntau, nsaha
      common /finess/ nemerg, ntau, nsaha
!
!     /tnew/
!
!i    newtmp         -
!r    amubar         - sum of n(el)/n(h) * atomic number
!r    arabum         - sum of n(el)/n(h)
!
      real amubar,arabum
      integer newtmp
      common /tnew  / newtmp, amubar, arabum
!
!     /mapem/          any use?
!
!i    ifcor(melmns)  -
!i    indexx(melmns) -
!
      integer ifcor,indexx
      common /mapem/ ifcor(melmns),indexx(melmns)
!
!     /mapex/          any use?
!
!i    ifcox(melmns)  -
!
      integer ifcox
      common /mapex/ ifcox(melmns)
!
!     /disk/           value of mu
!
!r    diskmu         - value of mu
!
      real diskmu
      common /disk/ diskmu
!
!     /spcbrd/         broadened synthetic spectrum
!
!r    xbroad(msynth) - array containing broadened spectrum wavelengths
!r    ybroad(msynth) - array containing broadened spectrum fluxes
!i    nbroad         - number of points in broadened spectrum
!i    cbroad         - title of broadened spectrum
!
      character*(mtitle) cbroad
      real xbroad, ybroad
      integer nbroad
      common /spcbrd/ xbroad(msynth), ybroad(msynth), nbroad, cbroad
!
!     /spcsyn/         synthetic spectrum
!
!c    csynth         - title of synthetic spectrum
!r    xsynth(msynth) - array containing synthetic spectrum wavelengths
!r    ysynth(msynth) - array containing synthetic spectrum fluxes
!r    zsynth(msynth) - array containing optical depth
!i    nsynth         - number of points in synthetic spectrum
!
      character*(mtitle) csynth
      real xsynth, ysynth, zsynth
      integer nsynth
      common /spcsyn/ xsynth(msynth), ysynth(msynth), zsynth(msynth),   &
     &                nsynth, csynth
!
!     /spcdip/         observational spectrum
!
!c    cdipso         - title of observational spectrum
!r    xdipso(mdipso) - array containing observed spectrum wavelengths
!r    ydipso(mdipso) - array containing observed spectrum fluxes
!i    ndipso         - number of points in observed spectrum
!r    radvel         - radial velocity of observed spectrum
!
      character*(mtitle) cdipso
      real xdipso, ydipso
      integer ndipso
      real radvel
      common /spcdip/ xdipso(mdipso), ydipso(mdipso),                   &
     &                ndipso, radvel, cdipso
!
!     /sins/           synthesis grid range
!
!d    xwavlo         - lower wavelength limit
!d    xwavhi         - upper wavelength limit
!d    xdelta         - step size
!
      double precision   xwavlo, xwavhi, xdelta
      common /sins  / xwavlo, xwavhi, xdelta
!
!     /syns/
!
!r    vesini         - rotational velocity
!r    gauwth         - instrumental broadening width
!r    xinstp         - delta lambda of empirical instrumental profile
!r    yinstp         - response of empirical instrumental profile
!i    ninstp         - number of points in empirical instrumental profile
!
      real vesini,gauwth
      real xinstp, yinstp
      integer ninstp
      common /syns  / vesini, gauwth,                                   &
     &                xinstp(minstp),yinstp(minstp),ninstp
!
!     /chicom/
!
!i    nfpram         - number of free parameters
!c    cfpram         - element name
!r    afpram(mfpram) - element starting abundance
!r    ywfunc(msynth) - weighting function
!i    nwfunc         - number of points in weighting function
!i    lineno         - line number (*kcs change 19-06-95*)
!
      character*(melmln) cfpram(mfpram)
      real afpram, ywfunc
      integer nfpram, nwfunc, lineno
      common /chicom/ ywfunc(msynth), nwfunc, lineno,                   &
     &                nfpram, afpram(mfpram), cfpram
!
!     /smear/          microturbulence
!
!r    vmicro         - microturbulent velocity (km/s)
!
      real vmicro
      common /smear/ vmicro
!
!     /xmol/           any use?
!
!r    fatom(mlayer,20,2) -
!r    fmol(mlayer,20)    -
!
      real fatom,fmol
      common /xmol/ fatom(mlayer,20,2),fmol(mlayer,20)
!
!     /state/
!
!r    arho(mlayer)   - density
!
      real arho
      common /state/ arho(mlayer)
!
!     /vgrad/
!
!r    gradv(mlayer)  -
!i    nograd         -
!
      real gradv
      integer nograd
      common /vgrad/ gradv(mlayer),nograd
!c
!c     /taunox/
!c
!ci    noxtau         -
!c
!      integer noxtau
!      common /taunox/ noxtau
!c
!
!     /result/         results from exact calculations
!
!r    resabs(mlines) - abundance results
!r    resews(mlines) - equivalent width results
!r    abvars(mlines) - abundance variances from errors
!r    covars(mcovar,mcovar) - abundance covariances from errors
!i    nreslt         - flag to indicate which type of results available
!l    labvar         - logical flag to indicate that errors are available
!
      integer mcovar
      parameter (mcovar = mlines)
!
      real resabs,resews,abvars,covars
      integer nreslt
      logical labvar
      common /result/resabs(mlines),resews(mlines),abvars(mlines),      &
     &                covars(mcovar,mcovar),nreslt,labvar
!
!     /prohe1/ and /pro447/ - he1 stark profile data from
!     Barnard, Cooper and Smith JQSRT 14, 1025, 1974 (for 4471)
!     and Shamey, unpublished PhD thesis, 1969 (for other lines)
!
      integer mnehe1, mlnhe1, mtehe1, mdlhe1, mne447, mte447, mdl447
!
      parameter (mnehe1 =  8)
      parameter (mlnhe1 =  3)
      parameter (mtehe1 =  4)
      parameter (mdlhe1 = 50)
      parameter (mne447 =  7)
      parameter (mte447 =  4)
      parameter (mdl447 = 80)
!
      real prfhe1,dlmhe1,xnehe1
      integer nwlam
      real prf447,dlm447,xne447
      real ti0,ti1,ti2
      integer jt
!
      common/prohe1/prfhe1(mdlhe1,mtehe1,mnehe1,mlnhe1),                &
     &              dlmhe1(mdlhe1,mnehe1,mlnhe1),                       &
     &              xnehe1(mnehe1), nwlam(mnehe1,mtehe1)
!
      common/pro447/prf447(mdl447,mte447,mne447),                       &
     &              dlm447(mdl447,mne447),xne447(mne447)
!
      common/proint/ti0(mlayer),ti1(mlayer),ti2(mlayer),jt(mlayer)
!
!     end 'commons.inc'
!-----------------------------------------------------------------------
!
      integer idno,i,j,k
!      real u(5),theta,eplow
	dimension u(5)
!
!     partfn.inc contains data for partition functions
!
!     note the use of 0 index to allow for ground state u values
!
      real uarray(92,5,0:22)
!-----------------------------------------------------------------------
!      include 'partfn.inc'
!      begin 'partfn.inc'
!     H  I    MMD            
      DATA (UARRAY( 1,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69316, &
     & 0.69327,0.69382,0.69572,0.70096,0.71285,0.73609,0.77620,0.83834, &
     & 0.92595,1.03958,1.17666,1.33230,1.50061,1.67583,   1235.9091    /
!     H  II   U=1 assumed
      DATA (UARRAY( 1,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     H  III  U=1 assumed
      DATA (UARRAY( 1,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     H  IV   U=1 assumed
      DATA (UARRAY( 1,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     H  V    U=1 assumed
      DATA (UARRAY( 1,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     He I    D&F            
      DATA (UARRAY( 2,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00017,0.00158,0.00693,0.02138,0.05521,0.11795,0.22423,0.37256, &
     & 0.57051,0.79720,1.04100,1.29148,1.54076,1.78359,   2234.6365    /
!     He II   D&F            
      DATA (UARRAY( 2,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315, &
     & 0.69365,0.69674,0.70715,0.73591,0.79924,0.91804,1.10345,1.34767, &
     & 1.62966,1.94367,2.25628,2.55889,2.84664,3.11725,   4945.7271    /
!     He III  U=1 assumed
      DATA (UARRAY( 2,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     He IV   U=1 assumed
      DATA (UARRAY( 2,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     He V    U=1 assumed
      DATA (UARRAY( 2,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     Li I    MMD            
      DATA (UARRAY( 3,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69320,0.69362,0.69519,0.69895, &
     & 0.70594,0.71707,0.73337,0.75623,0.78754,0.82965,0.88492,0.95525, &
     & 1.04153,1.14332,1.25894,1.38575,1.52064,1.66049,    490.0000    /
!     Li II   D&F            
      DATA (UARRAY( 3,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00011, &
     & 0.00368,0.02090,0.08201,0.23962,0.52895,0.94495,1.42905,1.91544, &
     & 2.37479,2.81410,3.20943,3.56732,3.89271,4.18951,   6874.4546    /
!     Li III  D&F            
      DATA (UARRAY( 3,3,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69316, &
     & 0.69490,0.70416,0.73799,0.82785,1.00799,1.30560,1.66404,1.94163, &
     & 2.15867,2.33690,2.48812,2.61944,2.73551,2.83949,  11131.9092    /
!     Li IV   U=1 assumed
      DATA (UARRAY( 3,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     Li V    U=1 assumed
      DATA (UARRAY( 3,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     Be I    MMD            
      DATA (UARRAY( 4,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00004,0.00080,0.00514,0.01775,0.04272, &
     & 0.08194,0.13515,0.20094,0.27792,0.36523,0.46261,0.57004,0.68729, &
     & 0.81366,0.94781,1.08796,1.23205,1.37798,1.52388,    847.2727    /
!     Be II   MMD            
      DATA (UARRAY( 4,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69343,0.69605,0.70472,0.72209,0.74852, &
     & 0.78273,0.82296,0.86775,0.91633,0.96866,1.02527,1.08694,1.15445, &
     & 1.22837,1.30882,1.39553,1.48779,1.58461,1.68483,   1655.0909    /
!     Be III  AEL-Kurucz     
      DATA (UARRAY( 4,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00002, &
     & 0.00997,0.02010,0.13220,0.23427,0.64666,0.93982,1.53561,1.90816, &
     & 2.45419,2.80665,3.25257,3.56078,3.92038,4.18432,  13990.2725    /
!     Be IV   AEL-Kurucz     
      DATA (UARRAY( 4,4,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315, &
     & 0.69340,0.69365,0.69640,0.69917,0.71250,0.72575,0.76561,0.80412, &
     & 0.88616,0.96224,1.08945,1.20261,1.36109,1.49786,  19792.0918    /
!     Be V    U=1 assumed
      DATA (UARRAY( 4,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   1000.0000    /
!     B  I    MMD            
      DATA (UARRAY( 5,1,K),K=0,22) /                                    &
     & 0.69315,1.77246,1.78209,1.78530,1.78692,1.78792,1.78874,1.78977, &
     & 1.79147,1.79433,1.79895,1.80606,1.81659,1.83164,1.85245,1.88028, &
     & 1.91621,1.96101,2.01500,2.07801,2.14939,2.22809,    754.1818    /
!     B  II   MMD            
      DATA (UARRAY( 5,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00007,0.00355,0.02496,0.07890,0.16615,0.27702, &
     & 0.40008,0.52688,0.65268,0.77555,0.89528,1.01251,1.12813,1.24297, &
     & 1.35757,1.47217,1.58670,1.70089,1.81429,1.92641,   2286.2727    /
!     B  III  MMD            
      DATA (UARRAY( 5,3,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69327,0.69672,0.71222,0.74466,0.79176,0.84819, &
     & 0.90893,0.97054,1.03113,1.08994,1.14697,1.20257,1.25728,1.31162, &
     & 1.36603,1.42084,1.47624,1.53230,1.58898,1.64617,   3447.2727    /
!     B  IV   AEL            
      DATA (UARRAY( 5,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00002, &
     & 0.00012,0.00054,0.00185,0.00512,0.01204,0.02484,0.04613,0.07851, &
     & 0.12412,0.18425,0.25906,0.34756,0.44782,0.55736,  23578.7266    /
!     B  V    AEL            
      DATA (UARRAY( 5,5,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315, &
     & 0.69318,0.69328,0.69358,0.69425,0.69558,0.69791,0.70164,0.70721, &
     & 0.71501,0.72542,0.73874,0.75517,0.77482,0.79772,  30929.0918    /
!     C  I    MMD            
      DATA (UARRAY( 6,1,K),K=0,22) /                                    &
     & 0.00000,2.15589,2.17693,2.18813,2.20238,2.22058,2.24122,2.26290, &
     & 2.28481,2.30670,2.32893,2.35253,2.37928,2.41159,2.45228,2.50413, &
     & 2.56933,2.64907,2.74325,2.85054,2.96868,3.09493,   1023.2727    /
!     C  II   MMD            
      DATA (UARRAY( 6,2,K),K=0,22) /                                    &
     & 0.69315,1.76450,1.77808,1.78281,1.78678,1.79387,1.80660,1.82589, &
     & 1.85151,1.88266,1.91837,1.95783,2.00045,2.04591,2.09405,2.14484, &
     & 2.19827,2.25430,2.31284,2.37371,2.43665,2.50132,   2216.0000    /
!     C  III  D&F            
      DATA (UARRAY( 6,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00156,0.05334,0.17175,0.27761,0.41897,0.58962, &
     & 0.75360,0.91086,1.06612,1.23073,1.41861,1.63595,1.89877,2.18103, &
     & 2.49005,2.80082,3.10354,3.39397,3.66969,3.92943,   4351.9092    /
!     C  IV   D&F            
      DATA (UARRAY( 6,4,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69424,0.71734,0.77164,0.82314,0.88740,0.96544, &
     & 1.04087,1.11273,1.18286,1.25601,1.33981,1.43372,1.51957,1.59863, &
     & 1.67189,1.74015,1.80405,1.86411,1.92076,1.97438,   5862.9087    /
!     C  V    AEL            
      DATA (UARRAY( 6,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00001, &
     & 0.00008,0.00032,0.00100,0.00255,0.00562,0.01098,0.01950,0.03207, &
     & 0.04949,0.07239,0.10119,0.13602,0.17677,0.22308,  35643.6367    /
!     N  I    MMD            
      DATA (UARRAY( 7,1,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38637,1.38866,1.40008,1.42620,1.46712,1.51939, &
     & 1.57866,1.64136,1.70529,1.76968,1.83497,1.90246,1.97385,2.05085, &
     & 2.13479,2.22639,2.32563,2.43183,2.54379,2.66000,   1320.9091    /
!     N  II   MMD            
      DATA (UARRAY( 7,2,K),K=0,22) /                                    &
     & 0.00000,2.15018,2.18298,2.21800,2.25771,2.29761,2.33574,2.37146, &
     & 2.40461,2.43522,2.46343,2.48938,2.51327,2.53528,2.55557,2.57432, &
     & 2.59165,2.60772,2.62265,2.63653,2.64947,2.66155,   2690.2727    /
!     N  III  MMD            
      DATA (UARRAY( 7,3,K),K=0,22) /                                    &
     & 0.69315,1.75334,1.77260,1.78235,1.79929,1.82915,1.87138,1.92293, &
     & 1.98049,2.04142,2.10394,2.16703,2.23025,2.29353,2.35700,2.42087, &
     & 2.48536,2.55057,2.61657,2.68329,2.75060,2.81830,   4311.4546    /
!     N  IV   D&F            
      DATA (UARRAY( 7,4,K),K=0,22) /                                    &
     & 0.00000,0.00001,0.00929,0.12079,0.32824,0.49996,0.68273,0.87862, &
     & 1.05550,1.21704,1.37194,1.53266,1.71302,1.92803,2.17959,2.46121, &
     & 2.75504,2.99305,3.18516,3.34624,3.48494,3.60673,   7042.9092    /
!     N  V    D&F            
      DATA (UARRAY( 7,5,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69756,0.74263,0.82634,0.90358,0.98802,1.07704, &
     & 1.15731,1.23001,1.29754,1.36448,1.43686,1.52249,1.60780,1.68641, &
     & 1.75928,1.82720,1.89081,1.95060,2.00702,2.06043,   8899.0908    /
!     O  I    MMD            
      DATA (UARRAY( 8,1,K),K=0,22) /                                    &
     & 1.60944,2.11187,2.15328,2.16882,2.18052,2.19338,2.20815,2.22434, &
     & 2.24130,2.25854,2.27590,2.29357,2.31209,2.33234,2.35544,2.38263, &
     & 2.41512,2.45394,2.49981,2.55303,2.61351,2.68076,   1237.6364    /
!     O  II   AEL-NIST       
      DATA (UARRAY( 8,2,K),K=0,22) /                                    &
     & 1.38629,1.38631,1.39236,1.43305,1.51502,1.61887,1.72661,1.82883, &
     & 1.92225,2.00675,2.08374,2.15527,2.22369,2.29138,2.36049,2.43284, &
     & 2.50969,2.59173,2.67906,2.77127,2.86759,2.96701,   3191.6365    /
!     O  III  MMD            
      DATA (UARRAY( 8,3,K),K=0,22) /                                    &
     & 0.00000,2.13943,2.19797,2.25819,2.31628,2.37029,2.42093,2.46934, &
     & 2.51639,2.56267,2.60853,2.65419,2.69979,2.74542,2.79113,2.83689, &
     & 2.88269,2.92844,2.97405,3.01941,3.06442,3.10895,   4989.6362    /
!     O  IV   D&F            
      DATA (UARRAY( 8,4,K),K=0,22) /                                    &
     & 0.69315,1.73978,1.76698,1.79762,1.85260,1.90471,1.96849,2.04827, &
     & 2.12956,2.21004,2.28891,2.36690,2.44539,2.52410,2.59707,2.66506, &
     & 2.72873,2.78859,2.84506,2.89852,2.94926,2.99755,   7037.5454    /
!     O  V    D&F            
      DATA (UARRAY( 8,5,K),K=0,22) /                                    &
     & 0.00000,0.00010,0.02892,0.21847,0.49382,0.70949,0.91767,1.12174, &
     & 1.30074,1.46350,1.62595,1.80253,1.95914,2.09452,2.21373,2.32023, &
     & 2.41648,2.50426,2.58496,2.65964,2.72912,2.79408,  10354.5459    /
!     F  I    D&F            
      DATA (UARRAY( 9,1,K),K=0,22) /                                    &
     & 1.38629,1.68371,1.73424,1.75069,1.75868,1.76661,1.77179,1.77464, &
     & 1.77679,1.77876,1.78107,1.78487,1.79158,1.80308,1.82234,1.85086, &
     & 1.89273,1.94772,2.01648,2.09887,2.19358,2.29875,   1583.4546    /
!     F  II   D&F            
      DATA (UARRAY( 9,2,K),K=0,22) /                                    &
     & 1.60944,2.12483,2.16519,2.19747,2.23227,2.26589,2.29826,2.32781, &
     & 2.35452,2.37922,2.40417,2.43322,2.47290,2.53028,2.61417,2.72901, &
     & 2.87922,3.05549,3.25172,3.46062,3.67711,3.89170,   3180.0000    /
!     F  III  D&F            
      DATA (UARRAY( 9,3,K),K=0,22) /                                    &
     & 1.38629,1.38675,1.42156,1.54566,1.70076,1.83500,1.95813,2.06836, &
     & 2.16528,2.25738,2.36055,2.49971,2.70048,2.97709,3.32474,3.70217, &
     & 4.08846,4.47295,4.79757,5.04223,5.23865,5.40274,   5695.0908    /
!     F  IV   D&F            
      DATA (UARRAY( 9,4,K),K=0,22) /                                    &
     & 0.00000,2.12867,2.21779,2.29770,2.36910,2.43573,2.49831,2.55859, &
     & 2.61813,2.67931,2.74869,2.83767,2.96065,3.13922,3.31582,3.46586, &
     & 3.59629,3.71166,3.81509,3.90882,3.99450,4.07342,   7921.6362    /
!     F  V    D&F            
      DATA (UARRAY( 9,5,K),K=0,22) /                                    &
     & 0.69315,1.72407,1.76283,1.81708,1.90397,1.98391,2.07202,2.17167, &
     & 2.26821,2.36070,2.45124,2.54305,2.62917,2.70846,2.78192,2.85035, &
     & 2.91440,2.97459,3.03136,3.08508,3.13606,3.18457,  10385.4541    /
!     Ne I    D&F            
      DATA (UARRAY(10,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00100,0.00487,0.01716,0.05285,0.13154,0.26958,0.48024,0.74622, &
     & 1.05865,1.38012,1.69769,2.00358,2.29378,2.56629,   1959.9091    /
!     Ne II   D&F            
      DATA (UARRAY(10,2,K),K=0,22) /                                    &
     & 1.38629,1.70104,1.74401,1.75809,1.76442,1.77070,1.77530,1.77781, &
     & 1.78064,1.78751,1.80848,1.86224,1.97263,2.17032,2.44285,2.77723, &
     & 3.12845,3.49153,3.84158,4.16889,4.47280,4.75391,   3733.6365    /
!     Ne III  D&F            
      DATA (UARRAY(10,3,K),K=0,22) /                                    &
     & 1.60944,2.12316,2.18161,2.23547,2.28342,2.32917,2.36941,2.40112, &
     & 2.42905,2.45637,2.49018,2.54262,2.63314,2.78091,2.98449,3.25663, &
     & 3.55730,3.86965,4.11894,4.31832,4.48448,4.62692,   5772.7271    /
!     Ne IV   D&F            
      DATA (UARRAY(10,4,K),K=0,22) /                                    &
     & 1.38629,1.38950,1.47996,1.66772,1.85178,2.00716,2.13849,2.24738, &
     & 2.34144,2.42725,2.51361,2.60992,2.70195,2.78622,2.86394,2.93605, &
     & 3.00330,3.06632,3.12560,3.18157,3.23456,3.28489,   8828.1816    /
!     Ne V    AEL-Kurucz     
      DATA (UARRAY(10,5,K),K=0,22) /                                    &
     & 0.00000,2.11743,2.23774,2.33520,2.41925,2.49706,2.57268,2.64299, &
     & 2.71265,2.77780,2.83924,2.89713,2.95211,3.00422,3.05352,3.10050, &
     & 3.14537,3.18832,3.22969,3.26942,3.30745,3.34409,  11473.6367    /
!     Na I    AEL-NIST       
      DATA (UARRAY(11,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69323,0.69365,0.69491, &
     & 0.69780,0.70346,0.71368,0.73134,0.76054,0.80636,0.87383,0.96654, &
     & 1.08538,1.22811,1.39000,1.56508,1.74741,1.93193,    467.0909    /
!     Na II   AEL-NIST       
      DATA (UARRAY(11,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00001,0.00009, &
     & 0.00056,0.00245,0.00813,0.02197,0.05041,0.10111,0.18104,0.29404, &
     & 0.43920,0.61111,0.80180,1.00292,1.20735,1.40975,   4299.0908    /
!     Na III  AEL-NIST       
      DATA (UARRAY(11,3,K),K=0,22) /                                    &
     & 1.38629,1.70092,1.74394,1.75934,1.76724,1.77204,1.77530,1.77774, &
     & 1.77993,1.78251,1.78646,1.79310,1.80411,1.82131,1.84652,1.88119, &
     & 1.92629,1.98204,2.04801,2.12312,2.20589,2.29463,   6513.6362    /
!     Na IV   AEL-NIST       
      DATA (UARRAY(11,4,K),K=0,22) /                                    &
     & 1.60944,2.11943,2.20346,2.27639,2.33559,2.38257,2.42053,2.45227, &
     & 2.47993,2.50513,2.52922,2.55338,2.57868,2.60609,2.63643,2.67030, &
     & 2.70810,2.74996,2.79582,2.84540,2.89827,2.95393,   8991.8184    /
!     Na V    AEL-NIST       
      DATA (UARRAY(11,5,K),K=0,22) /                                    &
     & 1.38629,1.39652,1.55209,1.78418,1.98464,2.14373,2.27220,2.37958, &
     & 2.47218,2.55403,2.62783,2.69552,2.75857,2.81813,2.87510,2.93019, &
     & 2.98394,3.03674,3.08885,3.14042,3.19155,3.24224,  12581.8184    /
!     Mg I    AEL-NIST       
      DATA (UARRAY(12,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00011,0.00104,0.00473,0.01393, &
     & 0.03143,0.05950,0.10012,0.15534,0.22746,0.31869,0.43040,0.56239, &
     & 0.71258,0.87734,1.05217,1.23253,1.41439,1.59451,    694.9091    /
!     Mg II   AEL-NIST       
      DATA (UARRAY(12,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69316,0.69339,0.69477,0.69883,0.70712, &
     & 0.72071,0.74030,0.76644,0.79987,0.84154,0.89258,0.95396,1.02630, &
     & 1.10957,1.20307,1.30547,1.41501,1.52972,1.64764,   1366.4546    /
!     Mg III  AEL-NIST       
      DATA (UARRAY(12,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00002,0.00017, &
     & 0.00091,0.00354,0.01071,0.02677,0.05752,0.10924,0.18701,0.29306, &
     & 0.42588,0.58074,0.75116,0.93048,1.11294,1.29412,   7283.6362    /
!     Mg IV   AEL-NIST       
      DATA (UARRAY(12,4,K),K=0,22) /                                    &
     & 1.38629,1.69535,1.74083,1.75718,1.76560,1.77075,1.77435,1.77729, &
     & 1.78027,1.78421,1.79046,1.80088,1.81780,1.84374,1.88104,1.93142, &
     & 1.99566,2.07343,2.16339,2.26346,2.37118,2.48406,   9937.2725    /
!     Mg V    AEL-NIST       
      DATA (UARRAY(12,5,K),K=0,22) /                                    &
     & 1.60944,2.11627,2.22620,2.31055,2.37320,2.42090,2.45916,2.49163, &
     & 2.52055,2.54731,2.57288,2.59801,2.62334,2.64947,2.67694,2.70621, &
     & 2.73764,2.77148,2.80784,2.84671,2.88800,2.93149,  12842.7275    /
!     Al I    AEL-NIST       
      DATA (UARRAY(13,1,K),K=0,22) /                                    &
     & 0.69315,1.60422,1.69544,1.72699,1.74298,1.75263,1.75911,1.76381, &
     & 1.76751,1.77085,1.77447,1.77914,1.78582,1.79566,1.80989,1.82977, &
     & 1.85642,1.89072,1.93315,1.98378,2.04228,2.10793,    544.0000    /
!     Al II   AEL-NIST       
      DATA (UARRAY(13,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00024,0.00339,0.01638,0.04655,0.09732, &
     & 0.16782,0.25501,0.35596,0.46910,0.59415,0.73147,0.88108,1.04210, &
     & 1.21255,1.38966,1.57030,1.75148,1.93063,2.10573,   1711.1818    /
!     Al III  AEL-NIST       
      DATA (UARRAY(13,3,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69329,0.69482,0.70063,0.71341,0.73441, &
     & 0.76360,0.80032,0.84398,0.89431,0.95138,1.01550,1.08694,1.16576, &
     & 1.25165,1.34394,1.44162,1.54347,1.64818,1.75446,   2585.4546    /
!     Al IV   AEL-NIST       
      DATA (UARRAY(13,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00003,0.00026, &
     & 0.00135,0.00501,0.01452,0.03499,0.07281,0.13435,0.22400,0.34263, &
     & 0.48719,0.65179,0.82937,1.01328,1.19811,1.37994,  10908.1816    /
!     Al V    AEL-NIST       
      DATA (UARRAY(13,5,K),K=0,22) /                                    &
     & 1.38629,1.68701,1.73612,1.75392,1.76313,1.76890,1.77320,1.77706, &
     & 1.78130,1.78695,1.79546,1.80879,1.82932,1.85947,1.90134,1.95628, &
     & 2.02467,2.10580,2.19806,2.29928,2.40702,2.51891,  13977.2725    /
!     Si I    AEL-NIST       
      DATA (UARRAY(14,1,K),K=0,22) /                                    &
     & 0.00000,1.92047,2.05673,2.11216,2.15312,2.18961,2.22305,2.25360, &
     & 2.28153,2.30733,2.33193,2.35673,2.38364,2.41504,2.45352,2.50155, &
     & 2.56110,2.63326,2.71810,2.81466,2.92117,3.03540,    740.8182    /
!     Si II   AEL-NIST       
      DATA (UARRAY(14,2,K),K=0,22) /                                    &
     & 0.69315,1.61512,1.70120,1.73090,1.74600,1.75555,1.76333,1.77174, &
     & 1.78253,1.79698,1.81600,1.84025,1.87020,1.90622,1.94851,1.99709, &
     & 2.05183,2.11235,2.17810,2.24838,2.32238,2.39924,   1485.4546    /
!     Si III  AEL-NIST       
      DATA (UARRAY(14,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00003,0.00212,0.01713,0.05939,0.13397,0.23585, &
     & 0.35620,0.48751,0.62533,0.76804,0.91566,1.06873,1.22750,1.39146, &
     & 1.55940,1.72954,1.89992,2.06865,2.23411,2.39500,   3041.8181    /
!     Si IV   AEL-NIST       
      DATA (UARRAY(14,4,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69316,0.69384,0.69878,0.71283,0.73829,0.77466, &
     & 0.82007,0.87261,0.93089,0.99418,1.06229,1.13528,1.21327,1.29622, &
     & 1.38384,1.47559,1.57071,1.66829,1.76738,1.86704,   4103.7271    /
!     Si V    AEL-NIST       
      DATA (UARRAY(14,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00004,0.00036, &
     & 0.00180,0.00643,0.01806,0.04232,0.08590,0.15496,0.25311,0.38005, &
     & 0.53169,0.70144,0.88213,1.06730,1.25195,1.43254,  15160.9092    /
!     P  I    AEL-NIST       
      DATA (UARRAY(15,1,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38676,1.39457,1.42117,1.46875,1.53169,1.60276, &
     & 1.67636,1.74912,1.81970,1.88830,1.95620,2.02531,2.09774,2.17542, &
     & 2.25977,2.35150,2.45057,2.55625,2.66732,2.78230,    953.0909    /
!     P  II   AEL-NIST       
      DATA (UARRAY(15,2,K),K=0,22) /                                    &
     & 0.00000,1.95492,2.09083,2.16887,2.23090,2.28194,2.32478,2.36189, &
     & 2.39545,2.42736,2.45934,2.49297,2.52966,2.57062,2.61681,2.66883, &
     & 2.72693,2.79096,2.86046,2.93470,3.01279,3.09379,   1792.7273    /
!     P  III  AEL-NIST       
      DATA (UARRAY(15,3,K),K=0,22) /                                    &
     & 0.69315,1.60598,1.69637,1.72772,1.74469,1.75880,1.77576,1.79864, &
     & 1.82872,1.86625,1.91094,1.96242,2.02028,2.08413,2.15356,2.22803, &
     & 2.30693,2.38948,2.47488,2.56227,2.65082,2.73975,   2741.4546    /
!     P  IV   AEL-NIST       
      DATA (UARRAY(15,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00024,0.00811,0.04655,0.13024,0.25289,0.39911, &
     & 0.55581,0.71539,0.87476,1.03354,1.19249,1.35250,1.51408,1.67706, &
     & 1.84074,2.00400,2.16560,2.32432,2.47911,2.62914,   4674.5454    /
!     P  V    AEL-NIST       
      DATA (UARRAY(15,5,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69321,0.69531,0.70627,0.73169,0.77178,0.82361, &
     & 0.88373,0.94944,1.01914,1.09218,1.16852,1.24849,1.33241,1.42048, &
     & 1.51258,1.60833,1.70707,1.80800,1.91021,2.01282,   5910.9087    /
!     S  I    AEL-NIST       
      DATA (UARRAY(16,1,K),K=0,22) /                                    &
     & 1.60944,1.95360,2.06256,2.10945,2.14344,2.17421,2.20334,2.23079, &
     & 2.25643,2.28031,2.30270,2.32408,2.34518,2.36687,2.39011,2.41591, &
     & 2.44518,2.47872,2.51706,2.56050,2.60907,2.66255,    941.5455    /
!     S  II   AEL-NIST       
      DATA (UARRAY(16,2,K),K=0,22) /                                    &
     & 1.38629,1.38640,1.40288,1.47555,1.58989,1.71457,1.83242,1.93794, &
     & 2.03125,2.11468,2.19133,2.26442,2.33693,2.41141,2.48978,2.57320, &
     & 2.66212,2.75633,2.85509,2.95736,3.06195,3.16765,   2127.2727    /
!     S  III  AEL-NIST       
      DATA (UARRAY(16,3,K),K=0,22) /                                    &
     & 0.00000,1.95747,2.12038,2.21949,2.29207,2.34879,2.39652,2.43993, &
     & 2.48220,2.52540,2.57079,2.61897,2.67012,2.72403,2.78032,2.83844, &
     & 2.89783,2.95791,3.01815,3.07811,3.13738,3.19565,   3181.8181    /
!     S  IV   AEL-NIST       
      DATA (UARRAY(16,4,K),K=0,22) /                                    &
     & 0.69315,1.59114,1.68852,1.72305,1.74531,1.77020,1.80401,1.84846, &
     & 1.90300,1.96632,2.03709,2.11416,2.19656,2.28343,2.37389,2.46709, &
     & 2.56214,2.65818,2.75440,2.85006,2.94455,3.03734,   4300.0000    /
!     S  V    AEL-NIST       
      DATA (UARRAY(16,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00099,0.02067,0.09295,0.22254,0.38830,0.56848, &
     & 0.75018,0.92785,1.10003,1.26708,1.42990,1.58935,1.74593,1.89979, &
     & 2.05077,2.19849,2.34251,2.48238,2.61771,2.74821,   6607.2729    /
!     Cl I    D&F            
      DATA (UARRAY(17,1,K),K=0,22) /                                    &
     & 1.38629,1.54406,1.64273,1.68284,1.70239,1.72157,1.73571,1.74335, &
     & 1.74967,1.75639,1.76733,1.78903,1.83045,1.90724,2.02455,2.19046, &
     & 2.39168,2.62630,2.87774,3.13159,3.38094,3.62122,   1182.7273    /
!     Cl II   D&F            
      DATA (UARRAY(17,2,K),K=0,22) /                                    &
     & 1.60944,2.00216,2.10431,2.16708,2.21974,2.26977,2.31177,2.34600, &
     & 2.37606,2.40433,2.43442,2.47250,2.52714,2.60705,2.71842,2.86824, &
     & 3.05083,3.25598,3.47594,3.70771,3.93739,4.16056,   2163.6365    /
!     Cl III  D&F            
      DATA (UARRAY(17,3,K),K=0,22) /                                    &
     & 1.38629,1.38821,1.45686,1.62239,1.79364,1.93980,2.06531,2.17073, &
     & 2.26339,2.35355,2.45778,2.60060,2.80664,3.07946,3.42081,3.78331, &
     & 4.16947,4.53782,4.88440,5.20799,5.50875,5.78753,   3627.2727    /
!     Cl IV   D&F            
      DATA (UARRAY(17,4,K),K=0,22) /                                    &
     & 0.00000,1.95171,2.14583,2.25173,2.32124,2.38622,2.43982,2.48449, &
     & 2.53036,2.58252,2.64876,2.74101,2.87277,3.05719,3.29699,3.53893, &
     & 3.73359,3.89646,4.03648,4.15928,4.26864,4.36721,   4860.0000    /
!     Cl V    D&F            
      DATA (UARRAY(17,5,K),K=0,22) /                                    &
     & 0.69315,1.57321,1.67902,1.72068,1.75994,1.79772,1.84604,1.91286, &
     & 1.98848,2.07035,2.15788,2.25312,2.35809,2.45673,2.54651,2.62889, &
     & 2.70499,2.77571,2.84176,2.90372,2.96206,3.01718,   6154.5454    /
!     Ar I    D&F            
      DATA (UARRAY(18,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00100,0.00548,0.02283,0.06789,0.16277,0.33526,0.57781,0.88481, &
     & 1.21404,1.55770,1.89516,2.21313,2.51007,2.78574,   1432.2727    /
!     Ar II   D&F            
      DATA (UARRAY(18,2,K),K=0,22) /                                    &
     & 1.38629,1.58522,1.67277,1.70514,1.72144,1.73747,1.74828,1.75455, &
     & 1.76088,1.77156,1.79600,1.85260,1.96533,2.14790,2.41191,2.71941, &
     & 3.06880,3.41501,3.74840,4.06411,4.36013,4.63601,   2510.9092    /
!     Ar III  D&F            
      DATA (UARRAY(18,3,K),K=0,22) /                                    &
     & 1.60944,2.01725,2.13904,2.22001,2.28087,2.33824,2.38392,2.41903, &
     & 2.45127,2.48541,2.52821,2.59022,2.68472,2.82074,3.01056,3.23690, &
     & 3.50519,3.78461,4.06529,4.34014,4.60464,4.85615,   3718.1819    /
!     Ar IV   D&F            
      DATA (UARRAY(18,4,K),K=0,22) /                                    &
     & 1.38629,1.39563,1.54038,1.76148,1.94811,2.10533,2.23236,2.33371, &
     & 2.42327,2.50958,2.60145,2.70945,2.84314,3.00811,3.21083,3.43656, &
     & 3.68333,3.92328,4.11665,4.27863,4.41798,4.54027,   5437.2729    /
!     Ar V    D&F            
      DATA (UARRAY(18,5,K),K=0,22) /                                    &
     & 0.00000,1.94177,2.16534,2.27318,2.34097,2.40446,2.45975,2.51145, &
     & 2.56944,2.63847,2.72816,2.85087,3.01893,3.25047,3.52798,3.84681, &
     & 4.12356,4.34008,4.51797,4.66894,4.80008,4.91600,   6821.8184    /
!     K  I    AEL-NIST       
      DATA (UARRAY(19,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69317,0.69337,0.69425,0.69663, &
     & 0.70153,0.71008,0.72367,0.74410,0.77365,0.81491,0.87042,0.94205, &
     & 1.03055,1.13526,1.25420,1.38447,1.52277,1.66583,    394.4546    /
!     K  II   AEL-NIST       
      DATA (UARRAY(19,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00004,0.00031, &
     & 0.00142,0.00470,0.01229,0.02700,0.05187,0.08962,0.14206,0.20968, &
     & 0.29160,0.38584,0.48970,0.60027,0.71479,0.83088,   2891.8181    /
!     K  III  AEL-NIST       
      DATA (UARRAY(19,3,K),K=0,22) /                                    &
     & 1.38629,1.59923,1.68228,1.71561,1.73343,1.74455,1.75234,1.75860, &
     & 1.76477,1.77227,1.78255,1.79704,1.81694,1.84311,1.87598,1.91554, &
     & 1.96139,2.01279,2.06886,2.12857,2.19095,2.25506,   4181.8184    /
!     K  IV   AEL-NIST       
      DATA (UARRAY(19,4,K),K=0,22) /                                    &
     & 1.60944,2.02305,2.16975,2.26791,2.33794,2.39019,2.43156,2.46675, &
     & 2.49907,2.53082,2.56357,2.59830,2.63547,2.67519,2.71730,2.76144, &
     & 2.80718,2.85403,2.90155,2.94929,2.99688,3.04401,   5538.1816    /
!     K  V    AEL-NIST       
      DATA (UARRAY(19,5,K),K=0,22) /                                    &
     & 1.38629,1.41121,1.63333,1.88864,2.08608,2.23561,2.35528,2.45759, &
     & 2.54998,2.63662,2.71966,2.80003,2.87804,2.95365,3.02673,3.09710, &
     & 3.16463,3.22926,3.29094,3.34971,3.40561,3.45875,   7514.5454    /
!     Ca I    AEL-NIST       
      DATA (UARRAY(20,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00002,0.00049,0.00378,0.01509,0.04113, &
     & 0.08771,0.15831,0.25406,0.37465,0.51896,0.68513,0.87021,1.07012, &
     & 1.28004,1.49501,1.71056,1.92304,2.12976,2.32892,    555.5455    /
!     Ca II   AEL-NIST       
      DATA (UARRAY(20,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69369,0.70452,0.74460,0.81836,0.91586,1.02449, &
     & 1.13492,1.24188,1.34317,1.43841,1.52826,1.61384,1.69637,1.77702, &
     & 1.85671,1.93614,2.01574,2.09572,2.17608,2.25668,   1078.9091    /
!     Ca III  AEL-NIST       
      DATA (UARRAY(20,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00010,0.00089,0.00452, &
     & 0.01549,0.04094,0.08973,0.17026,0.28755,0.44113,0.62504,0.83005, &
     & 1.04649,1.26613,1.48293,1.69293,1.89382,2.08444,   4655.4546    /
!     Ca IV   AEL-NIST       
      DATA (UARRAY(20,4,K),K=0,22) /                                    &
     & 1.38629,1.60122,1.68361,1.71658,1.73428,1.74604,1.75704,1.77238, &
     & 1.79735,1.83646,1.89246,1.96574,2.05457,2.15575,2.26549,2.38010, &
     & 2.49644,2.61207,2.72524,2.83475,2.93989,3.04026,   6104.5454    /
!     Ca V    AEL-NIST       
      DATA (UARRAY(20,5,K),K=0,22) /                                    &
     & 1.60944,2.02806,2.19752,2.30290,2.37410,2.42624,2.46790,2.50419, &
     & 2.53830,2.57226,2.60733,2.64422,2.68322,2.72430,2.76724,2.81172, &
     & 2.85732,2.90365,2.95032,2.99698,3.04332,3.08909,   7675.4546    /
!     Sc I    AEL-NIST       
      DATA (UARRAY(21,1,K),K=0,22) /                                    &
     & 1.38629,2.07848,2.18541,2.22363,2.24604,2.26974,2.30587,2.36111, &
     & 2.43724,2.53212,2.64156,2.76096,2.88625,3.01425,3.14259,3.26960, &
     & 3.39414,3.51540,3.63286,3.74620,3.85522,3.95985,    594.5455    /
!     Sc II   AEL-NIST       
      DATA (UARRAY(21,2,K),K=0,22) /                                    &
     & 1.09861,2.60135,2.78063,2.94704,3.08809,3.20599,3.30642,3.39416, &
     & 3.47258,3.54406,3.61033,3.67275,3.73237,3.79008,3.84655,3.90229, &
     & 3.95766,4.01288,4.06803,4.12313,4.17813,4.23291,   1163.6364    /
!     Sc III  AEL-NIST       
      DATA (UARRAY(21,3,K),K=0,22) /                                    &
     & 1.38629,2.22868,2.26521,2.27841,2.28721,2.29542,2.30401,2.31325, &
     & 2.32325,2.33410,2.34597,2.35908,2.37373,2.39024,2.40895,2.43015, &
     & 2.45409,2.48090,2.51062,2.54321,2.57851,2.61630,   2250.0000    /
!     Sc IV   AEL-NIST       
      DATA (UARRAY(21,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00006,0.00098,0.00607,0.02255, &
     & 0.06039,0.12928,0.23506,0.37770,0.55161,0.74804,0.95768,1.17248, &
     & 1.38635,1.59514,1.79624,1.98820,2.17033,2.34249,   6700.0000    /
!     Sc V    AEL-NIST       
      DATA (UARRAY(21,5,K),K=0,22) /                                    &
     & 1.38629,1.59898,1.68212,1.71550,1.73356,1.74562,1.75590,1.76708, &
     & 1.78109,1.79933,1.82260,1.85125,1.88513,1.92382,1.96663,2.01280, &
     & 2.06149,2.11195,2.16347,2.21545,2.26736,2.31881,   8336.3633    /
!     Ti I    AEL-NIST       
      DATA (UARRAY(22,1,K),K=0,22) /                                    &
     & 1.60944,2.59404,2.80378,2.89244,2.96679,3.05014,3.14664,3.25530, &
     & 3.37349,3.49820,3.62657,3.75621,3.88523,4.01225,4.13627,4.25662, &
     & 4.37285,4.48468,4.59198,4.69471,4.79289,4.88663,    620.0000    /
!     Ti II   AEL-NIST       
      DATA (UARRAY(22,2,K),K=0,22) /                                    &
     & 1.38629,3.39664,3.70584,3.87313,4.00811,4.12711,4.23375,4.33021, &
     & 4.41843,4.50003,4.57632,4.64829,4.71670,4.78209,4.84489,4.90539, &
     & 4.96378,5.02023,5.07482,5.12764,5.17874,5.22816,   1233.6364    /
!     Ti III  AEL-NIST       
      DATA (UARRAY(22,3,K),K=0,22) /                                    &
     & 1.60944,2.91323,3.02556,3.12594,3.21545,3.29255,3.35954,3.41906, &
     & 3.47327,3.52385,3.57218,3.61941,3.66650,3.71427,3.76336,3.81419, &
     & 3.86704,3.92196,3.97886,4.03755,4.09770,4.15897,   2497.2727    /
!     Ti IV   AEL-NIST       
      DATA (UARRAY(22,4,K),K=0,22) /                                    &
     & 1.38629,2.22107,2.26124,2.27490,2.28190,2.28652,2.29043,2.29446, &
     & 2.29915,2.30485,2.31190,2.32065,2.33145,2.34469,2.36073,2.37988, &
     & 2.40237,2.42833,2.45775,2.49054,2.52650,2.56534,   3932.7273    /
!     Ti V    AEL-NIST       
      DATA (UARRAY(22,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00001,0.00049,0.00495,0.02326,0.06964, &
     & 0.15573,0.28424,0.44817,0.63520,0.83291,1.03163,1.22508,1.40967, &
     & 1.58369,1.74662,1.89857,2.04005,2.17172,2.29431,   9036.3633    /
!     V  I    AEL-NIST       
      DATA (UARRAY(23,1,K),K=0,22) /                                    &
     & 1.38629,2.70531,3.08391,3.29986,3.44706,3.56010,3.65805,3.75158, &
     & 3.84618,3.94407,4.04558,4.15007,4.25649,4.36372,4.47072,4.57657, &
     & 4.68056,4.78211,4.88079,4.97631,5.06847,5.15718,    612.7273    /
!     V  II   AEL-NIST       
      DATA (UARRAY(23,2,K),K=0,22) /                                    &
     & 0.00000,3.07656,3.39846,3.62149,3.81455,3.99510,4.16354,4.31889, &
     & 4.46148,4.59251,4.71344,4.82562,4.93019,5.02811,5.12012,5.20682, &
     & 5.28869,5.36614,5.43949,5.50906,5.57508,5.63780,   1331.8182    /
!     V  III  AEL-NIST       
      DATA (UARRAY(23,3,K),K=0,22) /                                    &
     & 1.38629,3.15909,3.30531,3.46374,3.61942,3.75767,3.87831,3.98481, &
     & 4.08071,4.16887,4.25146,4.33013,4.40607,4.48013,4.55289,4.62466, &
     & 4.69560,4.76575,4.83504,4.90335,4.97055,5.03648,   2664.5454    /
!     V  IV   AEL-NIST       
      DATA (UARRAY(23,4,K),K=0,22) /                                    &
     & 1.60944,2.91903,3.07685,3.20316,3.29865,3.37095,3.42774,3.47461, &
     & 3.51535,3.55256,3.58807,3.62313,3.65865,3.69520,3.73313,3.77260, &
     & 3.81360,3.85602,3.89966,3.94428,3.98960,4.03537,   4246.3638    /
!     V  V    AEL-NIST       
      DATA (UARRAY(23,5,K),K=0,22) /                                    &
     & 1.38629,2.21440,2.25779,2.27257,2.28004,2.28468,2.28816,2.29139, &
     & 2.29503,2.29968,2.30595,2.31454,2.32617,2.34157,2.36139,2.38616, &
     & 2.41621,2.45164,2.49235,2.53800,2.58810,2.64206,   5930.0005    /
!     Cr I    AEL-NIST       
      DATA (UARRAY(24,1,K),K=0,22) /                                    &
     & 1.94591,1.94591,1.94629,1.95435,1.98524,2.04361,2.12375,2.21907, &
     & 2.32599,2.44331,2.57049,2.70663,2.85012,2.99878,3.15023,3.30222, &
     & 3.45276,3.60030,3.74367,3.88208,4.01503,4.14226,    614.9091    /
!     Cr II   AEL-NIST       
      DATA (UARRAY(24,2,K),K=0,22) /                                    &
     & 1.79176,1.79180,1.80617,1.90292,2.11116,2.39572,2.70837,3.01603, &
     & 3.30323,3.56547,3.80340,4.01955,4.21684,4.39794,4.56513,4.72026, &
     & 4.86482,5.00000,5.12678,5.24595,5.35819,5.46406,   1499.0909    /
!     Cr III  AEL-NIST       
      DATA (UARRAY(24,3,K),K=0,22) /                                    &
     & 0.00000,3.04470,3.16708,3.32901,3.52806,3.72408,3.90249,4.06220, &
     & 4.20597,4.33687,4.45740,4.56937,4.67402,4.77226,4.86471,4.95186, &
     & 5.03408,5.11170,5.18502,5.25430,5.31980,5.38175,   2813.6365    /
!     Cr IV   AEL-NIST       
      DATA (UARRAY(24,4,K),K=0,22) /                                    &
     & 1.38629,3.17459,3.39486,3.60979,3.78436,3.92027,4.02808,4.11688, &
     & 4.19320,4.26153,4.32487,4.38517,4.44364,4.50099,4.55757,4.61350, &
     & 4.66876,4.72325,4.77684,4.82939,4.88078,4.93089,   4463.6362    /
!     Cr V    AEL-NIST       
      DATA (UARRAY(24,5,K),K=0,22) /                                    &
     & 1.60944,2.93229,3.13182,3.27065,3.36592,3.43389,3.48500,3.52571, &
     & 3.56013,3.59098,3.62005,3.64850,3.67703,3.70605,3.73573,3.76609, &
     & 3.79707,3.82853,3.86032,3.89226,3.92420,3.95598,   6381.8184    /
!     Mn I    AEL-NIST       
      DATA (UARRAY(25,1,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79176,1.79179,1.79249,1.79677,1.81023,1.83929, &
     & 1.88905,1.96192,2.05742,2.17277,2.30388,2.44627,2.59570,2.74854, &
     & 2.90187,3.05347,3.20173,3.34554,3.48419,3.61728,    675.6364    /
!     Mn II   AEL-NIST       
      DATA (UARRAY(25,2,K),K=0,22) /                                    &
     & 1.94591,1.94596,1.95403,2.00082,2.10159,2.25352,2.44782,2.67019, &
     & 2.90471,3.13873,3.36438,3.57776,3.77753,3.96368,4.13691,4.29815, &
     & 4.44839,4.58860,4.71966,4.84236,4.95742,5.06549,   1421.4546    /
!     Mn III  AEL-NIST       
      DATA (UARRAY(25,3,K),K=0,22) /                                    &
     & 1.79176,1.79178,1.80729,1.93731,2.22033,2.56695,2.90340,3.20536, &
     & 3.47175,3.70774,3.91896,4.10998,4.28427,4.44440,4.59227,4.72936, &
     & 4.85686,4.97571,5.08673,5.19061,5.28796,5.37933,   3062.7273    /
!     Mn IV   AEL-NIST       
      DATA (UARRAY(25,4,K),K=0,22) /                                    &
     & 0.00000,3.06005,3.25044,3.49262,3.71370,3.89215,4.03474,4.15161, &
     & 4.25091,4.33841,4.41791,4.49187,4.56177,4.62849,4.69251,4.75408, &
     & 4.81332,4.87029,4.92503,4.97756,5.02791,5.07612,   4672.7271    /
!     Mn V    AEL-NIST       
      DATA (UARRAY(25,5,K),K=0,22) /                                    &
     & 1.38629,3.19618,3.49418,3.73916,3.91635,4.04563,4.14378,4.22192, &
     & 4.28725,4.34448,4.39663,4.44559,4.49250,4.53799,4.58239,4.62583, &
     & 4.66833,4.70985,4.75033,4.78972,4.82795,4.86500,   6636.3638    /
!     Fe I    AEL-NIST       
      DATA (UARRAY(26,1,K),K=0,22) /                                    &
     & 2.19722,2.63110,2.87325,2.98770,3.07326,3.15442,3.23776,3.32575, &
     & 3.41952,3.51929,3.62455,3.73434,3.84747,3.96269,4.07880,4.19472, &
     & 4.30950,4.42238,4.53272,4.64004,4.74400,4.84436,    715.4545    /
!     Fe II   AEL-NIST       
      DATA (UARRAY(26,2,K),K=0,22) /                                    &
     & 2.30259,3.17871,3.52634,3.71233,3.85170,3.97884,4.10379,4.22805, &
     & 4.35076,4.47094,4.58798,4.70161,4.81176,4.91845,5.02169,5.12151, &
     & 5.21791,5.31088,5.40045,5.48661,5.56940,5.64887,   1470.9091    /
!     Fe III  AEL-NIST       
      DATA (UARRAY(26,3,K),K=0,22) /                                    &
     & 2.19722,3.01774,3.13025,3.24754,3.40829,3.58432,3.75623,3.91781, &
     & 4.06886,4.21098,4.34585,4.47479,4.59864,4.71791,4.83285,4.94357, &
     & 5.05012,5.15250,5.25074,5.34489,5.43500,5.52116,   2785.7273    /
!     Fe IV   AEL-NIST       
      DATA (UARRAY(26,4,K),K=0,22) /                                    &
     & 1.79176,1.79242,1.88074,2.23469,2.68780,3.08680,3.40916,3.67024, &
     & 3.88757,4.07433,4.23941,4.38858,4.52551,4.65249,4.77098,4.88196, &
     & 4.98611,5.08398,5.17599,5.26255,5.34400,5.42070,   4981.8184    /
!     Fe V    AEL-NIST       
      DATA (UARRAY(26,5,K),K=0,22) /                                    &
     & 0.00000,3.07426,3.37724,3.71335,3.98060,4.18136,4.33510,4.45746, &
     & 4.55925,4.64765,4.72734,4.80127,4.87122,4.93823,5.00283,5.06527, &
     & 5.12565,5.18398,5.24024,5.29441,5.34649,5.39648,   6863.6362    /
!     Co I    AEL-NIST       
      DATA (UARRAY(27,1,K),K=0,22) /                                    &
     & 2.30259,2.48626,2.78050,2.99951,3.16837,3.30332,3.41586,3.51430, &
     & 3.60440,3.69001,3.77362,3.85670,3.94004,4.02398,4.10849,4.19337, &
     & 4.27829,4.36286,4.44668,4.52937,4.61059,4.69008,    714.5455    /
!     Co II   AEL-NIST       
      DATA (UARRAY(27,2,K),K=0,22) /                                    &
     & 2.19722,2.63286,3.06455,3.34623,3.55147,3.71763,3.86114,3.98982, &
     & 4.10792,4.21806,4.32193,4.42071,4.51515,4.60577,4.69289,4.77671, &
     & 4.85738,4.93497,5.00955,5.08120,5.14998,5.21596,   1550.0000    /
!     Co III  AEL-NIST       
      DATA (UARRAY(27,3,K),K=0,22) /                                    &
     & 2.30259,2.99910,3.18605,3.33707,3.48836,3.63159,3.76519,3.89124, &
     & 4.01219,4.12965,4.24439,4.35651,4.46574,4.57165,4.67382,4.77191, &
     & 4.86569,4.95506,5.03999,5.12056,5.19688,5.26911,   3044.5454    /
!     Co IV   AEL-NIST       
      DATA (UARRAY(27,4,K),K=0,22) /                                    &
     & 2.19722,3.04334,3.20676,3.43323,3.66818,3.87275,4.04391,4.18885, &
     & 4.31549,4.43007,4.53695,4.63887,4.73733,4.83305,4.92620,5.01671, &
     & 5.10440,5.18905,5.27052,5.34869,5.42350,5.49497,   4663.6362    /
!     Co V    AEL-NIST       
      DATA (UARRAY(27,5,K),K=0,22) /                                    &
     & 1.79176,1.79653,2.03359,2.58431,3.09255,3.47768,3.76735,3.99275, &
     & 4.17579,4.33084,4.46712,4.59041,4.70429,4.81092,4.91154,5.00691, &
     & 5.09747,5.18350,5.26521,5.34278,5.41639,5.48622,   7227.2729    /
!     Ni I    AEL-NIST       
      DATA (UARRAY(28,1,K),K=0,22) /                                    &
     & 2.19722,2.70588,2.97473,3.13804,3.24468,3.31966,3.37615,3.42175, &
     & 3.46126,3.49795,3.53413,3.57146,3.61101,3.65346,3.69912,3.74802, &
     & 3.79996,3.85461,3.91150,3.97017,4.03009,4.09078,    693.9091    /
!     Ni II   AEL-NIST       
      DATA (UARRAY(28,2,K),K=0,22) /                                    &
     & 1.79176,1.95801,2.15293,2.37365,2.59292,2.78748,2.95576,3.10283, &
     & 3.23441,3.35536,3.46957,3.58003,3.68897,3.79797,3.90800,4.01945, &
     & 4.13228,4.24608,4.36026,4.47410,4.58688,4.69795,   1650.0000    /
!     Ni III  AEL-NIST       
      DATA (UARRAY(28,3,K),K=0,22) /                                    &
     & 2.19722,2.68156,2.86902,2.99008,3.09326,3.18837,3.28252,3.38045, &
     & 3.48442,3.59474,3.71048,3.83004,3.95161,4.07349,4.19421,4.31259, &
     & 4.42778,4.53914,4.64631,4.74906,4.84731,4.94106,   3196.3635    /
!     Ni IV   AEL-NIST       
      DATA (UARRAY(28,4,K),K=0,22) /                                    &
     & 2.30259,3.04534,3.28133,3.49469,3.67457,3.81955,3.93911,4.04307, &
     & 4.13899,4.23195,4.32478,4.41865,4.51357,4.60889,4.70372,4.79714, &
     & 4.88835,4.97673,5.06184,5.14342,5.22132,5.29552,   4990.9092    /
!     Ni V    AEL-NIST       
      DATA (UARRAY(28,5,K),K=0,22) /                                    &
     & 2.19722,3.06109,3.31713,3.62714,3.89056,4.09528,4.25521,4.38457, &
     & 4.49423,4.59186,4.68258,4.76961,4.85472,4.93872,5.02177,5.10371, &
     & 5.18420,5.26288,5.33941,5.41352,5.48501,5.55375,   6863.6362    /
!     Cu I    AEL-NIST       
      DATA (UARRAY(29,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69318,0.69481,0.70499,0.73162,0.77676,0.83721, &
     & 0.90789,0.98444,1.06408,1.14560,1.22888,1.31442,1.40291,1.49494, &
     & 1.59080,1.69042,1.79339,1.89902,2.00647,2.11483,    702.1818    /
!     Cu II   AEL-NIST       
      DATA (UARRAY(29,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00246,0.04668,0.19317,0.42253,0.67653,0.91994, &
     & 1.14252,1.34561,1.53413,1.71328,1.88730,2.05907,2.23009,2.40067, &
     & 2.57031,2.73804,2.90272,3.06331,3.21892,3.36891,   1844.5454    /
!     Cu III  AEL-NIST       
      DATA (UARRAY(29,3,K),K=0,22) /                                    &
     & 1.79176,2.03368,2.14745,2.19477,2.22501,2.25801,2.30533,2.37250, &
     & 2.46038,2.56698,2.68897,2.82271,2.96476,3.11203,3.26188,3.41206, &
     & 3.56077,3.70657,3.84838,3.98546,4.11731,4.24365,   3348.1819    /
!     Cu IV   AEL-NIST       
      DATA (UARRAY(29,4,K),K=0,22) /                                    &
     & 2.19722,2.72492,2.93238,3.07622,3.18582,3.27171,3.34404,3.41154, &
     & 3.48088,3.55655,3.64085,3.73423,3.83569,3.94338,4.05510,4.16869, &
     & 4.28225,4.39427,4.50359,4.60944,4.71129,4.80885,   5018.1816    /
!     Cu V    AEL-NIST       
      DATA (UARRAY(29,5,K),K=0,22) /                                    &
     & 2.30259,3.07869,3.38238,3.63339,3.82084,3.96064,4.06953,4.16010, &
     & 4.24134,4.31934,4.39779,4.47849,4.56186,4.64743,4.73429,4.82137, &
     & 4.90767,4.99234,5.07473,5.15436,5.23093,5.30427,   7263.6362    /
!     Zn I    D&F            
      DATA (UARRAY(30,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00100,0.00380, &
     & 0.01014,0.02254,0.04506,0.08385,0.14393,0.23449,0.35595,0.51725, &
     & 0.70416,0.90986,1.13435,1.36126,1.58526,1.80271,    853.7273    /
!     Zn II   D&F            
      DATA (UARRAY(30,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69330,0.69368,0.69407,0.69605,0.70132, &
     & 0.71188,0.73044,0.76014,0.80601,0.87357,0.97007,1.09914,1.25953, &
     & 1.44352,1.65373,1.87280,2.09448,2.31396,2.52759,   1632.7273    /
!     Zn III  D&F            
      DATA (UARRAY(30,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00606,0.02609,0.04573,0.09655,0.20454, &
     & 0.35390,0.55020,0.79602,1.10656,1.48111,1.89274,2.34436,2.70724, &
     & 2.97292,3.18263,3.35590,3.50353,3.63214,3.74608,   3609.0908    /
!     Zn IV   G   state
      DATA (UARRAY(30,4,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,   5400.0000    /
!     Zn V    U=1 assumed
      DATA (UARRAY(30,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   7509.0908    /
!     Ga I    AEL            
      DATA (UARRAY(31,1,K),K=0,22) /                                    &
     & 0.69315,0.89709,1.20755,1.36978,1.46319,1.52318,1.56481,1.59536, &
     & 1.61880,1.63754,1.65320,1.66705,1.68019,1.69359,1.70821,1.72490, &
     & 1.74440,1.76734,1.79414,1.82508,1.86021,1.89943,    545.4545    /
!     Ga II   AEL            
      DATA (UARRAY(31,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00004,0.00082,0.00528,0.01832,0.04444, &
     & 0.08607,0.14354,0.21588,0.30170,0.39969,0.50863,0.62722,0.75399, &
     & 0.88715,1.02480,1.16496,1.30580,1.44571,1.58336,   1864.5454    /
!     Ga III  AEL            
      DATA (UARRAY(31,3,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69318,0.69373,0.69638,0.70327,0.71613, &
     & 0.73596,0.76323,0.79821,0.84112,0.89214,0.95132,1.01845,1.09302, &
     & 1.17419,1.26087,1.35182,1.44576,1.54145,1.63778,   2790.9092    /
!     Ga IV   AEL            
      DATA (UARRAY(31,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00007,0.00172,0.01154,0.04175,0.10482, &
     & 0.20698,0.34516,0.50931,0.68733,0.86883,1.04657,1.21622,1.37561, &
     & 1.52396,1.66133,1.78822,1.90531,2.01342,2.11331,   5818.1816    /
!     Ga V    G   state
      DATA (UARRAY(31,5,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,   7909.0908    /
!     Ge I    AEL-NIST       
      DATA (UARRAY(32,1,K),K=0,22) /                                    &
     & 0.00000,0.82167,1.36911,1.62014,1.77003,1.87539,1.95628,2.02152, &
     & 2.07587,2.12233,2.16319,2.20044,2.23607,2.27206,2.31038,2.35288, &
     & 2.40112,2.45621,2.51877,2.58883,2.66595,2.74925,    716.3636    /
!     Ge II   AEL-NIST       
      DATA (UARRAY(32,2,K),K=0,22) /                                    &
     & 1.38629,1.68305,1.99129,2.13481,2.21458,2.26505,2.30006,2.32647, &
     & 2.34834,2.36842,2.38877,2.41104,2.43659,2.46655,2.50180,2.54295, &
     & 2.59030,2.64385,2.70327,2.76802,2.83734,2.91037,   1448.1818    /
!     Ge III  AEL-NIST       
      DATA (UARRAY(32,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00052,0.00602,0.02606,0.06869,0.13586, &
     & 0.22425,0.32849,0.44347,0.56529,0.69123,0.81943,0.94850,1.07736, &
     & 1.20506,1.33080,1.45391,1.57381,1.69008,1.80242,   3110.0000    /
!     Ge IV   AEL-NIST       
      DATA (UARRAY(32,4,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38630,1.38650,1.38855,1.39576,1.41086,1.43478, &
     & 1.46704,1.50668,1.55283,1.60494,1.66274,1.72606,1.79468,1.86822, &
     & 1.94610,2.02760,2.11185,2.19796,2.28505,2.37231,   4155.4546    /
!     Ge V    AEL-NIST       
      DATA (UARRAY(32,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00003,0.00093,0.00730,0.02975,0.08238, &
     & 0.17694,0.31736,0.49834,0.70829,0.93419,1.16496,1.39266,1.61231, &
     & 1.82115,2.01793,2.20235,2.37465,2.53540,2.68531,   8500.0000    /
!     As I    MMD            
      DATA (UARRAY(33,1,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38671,1.39393,1.41895,1.46413,1.52437,1.59288, &
     & 1.66428,1.73523,1.80426,1.87131,1.93732,2.00384,2.07263,2.14533, &
     & 2.22325,2.30713,2.39712,2.49282,2.59339,2.69771,    891.8182    /
!     As II   MMD            
      DATA (UARRAY(33,2,K),K=0,22) /                                    &
     & 0.00000,1.02733,1.54267,1.78217,1.93181,2.03761,2.11751,2.18075, &
     & 2.23300,2.27827,2.31966,2.35980,2.40096,2.44505,2.49364,2.54780, &
     & 2.60813,2.67472,2.74722,2.82491,2.90687,2.99204,   1693.6364    /
!     As III  AEL-Kurucz     
      DATA (UARRAY(33,3,K),K=0,22) /                                    &
     & 0.69315,1.02045,1.32443,1.57872,1.78430,1.95468,2.09325,2.21478, &
     & 2.30962,2.39610,2.46610,2.53148,2.59026,2.64584,2.70179,2.75497, &
     & 2.81403,2.87011,2.93496,2.99626,3.06725,3.13353,   2576.3635    /
!     As IV   AEL            
      DATA (UARRAY(33,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00004,0.00244,0.01910,0.06519,0.14557,0.25427, &
     & 0.38120,0.51719,0.65566,0.79239,0.92486,1.05167,1.17211,1.28593, &
     & 1.39314,1.49390,1.58851,1.67728,1.76058,1.83877,   4557.2729    /
!     As V    AEL            
      DATA (UARRAY(33,5,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69316,0.69382,0.69863,0.71236,0.73726,0.77269, &
     & 0.81660,0.86673,0.92114,0.97834,1.03722,1.09695,1.15690,1.21657, &
     & 1.27559,1.33364,1.39049,1.44595,1.49989,1.55223,   5693.6362    /
!     Se I    AEL-Kurucz     
      DATA (UARRAY(34,1,K),K=0,22) /                                    &
     & 1.60944,1.63609,1.74515,1.84072,1.91228,1.97906,2.03023,2.07892, &
     & 2.11830,2.15618,2.18802,2.21888,2.24702,2.27438,2.30269,2.33022, &
     & 2.36181,2.39243,2.42966,2.46555,2.50959,2.55177,    886.3636    /
!     Se II   AEL-Kurucz     
      DATA (UARRAY(34,2,K),K=0,22) /                                    &
     & 1.38629,1.38641,1.40364,1.47710,1.60249,1.71396,1.82985,1.93368, &
     & 2.02353,2.10596,2.17844,2.24601,2.31155,2.37305,2.43668,2.49651, &
     & 2.55762,2.61521,2.66966,2.72130,2.77041,2.81721,   1954.5454    /
!     Se III  AEL-Kurucz     
      DATA (UARRAY(34,3,K),K=0,22) /                                    &
     & 0.00000,1.09507,1.61706,1.87473,2.00224,2.11530,2.18335,2.24706, &
     & 2.29455,2.33988,2.38416,2.42657,2.47696,2.52493,2.58174,2.63549, &
     & 2.68819,2.73825,2.78593,2.83144,2.87497,2.91669,   2909.0908    /
!     Se IV   AEL            
      DATA (UARRAY(34,4,K),K=0,22) /                                    &
     & 0.69315,1.02866,1.33127,1.46710,1.54189,1.58969,1.62462,1.65407, &
     & 1.68235,1.71205,1.74460,1.78064,1.82029,1.86330,1.90923,1.95753, &
     & 2.00761,2.05891,2.11089,2.16311,2.21516,2.26673,   3904.0000    /
!     Se V    AEL            
      DATA (UARRAY(34,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00019,0.00691,0.04165,0.12034,0.23880,0.38216, &
     & 0.53613,0.69100,0.84124,0.98410,1.11840,1.24383,1.36059,1.46908, &
     & 1.56983,1.66342,1.75040,1.83132,1.90670,1.97701,   6209.0908    /
!     Br I    D&F            
      DATA (UARRAY(35,1,K),K=0,22) /                                    &
     & 1.38629,1.38992,1.42800,1.47550,1.51688,1.55661,1.58517,1.60727, &
     & 1.62576,1.64231,1.65982,1.68397,1.72265,1.78449,1.88076,2.01221, &
     & 2.18548,2.38524,2.60287,2.82979,3.05864,3.28396,   1076.3636    /
!     Br II   D&F            
      DATA (UARRAY(35,2,K),K=0,22) /                                    &
     & 1.60944,1.67932,1.83598,1.95315,2.04202,2.12363,2.18725,2.23485, &
     & 2.27474,2.31170,2.35336,2.41280,2.50811,2.65335,2.85819,3.10422, &
     & 3.31086,3.48203,3.62814,3.75560,3.86863,3.97018,   1963.6364    /
!     Br III  D&F            
      DATA (UARRAY(35,3,K),K=0,22) /                                    &
     & 1.38629,1.38875,1.46393,1.63088,1.80009,1.94476,2.06736,2.16824, &
     & 2.25423,2.33297,2.41423,2.51134,2.63969,2.80855,3.00151,3.16320, &
     & 3.30235,3.42447,3.53330,3.63143,3.72079,3.80282,   3263.6365    /
!     Br IV   AEL            
      DATA (UARRAY(35,4,K),K=0,22) /                                    &
     & 0.00000,0.97236,1.55712,1.83181,1.99505,2.10372,2.18173,2.24156, &
     & 2.29073,2.33414,2.37507,2.41559,2.45695,2.49973,2.54406,2.58980, &
     & 2.63662,2.68415,2.73198,2.77973,2.82706,2.87372,   4300.0000    /
!     Br V    AEL            
      DATA (UARRAY(35,5,K),K=0,22) /                                    &
     & 0.69315,1.02817,1.33086,1.46682,1.54214,1.59175,1.63060,1.66624, &
     & 1.70252,1.74114,1.78249,1.82628,1.87187,1.91853,1.96556,2.01236, &
     & 2.05846,2.10348,2.14718,2.18939,2.23000,2.26896,   5427.2729    /
!     Kr I    AEL-NIST       
      DATA (UARRAY(36,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00008, &
     & 0.00054,0.00271,0.01000,0.02930,0.07135,0.14898,0.27256,0.44481, &
     & 0.65893,0.90183,1.15925,1.41949,1.67445,1.91932,   1272.3636    /
!     Kr II   AEL-NIST       
      DATA (UARRAY(36,2,K),K=0,22) /                                    &
     & 1.38629,1.40188,1.47120,1.53278,1.57732,1.60965,1.63387,1.65279, &
     & 1.66860,1.68348,1.69996,1.72102,1.74985,1.78951,1.84238,1.90979, &
     & 1.99171,2.08688,2.19306,2.30754,2.42750,2.55036,   2232.7273    /
!     Kr III  AEL-NIST       
      DATA (UARRAY(36,3,K),K=0,22) /                                    &
     & 1.60944,1.71167,1.89735,2.03581,2.13730,2.21350,2.27271,2.32088, &
     & 2.36276,2.40240,2.44328,2.48816,2.53896,2.59669,2.66147,2.73274, &
     & 2.80942,2.89019,2.97367,3.05857,3.14374,3.22826,   3354.5454    /
!     Kr IV   AEL-NIST       
      DATA (UARRAY(36,4,K),K=0,22) /                                    &
     & 1.38629,1.39754,1.54956,1.76918,1.95871,2.10738,2.22451,2.31995, &
     & 2.40112,2.47314,2.53944,2.60219,2.66270,2.72170,2.77953,2.83627, &
     & 2.89189,2.94630,2.99937,3.05100,3.10107,3.14953,   4772.7271    /
!     Kr V    AEL-NIST       
      DATA (UARRAY(36,5,K),K=0,22) /                                    &
     & 0.00000,1.10558,1.67318,1.93802,2.09466,2.20089,2.28162,2.34954, &
     & 2.41140,2.47072,2.52914,2.58722,2.64497,2.70213,2.75835,2.81329, &
     & 2.86667,2.91829,2.96799,3.01570,3.06138,3.10504,   5881.8179    /
!     Rb I    D&F            
      DATA (UARRAY(37,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69322,0.69341,0.69360,0.69415,0.69647, &
     & 0.70179,0.71090,0.72595,0.74919,0.78420,0.83500,0.90261,0.99293, &
     & 1.10184,1.22793,1.37270,1.52690,1.68693,1.84909,    379.6364    /
!     Rb II   D&F            
      DATA (UARRAY(37,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00131, &
     & 0.00830,0.03229,0.10317,0.26398,0.54261,0.91489,1.36745,1.81048, &
     & 2.25411,2.66366,3.03711,3.37840,3.69087,3.97743,   2500.0000    /
!     Rb III  D&F            
      DATA (UARRAY(37,3,K),K=0,22) /                                    &
     & 1.38629,1.41290,1.49615,1.55521,1.59177,1.62705,1.65483,1.67356, &
     & 1.69366,1.72400,1.78414,1.90132,2.09941,2.39972,2.75871,3.15835, &
     & 3.54983,3.93353,4.29607,4.60443,4.83978,5.03015,   3636.3635    /
!     Rb IV   G   state
      DATA (UARRAY(37,4,K),K=0,22) /                                    &
     & 1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944, &
     & 1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944, &
     & 1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,   4781.8184    /
!     Rb V    U=1 assumed
      DATA (UARRAY(37,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   6454.5454    /
!     Sr I    MMD            
      DATA (UARRAY(38,1,K),K=0,22) /                                    &
     & 0.00000,0.00019,0.06162,0.36423,0.76811,1.12061,1.39906,1.61806, &
     & 1.79404,1.93941,2.06304,2.17143,2.26956,2.36141,2.45020,2.53849, &
     & 2.62817,2.72047,2.81601,2.91485,3.01662,3.12068,    517.4545    /
!     Sr II   MMD            
      DATA (UARRAY(38,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69328,0.69754,0.71873,0.76633,0.83888,0.92855, &
     & 1.02685,1.12760,1.22734,1.32481,1.42025,1.51465,1.60929,1.70534, &
     & 1.80364,1.90461,2.00820,2.11404,2.22149,2.32977,   1002.4545    /
!     Sr III  FAK-Kurucz     
      DATA (UARRAY(38,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00050,0.00100, &
     & 0.01636,0.03150,0.15143,0.25851,0.63949,0.91472,1.45818,1.80835, &
     & 2.31747,2.65321,3.07659,3.37312,3.71794,3.97387,   3909.0908    /
!     Sr IV   D&F            
      DATA (UARRAY(38,4,K),K=0,22) /                                    &
     & 1.38629,1.41928,1.50806,1.56706,1.60276,1.63723,1.66467,1.68562, &
     & 1.71058,1.74910,1.81957,1.94664,2.15103,2.45178,2.80730,3.16000, &
     & 3.42022,3.62652,3.79745,3.94338,4.07071,4.18363,   5181.8184    /
!     Sr V    G   state
      DATA (UARRAY(38,5,K),K=0,22) /                                    &
     & 1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944, &
     & 1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,1.60944, &
     & 1.60944,1.60944,1.60944,1.60944,1.60944,1.60944,   6509.0908    /
!     Y  I    MMD            
      DATA (UARRAY(39,1,K),K=0,22) /                                    &
     & 1.38629,1.72451,1.96120,2.06346,2.12296,2.17254,2.22956,2.30278, &
     & 2.39425,2.50150,2.62002,2.74506,2.87265,2.99983,3.12456,3.24556, &
     & 3.36207,3.47371,3.58032,3.68190,3.77856,3.87045,    580.0000    /
!     Y  II   MMD            
      DATA (UARRAY(39,2,K),K=0,22) /                                    &
     & 0.00000,1.49229,2.16447,2.47010,2.67655,2.83911,2.97580,3.09490, &
     & 3.20100,3.29696,3.38471,3.46563,3.54076,3.61092,3.67673,3.73871, &
     & 3.79725,3.85270,3.90533,3.95536,4.00300,4.04842,   1111.8182    /
!     Y  III  MMD            
      DATA (UARRAY(39,3,K),K=0,22) /                                    &
     & 1.38629,2.00642,2.15742,2.22707,2.27133,2.30311,2.32797,2.34891, &
     & 2.36765,2.38532,2.40274,2.42060,2.43955,2.46016,2.48297,2.50841, &
     & 2.53681,2.56836,2.60310,2.64096,2.68174,2.72517,   1863.6364    /
!     Y  IV   G   state
      DATA (UARRAY(39,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5618.1816    /
!     Y  V    AEL            
      DATA (UARRAY(39,5,K),K=0,22) /                                    &
     & 1.38629,1.42729,1.52140,1.58410,1.62465,1.65348,1.67785,1.70382, &
     & 1.73639,1.77912,1.83366,1.89981,1.97596,2.05970,2.14838,2.23957, &
     & 2.33121,2.42174,2.51003,2.59532,2.67715,2.75527,   7000.0000    /
!     Zr I    D&F            
      DATA (UARRAY(40,1,K),K=0,22) /                                    &
     & 1.60944,1.99885,2.38962,2.64841,2.85830,3.03510,3.20833,3.37227, &
     & 3.52966,3.68161,3.82886,3.97258,4.11376,4.25306,4.39138,4.53069, &
     & 4.67076,4.81062,4.95042,5.09150,5.22889,5.34966,    621.8182    /
!     Zr II   MMD            
      DATA (UARRAY(40,2,K),K=0,22) /                                    &
     & 1.38629,2.63977,3.18110,3.52334,3.77199,3.96457,4.12049,4.25142, &
     & 4.36471,4.46513,4.55574,4.63858,4.71504,4.78610,4.85246,4.91468, &
     & 4.97317,5.02830,5.08036,5.12960,5.17624,5.22048,   1193.6364    /
!     Zr III  D&F            
      DATA (UARRAY(40,3,K),K=0,22) /                                    &
     & 1.60944,2.54458,2.85219,3.04300,3.17808,3.29707,3.39984,3.47982, &
     & 3.54952,3.61217,3.67037,3.72683,3.78225,3.83477,3.88466,3.93218, &
     & 3.97755,4.02094,4.06253,4.10246,4.14086,4.17784,   2089.0908    /
!     Zr IV   D&F            
      DATA (UARRAY(40,4,K),K=0,22) /                                    &
     & 1.38629,1.99775,2.13990,2.18740,2.21376,2.23944,2.25823,2.27346, &
     & 2.28848,2.30452,2.32419,2.35116,2.37986,2.40776,2.43491,2.46134, &
     & 2.48708,2.51219,2.53667,2.56057,2.58392,2.60673,   3121.8181    /
!     Zr V    G   state
      DATA (UARRAY(40,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   7409.0908    /
!     Nb I    AEL-Kurucz     
      DATA (UARRAY(41,1,K),K=0,22) /                                    &
     & 0.69315,2.32885,2.93858,3.24505,3.42890,3.58414,3.72079,3.84100, &
     & 3.96919,4.08284,4.20683,4.31716,4.43645,4.54304,4.65638,4.75821, &
     & 4.86358,4.95892,5.05618,5.14483,5.23346,5.31488,    625.4545    /
!     Nb II   AEL-Kurucz     
      DATA (UARRAY(41,2,K),K=0,22) /                                    &
     & 0.00000,2.55352,3.14467,3.51197,3.81480,4.04692,4.26077,4.43685, &
     & 4.60072,4.74149,4.87641,4.99528,5.10999,5.21289,5.31183,5.40186, &
     & 5.48614,5.56386,5.63613,5.70353,5.76653,5.82580,   1301.8182    /
!     Nb III  AEL-Kurucz     
      DATA (UARRAY(41,3,K),K=0,22) /                                    &
     & 1.38629,2.71243,2.99753,3.10772,3.17222,3.23280,3.28524,3.33507, &
     & 3.38526,3.43305,3.48219,3.52903,3.57587,3.62061,3.66409,3.70576, &
     & 3.74575,3.78421,3.82114,3.85675,3.89123,3.92457,   2276.3635    /
!     Nb IV   AEL            
      DATA (UARRAY(41,4,K),K=0,22) /                                    &
     & 1.60944,2.76028,3.13859,3.31134,3.41383,3.48475,3.53934,3.58444, &
     & 3.62344,3.65812,3.68951,3.71825,3.74474,3.76929,3.79212,3.81341, &
     & 3.83330,3.85192,3.86939,3.88580,3.90123,3.91578,   3481.8181    /
!     Nb V    AEL            
      DATA (UARRAY(41,5,K),K=0,22) /                                    &
     & 1.38629,1.99347,2.13739,2.19019,2.21795,2.23600,2.24987,2.26200, &
     & 2.27355,2.28507,2.29676,2.30869,2.32081,2.33305,2.34533,2.35757, &
     & 2.36969,2.38164,2.39337,2.40483,2.41601,2.42688,   4595.4546    /
!     Mo I    AEL-NIST       
      DATA (UARRAY(42,1,K),K=0,22) /                                    &
     & 1.94591,1.94591,1.94592,1.94676,1.95347,1.97572,2.02343,2.10333, &
     & 2.21689,2.36028,2.52616,2.70613,2.89270,3.08007,3.26420,3.44254, &
     & 3.61358,3.77657,3.93123,4.07760,4.21591,4.34649,    645.4545    /
!     Mo II   AEL-NIST       
      DATA (UARRAY(42,2,K),K=0,22) /                                    &
     & 1.79176,1.79178,1.80654,1.93343,2.22147,2.59016,2.95806,3.29104, &
     & 3.58293,3.83741,4.06029,4.25701,4.43207,4.58903,4.73071,4.85932, &
     & 4.97665,5.08416,5.18304,5.27429,5.35877,5.43718,   1468.1818    /
!     Mo III  AEL-NIST       
      DATA (UARRAY(42,3,K),K=0,22) /                                    &
     & 0.00000,2.60006,2.97445,3.27708,3.55024,3.78136,3.97576,4.14209, &
     & 4.28740,4.41668,4.53331,4.63958,4.73708,4.82698,4.91014,4.98727, &
     & 5.05896,5.12571,5.18795,5.24609,5.30046,5.35139,   2466.3635    /
!     Mo IV   AEL-NIST       
      DATA (UARRAY(42,4,K),K=0,22) /                                    &
     & 1.38629,2.86327,3.33931,3.64922,3.86032,4.01533,4.13827,4.24214, &
     & 4.33396,4.41750,4.49478,4.56689,4.63445,4.69784,4.75734,4.81320, &
     & 4.86563,4.91484,4.96104,5.00443,5.04520,5.08354,   4218.1816    /
!     Mo V    AEL-NIST       
      DATA (UARRAY(42,5,K),K=0,22) /                                    &
     & 1.60944,2.65849,3.03217,3.22714,3.34630,3.42855,3.49170,3.54499, &
     & 3.59362,3.64070,3.68807,3.73670,3.78699,3.83891,3.89218,3.94637, &
     & 4.00103,4.05569,4.10994,4.16342,4.21586,4.26702,   5563.6362    /
!     Tc I    AEL-Kurucz     
      DATA (UARRAY(43,1,K),K=0,22) /                                    &
     & 1.79176,1.79954,1.94934,2.20405,2.43326,2.61961,2.77911,2.91664, &
     & 3.04548,3.15962,3.27207,3.37314,3.47349,3.56468,3.65480,3.73747, &
     & 3.81871,3.89385,3.96714,4.03543,4.10180,4.16404,    661.8182    /
!     Tc II   AEL-Kurucz     
      DATA (UARRAY(43,2,K),K=0,22) /                                    &
     & 1.94591,2.00129,2.30582,2.57493,2.75494,2.90748,3.02096,3.12287, &
     & 3.21169,3.29327,3.36906,3.43951,3.50592,3.56819,3.62719,3.68290, &
     & 3.73569,3.78584,3.83359,3.87916,3.92273,3.96448,   1387.2727    /
!     Tc III  AEL-Kurucz     
      DATA (UARRAY(43,3,K),K=0,22) /                                    &
     & 5.22949,5.26216,5.29380,5.32447,5.35423,5.38542,5.42715,5.46948, &
     & 5.51751,5.56570,5.61714,5.66853,5.72168,5.77451,5.82779,5.88067, &
     & 5.93323,5.98531,6.03649,6.08697,6.13610,6.18294,   2818.1819    /
!     Tc IV   U=1 assumed
      DATA (UARRAY(43,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4181.8184    /
!     Tc V    U=1 assumed
      DATA (UARRAY(43,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5000.0000    /
!     Ru I    AEL-Kurucz     
      DATA (UARRAY(44,1,K),K=0,22) /                                    &
     & 2.39790,2.54888,2.74959,2.87253,3.03982,3.18311,3.33243,3.46233, &
     & 3.59693,3.71554,3.83515,3.94197,4.04934,4.14629,4.24397,4.33294, &
     & 4.42266,4.50499,4.58783,4.66434,4.74097,4.81214,    669.4545    /
!     Ru II   AEL-Kurucz     
      DATA (UARRAY(44,2,K),K=0,22) /                                    &
     & 2.30259,2.54063,2.84893,3.09849,3.34692,3.54504,3.74754,3.91587, &
     & 4.08233,4.22501,4.36489,4.48760,4.60825,4.71591,4.82037,4.91496, &
     & 5.00639,5.09016,5.16944,5.24290,5.31158,5.37584,   1523.6364    /
!     Ru III  AEL-Kurucz     
      DATA (UARRAY(44,3,K),K=0,22) /                                    &
     & 2.19722,2.74654,2.95410,3.03792,3.08603,3.12325,3.15316,3.18221, &
     & 3.20923,3.23554,3.26213,3.28803,3.31509,3.34144,3.36781,3.39350, &
     & 3.41838,3.44266,3.46652,3.48982,3.51244,3.53456,   2587.2727    /
!     Ru IV   G   state
      DATA (UARRAY(44,4,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,   4545.4546    /
!     Ru V    U=1 assumed
      DATA (UARRAY(44,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5454.5454    /
!     Rh I    AEL-Kurucz     
      DATA (UARRAY(45,1,K),K=0,22) /                                    &
     & 2.30259,2.33637,2.50719,2.69757,2.86357,3.00566,3.12940,3.23950, &
     & 3.33967,3.43071,3.51733,3.59704,3.67569,3.74860,3.82232,3.89098, &
     & 3.96118,4.02679,4.09404,4.15705,4.22121,4.28150,    678.1818    /
!     Rh II   AEL-Kurucz     
      DATA (UARRAY(45,2,K),K=0,22) /                                    &
     & 2.19722,2.31017,2.54337,2.73049,2.90771,3.05821,3.22086,3.36073, &
     & 3.50734,3.63518,3.76548,3.88074,3.99620,4.09969,4.19216,4.27680, &
     & 4.37620,4.46662,4.54362,4.61511,4.68121,4.74320,   1642.7273    /
!     Rh III  AEL-Kurucz     
      DATA (UARRAY(45,3,K),K=0,22) /                                    &
     & 2.30259,2.65229,2.98596,3.23673,3.44460,3.61661,3.76635,3.89655, &
     & 4.02109,4.13182,4.24380,4.34449,4.43815,4.52379,4.60260,4.67565, &
     & 4.74365,4.80732,4.86757,4.92439,4.97781,5.02851,   2822.7273    /
!     Rh IV   G   state
      DATA (UARRAY(45,4,K),K=0,22) /                                    &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,   4363.6362    /
!     Rh V    U=1 assumed
      DATA (UARRAY(45,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5909.0908    /
!     Pd I    AEL-Kurucz     
      DATA (UARRAY(46,1,K),K=0,22) /                                    &
     & 0.00000,0.00003,0.01700,0.14355,0.44600,0.67564,0.95031,1.16560, &
     & 1.36216,1.52637,1.67498,1.80436,1.92816,2.03834,2.15055,2.25145, &
     & 2.35820,2.45468,2.55825,2.65214,2.75231,2.84336,    757.2727    /
!     Pd II   AEL-Kurucz     
      DATA (UARRAY(46,2,K),K=0,22) /                                    &
     & 1.79176,1.82835,1.93821,2.02156,2.09625,2.16575,2.26754,2.35992, &
     & 2.48620,2.59832,2.73282,2.85136,2.98398,3.10106,3.22733,3.33944, &
     & 3.44882,3.54741,3.63688,3.71900,3.79499,3.86562,   1765.4546    /
!     Pd III  AEL-Kurucz     
      DATA (UARRAY(46,3,K),K=0,22) /                                    &
     & 2.19722,2.40292,2.68743,2.89149,3.03206,3.15528,3.26633,3.36626, &
     & 3.46854,3.56132,3.64584,3.72376,3.79606,3.86347,3.92662,3.98601, &
     & 4.04209,4.09518,4.14559,4.19358,4.23938,4.28318,   2992.7273    /
!     Pd IV   G   state
      DATA (UARRAY(46,4,K),K=0,22) /                                    &
     & 2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259, &
     & 2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259, &
     & 2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,   4818.1816    /
!     Pd V    U=1 assumed
      DATA (UARRAY(46,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5636.3638    /
!     Ag I    D&F            
      DATA (UARRAY(47,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69315,0.69392, &
     & 0.69594,0.70038,0.70886,0.72474,0.75245,0.79588,0.85966,0.94772, &
     & 1.06180,1.19691,1.34949,1.51584,1.69155,1.86682,    688.5455    /
!     Ag II   D&F            
      DATA (UARRAY(47,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00001,0.00773,0.03216,0.05600,0.11579,0.23005, &
     & 0.37905,0.55490,0.74957,0.95968,1.18264,1.41288,1.65289,1.89046, &
     & 2.12205,2.30996,2.46809,2.60459,2.72468,2.83188,   1952.7273    /
!     Ag III  D&F            
      DATA (UARRAY(47,3,K),K=0,22) /                                    &
     & 1.79176,1.87068,2.00201,2.07422,2.12268,2.16890,2.21693,2.27943, &
     & 2.36096,2.46550,2.59691,2.75133,2.88790,3.00804,3.11528,3.21212, &
     & 3.30041,3.38154,3.45657,3.52636,3.59160,3.65285,   3165.4546    /
!     Ag IV   G   state
      DATA (UARRAY(47,4,K),K=0,22) /                                    &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,   5090.9092    /
!     Ag V    U=1 assumed
      DATA (UARRAY(47,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   6181.8184    /
!     Cd I    D&F            
      DATA (UARRAY(48,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00100,0.00366, &
     & 0.01024,0.02358,0.04750,0.08879,0.15561,0.25444,0.39643,0.57066, &
     & 0.78337,1.01423,1.25365,1.49451,1.73100,1.92811,    817.3636    /
!     Cd II   D&F            
      DATA (UARRAY(48,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69330,0.69368,0.69407,0.69605,0.70185, &
     & 0.71464,0.73857,0.78057,0.84850,0.95093,1.09693,1.28638,1.51164, &
     & 1.75737,2.02480,2.29130,2.55135,2.80130,3.03905,   1536.7273    /
!     Cd III  D&F            
      DATA (UARRAY(48,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00345,0.01549,0.02738,0.05999,0.13829, &
     & 0.26074,0.43170,0.64781,0.90431,1.14764,1.34320,1.50670,1.64718, &
     & 1.77034,1.87998,1.97878,2.06869,2.15118,2.22738,   3406.3635    /
!     Cd IV   AEL            
      DATA (UARRAY(48,4,K),K=0,22) /                                    &
     & 1.79176,1.92298,2.05853,2.12591,2.16693,2.20113,2.23931,2.28633, &
     & 2.34272,2.40658,2.47515,2.54578,2.61634,2.68531,2.75171,2.81496, &
     & 2.87478,2.93111,2.98399,3.03355,3.07996,3.12343,   5363.6362    /
!     Cd V    U=1 assumed
      DATA (UARRAY(48,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   6545.4546    /
!     In I    D&F            
      DATA (UARRAY(49,1,K),K=0,22) /                                    &
     & 0.69315,0.69784,0.78569,0.92597,1.05202,1.16064,1.24065,1.30408, &
     & 1.35523,1.39762,1.43388,1.46718,1.50057,1.53755,1.58072,1.63465, &
     & 1.69979,1.77958,1.87407,1.98088,2.08884,2.18628,    525.9091    /
!     In II   D&F            
      DATA (UARRAY(49,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00114,0.00418,0.00721,0.01851,0.04486, &
     & 0.08883,0.15287,0.24679,0.37935,0.50654,0.61936,0.72074,0.81277, &
     & 0.89705,0.97477,1.04688,1.11414,1.17716,1.23645,   1714.5454    /
!     In III  D&F            
      DATA (UARRAY(49,3,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69375,0.69611,0.69846,0.70400,0.71765, &
     & 0.73905,0.77082,0.81888,0.89051,0.96309,1.03075,1.09412,1.15372, &
     & 1.20996,1.26320,1.31376,1.36188,1.40779,1.45168,   2548.1819    /
!     In IV   D&F            
      DATA (UARRAY(49,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00178,0.00834,0.01486,0.03473,0.08963, &
     & 0.18126,0.30870,0.44729,0.56898,0.67746,0.77532,0.86444,0.94627, &
     & 1.02191,1.09223,1.15792,1.21956,1.27763,1.33250,   4945.4546    /
!     In V    AEL            
      DATA (UARRAY(49,5,K),K=0,22) /                                    &
     & 1.79176,1.93401,2.06881,2.13396,2.17111,2.19506,2.21208,2.22551, &
     & 2.23740,2.24903,2.26122,2.27441,2.28877,2.30433,2.32097,2.33854, &
     & 2.35686,2.37572,2.39495,2.41439,2.43389,2.45332,   7000.0000    /
!     Sn I    AEL-Kurucz     
      DATA (UARRAY(50,1,K),K=0,22) /                                    &
     & 0.00000,0.07816,0.47570,0.84431,1.10470,1.31110,1.45750,1.58517, &
     & 1.68057,1.76765,1.83634,1.90061,1.95501,2.00660,2.05405,2.09935, &
     & 2.14416,2.18704,2.23155,2.27416,2.31991,2.36365,    667.4545    /
!     Sn II   AEL-Kurucz     
      DATA (UARRAY(50,2,K),K=0,22) /                                    &
     & 0.69315,0.71305,0.87591,1.05191,1.17204,1.27927,1.34484,1.40637, &
     & 1.45023,1.49224,1.53041,1.56719,1.60838,1.64795,1.69695,1.74367, &
     & 1.79248,1.83903,1.88348,1.92605,1.96689,2.00613,   1329.8182    /
!     Sn III  AEL-Kurucz     
      DATA (UARRAY(50,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.01242,0.02469,0.07927,0.13102, &
     & 0.23151,0.32283,0.44917,0.56133,0.67115,0.77010,0.86035,0.94314, &
     & 1.01939,1.09024,1.15656,1.21875,1.27717,1.33237,   2771.8181    /
!     Sn IV   AEL            
      DATA (UARRAY(50,4,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69337,0.69554,0.70303,0.71860,0.74330, &
     & 0.77693,0.81883,0.86830,0.92471,0.98746,1.05587,1.12912,1.20632, &
     & 1.28648,1.36862,1.45180,1.53520,1.61807,1.69982,   3703.0908    /
!     Sn V    AEL            
      DATA (UARRAY(50,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00002,0.00075,0.00603,0.02466,0.06807, &
     & 0.14512,0.25771,0.40027,0.56282,0.73495,0.90826,1.07702,1.23789, &
     & 1.38918,1.53034,1.66145,1.78294,1.89545,1.99966,   6570.9092    /
!     Sb I    AEL-Kurucz     
      DATA (UARRAY(51,1,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38689,1.39625,1.43365,1.46971,1.53620,1.59854, &
     & 1.67024,1.73715,1.80403,1.86671,1.92752,1.98484,2.04201,2.09609, &
     & 2.15231,2.20554,2.26255,2.31648,2.37535,2.43095,    785.3636    /
!     Sb II   AEL-Kurucz     
      DATA (UARRAY(51,2,K),K=0,22) /                                    &
     & 0.00000,0.16730,0.71062,1.10889,1.36469,1.56820,1.70311,1.82197, &
     & 1.90799,1.98719,2.05207,2.11300,2.17105,2.22592,2.28401,2.33892, &
     & 2.40107,2.45959,2.51850,2.57414,2.62648,2.67622,   1500.0000    /
!     Sb III  AEL-Kurucz     
      DATA (UARRAY(51,3,K),K=0,22) /                                    &
     & 0.69315,0.72532,0.92084,1.14103,1.23692,1.32442,1.39291,1.45701, &
     & 1.51328,1.56656,1.62455,1.67937,1.74233,1.80157,1.86821,1.93068, &
     & 1.99154,2.04890,2.10314,2.15459,2.20347,2.25007,   2300.0000    /
!     Sb IV   AEL            
      DATA (UARRAY(51,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00004,0.00240,0.01873,0.06409,0.14373,0.25213, &
     & 0.37939,0.51627,0.65595,0.79398,0.92768,1.05553,1.17678,1.29114, &
     & 1.39865,1.49950,1.59401,1.68253,1.76545,1.84316,   4018.1819    /
!     Sb V    AEL            
      DATA (UARRAY(51,5,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69316,0.69400,0.69960,0.71488,0.74188,0.77970, &
     & 0.82614,0.87885,0.93581,0.99540,1.05639,1.11783,1.17896,1.23926, &
     & 1.29829,1.35575,1.41145,1.46524,1.51705,1.56683,   5090.9092    /
!     Te I    AEL-Kurucz     
      DATA (UARRAY(52,1,K),K=0,22) /                                    &
     & 1.60944,1.60963,1.62190,1.66053,1.71760,1.77159,1.82885,1.88301, &
     & 1.93214,1.97897,2.02007,2.05954,2.09493,2.12911,2.16205,2.19394, &
     & 2.22706,2.25912,2.29482,2.32931,2.36799,2.40523,    819.0909    /
!     Te II   AEL-Kurucz     
      DATA (UARRAY(52,2,K),K=0,22) /                                    &
     & 1.38629,1.38650,1.40672,1.48330,1.60561,1.71449,1.82600,1.92631, &
     & 2.01418,2.09495,2.16934,2.23858,2.31105,2.37861,2.45658,2.52891, &
     & 2.61484,2.69398,2.78565,2.86962,2.96425,3.05070,   1690.9091    /
!     Te III  AEL-Kurucz     
      DATA (UARRAY(52,3,K),K=0,22) /                                    &
     & 0.00000,0.29516,0.94528,1.33741,1.59086,1.77207,1.90599,2.01347, &
     & 2.10590,2.19276,2.27404,2.36334,2.45125,2.54760,2.64723,2.75071, &
     & 2.86076,2.96790,3.08333,3.19002,3.28641,3.37432,   2818.1819    /
!     Te IV   AEL            
      DATA (UARRAY(52,4,K),K=0,22) /                                    &
     & 0.69315,0.73276,0.94333,1.12798,1.25520,1.34471,1.41206,1.46725, &
     & 1.51665,1.56416,1.61198,1.66119,1.71208,1.76450,1.81804,1.87223, &
     & 1.92655,1.98054,2.03383,2.08609,2.13710,2.18667,   3400.9092    /
!     Te V    AEL            
      DATA (UARRAY(52,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00017,0.00622,0.03823,0.11245,0.22663,0.36743, &
     & 0.52108,0.67769,0.83124,0.97844,1.11767,1.24832,1.37032,1.48393, &
     & 1.58959,1.68781,1.77914,1.86412,1.94326,2.01706,   5340.9092    /
!     I  I    D&F            
      DATA (UARRAY(53,1,K),K=0,22) /                                    &
     & 1.38629,1.38630,1.38788,1.39944,1.41794,1.43609,1.45713,1.47851, &
     & 1.49871,1.51781,1.53727,1.55984,1.58922,1.63161,1.69020,1.77369, &
     & 1.87962,2.00869,2.16150,2.32821,2.50353,2.68253,    950.3636    /
!     I  II   D&F            
      DATA (UARRAY(53,2,K),K=0,22) /                                    &
     & 1.60944,1.61208,1.65719,1.73778,1.82711,1.90912,1.98276,2.04566, &
     & 2.10057,2.15017,2.19758,2.24609,2.29313,2.33805,2.38105,2.42227, &
     & 2.46185,2.49994,2.53662,2.57201,2.60618,2.63923,   1735.4546    /
!     I  III  AEL-Kurucz     
      DATA (UARRAY(53,3,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.44316,1.63746,1.80008,1.93991,2.05720,2.16216, &
     & 2.26435,2.35707,2.46937,2.57032,2.69868,2.81242,2.94918,3.06947, &
     & 3.20417,3.32286,3.44870,3.56047,3.67427,3.77643,   2909.0908    /
!     I  IV   G   state
      DATA (UARRAY(53,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3818.1819    /
!     I  V    U=1 assumed
      DATA (UARRAY(53,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   6000.0000    /
!     Xe I    D&F            
      DATA (UARRAY(54,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00189,0.00879,0.03177,0.09346,0.22019,0.43200,0.73023,1.07469, &
     & 1.44051,1.81165,2.16200,2.48896,2.79231,3.07286,   1102.4546    /
!     Xe II   D&F            
      DATA (UARRAY(54,2,K),K=0,22) /                                    &
     & 1.38629,1.38649,1.39604,1.42252,1.45397,1.48447,1.51297,1.53825, &
     & 1.56318,1.59177,1.62982,1.68435,1.76173,1.86755,2.00451,2.16668, &
     & 2.34662,2.50410,2.64013,2.75984,2.86674,2.96332,   1927.2727    /
!     Xe III  D&F            
      DATA (UARRAY(54,3,K),K=0,22) /                                    &
     & 1.60944,1.61804,1.70052,1.81714,1.92162,2.01620,2.09625,2.16071, &
     & 2.21929,2.27776,2.34217,2.41811,2.50999,2.59942,2.68150,2.75736, &
     & 2.82786,2.89372,2.95551,3.01370,3.06869,3.12082,   2918.1816    /
!     Xe IV   D&F            
      DATA (UARRAY(54,4,K),K=0,22) /                                    &
     & 1.38629,1.40026,1.55226,1.75388,1.92356,2.06858,2.18821,2.28962, &
     & 2.38322,2.47295,2.55605,2.63278,2.70403,2.77054,2.83291,2.89161, &
     & 2.94706,2.99959,3.04950,3.09704,3.14242,3.18583,   4181.8184    /
!     Xe V    U=1 assumed
      DATA (UARRAY(54,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5181.8184    /
!     Cs I    D&F            
      DATA (UARRAY(55,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69315,0.69334,0.69471,0.69806, &
     & 0.70518,0.71806,0.73849,0.76888,0.81230,0.87023,0.94816,1.04278, &
     & 1.15770,1.28858,1.41759,1.53185,1.63439,1.72738,    353.9091    /
!     Cs II   D&F            
      DATA (UARRAY(55,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00086,0.00304, &
     & 0.01271,0.03626,0.08667,0.17574,0.31246,0.50394,0.73509,0.99493, &
     & 1.26203,1.53773,1.80439,2.05809,2.29754,2.52237,   2281.8181    /
!     Cs III  D&F            
      DATA (UARRAY(55,3,K),K=0,22) /                                    &
     & 1.38629,1.38723,1.40772,1.44540,1.48200,1.51730,1.54801,1.57214, &
     & 1.59255,1.61054,1.62702,1.64267,1.65822,1.67353,1.68861,1.70346, &
     & 1.71810,1.73252,1.74674,1.76077,1.77459,1.78823,   3181.8181    /
!     Cs IV   U=1 assumed
      DATA (UARRAY(55,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4181.8184    /
!     Cs V    U=1 assumed
      DATA (UARRAY(55,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5636.3638    /
!     Ba I    D&F            
      DATA (UARRAY(56,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00001,0.01196,0.04438,0.07578,0.15999,0.30150, &
     & 0.47450,0.66056,0.85084,1.03815,1.22193,1.40191,1.57598,1.75003, &
     & 1.91929,2.08456,2.24874,2.40868,2.56386,2.71394,    473.6364    /
!     Ba II   D&F            
      DATA (UARRAY(56,2,K),K=0,22) /                                    &
     & 0.69315,0.69442,0.76633,0.96529,1.17652,1.35082,1.49751,1.61757, &
     & 1.71724,1.80206,1.87626,1.94381,2.00795,2.07187,2.13796,2.20897, &
     & 2.28540,2.36955,2.46135,2.55956,2.66325,2.77117,    909.1818    /
!     Ba III  FAK-Kurucz     
      DATA (UARRAY(56,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00126,0.00271, &
     & 0.01542,0.03063,0.08559,0.15001,0.28920,0.44039,0.66686,0.89315, &
     & 1.15969,1.41458,1.67661,1.92492,2.16226,2.38827,   3272.7273    /
!     Ba IV   D&F            
      DATA (UARRAY(56,4,K),K=0,22) /                                    &
     & 1.38629,1.38787,1.41399,1.45578,1.49450,1.53177,1.56455,1.59502, &
     & 1.62943,1.67050,1.71151,1.75091,1.78881,1.82532,1.86055,1.89458, &
     & 1.92749,1.95935,1.99023,2.02019,2.04927,2.07753,   4454.5454    /
!     Ba V    U=1 assumed
      DATA (UARRAY(56,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5636.3638    /
!     La I    AEL-NIST       
      DATA (UARRAY(57,1,K),K=0,22) /                                    &
     & 1.38629,1.45996,1.72394,2.00235,2.25477,2.47555,2.67062,2.84691, &
     & 3.00979,3.16293,3.30866,3.44841,3.58296,3.71272,3.83786,3.95842, &
     & 4.07439,4.18576,4.29253,4.39475,4.49247,4.58579,    507.0000    /
!     La II   AEL-NIST       
      DATA (UARRAY(57,2,K),K=0,22) /                                    &
     & 1.60944,2.10929,2.69786,3.01771,3.22692,3.38708,3.52355,3.64700, &
     & 3.76181,3.86966,3.97130,4.06709,4.15735,4.24241,4.32258,4.39818, &
     & 4.46953,4.53693,4.60066,4.66097,4.71811,4.77230,   1005.4545    /
!     La III  AEL-Cowley     
      DATA (UARRAY(57,3,K),K=0,22) /                                    &
     & 1.38629,1.74336,2.04370,2.25449,2.41290,2.53274,2.62532,2.69877, &
     & 2.75854,2.80832,2.85057,2.88704,2.91893,2.94715,2.97234,2.99502, &
     & 3.01558,3.03431,3.05148,3.06727,3.08187,3.09540,   1743.3636    /
!     La IV   AEL-NIST       
      DATA (UARRAY(57,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00028,0.00328,0.01734,0.05695, &
     & 0.13685,0.26377,0.43262,0.63007,0.84121,1.05407,1.26081,1.45712, &
     & 1.64104,1.81209,1.97057,2.11718,2.25281,2.37838,   4540.9092    /
!     La V    AEL-NIST       
      DATA (UARRAY(57,5,K),K=0,22) /                                    &
     & 1.38629,1.38822,1.41686,1.46177,1.50386,1.53996,1.57302,1.60800, &
     & 1.64975,1.70168,1.76515,1.83959,1.92301,2.01280,2.10626,2.20101, &
     & 2.29516,2.38730,2.47649,2.56211,2.64384,2.72155,   5600.0000    /
!     Ce I    AEL-NIST       
      DATA (UARRAY(58,1,K),K=0,22) /                                    &
     & 2.19722,2.49091,2.86307,3.29297,3.69217,4.04009,4.34223,4.60832, &
     & 4.84647,5.06243,5.26008,5.44206,5.61023,5.76601,5.91054,6.04479, &
     & 6.16960,6.28577,6.39401,6.49498,6.58927,6.67744,    503.5454    /
!     Ce II   AEL-NIST       
      DATA (UARRAY(58,2,K),K=0,22) /                                    &
     & 2.07944,2.64666,3.67412,4.38347,4.87131,5.22695,5.49946,5.71681, &
     & 5.89582,6.04710,6.17758,6.29193,6.39343,6.48443,6.56667,6.64149, &
     & 6.70993,6.77282,6.83084,6.88456,6.93445,6.98091,    986.3636    /
!     Ce III  AEL-Cowley     
      DATA (UARRAY(58,3,K),K=0,22) /                                    &
     & 2.19722,2.71316,3.46714,3.92395,4.21627,4.41971,4.57051,4.68780, &
     & 4.78251,4.86129,4.92836,4.98652,5.03770,5.08325,5.12419,5.16126, &
     & 5.19504,5.22600,5.25449,5.28081,5.30523,5.32794,   1836.1818    /
!     Ce IV   AEL-NIST       
      DATA (UARRAY(58,4,K),K=0,22) /                                    &
     & 1.79176,2.20082,2.39110,2.46778,2.51106,2.54258,2.56991,2.59575, &
     & 2.62103,2.64601,2.67074,2.69527,2.71963,2.74387,2.76803,2.79217, &
     & 2.81631,2.84047,2.86464,2.88880,2.91294,2.93702,   3341.6365    /
!     Ce V    AEL-NIST       
      DATA (UARRAY(58,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00007,0.00066,0.00309,0.00950, &
     & 0.02229,0.04344,0.07406,0.11426,0.16323,0.21959,0.28161,0.34755, &
     & 0.41580,0.48501,0.55408,0.62216,0.68866,0.75316,   5959.0908    /
!     Pr I    AEL-NIST       
      DATA (UARRAY(59,1,K),K=0,22) /                                    &
     & 2.30259,2.32496,2.47850,2.69198,2.93829,3.21106,3.49703,3.78167, &
     & 4.05477,4.31102,4.54844,4.76698,4.96750,5.15129,5.31976,5.47435, &
     & 5.61638,5.74712,5.86768,5.97908,6.08224,6.17795,    496.7273    /
!     Pr II   AEL-NIST       
      DATA (UARRAY(59,2,K),K=0,22) /                                    &
     & 2.19722,2.82860,3.40471,3.93096,4.39411,4.78041,5.09672,5.35670, &
     & 5.57293,5.75521,5.91084,6.04523,6.16246,6.26562,6.35709,6.43877, &
     & 6.51213,6.57839,6.63852,6.69335,6.74353,6.78963,    959.0909    /
!     Pr III  AEL-Cowley     
      DATA (UARRAY(59,3,K),K=0,22) /                                    &
     & 2.30259,2.79969,3.31079,3.78442,4.23438,4.61808,4.93270,5.19065, &
     & 5.40492,5.58575,5.74064,5.87506,5.99302,6.09752,6.19084,6.27475, &
     & 6.35065,6.41966,6.48270,6.54053,6.59376,6.64294,   1965.8182    /
!     Pr IV   AEL-NIST       
      DATA (UARRAY(59,4,K),K=0,22) /                                    &
     & 2.19722,2.87228,3.36656,3.62039,3.78202,3.90085,3.99790,4.08282, &
     & 4.16015,4.23210,4.29980,4.36394,4.42495,4.48317,4.53885,4.59218, &
     & 4.64333,4.69244,4.73962,4.78496,4.82854,4.87045,   3543.6365    /
!     Pr V    AEL-NIST       
      DATA (UARRAY(59,5,K),K=0,22) /                                    &
     & 1.79176,2.24902,2.42260,2.48997,2.52568,2.54855,2.56576,2.58064, &
     & 2.59482,2.60906,2.62367,2.63872,2.65419,2.66997,2.68595,2.70203, &
     & 2.71811,2.73410,2.74993,2.76553,2.78087,2.79591,   5230.0000    /
!     Nd I    AEL-NIST       
      DATA (UARRAY(60,1,K),K=0,22) /                                    &
     & 2.19722,2.24598,2.46092,2.69166,2.91046,3.13422,3.37652,3.63683, &
     & 3.90581,4.17309,4.43096,4.67494,4.90301,5.11468,5.31035,5.49092, &
     & 5.65746,5.81114,5.95309,6.08438,6.20602,6.31890,    502.2727    /
!     Nd II   AEL-NIST       
      DATA (UARRAY(60,2,K),K=0,22) /                                    &
     & 2.07944,2.74032,3.35596,3.85317,4.29902,4.69953,5.05174,5.35861, &
     & 5.62609,5.86033,6.06664,6.24941,6.41220,6.55793,6.68901,6.80741, &
     & 6.91482,7.01262,7.10199,7.18394,7.25933,7.32887,    975.4545    /
!     Nd III  AEL-Cowley     
      DATA (UARRAY(60,3,K),K=0,22) /                                    &
     & 2.19722,2.86293,3.41418,3.88981,4.37189,4.81647,5.20144,5.52731, &
     & 5.80221,6.03520,6.23417,6.40554,6.55438,6.68470,6.79964,6.90171, &
     & 6.99291,7.07487,7.14890,7.21608,7.27732,7.33336,   2009.0909    /
!     Nd IV   AEL-NIST       
      DATA (UARRAY(60,4,K),K=0,22) /                                    &
     & 2.30259,3.03240,3.56707,3.91733,4.15528,4.32262,4.44516,4.53820, &
     & 4.61103,4.66948,4.71738,4.75732,4.79112,4.82008,4.84516,4.86710, &
     & 4.88645,4.90363,4.91900,4.93282,4.94531,4.95667,   3672.7273    /
!     Nd V    G   state
      DATA (UARRAY(60,5,K),K=0,22) /                                    &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,   5454.5454    /
!     Pm I    AEL-NIST       
      DATA (UARRAY(61,1,K),K=0,22) /                                    &
     & 1.79176,1.92898,2.26849,2.56653,2.81412,3.01933,3.19030,3.33504, &
     & 3.46123,3.57571,3.68405,3.79023,3.89667,4.00442,4.11348,4.22320, &
     & 4.33267,4.44086,4.54689,4.65002,4.74971,4.84560,    504.9091    /
!     Pm II   AEL-NIST       
      DATA (UARRAY(61,2,K),K=0,22) /                                    &
     & 1.60944,2.49893,3.18117,3.59555,3.86944,4.07133,4.23831,4.38935, &
     & 4.53271,4.67086,4.80366,4.93022,5.04971,5.16169,5.26603,5.36292, &
     & 5.45269,5.53581,5.61275,5.68402,5.75011,5.81147,    990.9091    /
!     Pm III  AEL-Cowley     
      DATA (UARRAY(61,3,K),K=0,22) /                                    &
     & 1.79176,2.79817,3.48990,3.99005,4.46004,4.89876,5.29251,5.63870, &
     & 5.94091,6.20462,6.43542,6.63835,6.81771,6.97709,7.11946,7.24730, &
     & 7.36265,7.46719,7.56234,7.64928,7.72901,7.80239,   2027.2727    /
!     Pm IV   AEL-NIST       
      DATA (UARRAY(61,4,K),K=0,22) /                                    &
     & 2.19722,3.13549,3.63147,3.88574,4.04042,4.14396,4.21790,4.27327, &
     & 4.31624,4.35054,4.37854,4.40182,4.42149,4.43831,4.45287,4.46560, &
     & 4.47681,4.48676,4.49566,4.50365,4.51088,4.51745,   3736.3635    /
!     Pm V    G   state
      DATA (UARRAY(61,5,K),K=0,22) /                                    &
     & 2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259, &
     & 2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,2.30259, &
     & 2.30259,2.30259,2.30259,2.30259,2.30259,2.30259,   5454.5454    /
!     Sm I    AEL-NIST       
      DATA (UARRAY(62,1,K),K=0,22) /                                    &
     & 0.00000,1.08464,1.79419,2.23370,2.53343,2.75397,2.93392,3.09823, &
     & 3.26142,3.42995,3.60477,3.78371,3.96344,4.14068,4.31282,4.47807, &
     & 4.63533,4.78409,4.92423,5.05592,5.17948,5.29534,    513.0909    /
!     Sm II   AEL-NIST       
      DATA (UARRAY(62,2,K),K=0,22) /                                    &
     & 0.69315,2.14542,2.97054,3.44946,3.79814,4.08803,4.34489,4.57954, &
     & 4.79631,4.99667,5.18113,5.35023,5.50471,5.64558,5.77396,5.89101, &
     & 5.99787,6.09559,6.18514,6.26739,6.34310,6.41297,   1006.3636    /
!     Sm III  AEL-Cowley     
      DATA (UARRAY(62,3,K),K=0,22) /                                    &
     & 0.00000,2.56763,3.14679,3.45310,3.80934,4.26304,4.74788,5.20610, &
     & 5.61532,5.97316,6.28455,6.55607,6.79398,7.00368,7.18964,7.35553, &
     & 7.50434,7.63853,7.76012,7.87078,7.97191,8.06467,   2127.2727    /
!     Sm IV   AEL-NIST       
      DATA (UARRAY(62,4,K),K=0,22) /                                    &
     & 1.79176,3.16853,3.80984,4.14250,4.35344,4.49957,4.60654,4.68806, &
     & 4.75214,4.80380,4.84630,4.88186,4.91204,4.93797,4.96049,4.98022, &
     & 4.99765,5.01316,5.02705,5.03956,5.05089,5.06119,   3763.6365    /
!     Sm V    G   state
      DATA (UARRAY(62,5,K),K=0,22) /                                    &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,2.19722, &
     & 2.19722,2.19722,2.19722,2.19722,2.19722,2.19722,   5454.5454    /
!     Eu I    AEL-NIST       
      DATA (UARRAY(63,1,K),K=0,22) /                                    &
     & 2.07944,2.07944,2.07944,2.07948,2.08029,2.08556,2.10261,2.13975, &
     & 2.20274,2.29307,2.40853,2.54480,2.69690,2.86016,3.03045,3.20429, &
     & 3.37884,3.55179,3.72137,3.88626,4.04555,4.19865,    515.5455    /
!     Eu II   AEL-NIST       
      DATA (UARRAY(63,2,K),K=0,22) /                                    &
     & 2.19722,2.26877,2.41500,2.52852,2.64578,2.77512,2.91194,3.05222, &
     & 3.19431,3.33747,3.48089,3.62339,3.76358,3.90014,4.03198,4.15829, &
     & 4.27861,4.39268,4.50051,4.60219,4.69796,4.78808,   1021.9091    /
!     Eu III  AEL-Cowley     
      DATA (UARRAY(63,3,K),K=0,22) /                                    &
     & 2.07944,2.07944,2.08032,2.10989,2.28092,2.69247,3.27108,3.87848, &
     & 4.44039,4.93667,5.36896,5.74507,6.07354,6.36199,6.61685,6.84339, &
     & 7.04593,7.22799,7.39245,7.54172,7.67777,7.80227,   2265.4546    /
!     Eu IV   AEL-NIST       
      DATA (UARRAY(63,4,K),K=0,22) /                                    &
     & 0.00000,2.90584,3.36862,3.56579,3.68698,3.77256,3.83708,3.88766, &
     & 3.92842,3.96198,3.99008,4.01395,4.03447,4.05229,4.06792,4.08174, &
     & 4.09403,4.10505,4.11496,4.12395,4.13212,4.13958,   3881.8181    /
!     Eu V    G   state
      DATA (UARRAY(63,5,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,   5454.5454    /
!     Gd I    AEL-NIST       
      DATA (UARRAY(64,1,K),K=0,22) /                                    &
     & 1.60944,2.50992,2.96090,3.19398,3.35587,3.49862,3.64064,3.78763, &
     & 3.94021,4.09688,4.25541,4.41351,4.56917,4.72081,4.86727,5.00776, &
     & 5.14186,5.26936,5.39025,5.50468,5.61284,5.71503,    559.0909    /
!     Gd II   AEL-NIST       
      DATA (UARRAY(64,2,K),K=0,22) /                                    &
     & 1.79176,3.01297,3.58299,3.97807,4.28383,4.53472,4.74939,4.93851, &
     & 5.10827,5.26235,5.40304,5.53193,5.65023,5.75896,5.85901,5.95121, &
     & 6.03628,6.11490,6.18768,6.25517,6.31786,6.37621,   1099.0909    /
!     Gd III  AEL-Cowley     
      DATA (UARRAY(64,3,K),K=0,22) /                                    &
     & 1.60944,3.21714,3.72512,4.01631,4.25553,4.50639,4.78903,5.09712, &
     & 5.41456,5.72661,6.02377,6.30125,6.55749,6.79274,7.00820,7.20545, &
     & 7.38618,7.55204,7.70454,7.84507,7.97488,8.09505,   1875.4546    /
!     Gd IV   AEL-NIST       
      DATA (UARRAY(64,4,K),K=0,22) /                                    &
     & 2.07944,2.07946,2.08586,2.12377,2.19301,2.27503,2.35648,2.43144, &
     & 2.49821,2.55691,2.60835,2.65348,2.69322,2.72838,2.75963,2.78755, &
     & 2.81262,2.83524,2.85573,2.87438,2.89140,2.90701,   4000.0000    /
!     Gd V    G   state
      DATA (UARRAY(64,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5454.5454    /
!     Tb I    AEL-NIST       
      DATA (UARRAY(65,1,K),K=0,22) /                                    &
     & 2.77259,3.41456,3.77164,4.00162,4.19972,4.38315,4.55646,4.72168, &
     & 4.88013,5.03266,5.17972,5.32135,5.45742,5.58772,5.71205,5.83031, &
     & 5.94249,6.04867,6.14901,6.24373,6.33309,6.41737,    533.0909    /
!     Tb II   AEL-NIST       
      DATA (UARRAY(65,2,K),K=0,22) /                                    &
     & 2.83321,3.07186,3.55822,3.99836,4.33548,4.59616,4.80602,4.98191, &
     & 5.13428,5.26953,5.39156,5.50284,5.60498,5.69913,5.78615,5.86675, &
     & 5.94151,6.01098,6.07561,6.13583,6.19202,6.24452,   1047.2727    /
!     Tb III  AEL-Cowley     
      DATA (UARRAY(65,3,K),K=0,22) /                                    &
     & 2.77259,2.92747,3.41856,3.95986,4.44774,4.88216,5.27926,5.64853, &
     & 5.99393,6.31712,6.61905,6.90056,7.16248,7.40575,7.63141,7.84059, &
     & 8.03444,8.21415,8.38085,8.53563,8.67952,8.81346,   1991.8182    /
!     Tb IV   AEL-NIST       
      DATA (UARRAY(65,4,K),K=0,22) /                                    &
     & 2.56495,3.11617,3.43486,3.57219,3.65234,3.71272,3.76674,3.81889, &
     & 3.87028,3.92083,3.97013,4.01780,4.06357,4.10728,4.14885,4.18827, &
     & 4.22559,4.26089,4.29425,4.32577,4.35557,4.38374,   3579.0908    /
!     Tb V    G   state
      DATA (UARRAY(65,5,K),K=0,22) /                                    &
     & 2.07944,2.07944,2.07944,2.07944,2.07944,2.07944,2.07944,2.07944, &
     & 2.07944,2.07944,2.07944,2.07944,2.07944,2.07944,2.07944,2.07944, &
     & 2.07944,2.07944,2.07944,2.07944,2.07944,2.07944,   5454.5454    /
!     Dy I    AEL-NIST       
      DATA (UARRAY(66,1,K),K=0,22) /                                    &
     & 2.83321,2.83323,2.83690,2.85900,2.90907,2.99100,3.10688,3.25671, &
     & 3.43711,3.64116,3.86001,4.08496,4.30893,4.52685,4.73557,4.93337, &
     & 5.11953,5.29400,5.45710,5.60941,5.75160,5.88439,    539.9091    /
!     Dy II   AEL-NIST       
      DATA (UARRAY(66,2,K),K=0,22) /                                    &
     & 2.89037,3.14719,3.35926,3.55466,3.78818,4.06635,4.36817,4.66925, &
     & 4.95406,5.21583,5.45303,5.66670,5.85885,6.03182,6.18786,6.32901, &
     & 6.45710,6.57371,6.68023,6.77782,6.86752,6.95020,   1060.9091    /
!     Dy III  AEL-Cowley     
      DATA (UARRAY(66,3,K),K=0,22) /                                    &
     & 2.83321,2.87825,3.10227,3.44312,3.88482,4.35886,4.81079,5.21976, &
     & 5.58251,5.90229,6.18419,6.43340,6.65463,6.85193,7.02874,7.18793, &
     & 7.33192,7.46270,7.58197,7.69116,7.79146,7.88391,   2072.7273    /
!     Dy IV   AEL-NIST       
      DATA (UARRAY(66,4,K),K=0,22) /                                    &
     & 2.77259,3.08911,3.50257,3.74109,3.88573,3.98124,4.04864,4.09862, &
     & 4.13712,4.16766,4.19247,4.21302,4.23032,4.24507,4.25781,4.26892, &
     & 4.27869,4.28734,4.29507,4.30201,4.30827,4.31395,   3763.6365    /
!     Dy V    G   state
      DATA (UARRAY(66,5,K),K=0,22) /                                    &
     & 2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495, &
     & 2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495, &
     & 2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,   5454.5454    /
!     Ho I    AEL-NIST       
      DATA (UARRAY(67,1,K),K=0,22) /                                    &
     & 2.77259,2.77259,2.77335,2.78287,2.81708,2.88935,3.00373,3.15566, &
     & 3.33572,3.53314,3.73809,3.94287,4.14215,4.33263,4.51254,4.68119, &
     & 4.83854,4.98499,5.12112,5.24765,5.36531,5.47482,    547.4545    /
!     Ho II   AEL-NIST       
      DATA (UARRAY(67,2,K),K=0,22) /                                    &
     & 2.83321,3.15243,3.31272,3.43676,3.56449,3.69784,3.83558,3.97516, &
     & 4.11342,4.24750,4.37530,4.49556,4.60774,4.71178,4.80797,4.89675, &
     & 4.97865,5.05424,5.12407,5.18867,5.24853,5.30409,   1072.7273    /
!     Ho III  AEL-Cowley     
      DATA (UARRAY(67,3,K),K=0,22) /                                    &
     & 2.77259,2.79230,2.95940,3.27033,3.69457,4.16101,4.61160,5.02086, &
     & 5.38319,5.70134,5.98061,6.22656,6.44420,6.63783,6.81100,6.96668, &
     & 7.10732,7.23494,7.35124,7.45765,7.55534,7.64535,   2076.3635    /
!     Ho IV   AEL-NIST       
      DATA (UARRAY(67,4,K),K=0,22) /                                    &
     & 2.83321,3.00264,3.40396,3.75627,4.02298,4.22062,4.36960,4.48475, &
     & 4.57594,4.64973,4.71057,4.76153,4.80480,4.84199,4.87429,4.90258, &
     & 4.92758,4.94981,4.96972,4.98764,5.00386,5.01861,   3863.6365    /
!     Ho V    G   state
      DATA (UARRAY(67,5,K),K=0,22) /                                    &
     & 2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259, &
     & 2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259, &
     & 2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,   5454.5454    /
!     Er I    AEL-NIST       
      DATA (UARRAY(68,1,K),K=0,22) /                                    &
     & 2.56495,2.56495,2.56625,2.58144,2.63281,2.73265,2.87782,3.05752, &
     & 3.26044,3.47726,3.70051,3.92432,4.14417,4.35684,4.56024,4.75315, &
     & 4.93506,5.10596,5.26611,5.41601,5.55623,5.68741,    555.2727    /
!     Er II   AEL-NIST       
      DATA (UARRAY(68,2,K),K=0,22) /                                    &
     & 2.63906,3.03065,3.17551,3.34503,3.57638,3.84777,4.13026,4.40434, &
     & 4.66053,4.89547,5.10894,5.30215,5.47685,5.63493,5.77823,5.90843, &
     & 6.02703,6.13538,6.23463,6.32580,6.40979,6.48736,   1084.5454    /
!     Er III  AEL-Cowley     
      DATA (UARRAY(68,3,K),K=0,22) /                                    &
     & 2.56495,2.58873,2.76494,3.06420,3.44832,3.85439,4.23906,4.58639, &
     & 4.89443,5.16628,5.40638,5.61911,5.80838,5.97755,6.12944,6.26644, &
     & 6.39052,6.50338,6.60640,6.70079,6.78755,6.86757,   2067.2727    /
!     Er IV   AEL-NIST       
      DATA (UARRAY(68,4,K),K=0,22) /                                    &
     & 2.77259,2.87361,3.18429,3.43847,3.61728,3.74483,3.83909,3.91114, &
     & 3.96784,4.01354,4.05112,4.08254,4.10919,4.13208,4.15194,4.16934, &
     & 4.18470,4.19836,4.21059,4.22161,4.23157,4.24063,   3881.8181    /
!     Er V    G   state
      DATA (UARRAY(68,5,K),K=0,22) /                                    &
     & 2.83321,2.83321,2.83321,2.83321,2.83321,2.83321,2.83321,2.83321, &
     & 2.83321,2.83321,2.83321,2.83321,2.83321,2.83321,2.83321,2.83321, &
     & 2.83321,2.83321,2.83321,2.83321,2.83321,2.83321,   5454.5454    /
!     Tm I    AEL-NIST       
      DATA (UARRAY(69,1,K),K=0,22) /                                    &
     & 2.07944,2.07944,2.07945,2.07989,2.08288,2.09303,2.11797,2.16728, &
     & 2.24935,2.36796,2.52083,2.70081,2.89862,3.10537,3.31391,3.51914, &
     & 3.71775,3.90780,4.08833,4.25899,4.41985,4.57123,    562.1818    /
!     Tm II   AEL-NIST       
      DATA (UARRAY(69,2,K),K=0,22) /                                    &
     & 2.19722,2.64816,2.71010,2.75077,2.81805,2.93848,3.12085,3.35159, &
     & 3.60654,3.86467,4.11288,4.34485,4.55846,4.75378,4.93191,5.09436, &
     & 5.24270,5.37841,5.50288,5.61733,5.72284,5.82036,   1095.4546    /
!     Tm III  AEL-Cowley     
      DATA (UARRAY(69,3,K),K=0,22) /                                    &
     & 2.07944,2.08144,2.12036,2.23566,2.48880,2.85577,3.24781,3.60910, &
     & 3.92340,4.19203,4.42131,4.61799,4.78794,4.93594,5.06582,5.18061, &
     & 5.28275,5.37417,5.45647,5.53093,5.59860,5.66037,   2152.7273    /
!     Tm IV   AEL-NIST       
      DATA (UARRAY(69,4,K),K=0,22) /                                    &
     & 2.56495,2.69507,3.01993,3.26140,3.42502,3.53974,3.62374,3.68760, &
     & 3.73767,3.77793,3.81098,3.83858,3.86197,3.88204,3.89944,3.91468, &
     & 3.92813,3.94009,3.95079,3.96042,3.96914,3.97706,   3881.8181    /
!     Tm V    G   state
      DATA (UARRAY(69,5,K),K=0,22) /                                    &
     & 2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259, &
     & 2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,2.77259, &
     & 2.77259,2.77259,2.77259,2.77259,2.77259,2.77259,   5454.5454    /
!     Yb I    AEL-NIST       
      DATA (UARRAY(70,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00008,0.00089,0.00485,0.01724, &
     & 0.04620,0.10116,0.19003,0.31624,0.47737,0.66622,0.87340,1.08983, &
     & 1.30812,1.52300,1.73103,1.93021,2.11952,2.29864,    568.5455    /
!     Yb II   AEL-NIST       
      DATA (UARRAY(70,2,K),K=0,22) /                                    &
     & 2.07944,2.07944,2.07944,2.07951,2.08069,2.08730,2.10789,2.15295, &
     & 2.23105,2.34533,2.49247,2.66441,2.85143,3.04474,3.23759,3.42546, &
     & 3.60562,3.77665,3.93799,4.08960,4.23179,4.36502,   1107.6364    /
!     Yb III  AEL-Cowley     
      DATA (UARRAY(70,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00071,0.03436,0.23194,0.64020,1.12188,1.56713, &
     & 1.94574,2.26100,2.52378,2.74461,2.93210,3.09293,3.23222,3.35393, &
     & 3.46114,3.55625,3.64118,3.71747,3.78636,3.84887,   2277.2727    /
!     Yb IV   AEL-NIST       
      DATA (UARRAY(70,4,K),K=0,22) /                                    &
     & 2.07944,2.09761,2.19034,2.27716,2.35013,2.43919,2.57397,2.76119, &
     & 2.98484,3.22216,3.45563,3.67549,3.87773,4.06159,4.22792,4.37823, &
     & 4.51420,4.63748,4.74956,4.85177,4.94528,5.03109,   3960.0000    /
!     Yb V    G   state
      DATA (UARRAY(70,5,K),K=0,22) /                                    &
     & 2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495, &
     & 2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,2.56495, &
     & 2.56495,2.56495,2.56495,2.56495,2.56495,2.56495,   5454.5454    /
!     Lu I    AEL-NIST       
      DATA (UARRAY(71,1,K),K=0,22) /                                    &
     & 1.38629,1.39076,1.46612,1.58964,1.70783,1.80921,1.89514,1.96920, &
     & 2.03512,2.09623,2.15540,2.21491,2.27654,2.34146,2.41040,2.48359, &
     & 2.56091,2.64195,2.72608,2.81256,2.90063,2.98952,    493.2727    /
!     Lu II   AEL-NIST       
      DATA (UARRAY(71,2,K),K=0,22) /                                    &
     & 0.00000,0.00001,0.01021,0.11144,0.34635,0.64392,0.94060,1.21299, &
     & 1.45738,1.67602,1.87224,2.04916,2.20941,2.35518,2.48829,2.61025, &
     & 2.72236,2.82572,2.92128,3.00986,3.09216,3.16882,   1263.6364    /
!     Lu III  AEL-Cowley     
      DATA (UARRAY(71,3,K),K=0,22) /                                    &
     & 0.69315,0.72331,0.98424,1.28022,1.50242,1.66336,1.78383,1.87783, &
     & 1.95395,2.01749,2.07180,2.11912,2.16093,2.19833,2.23209,2.26281, &
     & 2.29092,2.31680,2.34072,2.36293,2.38361,2.40294,   1905.4546    /
!     Lu IV   AEL-NIST       
      DATA (UARRAY(71,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00103,0.01893,0.10711,0.31832,0.63685, &
     & 0.99934,1.35501,1.68100,1.97162,2.22850,2.45553,2.65693,2.83649, &
     & 2.99748,3.14260,3.27410,3.39382,3.50331,3.60384,   4113.6362    /
!     Lu V    AEL-NIST       
      DATA (UARRAY(71,5,K),K=0,22) /                                    &
     & 2.07944,2.12430,2.24960,2.33839,2.39747,2.44348,2.49088,2.55046, &
     & 2.62738,2.72121,2.82788,2.94202,3.05869,3.17407,3.28561,3.39181, &
     & 3.49191,3.58566,3.67313,3.75457,3.83034,3.90083,   6072.7271    /
!     Hf I    AEL-Kurucz     
      DATA (UARRAY(72,1,K),K=0,22) /                                    &
     & 1.60944,1.61627,1.71408,1.88589,2.07943,2.24152,2.40723,2.54935, &
     & 2.70000,2.83092,2.97581,3.10234,3.24256,3.36552,3.49880,3.61640, &
     & 3.74082,3.85146,3.96623,4.06918,4.17412,4.26909,    636.3636    /
!     Hf II   AEL-Kurucz     
      DATA (UARRAY(72,2,K),K=0,22) /                                    &
     & 1.38629,1.47288,1.89177,2.29373,2.63167,2.88383,3.12116,3.31286, &
     & 3.50070,3.65879,3.81661,3.95289,4.08943,4.20956,4.32898,4.43566, &
     & 4.54072,4.63579,4.72924,4.81471,4.89718,4.97336,   1354.5454    /
!     Hf III  AEL-Kurucz     
      DATA (UARRAY(72,3,K),K=0,22) /                                    &
     & 1.66771,2.23410,2.59348,2.85729,3.06584,3.23831,3.36940,3.46781, &
     & 3.55388,3.62241,3.68644,3.74063,3.79203,3.83876,3.88287,3.92424, &
     & 3.96336,4.00062,4.03586,4.06980,4.10191,4.13303,   1909.0909    /
!     Hf IV   U=1 assumed
      DATA (UARRAY(72,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3027.2727    /
!     Hf V    U=1 assumed
      DATA (UARRAY(72,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4545.4546    /
!     Ta I    AEL-Kurucz     
      DATA (UARRAY(73,1,K),K=0,22) /                                    &
     & 1.38629,1.41313,1.60887,1.87492,2.16406,2.38817,2.63998,2.84102, &
     & 3.06575,3.24922,3.45202,3.62061,3.80241,3.95622,4.11731,4.25606, &
     & 4.39790,4.52212,4.64661,4.75734,4.86702,4.96585,    716.3636    /
!     Ta II   AEL-Kurucz     
      DATA (UARRAY(73,2,K),K=0,22) /                                    &
     & 1.09861,1.75608,2.46706,2.95913,3.35199,3.63331,3.89088,4.09552, &
     & 4.28267,4.44025,4.58444,4.71043,4.82553,4.92873,5.02257,5.10836, &
     & 5.18655,5.25907,5.32496,5.38678,5.44307,5.49636,   1472.7274    /
!     Ta III  FAK-Kurucz     
      DATA (UARRAY(73,3,K),K=0,22) /                                    &
     & 0.00000,0.48234,1.66163,2.18784,2.53083,2.78575,3.00902,3.20529, &
     & 3.37933,3.54253,3.68680,3.82861,3.95278,4.07855,4.19026,4.30025, &
     & 4.40124,4.49820,4.58929,4.67570,4.75845,4.83587,   2000.0000    /
!     Ta IV   U=1 assumed
      DATA (UARRAY(73,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3000.0000    /
!     Ta V    U=1 assumed
      DATA (UARRAY(73,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4090.9092    /
!     W  I    AEL-Kurucz     
      DATA (UARRAY(74,1,K),K=0,22) /                                    &
     & 0.00000,0.12812,0.79278,1.38380,1.80466,2.09997,2.35927,2.56501, &
     & 2.77806,2.95367,3.14677,3.30862,3.48698,3.63834,3.80283,3.94412, &
     & 4.09460,4.22543,4.36199,4.48217,4.60542,4.71514,    725.4545    /
!     W  II   AEL-Kurucz     
      DATA (UARRAY(74,2,K),K=0,22) /                                    &
     & 0.69315,1.26700,2.02873,2.57577,3.05926,3.38382,3.71114,3.95733, &
     & 4.19166,4.38138,4.56095,4.71313,4.85827,4.98499,5.10475,5.21169, &
     & 5.31246,5.40400,5.48975,5.56871,5.64261,5.71142,   1609.0909    /
!     W  III  FAK-Kurucz     
      DATA (UARRAY(74,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,2.22565,2.89962,3.29861,3.58304,3.83485,4.04231, &
     & 4.22646,4.38625,4.52917,4.65692,4.77193,4.87639,4.97110,5.05778, &
     & 5.13722,5.21033,5.27796,5.34015,5.39829,5.45179,   2181.8181    /
!     W  IV   U=1 assumed
      DATA (UARRAY(74,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3181.8181    /
!     W  V    U=1 assumed
      DATA (UARRAY(74,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4363.6362    /
!     Re I    AEL-Kurucz     
      DATA (UARRAY(75,1,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79178,1.79327,1.81719,1.84059,1.93805,2.02688, &
     & 2.19725,2.34279,2.54240,2.70872,2.90634,3.07128,3.25344,3.40747, &
     & 3.57174,3.71279,3.85941,3.98725,4.11783,4.23331,    715.4545    /
!     Re II   AEL-Kurucz     
      DATA (UARRAY(75,2,K),K=0,22) /                                    &
     & 1.94591,1.94591,1.94983,1.99746,2.20942,2.38428,2.70539,2.94806, &
     & 3.22541,3.44232,3.66034,3.83924,4.01305,4.16109,4.30368,4.42846, &
     & 4.54754,4.65393,4.75486,4.84653,4.93374,5.01396,   1509.0909    /
!     Re III  FAK-Kurucz     
      DATA (UARRAY(75,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,1.60326,2.57564,3.05918,3.38374,3.71117,3.95743, &
     & 4.19171,4.38138,4.56095,4.71313,4.85827,4.98499,5.10473,5.21165, &
     & 5.31245,5.40402,5.48976,5.56873,5.64262,5.71143,   2363.6365    /
!     Re IV   U=1 assumed 
      DATA (UARRAY(75,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3454.5454    /
!     Re V    U=1 assumed
      DATA (UARRAY(75,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4636.3638    /
!     Os I    AEL-Kurucz     
      DATA (UARRAY(76,1,K),K=0,22) /                                    &
     & 2.19722,2.20153,2.27202,2.41916,2.60252,2.76394,2.94580,3.10652, &
     & 3.27569,3.42755,3.58329,3.72525,3.86818,4.00013,4.13097,4.25302, &
     & 4.37255,4.48489,4.59368,4.69628,4.79457,4.88406,    790.9091    /
!     Os II   AEL-Kurucz     
      DATA (UARRAY(76,2,K),K=0,22) /                                    &
     & 2.30259,2.34823,2.57469,2.82263,3.06320,3.25698,3.43675,3.58908, &
     & 3.72281,3.84075,3.94382,4.03724,4.11977,4.19599,4.26380,4.32730, &
     & 4.38420,4.43804,4.48643,4.53259,4.57470,4.61511,   1545.4546    /
!     Os III  FAK-Kurucz     
      DATA (UARRAY(76,3,K),K=0,22) /                                    &
     & 0.76530,1.32463,1.68119,1.94349,2.15109,2.32291,2.57512,2.82959, &
     & 3.07050,3.30308,3.50357,3.69598,3.86011,4.01894,4.15596,4.28879, &
     & 4.40648,4.51801,4.61974,4.71500,4.80347,4.88626,   2272.7273    /
!     Os IV   U=1 assumed
      DATA (UARRAY(76,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3636.3635    /
!     Os V    U=1 assumed
      DATA (UARRAY(76,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4909.0908    /
!     Ir I    AEL-Kurucz     
      DATA (UARRAY(77,1,K),K=0,22) /                                    &
     & 2.30259,2.30974,2.39996,2.55209,2.72725,2.87662,3.03091,3.16577, &
     & 3.29974,3.41899,3.53914,3.64772,3.76022,3.86314,3.97176,4.07205, &
     & 4.17792,4.27630,4.37910,4.47512,4.57340,4.66548,    818.1818    /
!     Ir II   FAK-Kurucz     
      DATA (UARRAY(77,2,K),K=0,22) /                                    &
     & 2.19722,2.19722,2.19722,2.31103,2.48186,2.62772,2.76057,2.91429, &
     & 3.04749,3.19033,3.32205,3.45010,3.57833,3.69340,3.81712,3.92720, &
     & 4.03948,4.14530,4.24642,4.34728,4.43889,4.53418,   1545.4546    /
!     Ir III  FAK-Kurucz     
      DATA (UARRAY(77,3,K),K=0,22) /                                    &
     & 1.13008,2.01589,2.47813,2.79303,3.03216,3.22500,3.40084,3.55311, &
     & 3.68638,3.80433,3.90821,4.00157,4.08523,4.16132,4.23047,4.29379, &
     & 4.35209,4.40571,4.45560,4.50151,4.54481,4.58496,   2454.5454    /
!     Ir IV   U=1 assumed
      DATA (UARRAY(77,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3545.4546    /
!     Ir V    U=1 assumed
      DATA (UARRAY(77,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5181.8184    /
!     Pt I    AEL-Kurucz     
      DATA (UARRAY(78,1,K),K=0,22) /                                    &
     & 1.94591,2.34102,2.63268,2.76561,2.84488,2.91830,2.97707,3.03256, &
     & 3.08534,3.13547,3.18696,3.23593,3.28972,3.34077,3.39903,3.45409, &
     & 3.51734,3.57684,3.64390,3.70676,3.77621,3.84116,    818.1818    /
!     Pt II   AEL-Kurucz     
      DATA (UARRAY(78,2,K),K=0,22) /                                    &
     & 1.79176,1.82045,2.02715,2.28334,2.53648,2.73831,2.94319,3.11317, &
     & 3.28154,3.42561,3.56747,3.69170,3.81472,3.92426,4.03388,4.13267, &
     & 4.23237,4.32303,4.41477,4.49881,4.58353,4.66163,   1687.2727    /
!     Pt III  FAK-Kurucz     
      DATA (UARRAY(78,3,K),K=0,22) /                                    &
     & 2.19722,2.19722,2.19722,2.38430,2.56381,2.71594,2.88761,3.04136, &
     & 3.20252,3.34924,3.49808,3.63619,3.77320,3.90256,4.02853,4.14918, &
     & 4.26479,4.37663,4.48253,4.58573,4.68240,4.77658,   2545.4546    /
!     Pt IV   U=1 assumed
      DATA (UARRAY(78,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3727.2727    /
!     Pt V    U=1 assumed
      DATA (UARRAY(78,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5000.0000    /
!     Au I    AEL-Kurucz     
      DATA (UARRAY(79,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69430,0.70890,0.76382,0.81580,0.89976,0.97723, &
     & 1.06160,1.13940,1.21860,1.29198,1.36993,1.44224,1.52378,1.59916, &
     & 1.68722,1.76816,1.86256,1.94882,2.04726,2.13687,    838.1818    /
!     Au II   AEL-Kurucz     
      DATA (UARRAY(79,2,K),K=0,22) /                                    &
     & 0.00000,0.00007,0.02637,0.18604,0.51223,0.75756,1.04912,1.27459, &
     & 1.49930,1.68266,1.86484,2.01890,2.17419,2.30859,2.44711,2.56876, &
     & 2.69528,2.80758,2.92477,3.02966,3.13855,3.23673,   1863.6364    /
!     Au III  AEL-Kurucz     
      DATA (UARRAY(79,3,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79176,1.87441,1.98781,2.08966,2.24930,2.38692, &
     & 2.55141,2.69263,2.84376,2.97502,3.10994,3.22881,3.34779,3.45411, &
     & 3.55891,3.65376,3.74597,3.83038,3.91162,3.98675,   2727.2727    /
!     Au IV   U=1 assumed
      DATA (UARRAY(79,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4000.0000    /
!     Au V    U=1 assumed
      DATA (UARRAY(79,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5272.7271    /
!     Hg I    D&F            
      DATA (UARRAY(80,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00100, &
     & 0.00407,0.01077,0.02494,0.05031,0.09704,0.17159,0.27666,0.42311, &
     & 0.59711,0.79319,1.01045,1.23248,1.45311,1.66820,    948.1818    /
!     Hg II   D&F            
      DATA (UARRAY(80,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69471,0.69905,0.70338,0.71608,0.74210, &
     & 0.78115,0.83381,0.90223,0.98911,1.09849,1.23289,1.38876,1.57303, &
     & 1.77144,1.97930,2.19685,2.41395,2.62623,2.83129,   1704.6364    /
!     Hg III  D&F            
      DATA (UARRAY(80,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00047,0.03912,0.15146,0.25244,0.40258,0.62144, &
     & 0.85139,1.08171,1.30783,1.52308,1.70483,1.85857,1.99179,2.10933, &
     & 2.21450,2.30965,2.39653,2.47646,2.55048,2.61939,   3109.0908    /
!     Hg IV   D&F            
      DATA (UARRAY(80,4,K),K=0,22) /                                    &
     & 2.30259,2.76911,3.26235,3.54848,3.74198,3.90404,4.03966,4.15648, &
     & 4.26250,4.35930,4.44786,4.52908,4.60420,4.67407,4.73937,4.80067, &
     & 4.85843,4.91303,4.96481,5.01403,5.06095,5.10576,   4181.8184    /
!     Hg V    U=1 assumed
      DATA (UARRAY(80,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5545.4546    /
!     Tl I    AEL-Kurucz     
      DATA (UARRAY(81,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69323,0.69515,0.71173,0.72806,0.76434,0.79936, &
     & 0.84399,0.88672,0.93243,0.97616,1.02434,1.07034,1.12700,1.18068, &
     & 1.25194,1.31854,1.40696,1.48829,1.59149,1.68504,    555.0909    /
!     Tl II   AEL-Kurucz     
      DATA (UARRAY(81,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00100,0.00200,0.01045,0.01883, &
     & 0.04833,0.07700,0.13719,0.19398,0.29197,0.38123,0.51857,0.63936, &
     & 0.80569,0.94830,1.12669,1.27807,1.45422,1.60394,   1856.3636    /
!     Tl III  AEL-Kurucz     
      DATA (UARRAY(81,3,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.70087,0.70853,0.74123,0.77290, &
     & 0.83878,0.90060,0.99772,1.08625,1.21020,1.32049,1.46274,1.58727, &
     & 1.73782,1.86866,2.01825,2.14837,2.29161,2.41689,   2709.0908    /
!     Tl IV   AEL            
      DATA (UARRAY(81,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00008,0.00462,0.03500,0.11694,0.25537,0.43516, &
     & 0.63569,0.84056,1.03961,1.22758,1.40229,1.56327,1.71097,1.84627, &
     & 1.97022,2.08386,2.18822,2.28425,2.37279,2.45461,   4609.0908    /
!     Tl V    G   state
      DATA (UARRAY(81,5,K),K=0,22) /                                    &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,1.79176, &
     & 1.79176,1.79176,1.79176,1.79176,1.79176,1.79176,   5818.1816    /
!     Pb I    D&F            
      DATA (UARRAY(82,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00077,0.03314,0.09736,0.15771,0.26351,0.39091, &
     & 0.52086,0.64759,0.77099,0.89405,1.02088,1.15725,1.30752,1.47175, &
     & 1.64861,1.83999,2.03655,2.23490,2.43186,2.62486,    674.0909    /
!     Pb II   D&F            
      DATA (UARRAY(82,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69435,0.71356,0.75721,0.79903,0.84921,0.90909, &
     & 0.96753,1.02389,1.07929,1.13650,1.19875,1.27113,1.34365,1.45840, &
     & 1.57527,1.68533,1.78446,1.87465,1.95737,2.03377,   1366.1818    /
!     Pb III  D&F            
      DATA (UARRAY(82,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00185,0.00744,0.01301,0.03000,0.06643, &
     & 0.12260,0.19964,0.30219,0.43623,0.60989,0.82183,1.06995,1.34303, &
     & 1.63262,1.88630,2.08846,2.25655,2.40041,2.52616,   2902.7273    /
!     Pb IV   D&F            
      DATA (UARRAY(82,4,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69429,0.69849,0.70267,0.71314,0.73985, &
     & 0.78342,0.84869,0.94050,1.06509,1.22645,1.42787,1.61417,1.77115, &
     & 1.90679,2.02622,2.13289,2.22927,2.31718,2.39798,   3847.2727    /
!     Pb V    D&F            
      DATA (UARRAY(82,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00003,0.01612,0.06797,0.11726,0.22038,0.41598, &
     & 0.65464,0.91529,1.17731,1.42951,1.64947,1.82966,1.98229,2.11467, &
     & 2.23157,2.33621,2.43094,2.51746,2.59709,2.67085,   6254.5454    /
!     Bi I    MMD            
      DATA (UARRAY(83,1,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38630,1.38657,1.38867,1.39514,1.40790,1.42749, &
     & 1.45337,1.48444,1.51953,1.55759,1.59783,1.63970,1.68285,1.72708, &
     & 1.77228,1.81835,1.86522,1.91283,1.96107,2.00982,    662.4545    /
!     Bi II   MMD            
      DATA (UARRAY(83,2,K),K=0,22) /                                    &
     & 0.00000,0.00001,0.00692,0.06518,0.19624,0.36809,0.54670,0.71482, &
     & 0.86711,1.00378,1.12718,1.24016,1.34542,1.44524,1.54134,1.63495, &
     & 1.72684,1.81740,1.90678,1.99493,2.08168,2.16684,   1516.3636    /
!     Bi III  MMD            
      DATA (UARRAY(83,3,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69635,0.72017,0.77020,0.83527,0.90483,0.97369, &
     & 1.04057,1.10616,1.17183,1.23898,1.30873,1.38176,1.45831,1.53824, &
     & 1.62111,1.70628,1.79300,1.88050,1.96807,2.05505,   2323.6365    /
!     Bi IV   AEL            
      DATA (UARRAY(83,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00001,0.00076,0.00720,0.02848,0.07193,0.13935, &
     & 0.22780,0.33203,0.44654,0.56653,0.68829,0.80910,0.92710,1.04108, &
     & 1.15028,1.25433,1.35309,1.44657,1.53491,1.61830,   4118.1816    /
!     Bi V    AEL            
      DATA (UARRAY(83,5,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69362,0.69787,0.71248,0.74308,0.79133, &
     & 0.85533,0.93130,1.01513,1.10325,1.19291,1.28215,1.36963,1.45450, &
     & 1.53627,1.61466,1.68957,1.76099,1.82897,1.89365,   5090.9092    /
!     Po I    AEL-Kurucz     
      DATA (UARRAY(84,1,K),K=0,22) /                                    &
     & 1.60944,1.60944,1.60961,1.61124,1.61691,1.62255,1.63383,1.64499, &
     & 1.66174,1.67823,1.69983,1.72098,1.74780,1.77393,1.80690,1.83883, &
     & 1.87918,1.91796,1.96649,2.01276,2.06900,2.12224,    766.3636    /
!     Po II   FAK-Kurucz     
      DATA (UARRAY(84,2,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38629,1.48996,1.61553,1.72707,1.83945,1.94046, &
     & 2.02862,2.10963,2.18522,2.25549,2.32970,2.39879,2.47947,2.55413, &
     & 2.64333,2.72522,2.81988,2.90637,3.00320,3.09148,   1727.2727    /
!     Po III  FAK-Kurucz     
      DATA (UARRAY(84,3,K),K=0,22) /                                    &
     & 0.00000,0.34912,0.87269,1.21457,1.46888,1.67145,1.80368,1.92044, &
     & 2.00989,2.09198,2.16797,2.23858,2.31550,2.38693,2.47107,2.54867, &
     & 2.64121,2.72590,2.82525,2.91562,3.01800,3.11086,   2454.5454    /
!     Po IV   U=1 assumed
      DATA (UARRAY(84,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3454.5454    /
!     Po V    U=1 assumed
      DATA (UARRAY(84,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5545.4546    /
!     At I    FAK-Kurucz     
      DATA (UARRAY(85,1,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38629,1.39328,1.40977,1.42600,1.44421,1.46210, &
     & 1.48070,1.49896,1.51601,1.53277,1.55414,1.57505,1.61253,1.64865, &
     & 1.71818,1.78320,1.89435,1.99439,2.14194,2.27051,    845.4545    /
!     At II   FAK-Kurucz     
      DATA (UARRAY(85,2,K),K=0,22) /                                    &
     & 1.40884,1.53638,1.64949,1.75109,1.84331,1.92774,1.99857,2.06472, &
     & 2.11902,2.17053,2.22207,2.27109,2.31896,2.36463,2.40875,2.45100, &
     & 2.49114,2.52972,2.56688,2.60270,2.63763,2.67138,   1818.1818    /
!     At III  FAK-Kurucz     
      DATA (UARRAY(85,3,K),K=0,22) /                                    &
     & 0.89778,1.18482,1.40758,1.58965,1.74362,1.87701,1.99203,2.09518, &
     & 2.18869,2.27420,2.36931,2.45616,2.56417,2.66165,2.78130,2.88817, &
     & 3.01187,3.12193,3.24240,3.34990,3.46260,3.56387,   2636.3635    /
!     At IV   U=1 assumed
      DATA (UARRAY(85,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3727.2727    /
!     At V    U=1 assumed
      DATA (UARRAY(85,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4636.3638    /
!     Rn I    AEL-Kurucz     
      DATA (UARRAY(86,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00150,0.00301,0.01293,0.02278,0.06207,0.09993,0.19692,0.28548, &
     & 0.45036,0.59201,0.80047,0.97299,1.19183,1.37128,    976.9091    /
!     Rn II   FAK-Kurucz     
      DATA (UARRAY(86,2,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.38630,1.39576,1.41303,1.43058,1.44961,1.46865, &
     & 1.48804,1.50666,1.52442,1.54458,1.56674,1.59958,1.63795,1.70126, &
     & 1.77055,1.87856,1.98502,2.13520,2.27057,2.38978,   1909.0909    /
!     Rn III  FAK-Kurucz     
      DATA (UARRAY(86,3,K),K=0,22) /                                    &
     & 1.40886,1.53241,1.64235,1.74140,1.83151,1.91417,1.98481,2.04966, &
     & 2.10483,2.15537,2.20520,2.25340,2.30008,2.34508,2.38839,2.43010, &
     & 2.46993,2.50804,2.54474,2.58015,2.61447,2.64787,   2636.3635    /
!     Rn IV   U=1 assumed
      DATA (UARRAY(86,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4000.0000    /
!     Rn V    U=1 assumed
      DATA (UARRAY(86,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5000.0000    /
!     Fr I    FAK-Kurucz     
      DATA (UARRAY(87,1,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69315,0.69315,0.69340,0.69365,0.69639,0.69914, &
     & 0.71074,0.72221,0.75141,0.77978,0.83922,0.89535,0.99563,1.08677, &
     & 1.22748,1.35087,1.47578,1.58681,1.68664,1.77741,    363.6364    /
!     Fr II   FAK-Kurucz     
      DATA (UARRAY(87,2,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00150,0.00300,0.01292,0.02274,0.06204,0.09985,0.19680,0.28518, &
     & 0.45012,0.59167,0.80021,0.97267,1.19165,1.37118,   2000.0000    /
!     Fr III  FAK-Kurucz     
      DATA (UARRAY(87,3,K),K=0,22) /                                    &
     & 1.32641,1.36533,1.40278,1.43889,1.47373,1.50741,1.53568,1.56318, &
     & 1.58258,1.60161,1.61731,1.63276,1.64760,1.66222,1.67663,1.69083, &
     & 1.70475,1.71847,1.73210,1.74554,1.75881,1.77190,   3000.0000    /
!     Fr IV   U=1 assumed
      DATA (UARRAY(87,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   3909.0908    /
!     Fr V    U=1 assumed
      DATA (UARRAY(87,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5363.6362    /
!     Ra I    AEL-Kurucz     
      DATA (UARRAY(88,1,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00200,0.00402,0.02766,0.05087, &
     & 0.12455,0.19325,0.31682,0.42691,0.57669,0.70706,0.86134,0.99503, &
     & 1.14264,1.27130,1.40725,1.52699,1.65040,1.76025,    479.7273    /
!     Ra II   AEL-Kurucz     
      DATA (UARRAY(88,2,K),K=0,22) /                                    &
     & 0.69315,0.69315,0.69337,0.69924,0.73978,0.77886,0.86711,0.94821, &
     & 1.05221,1.14641,1.24744,1.33920,1.43439,1.52130,1.61353,1.69797, &
     & 1.78975,1.87381,1.96630,2.05096,2.14319,2.22763,    922.1818    /
!     Ra III  FAK-Kurucz     
      DATA (UARRAY(88,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00150,0.00300,0.01292,0.02274,0.06203,0.09984,0.19679,0.28517, &
     & 0.45012,0.59170,0.80023,0.97269,1.19151,1.37093,   3090.9092    /
!     Ra IV   U=1 assumed
      DATA (UARRAY(88,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4181.8184    /
!     Ra V    U=1 assumed
      DATA (UARRAY(88,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5272.7271    /
!     Ac I    AEL-Kurucz     
      DATA (UARRAY(89,1,K),K=0,22) /                                    &
     & 1.38629,1.39523,1.49614,1.62903,1.75347,1.86401,1.99602,2.11262, &
     & 2.25867,2.38612,2.53136,2.65816,2.79357,2.91281,3.03588,3.14546, &
     & 3.25692,3.35719,3.45851,3.55051,3.64283,3.72734,    627.2727    /
!     Ac II   AEL-Kurucz     
      DATA (UARRAY(89,2,K),K=0,22) /                                    &
     & 0.00000,0.01157,0.31023,0.82812,1.31453,1.64035,1.94641,2.18042, &
     & 2.40153,2.58249,2.75874,2.90854,3.05542,3.18346,3.30835,3.41936, &
     & 3.52695,3.62407,3.71722,3.80243,3.88373,3.95891,   1100.0000    /
!     Ac III  AEL-Kurucz     
      DATA (UARRAY(89,3,K),K=0,22) /                                    &
     & 0.69315,1.46732,1.80028,1.97951,2.08888,2.18751,2.26504,2.33700, &
     & 2.40243,2.46384,2.51890,2.57108,2.61849,2.66376,2.70438,2.74341, &
     & 2.77851,2.81241,2.84346,2.87357,2.90087,2.92745,   1818.1818    /
!     Ac IV   U=1 assumed
      DATA (UARRAY(89,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   4454.5454    /
!     Ac V    U=1 assumed
      DATA (UARRAY(89,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5636.3638    /
!     Th I    AEL-Kurucz     
      DATA (UARRAY(90,1,K),K=0,22) /                                    &
     & 1.60944,1.60944,1.60944,1.82475,2.08730,2.29512,2.57248,2.78942, &
     & 3.04553,3.24925,3.47987,3.66715,3.87337,4.04425,4.22803,4.38323, &
     & 4.54536,4.68484,4.82864,4.95434,5.08132,5.19398,    545.4545    /
!     Th II   FAK-Kurucz     
      DATA (UARRAY(90,2,K),K=0,22) /                                    &
     & 1.38629,1.38629,1.43092,1.78621,2.04783,2.25500,2.50987,2.73081, &
     & 2.96723,3.17762,3.39015,3.58564,3.77556,3.95533,4.12496,4.28909, &
     & 4.44026,4.58817,4.72290,4.85628,4.97677,5.09622,   1045.4546    /
!     Th III  FAK-Kurucz     
      DATA (UARRAY(90,3,K),K=0,22) /                                    &
     & 1.60944,1.60944,1.60944,1.82447,2.08719,2.29505,2.57243,2.78937, &
     & 3.04547,3.24917,3.47983,3.66715,3.87340,4.04432,4.22808,4.38325, &
     & 4.54541,4.68491,4.82873,4.95445,5.08141,5.19406,   1818.1818    /
!     Th IV   U=1 assumed
      DATA (UARRAY(90,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   2618.1819    /
!     Th V    U=1 assumed
      DATA (UARRAY(90,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5909.0908    /
!     Pa I    AEL-Kurucz     
      DATA (UARRAY(91,1,K),K=0,22) /                                    &
     & 2.48491,2.48491,2.57955,3.54946,4.03243,4.35675,4.64107,4.86220, &
     & 5.05150,5.21060,5.35636,5.48355,5.60946,5.72128,5.83741,5.94145, &
     & 6.05099,6.14970,6.25271,6.34610,6.44247,6.53037,    545.4545    /
!     Pa II   FAK-Kurucz     
      DATA (UARRAY(91,2,K),K=0,22) /                                    &
     & 2.19722,2.19722,2.57805,3.54889,4.03208,4.35647,4.64096,4.86220, &
     & 5.05150,5.21060,5.35636,5.48355,5.60946,5.72128,5.83741,5.94145, &
     & 6.05103,6.14979,6.25280,6.34619,6.44256,6.53045,   1090.9091    /
!     Pa III  FAK-Kurucz     
      DATA (UARRAY(91,3,K),K=0,22) /                                    &
     & 0.00000,0.00000,2.57661,3.54860,4.03206,4.35657,4.64098,4.86217, &
     & 5.05145,5.21055,5.35634,5.48355,5.60948,5.72132,5.83744,5.94147, &
     & 6.05102,6.14976,6.25279,6.34619,6.44255,6.53043,   1818.1818    /
!     Pa IV   U=1 assumed
      DATA (UARRAY(91,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   2636.3635    /
!     Pa V    U=1 assumed
      DATA (UARRAY(91,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5454.5454    /
!     U  I    AEL-Kurucz     
      DATA (UARRAY(92,1,K),K=0,22) /                                    &
     & 2.56495,2.56495,2.79627,3.04228,3.23957,3.40432,3.62669,3.80852, &
     & 4.02397,4.20114,4.40282,4.57058,4.75702,4.91410,5.08354,5.22839, &
     & 5.38007,5.51175,5.64608,5.76449,5.88282,5.98863,    545.4545    /
!     U  II   FAK-Kurucz     
      DATA (UARRAY(92,2,K),K=0,22) /                                    &
     & 2.30259,2.46891,2.79601,3.04208,3.23940,3.40413,3.62660,3.80852, &
     & 4.02397,4.20114,4.40282,4.57058,4.75702,4.91410,5.08354,5.22839, &
     & 5.38013,5.51187,5.64620,5.76460,5.88292,5.98872,   1090.9091    /
!     U  III  FAK-Kurucz     
      DATA (UARRAY(92,3,K),K=0,22) /                                    &
     & 2.56495,2.56495,2.79596,3.04208,3.23944,3.40419,3.62662,3.80847, &
     & 4.02392,4.20107,4.40279,4.57058,4.75705,4.91416,5.08358,5.22841, &
     & 5.38012,5.51182,5.64617,5.76460,5.88291,5.98870,   1818.1818    /
!     U  IV   U=1 assumed
      DATA (UARRAY(92,4,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   2636.3635    /
!     U  V    U=1 assumed
      DATA (UARRAY(92,5,K),K=0,22) /                                    &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,   5454.5454    /
!      end 'partfn.inc'
!-----------------------------------------------------------------------
!
!     uarray - contains partition function data given as log(u)
!              0 is ground state value
!              22 contains temperature step size = ip*2000./22.
!
!
!     array of g values
!
!      real g(92,5)
!-----------------------------------------------------------------------
!      include 'gvalues.inc'
!      begin 'gvalues.inc'
!     h 
      data (g( 1,k),k=1,5) /   2.0,   0.0,   0.0,   0.0,   0.0 /
!     he
      data (g( 2,k),k=1,5) /   4.0,   2.0,   0.0,   0.0,   0.0 /
!     li
      data (g( 3,k),k=1,5) /   2.0,   4.0,   2.0,   0.0,   0.0 /
!     be
      data (g( 4,k),k=1,5) /   4.0,   2.0,   4.0,   2.0,   0.0 /
!     b 
      data (g( 5,k),k=1,5) /   2.0,   4.0,   2.0,   4.0,   2.0 /
!     c 
      data (g( 6,k),k=1,5) /  12.0,   2.0,   4.0,   2.0,   4.0 /
!     n 
      data (g( 7,k),k=1,5) /  18.0,  12.0,   2.0,   4.0,   2.0 /
!     o 
      data (g( 8,k),k=1,5) /   8.0,  18.0,  12.0,   2.0,   4.0 /
!     f 
      data (g( 9,k),k=1,5) /  18.0,   8.0,  18.0,  12.0,   2.0 /
!     ne
      data (g(10,k),k=1,5) /  12.0,  18.0,   8.0,  18.0,  12.0 /
!     na
      data (g(11,k),k=1,5) /   2.0,  12.0,  18.0,   8.0,  18.0 /
!     mg
      data (g(12,k),k=1,5) /   4.0,   2.0,  12.0,  18.0,   8.0 /
!     al
      data (g(13,k),k=1,5) /   2.0,   4.0,   2.0,  12.0,  18.0 /
!     si
      data (g(14,k),k=1,5) /  12.0,   2.0,   4.0,   2.0,  12.0 /
!     p 
      data (g(15,k),k=1,5) /  18.0,  12.0,   2.0,   4.0,   2.0 /
!     s 
      data (g(16,k),k=1,5) /   8.0,  18.0,  12.0,   2.0,   4.0 /
!     cl
      data (g(17,k),k=1,5) /  18.0,   8.0,  18.0,  12.0,   2.0 /
!     ar
      data (g(18,k),k=1,5) /  12.0,  18.0,   8.0,  18.0,  12.0 /
!     k 
      data (g(19,k),k=1,5) /   2.0,  12.0,  18.0,   8.0,  18.0 /
!     ca
      data (g(20,k),k=1,5) /   4.0,   2.0,  12.0,  18.0,   8.0 /
!     sc
      data (g(21,k),k=1,5) /  30.0,  20.0,   2.0,  12.0,  18.0 /
!     ti
      data (g(22,k),k=1,5) /  56.0,  42.0,  20.0,   2.0,  12.0 /
!     v 
      data (g(23,k),k=1,5) /  50.0,  56.0,  42.0,  20.0,   2.0 /
!     cr
      data (g(24,k),k=1,5) /  12.0,  50.0,  56.0,  42.0,  20.0 /
!     mn
      data (g(25,k),k=1,5) /  14.0,  12.0,  50.0,  56.0,  42.0 /
!     fe
      data (g(26,k),k=1,5) /  60.0,  50.0,  12.0,  50.0,  56.0 /
!     co
      data (g(27,k),k=1,5) /  42.0,  56.0,  50.0,  12.0,  50.0 /
!     ni
      data (g(28,k),k=1,5) /  20.0,  42.0,  56.0,  50.0,  12.0 /
!     cu
      data (g(29,k),k=1,5) /   2.0,  20.0,  42.0,  56.0,  50.0 /
!     zn
      data (g(30,k),k=1,5) /   4.0,   2.0,  20.0,   0.0,   0.0 /
!     ga
      data (g(31,k),k=1,5) /   2.0,   4.0,   2.0,  20.0,   0.0 /
!     ge
      data (g(32,k),k=1,5) /  12.0,   2.0,   4.0,   2.0,  20.0 /
!     as
      data (g(33,k),k=1,5) /  18.0,  12.0,   2.0,   4.0,   2.0 /
!     se
      data (g(34,k),k=1,5) /   8.0,  18.0,  12.0,   2.0,   4.0 /
!     br
      data (g(35,k),k=1,5) /  18.0,   8.0,  18.0,  12.0,   2.0 /
!     kr
      data (g(36,k),k=1,5) /  12.0,  18.0,   8.0,  18.0,  12.0 /
!     rb
      data (g(37,k),k=1,5) /   2.0,  12.0,  18.0,   0.0,   0.0 /
!     sr
      data (g(38,k),k=1,5) /   4.0,   2.0,  12.0,  18.0,   0.0 /
!     y 
      data (g(39,k),k=1,5) /   2.0,  20.0,   2.0,  12.0,  18.0 /
!     zr
      data (g(40,k),k=1,5) /  56.0,  42.0,  20.0,   2.0,  12.0 /
!     cb
      data (g(41,k),k=1,5) /  50.0,  56.0,  42.0,  20.0,   2.0 /
!     mo
      data (g(42,k),k=1,5) /  12.0,  50.0,  56.0,  42.0,  20.0 /
!     tc
      data (g(43,k),k=1,5) /  14.0,  12.0,   0.0,   0.0,   0.0 /
!     ru
      data (g(44,k),k=1,5) /  56.0,  50.0,  12.0,   0.0,   0.0 /
!     rh
      data (g(45,k),k=1,5) /  42.0,  56.0,  50.0,   0.0,   0.0 /
!     pd
      data (g(46,k),k=1,5) /  20.0,  42.0,  56.0,   0.0,   0.0 /
!     ag
      data (g(47,k),k=1,5) /   2.0,  20.0,  42.0,   0.0,   0.0 /
!     cd
      data (g(48,k),k=1,5) /   4.0,   2.0,  20.0,   0.0,   0.0 /
!     in
      data (g(49,k),k=1,5) /   2.0,   4.0,   2.0,  20.0,   0.0 /
!     sn
      data (g(50,k),k=1,5) /  12.0,   2.0,   4.0,   2.0,  20.0 /
!     sb
      data (g(51,k),k=1,5) /  18.0,  12.0,   2.0,   4.0,   2.0 /
!     te
      data (g(52,k),k=1,5) /   8.0,  18.0,  12.0,   2.0,   4.0 /
!     i 
      data (g(53,k),k=1,5) /  18.0,   8.0,  18.0,   0.0,   2.0 /
!     xe
      data (g(54,k),k=1,5) /  12.0,  18.0,   8.0,   0.0,   0.0 /
!     cs
      data (g(55,k),k=1,5) /   2.0,  12.0,   0.0,   0.0,   0.0 /
!     ba
      data (g(56,k),k=1,5) /   4.0,   2.0,  12.0,   0.0,   0.0 /
!     la
      data (g(57,k),k=1,5) /  42.0,  20.0,   2.0,  12.0,  18.0 /
!     ce
      data (g(58,k),k=1,5) /  88.0,  66.0,  28.0,   2.0,  12.0 /
!     pr
      data (g(59,k),k=1,5) / 130.0, 104.0,  66.0,  28.0,   2.0 /
!     nd
      data (g(60,k),k=1,5) / 156.0, 130.0, 104.0,  66.0,   0.0 /
!     pm
      data (g(61,k),k=1,5) / 154.0, 132.0, 130.0, 104.0,   0.0 /
!     sm
      data (g(62,k),k=1,5) / 112.0,  98.0, 132.0, 130.0,   0.0 /
!     eu
      data (g(63,k),k=1,5) /  90.0,  16.0,  98.0, 132.0,   0.0 /
!     gd
      data (g(64,k),k=1,5) / 100.0,  90.0,  16.0,  98.0,   0.0 /
!     tb
      data (g(65,k),k=1,5) / 154.0, 132.0,  98.0,  16.0,   0.0 /
!     dy
      data (g(66,k),k=1,5) / 156.0, 130.0, 132.0,  98.0,   0.0 /
!     ho
      data (g(67,k),k=1,5) / 130.0, 104.0, 130.0, 132.0,   0.0 /
!     er
      data (g(68,k),k=1,5) /  88.0,  66.0, 104.0, 130.0,   0.0 /
!     tm
      data (g(69,k),k=1,5) /  42.0,  28.0,  66.0, 104.0,   0.0 /
!     yb
      data (g(70,k),k=1,5) /   4.0,   2.0,  28.0,  66.0,   0.0 /
!     lu
      data (g(71,k),k=1,5) /   2.0,   4.0,   2.0,  28.0,  66.0 /
!     hf
      data (g(72,k),k=1,5) /  20.0,  42.0,   0.0,   0.0,   0.0 /
!     ta
      data (g(73,k),k=1,5) /  70.0,   0.0,   0.0,   0.0,   0.0 /
!     w 
      data (g(74,k),k=1,5) /  60.0,   0.0,   0.0,   0.0,   0.0 /
!     re
      data (g(75,k),k=1,5) /  14.0,   0.0,   0.0,   0.0,   0.0 /
!     os
      data (g(76,k),k=1,5) /  60.0,   0.0,   0.0,   0.0,   0.0 /
!     ir
      data (g(77,k),k=1,5) /  42.0,   0.0,   0.0,   0.0,   0.0 /
!     pt
      data (g(78,k),k=1,5) /  20.0,  42.0,   0.0,   0.0,   0.0 /
!     au
      data (g(79,k),k=1,5) /   2.0,  20.0,   0.0,   0.0,   0.0 /
!     hg
      data (g(80,k),k=1,5) /   4.0,   2.0,  20.0,   0.0,   0.0 /
!     tl
      data (g(81,k),k=1,5) /   2.0,   4.0,   2.0,  20.0,   0.0 /
!     pb
      data (g(82,k),k=1,5) /  12.0,   2.0,   4.0,   2.0,  20.0 /
!     bi
      data (g(83,k),k=1,5) /  18.0,  12.0,   2.0,   4.0,   2.0 /
!     po
      data (g(84,k),k=1,5) /   8.0,   0.0,   0.0,   0.0,   0.0 /
!     at
      data (g(85,k),k=1,5) /   0.0,   0.0,   0.0,   0.0,   0.0 /
!     rn
      data (g(86,k),k=1,5) /  12.0,   0.0,   0.0,   0.0,   0.0 /
!     fr
      data (g(87,k),k=1,5) /   0.0,   0.0,   0.0,   0.0,   0.0 /
!     ra
      data (g(88,k),k=1,5) /   4.0,   2.0,   0.0,   0.0,   0.0 /
!     ac
      data (g(89,k),k=1,5) /   2.0,   4.0,   0.0,   0.0,   0.0 /
!     th
      data (g(90,k),k=1,5) /  56.0,   0.0,   0.0,   0.0,   0.0 /
!     pa
      data (g(91,k),k=1,5) /  66.0,   0.0,   0.0,   0.0,   0.0 /
!     u 
      data (g(92,k),k=1,5) / 104.0,   0.0,   0.0,   0.0,   0.0 /
!      end 'gvalues.inc'
!-----------------------------------------------------------------------
!
!     array of cutoff values
!
!      real dele(92,5)
!-----------------------------------------------------------------------
!      include 'cutoff.inc'
!      begin 'cutoff.inc'
!     h 
      data (dele( 1,k),k=1,5) /-0.087 , 0.    , 0.    , 0.    , 0.    /
!     he
      data (dele( 2,k),k=1,5) / 0.10  , 0.10  , 0.    , 0.    , 0.    /
!     li
      data (dele( 3,k),k=1,5) /-0.1233, 0.10  , 0.10  , 0.    , 0.    /
!     be
      data (dele( 4,k),k=1,5) /-0.242 ,-0.968 , 0.    , 0.    , 0.    /
!     b 
      data (dele( 5,k),k=1,5) /-0.1233,-1.287 , 0.    , 0.    , 0.    /
!     c 
      data (dele( 6,k),k=1,5) /-0.1233,-0.9668, 0.10  , 0.10  , 0.    /
!     n 
      data (dele( 7,k),k=1,5) /-0.322 ,-0.1   , 0.    , 0.10  , 0.10  /
!     o 
      data (dele( 8,k),k=1,5) /-0.242 ,-1.798 , 0.    , 0.10  , 0.10  /
!     f 
      data (dele( 9,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.10  /
!     ne
      data (dele(10,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.    /
!     na
      data (dele(11,k),k=1,5) /-0.1882, 0.10  , 0.10  , 0.10  , 0.    /
!     mg
      data (dele(12,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     al
      data (dele(13,k),k=1,5) /-0.1882, 0.10  , 0.10  , 0.10  , 0.    /
!     si
      data (dele(14,k),k=1,5) /-0.242 ,-0.603 , 0.    ,-0.1   , 0.    /
!     p 
      data (dele(15,k),k=1,5) /-0.322 ,-1.79  , 0.    , 0.    , 0.    /
!     s 
      data (dele(16,k),k=1,5) /-0.1   ,-0.1   ,-0.1   , 0.10  , 0.10  /
!     cl
      data (dele(17,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.10  /
!     ar
      data (dele(18,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.10  /
!     k 
      data (dele(19,k),k=1,5) /-0.1233, 0.10  , 0.10  , 0.10  , 0.10  /
!     ca
      data (dele(20,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.10  /
!     sc
      data (dele(21,k),k=1,5) /-0.45  ,-1.287 , 0.    , 0.    , 0.    /
!     ti
      data (dele(22,k),k=1,5) / 0.10  ,-2.0   , 0.    , 0.10  , 0.    /
!     v 
      data (dele(23,k),k=1,5) /-0.50  ,-4.0   , 0.    , 0.    , 0.    /
!     cr
      data (dele(24,k),k=1,5) /-0.1   ,-3.0   , 0.    , 0.10  , 0.10  /
!     mn
      data (dele(25,k),k=1,5) /-0.1   ,-0.75  , 0.    , 0.    , 0.    /
!     fe
      data (dele(26,k),k=1,5) /-1.0   ,-3.0   , 0.10  , 0.10  , 0.10  /
!     co
      data (dele(27,k),k=1,5) / 0.10  ,-2.0   , 0.10  , 0.    , 0.    /
!     ni
      data (dele(28,k),k=1,5) / 0.10  ,-1.0   , 0.    , 0.    , 0.    /
!     cu
      data (dele(29,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     zn
      data (dele(30,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     ga
      data (dele(31,k),k=1,5) /-0.1234,-1.288 , 0.    , 0.    , 0.    /
!     ge
      data (dele(32,k),k=1,5) /-0.242 ,-0.753 , 0.    , 0.    , 0.    /
!     as
      data (dele(33,k),k=1,5) /-0.322 ,-1.287 , 0.    , 0.    , 0.    /
!     se
      data (dele(34,k),k=1,5) /-0.322 ,-1.287 , 0.    , 0.    , 0.    /
!     br
      data (dele(35,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     kr
      data (dele(36,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     rb
      data (dele(37,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     sr
      data (dele(38,k),k=1,5) /-0.103 ,-0.493 , 0.    , 0.10  , 0.    /
!     y 
      data (dele(39,k),k=1,5) /-0.50  ,-1.798 , 0.    , 0.    , 0.    /
!     zr
      data (dele(40,k),k=1,5) / 0.10  , 0.    , 0.10  , 0.10  , 0.    /
!     cb
      data (dele(41,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     mo
      data (dele(42,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.10  /
!     tc
      data (dele(43,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ru
      data (dele(44,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     rh
      data (dele(45,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     pd
      data (dele(46,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ag
      data (dele(47,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     cd
      data (dele(48,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     in
      data (dele(49,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.    /
!     sn
      data (dele(50,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     sb
      data (dele(51,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     te
      data (dele(52,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     i 
      data (dele(53,k),k=1,5) / 0.10  , 0.10  , 0.    , 0.    , 0.    /
!     xe
      data (dele(54,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.    /
!     cs
      data (dele(55,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.    , 0.    /
!     ba
      data (dele(56,k),k=1,5) / 0.10  , 0.10  , 0.    , 0.10  , 0.    /
!     la
      data (dele(57,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ce
      data (dele(58,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     pr
      data (dele(59,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     nd
      data (dele(60,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     pm
      data (dele(61,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     sm
      data (dele(62,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     eu
      data (dele(63,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     gd
      data (dele(64,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     tb
      data (dele(65,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     dy
      data (dele(66,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ho
      data (dele(67,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     er
      data (dele(68,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     tm
      data (dele(69,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     yb
      data (dele(70,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     lu
      data (dele(71,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     hf
      data (dele(72,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ta
      data (dele(73,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     w 
      data (dele(74,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     re
      data (dele(75,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     os
      data (dele(76,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ir
      data (dele(77,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     pt
      data (dele(78,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     au
      data (dele(79,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     hg
      data (dele(80,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.    /
!     tl
      data (dele(81,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     pb
      data (dele(82,k),k=1,5) / 0.10  , 0.10  , 0.10  , 0.10  , 0.10  /
!     bi
      data (dele(83,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     po
      data (dele(84,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     at
      data (dele(85,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     rn
      data (dele(86,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     fr
      data (dele(87,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ra
      data (dele(88,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     ac
      data (dele(89,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     th
      data (dele(90,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     pa
      data (dele(91,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!     u 
      data (dele(92,k),k=1,5) / 0.    , 0.    , 0.    , 0.    , 0.    /
!      end 'cutoff.inc'
!-----------------------------------------------------------------------
!
!	include 'atmdat.f'
!
!     first ionisation potentials
!
      data aipot1/                                                      &
     &    13.595,   24.581,    5.390,    9.320,    8.296,   11.256,     &
     &    14.530,   13.614,   17.418,   21.559,    5.138,    7.644,     &
     &     5.984,    8.149,   10.484,   10.357,   13.010,   15.755,     &
     &     4.339,    6.111,    6.540,    6.820,    6.740,    6.764,     &
     &     7.432,    7.870,    7.860,    7.633,    7.724,    9.391,     &
     &     6.000,    7.880,    9.810,    9.750,   11.840,   13.996,     &
     &     4.176,    5.692,    6.380,    6.840,    6.880,    7.100,     &
     &     7.280,    7.364,    7.460,    8.330,    7.574,    8.991,     &
     &     5.785,    7.342,    8.639,    9.010,   10.454,   12.127,     &
     &     3.893,    5.210,    5.577,    5.539,    5.464,    5.525,     &
     &     5.554,    5.644,    5.671,    6.150,    5.864,    5.939,     &
     &     6.022,    6.108,    6.184,    6.254,    5.426,    7.000,     &
     &     7.880,    7.980,    7.870,    8.700,    9.000,    9.000,     &
     &     9.220,   10.430,    6.106,    7.415,    7.287,    8.430,     &
     &     9.300,   10.746,    4.000,    5.277,    6.9  ,    6.   ,     &
     &     6.   ,    6.                                            /
!
!     second ionisation potentials
!
      data aipot2/                                                      &
     &     0.   ,   54.403,   75.619,   18.206,   25.149,   24.376,     &
     &    29.593,   35.108,   34.980,   41.070,   47.290,   15.031,     &
     &    18.823,   16.340,   19.720,   23.400,   23.800,   27.620,     &
     &    31.810,   11.868,   12.800,   13.570,   14.650,   16.490,     &
     &    15.636,   16.180,   17.050,   18.150,   20.290,   17.960,     &
     &    20.510,   15.930,   18.630,   21.500,   21.600,   24.560,     &
     &    27.500,   11.027,   12.230,   13.130,   14.320,   16.150,     &
     &    15.260,   16.760,   18.070,   19.420,   21.480,   16.904,     &
     &    18.860,   14.628,   16.500,   18.600,   19.090,   21.200,     &
     &    25.100,   10.001,   11.060,   10.85 ,   10.55 ,   10.73 ,     &
     &    10.90 ,   11.07 ,   11.241,   12.09 ,   11.52 ,   11.67 ,     &
     &    11.80 ,   11.93 ,   12.05 ,   12.184,   13.9  ,   14.900,     &
     &    16.200,   17.700,   16.600,   17.000,   17.000,   18.560,     &
     &    20.500,   18.751,   20.420,   15.028,   16.680,   19.000,     &
     &    20.000,   21.000,   22.000,   10.144,   12.1  ,   11.5  ,     &
     &    12.   ,   12.                                            /
!
!     third ionisation potentials
!
      data aipot3/                                                      &
     &     0.   ,    0.   ,  122.451,  153.893,   37.920,   47.871,     &
     &    47.426,   54.886,   62.646,   63.500,   71.650,   80.120,     &
     &    28.440,   33.460,   30.156,   35.000,   39.900,   40.900,     &
     &    46.000,   51.210,   24.750,   27.470,   29.310,   30.950,     &
     &    33.690,   30.643,   33.490,   35.160,   36.830,   39.700,     &
     &    30.700,   34.210,   28.340,   32.000,   35.900,   36.900,     &
     &    40.000,   43.000,   20.500,   22.980,   25.040,   27.130,     &
     &    31.000,   28.460,   31.050,   32.920,   34.820,   37.470,     &
     &    28.030,   30.490,   25.300,   31.000,   32.000,   32.100,     &
     &    35.000,   36.000,   19.177,   20.198,   21.624,   22.1  ,     &
     &    22.3  ,   23.4  ,   24.92 ,   20.63 ,   21.91 ,   22.8  ,     &
     &    22.84 ,   22.74 ,   23.68 ,   25.05 ,   20.960,   21.000,     &
     &    22.000,   24.000,   26.000,   25.000,   27.000,   28.000,     &
     &    30.000,   34.200,   29.800,   31.930,   25.560,   27.000,     &
     &    29.000,   29.000,   33.000,   34.000,   20.   ,   20.0  ,     &
     &    20.   ,   20.                                            /
!
!     fourth ionization potentials
!
      data aipot4/                                                      &
     &     0.   ,    0.   ,    0.   ,  217.713,  259.366,   64.492,     &
     &    77.472,   77.413,   87.138,   97.11 ,   98.91 ,  109.31 ,     &
     &   119.99 ,   45.141,   51.42 ,   47.30 ,   53.46 ,  59.81  ,     &
     &    60.92 ,   67.15 ,   73.7  ,   43.26 ,   46.71 ,  49.1   ,     &
     &    51.4  ,   54.8  ,   51.3  ,   54.9  ,   55.2  ,   59.4  ,     &
     &    64.   ,   45.71 ,   50.13 ,   42.944,   47.3  ,   52.5  ,     &
     &    52.6  ,   57.   ,   61.8  ,   34.34 ,   38.3  ,   46.4  ,     &
     &    46.   ,   50.   ,   48.   ,   53.   ,   56.   ,   59.   ,     &
     &    54.4  ,   40.734,   44.2  ,   37.41 ,   42.   ,   46.   ,     &
     &    46.   ,   49.   ,   49.95 ,   36.758,   38.98 ,   40.4  ,     &
     &    41.1  ,   41.4  ,   42.7  ,   44.0  ,   39.37 ,   41.4  ,     &
     &    42.5  ,   42.7  ,   42.7  ,   43.56 ,   45.25 ,   33.3  ,     &
     &    33.   ,   35.   ,   38.   ,   40.   ,   39.   ,   41.   ,     &
     &    44.   ,   46.   ,   50.7  ,   42.32 ,   45.3  ,   38.   ,     &
     &    41.   ,   44.   ,   43.   ,   46.   ,   49.   ,   28.8  ,     &
     &    29.   ,   29.                                            /
!
!     fifth ionization potentials
!
      data aipot5/                                                      &
     &     0.   ,    0.   ,    0.   ,    0.   ,  340.22 ,  392.08 ,     &
     &    97.89 ,  113.90 ,  114.24 ,  126.21 ,  138.40 ,  141.27 ,     &
     &   153.75 ,  166.77 ,   65.02 ,   72.68 ,   67.7  ,   75.04 ,     &
     &    82.66 ,   84.43 ,   91.7  ,   99.4  ,   65.23 ,   70.2  ,     &
     &    73.0  ,   75.5  ,   79.5  ,   75.5  ,   79.9  ,   82.6  ,     &
     &    87.   ,   93.5  ,   62.63 ,   68.3  ,   59.7  ,   64.7  ,     &
     &    71.0  ,   71.6  ,   77.0  ,   81.5  ,   50.55 ,   61.2  ,     &
     &    55.   ,   60.   ,   65.   ,   62.   ,   68.   ,   72.   ,     &
     &    77.   ,   72.28 ,   56.   ,   58.75 ,   66.   ,   57.   ,     &
     &    62.   ,   62.   ,   61.6  ,   65.55 ,   57.53 ,   60.   ,     &
     &    60.   ,   60.   ,   60.   ,   60.   ,   60.   ,   60.   ,     &
     &    60.   ,   60.   ,   60.   ,   60.   ,   66.8  ,   50.   ,     &
     &    45.   ,   48.   ,   51.   ,   54.   ,   57.   ,   55.   ,     &
     &    58.   ,   61.   ,   64.   ,   68.8  ,   56.0  ,   61.   ,     &
     &    51.   ,   55.   ,   59.   ,   58.   ,   62.   ,   65.   ,     &
     &    60.   ,   60.                                            /
!-----------------------------------------------------------------------
!     interpolate partition function for theta
!
!     uses linear interpolation
!
      do j=1,5
!
!     locate position in table and adjust for bounds if necessary
!
!        i=int(5040./(theta*uarray(idno,j,22)))
!        i=max(0,min(i,20))
!
        i=max(0,min(int(5040./(theta*uarray(idno,j,22))),20))
!
        u(j)=exp((5040./theta/uarray(idno,j,22)-real(i))                &
     &       *(uarray(idno,j,i+1)-uarray(idno,j,i))                     &
     &       +uarray(idno,j,i))
      enddo
!
!     add braket function if required
!
      if (dele(idno,1) .gt. 0.0001) then
        u(1) = u(1) + g(idno,1)*exp(-aipot1(idno)*alge10*theta)         &
     &                *braket(1.d0,eplow,dele(idno,1),theta)
      endif
!
      if (dele(idno,2) .gt. 0.0001) then
        u(2) = u(2) + g(idno,2)*exp(-aipot2(idno)*alge10*theta)         &
     &                *braket(2.d0,2*eplow,dele(idno,2),theta)
      endif
!
      if (dele(idno,3) .gt. 0.0001) then
        u(3) = u(3) + g(idno,3)*exp(-aipot3(idno)*alge10*theta)         &
     &                *braket(3.d0,3*eplow,dele(idno,3),theta)
      endif
!
      if (dele(idno,4) .gt. 0.0001) then
        u(4) = u(4) + g(idno,4)*exp(-aipot4(idno)*alge10*theta)         &
     &                *braket(4.d0,4*eplow,dele(idno,4),theta)
      endif
!
      if (dele(idno,5) .gt. 0.0001) then
        u(5) = u(5) + g(idno,5)*exp(-aipot5(idno)*alge10*theta)         &
     &                *braket(5.d0,5*eplow,dele(idno,5),theta)
      endif
!
      return
      end
!
!***********************************************************************
      function braket (a,b,c,d)
        IMPLICIT REAL*8 (A-H, O-Z)
!
!     this function is used as an aid to calculation of the
!     hydrogenic portion of the partition function. use values of
!
!     a = z (z=1.0 for neutral, z=2.0 for once ionized, etc.)
!     b = eplow (in ev) *z
!     c = dele of explicit part of partition function nmax+0.5, in ev
!     d = theta
!
!     note that c=13.595*z*z/(nmax**2)
!
!     see kurucz (1970) smithsonian ap obs sp. rep. 309, pp 63-65
!
!-----------------------------------------------------------------------
!      include 'params.inc'
!      begin 'params.inc'
!
!i    melmln         - size of element abreviations
!
      integer    melmln
      parameter (melmln =     2)
!
!i    melmns         - number of elements
!
      integer    melmns
      parameter (melmns =    92)
!
!i    mlayer         - maximum number of model depths
!
      integer    mlayer
      parameter (mlayer =   200)
!
!i    mlines         - maximum number of lines in line data buffer
!
      integer    mlines
      parameter (mlines =  1000)
!
!i    msynth         - size of "synth" arrays
!
      integer    msynth
      parameter (msynth = 10000)
!
!i    mbroad         - size of "broad" arrays
!
      integer    mbroad
      parameter (mbroad = 10000)
!
!i    mtitle         - size of character strings for titles and filenames
!
      integer    mtitle
      parameter (mtitle =    80)
!
!i    mdipso         - size of "dipso" arrays
!
      integer    mdipso
      parameter (mdipso = 100000)
!
!i    mfpram         - maximum number of free parameters for chifit
!
      integer    mfpram
      parameter (mfpram =     6)
!
!i    mchmax         -
!
      integer    mchmax
      parameter (mchmax =    50)
!
!i    mdskin         - file input unit
!
      integer    mdskin
      parameter (mdskin =     2)
!
!i    mdskot         - file output unit
!
      integer    mdskot
      parameter (mdskot =     3)
!
!i    mchlin         - default input unit
!
      integer    mchlin
      parameter (mchlin =     5)
!
!i    mchlot         - default output unit
!
      integer    mchlot
      parameter (mchlot =     6)
!
!i    mchses         - unit for session log file
!
      integer    mchses
      parameter (mchses =    10)
!
!r    alge10         - value of log(10)
!
      real alge10
      parameter (alge10 = 2.302585093)
!
!r    around         - value by which to round-up errors
!
      real around
      parameter (around = 0.01)
!
!i    minstp         - maximum number of points in instrumental profile
!
      integer minstp
      parameter (minstp = 1000)
!      end 'params.inc'
!-----------------------------------------------------------------------
!
      x=(13.595*a*a/b)**1.5
      y=(13.595*a*a/c)**1.5
      v=b*alge10*d
      w=c*alge10*d
      braket=x*(0.3333333333+v-0.5*v*v-0.05555555556*v*v*v)             &
     &      -y*(0.3333333333+w-0.5*w*w-0.05555555556*w*w*w)
      return
      end