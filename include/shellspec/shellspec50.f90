!
	program shellrun
	call shellspec()
	call exit(0)
	end program shellrun
!
	subroutine shellspec()
!	program shellspec
!
!	is designed to calculate lightcurves, spectra and images of 
!	interacting binaries and extrasolar planets immersed in a moving
!       circumstellar matter (CM). 
!       It solves simple radiative transfer along the line of sight 
!       in 3D moving media. Roche model can be used as a boundary 
!       condition for the radiative tranfer. The scattered light from 
!       the two stars can be taken into account assuming that CM is 
!       optically thin. Shadows can be taken into account as well.
!	The assumptions include LTE and optional known state quantities
!       and velocity field in 3D. These can be taken from 
!       the 3D hydrodynamical simulations. Alternatively,
!	optional (non)transparent objects such as:
!	a central STAR, COMPANION star, SPOT, STREAM, RING,
!       DISK, ENVELOPE, NEBULA, FLOW, JET, UFO, SHELL or an empty space 
!       may be defined in 3D and their composite synthetic spectrum 
!       calculated. The stars may have either the Roche or spherical 
!       geometry, optional velocity or rotation. They are subject
!       to the gravity darkening, limb darkening, and irradiation effect
!	including the heating, reflection and day-night heat 
!       redistribution. They may be ascribed a precalculated spectrum.
!	Synthetic light curves or trailing spectrograms can be produced
!	by changing your view points on the 3D object.
!	Opacities: lines, HI (bound-free, free-free), 
!               H- (bound-free, free-free),
!		Thomson scattering, Rayleigh scattering, 
!               Mie absorption and scattering on dust. 
!	Emissivities: thermal from the gas, 
!               Thomson and Rayleigh scattering from the stars, 
!               thermal from the dust and Mie scattering from stars 
!               on dust
!
!       INPUT:
!		shellspec.in -(9) main input (geometry, objects...)
!		line.dat - (8) atomic data for the sp. lines
!                       (optional if iline=1)
!               shellspec.mod - (10) input 3D model of the shell
!			(optional if imodel=2)
!		abundances  - (7) abundances 
!			(optional if ichemc=1 or ielnd=1)
!		lambda - (15) list of lambda [A] to calculate
!			(optional if loglam=2)
!		phases - (15) list of orbital phases to calculate
!			(optional if nphase=0)
!		starspec1 - (12) spectrum of the primary star
!			(optional if lunt1>0)
!		starspec2 - (13) spectrum of the secondary
!			(optional if lunt2>0)
!		starspec3 - (14) spectrum of the third body 
!			(optional if lunt3>0)
!		albedo1 - (12) monochromatic albedo of the primary star
!			(optional if ialbst=1 and irrst=1)
!		albedo2 - (13) monochromatic albedo of the secondary
!			(optional if ialbcp=1 and irrcp=1)
!               dust_opac - (12) dust opacities (optional if imie>0)
!		mie_phase - (13) dust phase func. (optional if imiepf=1)
!               gas_opac - (16) molec x-sections (optional if iopac=1)
!               chem_eq_tab - (16) molec popul tab (optional if iopac=1)
!               wind_prof -(17) vertic. density wind prof (if idennb=1)
!
!       OUTPUT:
!		shellspec.out - (2) more detailed output
!		2Dimage_xxx - (100+iang) 2D images at some frequency
!		shellspectrum - (4) spectrum of the shell
!		lightcurve - (11) light curve or trailed spectrogram
!		errors - (3) error messages
!
!	COMMENT: CGS units are used if not specified otherwise
!		 Nominal values of solar mass and radius are used
!		 Cartesian 3D coordinates are used
!		 real(kind=8) i.e. double precision is used and 
!		 use e.g. 1.d40 and not 1.e40
!		 ?? - BE CAUTIOUS
!
!	Non original routines: 
!	 	pfdwor  -from UCLSYN	
!		voigt, state0, gaunt, gfree   -from SYNSPEC
!		hunt,locate    -from Numerical recipes
!
!	Author/contact: Jan Budaj
!			http://www.ta3.sk/~budaj
!			budaj@ta3.sk
!
!	Reference:
!		Budaj J., Richards M.T., 2004, 
!		Contrib. Astron. Obs. Skalnate Pleso, 34, 167
!
!	Things to be improved - added:
!	-introduce sets with continuum opacity at the edges
!	-take the rotation of star1 into account in the scattered 
!	 light (interpolate in precalculated rotationally broadened  
!	 intensity for several inclinations),
!       -irradiation (heating) of the companion in case icomp=1
!	 (reflection is already included)
!	-gravity reddening of the primary star
!	-check H- opacity
!	-fix bug at the Roche surface at L1 point for f=1.
!       -elnd has slow or no convergence for T approx< 90 K
!	-make the code both F77 and F90: place & on 73-rd character
!	 and on 6-th character, use 'space' instead of 'tab', use !/c
!        some new variables have 7 characters long names
!	-broadening of H lines
!       -save memory by declaring velocities & temperatures real*4
!
!		ndimf1,2,3 -main space xyz dimension parameters
!                     for the body frozen grid
!               ndim1,2,3 - main space xyz dimension parameters
!                      for the line of sight grid
!               mion - max. number of ions dimension (fixed=9)
!		matom - highest possible atom number 
!		mfreq - max. number of frequences
!		mline - max. number of sp. lines considered
!		mphase - max. number of grid rotations
!		mstarx - max. number of frequences for albedo or 
!			for star1, star2, star3 spectra
!		mspecx - max. number of precalculated spectra/files
!			for star1 for mspecx different temperatures
!               npfang -number of angles for the Mie phase function
!		mtemp,mdens,melm -number of temperatures, densities,
!			and molecules in the chemistry table
!
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
	include 'param.inc'
	integer, parameter:: mtemp=57,mdens=73,melm=11
        character(4):: dyp(matom)
        character(8):: image1
        character(3):: image2
        character(11):: image12
	real(i8):: ar(ndim1),at(ndim2),az(ndim3),azn(ndim3)
        real(i8):: far(ndimf1),fat(ndimf2),faz(ndimf3)
!	line of sight quantities are in 1D fields
        real(i8):: atemp(ndim3),adens(ndim3)
        real(i8):: ane(ndim3),avr(ndim3)
        real(i8):: avt(ndim3),avz(ndim3)
        real(i8):: avtrb(ndim3)
        real(i8):: adustd(ndim3),adustt(ndim3)
	integer(i1):: lshade(ndim3)
!	most memory consuming 3D fields and the body frozen grid
        real(i8):: ftemp(ndimf1,ndimf2,ndimf3)
        real(i8):: fdens(ndimf1,ndimf2,ndimf3),fne(ndimf1,ndimf2,ndimf3)
        real(i8):: fvr(ndimf1,ndimf2,ndimf3),fvt(ndimf1,ndimf2,ndimf3)
        real(i8):: fvz(ndimf1,ndimf2,ndimf3),fvtrb(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustd(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustt(ndimf1,ndimf2,ndimf3)
	integer(i1):: kshade(ndimf1,ndimf2,ndimf3)
        real(i8):: aint(ndim1,ndim2,mfreq)
!	line of sight with 1D fields up to untransparent point
        real(i8):: vr(ndim3),vt(ndim3),vz(ndim3),vtrb(ndim3)
        real(i8):: dustd(ndim3),dustt(ndim3)
        real(i8):: temp(ndim3),dens(ndim3),ekonc(ndim3),akonc(ndim3)
	integer(i1):: mshade(ndim3)
        real(i8):: alam(mfreq),flux(mfreq,mphase),area(ndim1,ndim2)
        real(i8):: fluxn(mfreq,mphase)
        real(i8):: alpha(mphase)
        real(i8):: sta(6,mion),d(3,matom),xi(8,matom)
        real(i8):: rr(ndim3,mion),stavs(ndim3,mion)
        integer:: necod(matom)
        real(i8):: denxnb(mstarx),denznb(mstarx)
!	H+He
        real(i8):: hi(ndim3),hneg(ndim3),h2(ndim3),hipf(ndim3,2)
        real(i8):: hei(ndim3)
        real(i8):: ophbf(ndim3),ophbf1(ndim3),ophbf2(ndim3)
        real(i8):: ophff(ndim3),ophff1(ndim3),ophff2(ndim3)
        real(i8):: ophrs(ndim3),ophrs1(ndim3),ophrs2(ndim3)
        real(i8):: ophn(ndim3),ophn1(ndim3),ophn2(ndim3)
!	dust+gas
        real(i8):: opmiex(mstarx,mspecx)
        real(i8):: opmies(mstarx,mspecx),opmiea(mstarx,mspecx)
	real(i8):: dtlow(mspecx),dthig(mspecx),drmf(mspecx)
	integer:: nmie(mspecx)
	real(i8):: pmiex(mstarx),pmies(mstarx),pmiea(mstarx)
	real(i8):: opmis(mspecx),opmia(mspecx)
        real(i8):: pfmiex(mstarx),pfang(npfang),pfmie(mstarx,npfang)
	real(i8):: pfmief(npfang)
        real(i8):: xsecx(mstarx),xsecy(mstarx,mspecx)
        real(i8):: txsec(mspecx),amixf(ndim3)
	real(i8):: popt(mtemp),popd(mdens),pop(mtemp,mdens,melm)
!	stellar spectra+albedos
        real(i8):: xstar1(mstarx,mspecx),star1(mstarx,mspecx)
        real(i8):: tspec1(mspecx),tspec2(mspecx)
        integer:: nstar1(mspecx),nstar2(mspecx)
        real(i8):: wstar1(mstarx),fstar1(mstarx)
	real(i8):: xstar2(mstarx,mspecx),star2(mstarx,mspecx)
	real(i8):: wstar2(mstarx),fstar2(mstarx) 
        real(i8):: xstar3(mstarx),star3(mstarx)
        real(i8):: alb1x(mstarx),alb1y(mstarx)
        real(i8):: alb2x(mstarx),alb2y(mstarx)
!        
        integer:: iat(mline),iion(mline),nion(mline)
        real(i8):: elo(mline),eup(mline),glo(mline),bij(mline)
        real(i8):: gr0(mline),gs0(mline),gw0(mline)
        real(i8):: gr(mline),gs(mline),gw(mline)
        real(i8):: aa(ndim3,mline),dop(ndim3,mline),popul(ndim3,mline)
        real(i8):: zip(mion-1,mline),zipp(mion-1)
        real(i8):: abund(mline),hmc(mline),wlab(mline)
	real(i8):: bolt
        real(i8):: bol,clight,hjed,plank,rydb,rydbh,pi,elmass,stefb
        real(i8):: everg,evcm,ergcm,rsol,grav,gmsol,emsol,pc,au
        real(i8):: dcut1,dcut2,dcut3,dcutn,denvac
        real(i8):: alam1,alamn,alams,cutoff,eps,offset
        real(i8):: phase1,phasen,dinc,dd,vgamma,rv,ebv
	real(i8):: xunt1,yunt1,xunt2,yunt2,xunt3,yunt3
        real(i8):: rmdfx1,rmdfx2,rmdfy1,rmdfy2,rmdfz1,rmdfz2
        real(i8):: stepfx,stepfy,stepfz,gainfx,gainfy,gainfz
        real(i8):: rmdfx3,rmdfx4,rmdfy3,rmdfy4,rmdfz3,rmdfz4
        real(i8):: stepfxb,stepfyb,stepfzb,gainfxb,gainfyb,gainfzb       
        real(i8):: rmdx1,rmdx2,rmdy1,rmdy2,rmdz1,rmdz2
        real(i8):: stepx,stepy,stepz,gainx,gainy,gainz
        real(i8):: rmdx3,rmdx4,rmdy3,rmdy4,rmdz3,rmdz4
        real(i8):: stepxb,stepyb,stepzb,gainxb,gainyb,gainzb       
        real(i8):: rstar,tstar,emstar
        real(i8):: xstar,ystar,zstar,vrotst,drotst,hst
        real(i8):: vxst,vyst,vzst,dlst,dlst2,dgst,ffst,albst,htst,htsta
	real(i8):: xspst,yspst,zspst,aspst,tspst
	real(i8):: rcp,tempcp,qq,vrxcp,vrycp,vrzcp,vrotcp
	real(i8):: xcp,ycp,zcp,vxcp,vycp,vzcp,dlcp,dlcp2,dgcp,ffcp
	real(i8):: albcp,htcp,htcpa
	real(i8):: emen,qqen,aen,ffen,hen
	real(i8):: tempen,densen,aneen,vtrben,dstden,dstten
	real(i8):: vrxsp,vrysp,vrzsp,vrotsp,rsp
	real(i8):: xsp,ysp,zsp,vxsp,vysp,vzsp
	real(i8):: tempsp,denssp,anesp,vtrbsp,dstdsp,dsttsp
	real(i8):: v1sm,v2sm,r1sm,r2sm,x1sm,y1sm,z1sm,x2sm,y2sm,z2sm
	real(i8):: vxsm,vysm,vzsm,xsm,ysm,zsm,psm
	real(i8):: tempsm,denssm,anesm,vtrbsm,edensm,dstdsm,dsttsm
	real(i8):: rrg,emrg,b1rg,b2rg,a1rg,a2rg,dr1rg,dr2rg,xrg,yrg,zrg
	real(i8):: xpolrg,ypolrg,zpolrg,vxrg,vyrg,vzrg
	real(i8):: temprg,densrg,anerg,vtrbrg
	real(i8):: edenrg,dstdrg,ede2rg,dst2rg,dsttrg
	real(i8):: adisc,rindc,routdc,emdc,rdc,xdc,ydc,zdc
	real(i8):: xdisc,ydisc,zdisc,vxdc,vydc,vzdc
	real(i8):: tempdc,densdc,anedc,vtrbdc,edendc,etmpdc
	real(i8):: dstddc,dsttdc
	real(i8):: aneb,rinnb,routnb,emnb,rnb,hinvnb,tinvnb,hwindnb
	real(i8):: hvelnb,vnb,evelnb
	real(i8):: xneb,yneb,zneb,vxnb,vynb,vznb,hcnb,hshdnb
	real(i8):: tempnb,densnb,anenb,vtrbnb,edennb,etmpnb
	real(i8):: dstdnb,dsttnb
	real(i8):: v1fw,v2fw,r1fw,r2fw,x1fw,y1fw,z1fw,x2fw,y2fw,z2fw
	real(i8):: vxfw,vyfw,vzfw,xfw,yfw,zfw,pfw
	real(i8):: tempfw,densfw,anefw,vtrbfw,edenfw,dstdfw,dsttfw
	real(i8):: ajet,rinjt,routjt
	real(i8):: xjt,yjt,zjt,xjet,yjet,zjet
	real(i8):: vjt,eveljt,rcjt,vxjt,vyjt,vzjt
	real(i8):: tempjt,densjt,anejt,vtrbjt,dstdjt,dsttjt,etmpjt,asymjt
	real(i8):: aufo,rinuf,routuf,emuf,ruf,xuf,yuf,zuf
	real(i8):: xufo,yufo,zufo,vxuf,vyuf,vzuf
	real(i8):: tempuf,densuf,aneuf,vtrbuf,edenuf,etmpuf
	real(i8):: dstduf,dsttuf
	real(i8):: rinsh,routsh,vsh,evelsh,rcsh,vxsh,vysh,vzsh
	real(i8):: tempsh,denssh,anesh,vtrbsh,dstdsh,dsttsh,etmpsh
	real(i8):: temp0,dens0,ane0,v0
	real(i8):: alphas,abhyd,abhel,wm,pom,opdepm,aint0,aintb,amixr
	real(i8):: freq,plocha,der,cont,fluxl,fluxs,refw,dvel,alamg
	real(i8):: vxstr,vystr,vzstr,xcpr,ycpr,zcpr,vxcpr,vycpr,vzcpr
	real(i8):: darkl,cdelta,vdif,opdep,amag,aflx
	integer:: loglam,imodel,irotat,ipart,ichemc,ielnd
	integer:: ithom,irayl,imie,imiepf,ihyd,iopac,iline,ionu,ior,iot
	integer:: nphase,iext,lunt1,lunt2,lunt3
	integer:: istar,icomp,ienv,ispot,ism,iring,idisc
	integer:: inebl,iflow,ijet,iufo,ishell
	integer:: idifst,irrst,ialbst,ispst,irrcp,ialbcp,itrg,itdc
	integer:: iinvnb,idennb,ivelnb,ishdnb,itnb,ituf,ivjt
	integer:: nbodf1,nbodf1a,nbodf1b,nfreq
	integer:: nbodf2,nbodf2a,nbodf2b
	integer:: nbodf3,nbodf3a,nbodf3b
	integer:: nbod1,nbod1a,nbod1b
	integer:: nbod2,nbod2a,nbod2b
	integer:: nbod3,nbod3a,nbod3b
	integer:: nalb1,nalb2,ndennb,iang,ii,jj,iprint,ivac1,kk,ivacn
	integer:: nu,nivac,ilin,ll,idst,ij,ipfang,i,j,k,iunt
	integer:: nspec1,nspec2,nstar3,ntxsec,nxsec,nline,ndust,npfmie
	integer:: jstar1,jstar2
!-----------------------------------------------------------------------
        write(6,'(A)')
        write(6,'(A)')'                       Shellspec'
        write(6,'(A)')
        write(6,'(A)')'       INPUT: '
        write(6,'(A)')' shellspec.in - main input (geometry,objects...)'
        write(6,'(A)')' line.dat - atomic data for the sp. lines'
        write(6,'(A)')'                          (optional if iline=1)'        
        write(6,'(A)')' shellspec.mod - 3D model of the shell '
        write(6,'(A)')'                          (optional if imodel=2)'
        write(6,'(A)')' abundances  - element abundances'
        write(6,'(A)')'                          (optional if ichemc=1)'
        write(6,'(A)')' lambda - list of lambda  (optional if loglam=2)'
        write(6,'(A)')' phases - orbital phases (optional if nphase=0)'  
        write(6,'(A)')' starspec1 - star spectrum (optional if lunt1>0)'
        write(6,'(A)')' starspec2 - star spectrum (optional if lunt2>0)'
        write(6,'(A)')' starspec3 - star spectrum (optional if lunt3>0)'
        write(6,'(A)')' albedo1 - albedo (optional if ialbst,irrst=1)'
        write(6,'(A)')' albedo2 - albedo (optional if ialbcp,irrcp=1)'        
        write(6,'(A)')' dust_opac - dust opacity (optional if imie>0)'                
        write(6,'(A)')' mie_phase - dust phase f.(optional if imiepf=1)'
        write(6,'(A)')' gas_opac - molec x-sectio (optional if iopac=1)'
        write(6,'(A)')' chem_eq_tab - molec popul (optional if iopac=1)'
	write(6,'(A)')' wind_prof - wind density (optional if idennb=1)'
        write(6,'(A)')'       OUTPUT: '
        write(6,'(A)')' shellspectrum - spectrum of the shell'
        write(6,'(A)')' lightcurve - LC or trailed spectrogram'        
        write(6,'(A)')' shellspec.out - more detail output'
        write(6,'(A)')' 2Dimage_xxx - 2D images at some frequency'
        write(6,'(A)')' errors - error messages'
	write(6,*)
	open(2,file='shellspec.out',status='unknown')
	open(3,file='errors',status='unknown')
	open(4,file='shellspectrum',status='unknown')
	open(11,file='lightcurve',status='unknown')
!	constants and units from NIST(2002)
	bol=1.3806503d-16
	clight=2.99792458d10
	hjed=1.66053873d-24
	plank=6.62606876d-27
	rydb=109737.31568549d0
	rydbh=109677.585d0
	pi=3.1415926535897931d0
!	pi=dacos(-1.d0) seems to have the same precission
	elmass=9.10938188d-28
	stefb=5.670400d-5
!	1EV=X ERG, X CM^-1; 1ERG= X CM^-1 
	everg=1.602176462d-12
	evcm=8.06554477d3
	ergcm=5.03411762d15
!	Sun
!	Lang 1992
!	rsol=6.9599d10
!	Cox 2000
!	rsol=6.95508d10
!	Nominal solar values, G*Msol is used instead of Msol
!	Prsa, Harmanec et al. 2016
	rsol=6.957d10
	grav=6.67408d-8
	gmsol=1.3271244d26
	emsol=gmsol/grav
	pc=3.0857d18
	au=1.4960d13
!-----------------------------------------------------------------------
!       Definition of the input quantities:                             
!       alam1, alamn, alams -start, end and step of wavelength in [A]    
!         Note that if loglam=0 the gas continuum opacity is calculated
!         only at alam1 & alamn and then interpolated.
!         So keep the iterval <alam1, alamn> short enough in that case.
!       loglam=0 equidistant step in lambda. Good for short interval.
!         H continuum opacity is calculated only at alam1 & alamn
!         and interpolated for given lambda to speed up calculations.
!       loglam=1 equidistant step in log(lambda), for long intervals.
!         Number of steps will be the same as for loglam=0.
!         H continum opacity is calculated at each lambda.
!       loglam=2 lambda's are read from file lambda. It sets nfreq. 
!         Good for comparison with the observations.
!         H continum opacity is calculated at each lambda.
!         (Dust opacity is interpolated to lambda independ. of loglam)	
!       cutoff - extension of the <alam1,alamn> interval in [A] when
!               reading the gas_opac table. Assuming that 
!               broadening by the velocity field dominates:
!               cutoff>maximal radial velocity/c*lambda 
!       imodel=1  calculate your own input shell model                  
!       imodel=2  read input shell model from `shellspec.mod'.          
!               You can ignore most of input below defining geometry,   
!               the velocity field and state quantities of objects but
!               you must still input the data for the scattering:           
!                 rstar,tstar,vxst,vyst,vzst                              
!               for the coordinate rotation:                        
!                 temp0,ane0,xcp,ycp,zcp 
!               and for the limb darkening:
!                 istar,rstar,tstar,dlst,dlst2,
!                 icomp,rcp,tempcp,dlcp,dlcp2,xcp,qq
!               and switches: lunt1,lunt2,lunt3,ithom,irayl,
!                 imie,imiepf,ihyd,iopac,iline,eps
!       irotat -option of interpolation from the body frozen grid
!               to the line of sight grid during the coord. rotation
!               0=linear interpolation, good for continuous fields,
!                 otherwise the result may depend on discontinuities
!                 or background (temp0,ane0,...) 
!                 does not support shadows (sets lshade=3)
!               1=nearest neighbour approximation, may be less smooth
!                 but can handle discontinuities and shadows
!       ipart  -option of partition functions 
!               [1-built in Dworetsky & Smalley, 2-Irwin]
!               (only ipart=1 is implemented so far)
!       ichemc -option of abundances, if ielnd=1 then ichemc=1                                     
!               [0-default solar, 1-read from file `abundances']       
!       ielnd=1 electron number densities provided in the input model
!               are ignored and code calculates el.num.dens. 
!               assuming LTE, from known temperature, density and 
!               chemical composition. File 'abundances' is read and must
!               contain 3.column which specifies which elements are
!               considered in Ne calculations, this sets ichemc=1
!       ielnd=0 electron number densities are known apriori and are
!               specified in the input model
!       ithom=0 Thomson scattering is off
!       ithom=1 Thomson scattering from stars is on 
!               (assumes optically thin environment)-check also shadows
!       irayl=0 Rayleigh scattering on neutral hydrogen is off.
!               If Lyman lines are treated explicitely in the linelist
!               set irayl=0 not to count the contribution twice
!       irayl=1 Rayleigh scattering from stars on neutral hydrogen is on 
!               (assuming optically thin environment)-check shadows
!       imie=0  Mie scattering and absorption on dust is off
!       imie=1  Mie scattering+absorption opacity is on.
!               Several species or input files can be included.
!               dust_opac file with tables must be provided.
!               Mie thermal and scattering emissivity on dust is on.
!               It is scattering of light from the stars assuming
!               optically thin medium.
!               Scattering emission can be isotropic or  
!               non-isotropic (see imiepf).
!       imie=2  Mie scattering+absorption opacity is on
!               Mie thermal emissivity is on, but
!               Mie scattering emissivity is off
!       imie=3  Mie scattering+absorption opacity is on
!               Mie thermal emissivity is on
!               Mie scattering emissivity is on but is isotropic and
!               assumes J=B(T) i.e. it is not scattered light from stars
!               Check for shadows.
!       imiepf  angular dependence of the scattered light from stars,
!               has an effect only if imie=1 
!       imiepf=1 angular dependent scattering emissivity, 
!               reads extra table with phase functions (mie_phase),
!               otherwise it is isotropic
!               In case there are several species in dust_opac 
!               this will redistribute the total scattering opacity. 
!       ihyd=1  hydrogen bound-free and free-free opacity is turned on
!               assuming only atomic H (no molecules)
!       iopac=1 additional tabulated gas true opacity is added
!               reads extra table with gas opacities (no scattering)
!       iline=0 No line opacity   
!       iline=1 line opacity is included. Spectral line parameters must
!               be specified in the file 'line.dat'
!       eps -artificial number <0.,1.> for test purpose which splits    
!               the line opacity (emissivity) into the true             
!               absorbtion (eps->1.) and coherent scattering (eps->0.). 
!               In LTE eps=1. ( S=eps*B+(1-eps)*J )  
!               If ithom=irayl=0 set also eps=1. for consistency                   
!       ionu, ior, iot -sequential indexes of frequency, x, and y point 
!               for which you want a more detailed output along the line
!               of sight (specified by x,y) 
!       offset -vertical shift applied to the normalized spectra output
!               to plot many spectra from different rotation phases        
!       phase1, phasen - start, end of the phase interval you want 
!               to cover [deg] (e.g. if xcp>0,ycp=zcp=0, dinc=90 then
!               phase1=-90 will start from the primary eclipse)
!       nphase -number of rotations (different view points) within 
!               the interval above
!               if nphase=0 then phase1 and phasen are ignored, and it
!               reads one column from the file `phases' with phases. 
!               These are values <0,1> and count from the x axis 
!               so that phase=0.0 or 1.0 is primary eclipse 
!               if xcp>0,ycp=zcp=0, dinc=90
!       dinc   -angle between rotation axis of the model and the line
!               of sight [deg], dinc=90.0 is edge on.
!       dd   -distance from the Earth in [pc]
!       vgamma -gamma velocity. It is applied only to the spectra and
!               lightcurves in the end 
!               (independently on the phase/view point)
!               2D images correspond to the original=final lambda(ionu)
!       iext=0  no reddening/extinction   
!       iext=1  reddening according to Cardelli, Clayton & Mathis 1989
!               it is applied to 2D images, spectra & lightcurves
!               covers: 0.1<=lambda[mic]<=10/3
!       rv      =A(V)/E(B-V)
!       ebv     =E(B-V) [mag]
!----------------intrinsic spectra specifications:
!       lunt1=0 all objects with density from <dcut1,dcut2> interval are 
!               nontransparent blackbodies with the same temperatures as
!               in the case of transparency. 
!       lunt1>0 all objects with density within <dcut1,dcut2> are 
!               nontransparent and have an intrinsic intensity spectrum.
!               The spectrum is read from file `starspec1'.
!       lunt1=1 the x,y column input required with wavelength [A] and
!               H_lambda flux [erg/cm^2/s/A] (as an output of SYNSPEC)
!       lunt1=2 the x,y column input required with wavelength [A] and
!               I_nu intensity [erg/cm^2/s/Hz/sterad] 
!       lunt1=3 the 4 column input required with idummy,frequency [Hz],
!               dummy, F_nu flux [erg/cm^2/s/Hz]
!               (output of coolTlusty, unit 21, first 2 rows are dummy) 
!       xunt1  -multiplication factor applied to starspec1 x-column
!               if it is not in the correct-required units
!               (otherwise set it =1.)
!       yunt1  -multiplication factor applied to starspec1 y-column
!               if it is not in the correct-required units
!               (otherwise set it =1.)
!       lunt2,xunt2,yunt2 -the same meaning as above except that these
!               deal with density interval <dcut2,dcut3> and
!               the spectrum is read from file `starspec2'. 
!       lunt3,xunt3,yunt3 -the same meaning as above except that these
!               deal with density interval <dcut3,dcutn> and
!               the spectrum is read from file `starspec3'.          
!-----------definitions of grids:
!       There are two main grids: 
!         -body frozen grid - defines your model 'frozen' at a moment
!         -line of sight grid - defines grid for rad.transfer calcul.
!       Each grid can be composed of two subgrids: A,B
!       If subgrids overlap then the denser subgrid has priority.
!       Subgrids A,B are merged into one grid.
!       rmdfx1<rmdfx2, rmdfy1<rmdfy2, rmdfz1<rmdfz2 - define 
!         the box A in the body frozen frame (if imodel=1) [R_sol]
!       rmdfx3<rmdfx4, rmdfy3<rmdfy4, rmdfz3<rmdfz4 - define 
!         the box B in the body frozen frame (if imodel=1) [R_sol]
!       stepfxyz>0.,stepfxyzb>0. -is a mean distance between the x,y,z 
!         grid points of box A,B respectively [R_sol]
!         It determines the number of grid points:
!       nbodf1, nbodf2, nbodf3 -number of grid points in x, y, z 
!         direction in body frozen coordinates of the model.
!         (Points are overridden 
!         by the values from `shellspec.mod' if imodel=2)
!       gainfxyz, gainfxyzb -grid step multiplication factors
!         of the body frozen grids A,B to allow for logarithmic grid 
!         [gainfx=(x_{i+1}-x_{i})/(x_{i}-x_{i-1})]
!         e.g. gainfx=1. for equidistant step
!         gainfx>1. step increases symetrically from the middle to 
!         the left and to the right
!       rmdx1<rmdx2, rmdy1<rmdy2, rmdz1<rmdz2 - define the box A
!         of the observer's line of sight frame all in R_sol. 
!       rmdx3<rmdx4, rmdy3<rmdy4, rmdz3<rmdz4 - define the box B
!         of the observer's line of sight frame all in R_sol. 
!         Observer looks along the opposite z-direction.
!       stepxyz>0.,stepxyzb>0. -is a mean distance between the x,y,z 
!         grid points in box A,B respectively [R_sol]
!         They determine the number of grid points:
!       nbod1, nbod2, nbod3 -number of grid points in x, y, z 
!         in the line of sight observer's frame.
!       gainxyz, gainxyzb -grid step multiplication factors
!         (common ratio of the geometric sequence) of the line of sight
!         grids A,B respectively [gainx=(x_{i+1}-x_{i})/(x_{i}-x_{i-1})]
!         e.g. gainx=1. for equidistant step
!         gainx>1. step increases symetrically from the middle to 
!         the left and to the right
!       To turn off the B grids simply set their dimension to zero e.g. 
!         rmdx3=rmdx4, rmdy3=rmdy4, rmdz3=rmdz4  
!-----------------object definitions:
!       istar,icomp,ispot,ism,iring,idisc,ienv,
!       inebl,iflow,ijet,iufo,ishell
!       These are main on/off switches for the objects.
!       They are ordered according to priority.
!       Priority is determined in smod1.
!       It is important in case objects happen to overlap.
!-------primary star (central object)----------
!       istar=0 accompanied by rstar=0 will switch off the primary
!       istar=1 central object is a nontransparent uniformly rotating
!         sphere. Its density is set to <dcut1,dcut2>. It can be either 
!         black body with T=tstar if lunt1=0 or has its intrinsic 
!         intensity spectrum if lunt1>0. In case of scattering or 
!         reflection of its light by other objects its rotation is 
!         ignored.
!         Code ignores: dgst,ffst,qq                 
!       istar=2 central object is a detached component of a binary.       
!         It has a Roche shape defined by ffst<=1, synchronous rotation,
!         is nonstrasparent with density within <dcut1,dcut2>.  
!         It can be either black body with T=tstar if lunt1=0 or 
!         has its intrinsic intensity spectrum if lunt1>0. 
!         You must also set: xcp>0,qq>0
!         Code also calculates/ignores: xstar,ystar,zstar,vrotst
!         ,drotst,hst,rstar
!       istar=3 central object is a figure 8 contact system. It has
!         a Roche shape defined by 1<ffst<=2, synchronous rotation,
!         is nonstrasparent with density within <dcut1,dcut2>.  
!         It can be either black body with T=tstar if lunt1=0 or 
!         has its intrinsic intensity spectrum if lunt1>0. 
!         You must also set: xcp>0,qq>0
!         Code also calculates/ignores: 
!           xstar,ystar,zstar,vrotst,drotst,hst,rstar,icomp
!       if istar>1 or icomp>1 or (istar>0 and icomp>0 and vxst>clight)
!         then code calculates (from emstar,xcp,qq):
!         ycp,zcp,vxst,vyst,vzst,vxcp,vycp,vzcp
!         assuming circular orbit.
!       rstar -radius of the central star in [R_sol] 
!         if istar>1 (Roche Geometry) this value will be used for
!         scattering in the circumstellar matter and irradiation effect 
!         on the companion which use spherical approximation
!       tstar -effective temperature of the central star in [K]
!         without gravity darkening and irradiation. This value will 
!         be used for scattering in the circumstellar matter (in case 
!         of black body) and irradiation effect on the companion
!         if istar=2 it is the temperature at the rotation pole 
!         if istar=3 it is the temperature at the rotation pole of
!         the more massive star
!       emstar -mass of the central star in [M_sol]
!       xstar,ystar,zstar -define unit aiming vector of the rotational        
!         axis of the central star         
!       vrotst -equatorial rotation velocity of the central star [km/s]
!         in case istar=1 corresponding to the equatorial angular vel.  
!       idifst -on/off differential rotation only for istar=1
!       idifst=0 no differential rotation
!       idifst=1 smooth differential rotation
!         omega(phi)=omega_eq-(omega_eq-omega_pol)*sin(phi)**2
!       idifst=2 step function differential rotation
!         omega(phi)=omega_eq  for z/rstar<hst
!         omega(phi)=omega_pol for z/rstar>hst
!       drotst - the ratio of angular velocity at the rotation pole to
!         the angular vel. at the equator: drotst=omega_pol/omega_eq.
!       hst -break in the step function =z/rstar for idifst=2
!       vxst, vyst, vzst -net velocity components                 
!         of the center of mass of the central star [km/s]
!         (if vxst>clight and istar>0 and icomp>0 then see istar)
!       dlst -limb darkening coefficient of the central star
!       dlst2 -second limb darkening coefficient 
!            I=1-dlst*(1-mu)-dlst2*(1-mu)**2
!       dgst -gravity darkening coefficient (beta) of the central star
!         (0.25 for radiative, 0.08 for convective atmospheres)
!         It is dummy if istar=1.
!       ffst<=1 -Roche lobe fill-in factor of the primary. Its is
!         the distance of the inner substellar point of the primary
!         (between the stars) from the center of the primary relative
!         to the distance to L1, the Roche lobe is reproduced if ff=1
!       1<ffst<=2 -Roche lobe fill-out factor of the contact system
!               ffst=(C1-C)/(C1-C2)+1
!               It is dummy if istar=1.
!       irrst=0 -irradiation and reflection effect is off
!           (ialbst,albst,htst,htsta have no meaning in this case)
!       irrst=1 -irradiation of the object from the companion is on.
!           Irradiation (heating) applies only if istar=1,2.
!           Reflection of the sp. of companion applies if istar=1,2
!           (rcp,tempcp>0 ... are presumed).
!       ialbst=1 monochomatic albedo is red from file=albedo1 
!           (if irrst=1). It should be compatible with Bond albedo.
!       albst  -Bond albedo <0,1>
!       htst   -heat redistribution parameter in case of irradiation, 
!           fraction of the heat absorbed on the day side which is 
!           redistributed over the day-night sides, <0,1>, 
!           0-nothing is redistributed and nothing goes to the night, 
!           1-all the energy (which is not reflected) impinging on 
!           the planet is evenly distributed over the day-night sides.
!           It is analoguous to the so called Pn parameter of A.Burrows
!           (a fraction of the irradiating energy impinging on 
!           the day side which is transfered to and irradiated from 
!           the night side), Pn=(1-albst)*htst/2
!       htsta  -degree of the inhomegenity of the heat transport, <0,1>.
!           1-homegeneous, 0-cosine dependence
!           T**4=T0**4(htsta+4(1-htsta)/pi*cos_latitude)
!       ispst=1/0 will turn on/off a spot on the star if istar=1
!           (it has the shape of a circle)
!       xspst,yspst,zspst -define unit aiming vector of the location
!           of the spot center on the surface 
!       aspst -angular radius of the spot in [deg]
!       tspst -ratio of the spot temperature to the ambient temperature
!           (i.e. temperature accounting for the reflection effect...)
!-------
!       temp*,dens*,ane* - state quantities in various objects   
!               temperature, density, electron number density [K,CGS]
!       vtrb* - microturbulent velocity in various objects [km/s],
!               it does not apply to nontransparent objects
!       dstd* - density of dust in various objects [g/cm^3]
!               you must also set dens*>0. to have an effect
!       dstt* - temperature of dust in various objects [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       ishd* - define shadows for scattering (in some objects)
!             =0 no scattering
!             =1 scattering from central star only
!             =2 scattering from secondary star only
!             =3 scattering from both stars
!-------companion or secondary star
!       icomp=0 secondary off
!       icomp=1 secondary on, it is a uniformly rotating nontransparent
!               sphere. It may be a blackbody with T=tempcp if lunt2=0
!               or has its own spectrum if lunt2>0. Its density is set 
!               to <dcut2,dcut3>. Code ignores: dgcp,ffcp,qq
!       icomp=2 secondary is a detached component of a binary.
!         It has a Roche shape defined by ffcp<=1, synchronous rotation,
!         is nonstrasparent with density within <dcut2,dcut3>.  
!         It can be either black body with T=tempcp if lunt2=0 or 
!         has its intrinsic intensity spectrum if lunt2>0. 
!         You must set: xcp>0,qq>0,emstar>0
!         Code also calculates/ignores: vrxcp,vrycp,vrzcp,vrotcp,rcp
!       rcp  -radius of the spherical companion [R_sol],
!         if icomp=2 this input is used only for the scattering 
!         and irradiation from the object otherwise it is superfluous 
!       tempcp -see primary star above, this value is used for 
!         the scattering on the circumstellar material and irradiation 
!         of the primary [K]
!       qq -mass ratio (companion/star), important only for Roche geom.
!           if istar>1 or icomp>1
!       vrxcp, vrycp, vrzcp -define unit aiming vector of the rotational      
!         axis of the secondary star (companion)
!       vrotcp -equatorial rotation velocity of the companion [km/s]
!       xcp,ycp,zcp -location of the center (of mass) of 
!               the companion [R_sol]
!       vxcp,vycp,vzcp -components of the velocity vector of the center 
!               (of mass) of the companion [km/s]
!       dlcp -limb darkening coefficient of the secondary star
!       dlcp2 -second limb darkening coefficient (the same as dlst2)
!       dgcp -gravity darkening coefficient (beta) of the secondary
!       ffcp<=1 -Roche lobe filling factor of the secondary is 
!         the distance of the inner substellar point of the secondary
!         from the center of the secondary relative to 1-L1, 
!         the Roche lobe is reproduced if ffcp=1 
!       irrcp=0 -irradiation and reflection effect is off
!           (ialbcp,albcp,htcp,htcpa have no meaning in this case)
!       irrcp=1 -irradiation of the secondary from the primary is on.
!           Irradiation (heating) applies only if icomp=2.
!           Reflection (of the spectrum of primary) applies if icomp=1,2
!           (istar=1,2 and rstar,tstar>0 are presumed)
!       ialbcp=1  monochomatic albedo is red from file=albedo2
!           (if irrcp=1). It should be compatible with the Bond albedo.
!       albcp  -Bond albedo <0,1>
!       htcp   -heat transport parameter in case of the irradiation. 
!           The same as htst, <0,1>.
!       htcpa  -degree of the inhomegenity of the heat transport, <0,1>. 
!           1-homegeneous, 0-cosine dependence, the same as htsta.
!-------spot or third star
!       ispot=0 spot is off
!       ispot=1 spot is on, it is a uniformly rotating sphere
!       vrxsp, vrysp, vrzsp -define unit aiming vector of the rotational      
!               axis of the spot
!       vrotsp -equatorial rotation velocity of the spot [km/s]
!       rsp  -radius of the spherical spot [R_sol]
!       xsp,ysp,zsp -location of the center of the spot [R_sol]
!       vxsp,vysp,vzsp -components of the velocity vector of the center 
!               of the spot [km/s]
!       tempsp -constant temperature [K]
!       denssp -gas density [g/cm^3]
!       anesp -electron number density [cm^-3]
!       dstdsp -dust density [g/cm^3]
!               you must also set denssp>0. to have an effect
!       dsttsp -dust temperature [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbsp -microturbulence [km/s]
!-------stream
!       ism=0/1  -stream off/on
!       v1sm   -stream velocity at the beginnig of stream [km/s]
!       v2sm   -stream velocity at the end of stream [km/s]
!         velocity is directed from beginning to end 
!       r1sm    -radius of the stream at the beginning [R_sol]
!       r2sm    -radius of the stream at the end [R_sol]
!         notice that although the radius changes the streamlines 
!         are parallel (contrary to jet) 
!       x1sm,y1sm,z1sm -position of the beginning of the stream [R_sol]
!       x2sm,y2sm,z2sm -position of the end of the stream [R_sol]
!       vxsm, vysm, vzsm -net velocity [km/s]
!         you can use it also to mimic orbital drag or if the center
!         of rotation is not at the center of coordinates
!       xsm,ysm,zsm -rotational vector of stream
!       psm -rotational period of stream in days
!       tempsm -temperature [K], constant along the stream 
!       denssm - gas density at the beginning [g/cm^3] and scales along
!         the stream to satisfy the continuity equation: 
!         density=denssm*v1sm*r1sm**2/(vsm*rsm**2)*exp(t/rsol*edensm)
!         where t is distance along the stream and exp term allows
!         e.g. for a dust destruction
!       anesm - electron number density at the beginning [cm^-3],
!         similar to the density but if ielnd=1 then it is overriden by
!         the calculation from the state quantities
!       edensm -density dependence exponent to enable the modeling 
!         of additional phenomena
!       dstdsm -dust density [g/cm^3], it changes along the stream like 
!         the gas density
!               you must also set denssm>0. to have an effect
!       dsttsm -dust temperature [K], constant along the stream
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbsm -microturbulence velocity [km/s]
!-------ring
!       iring>0 ring is on
!       rrg -radius of the ring [R_sol]
!       emrg -mass in its center to calculate velocities [Msol]
!       b1rg, b2rg -specifies the arc from-to in [deg], b1rg><b2rg
!         The location of the zero angle is not simple to explain so 
!         test it first or consult subroutine trans. In many cases
!         it will be along the x axis.
!       a1rg,a2rg -vertical half width of the ring at the beginning
!         and end of the arc in [R_sol]
!       dr1rg, dr2rg -horizontal half thickness the ring at 
!         the beginning and end of the arc in [R_sol]
!         The crosssection, C, of the ring may vary along the arc and
!         is C1=4*a1rg*dr1rg at the beginning.
!       xrg, yrg, zrg -location of the center in [R_sol]
!       xpolrg, ypolrg, zpolrg -orientation of the polar axis
!       vxrg, vyrg, vzrg -net overall space velocity [km/s]
!       edenrg, ede2rg -density dependence exponent to enable 
!         the modeling of additional phenomena. Density, dust density 
!         and electron number density change along the ring (arc) 
!         to safisfy continuity equation+additional phenomenon 
!         e.g. destruction (lifetime) of dust grains along the arc. 
!       dstdrg, dst2rg -dust density at the beginning (b1rg) [g/cm^3].
!               you must also set densrg>0. to have an effect
!       itrg=1 then
!         gas density=densrg*C1/C*[|t-b1rg|/pi+1]**edenrg
!         electron num. density=anerg*C1/C*[|t-b1rg|/pi+1]**edenrg
!         dust density=dstdrg*C1/C*[|t-b1rg|/pi+1]**edenrg+
!                 dst2rg*C1/C*[|t-b1rg|/pi+1]**ede2rg
!       itrg>or< 1 then
!         gas density=densrg*C1/C*dexp[|t-b1rg|/pi]**edenrg
!         electron num. density=anerg*C1/C*dexp[|t-b1rg|/pi]**edenrg
!         dust density=dstdrg*C1/C*dexp[|t-b1rg|/pi]**edenrg+
!                dst2rg*C1/C*dexp[|t-b1rg|/pi]**ede2rg
!         where t-is angle along the arc.
!       densrg -gas density at b1rg [g/cm^3]
!       anerg -electron number density at b1rg [cm^-3]
!       temprg -constant gas temperature [K]
!       dsttrg -constant dust temperature [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbrg -microturbulence [km/s]
!-------disk (accretion disk around some object)
!       idisc=0 switch off the disc		
!       idisc=1 disc has the shape of a rotating wedge 
!               limited by inner and outer radii (spherical surfaces)
!               tip of the wedge is at the center of the object
!       idisc=2 disc has the shape of a slab
!               limited by inner and outer radii (spherical surfaces)
!       idisc=3 disc has the shape of a rotating ellipsoid
!               limited by inner spherical and outer ellipsoidal surface
!       adisc -angular halfwidth of the disc wedge [deg] 
!              (if idisc=1)
!             -half of the thickness of the disc slab [R_sol] 
!              (if idisc=2)
!             -semiaxis of the ellipsoid along the rotational axis 
!              [R_sol] (if idisc=3)	
!       rindc -inner radius of the disc [R_sol]
!       routdc -outer radius of the disc [R_sol] or
!              -semiaxis of the ellipsoid perpendicular to the rotation 
!               axis, if idisc=3, [R_sol]
!       emdc -mass of the object in the disk center [M_sol]
!         it determines its Keplerian velocity
!       rdc  -radius of the object in the disk center [R_sol]
!         it determines its temperature structure if itdc=2
!       xdc,ydc,zdc -location of the disk center in [R_sol]
!       xdisc,ydisc,zdisc -components of the unit aiming vector of 
!               the rotational axis of the Keplerian disc around emstar
!       vxdc, vydc, vzdc -net velocity components                 
!               of the center of the disc [km/s]
!       densdc -gas density at rindc [g/cm^3]
!       anedc  -electron num. density at rindc [cm^-3] 
!       tempdc -characteristic gas temperature [K], see below
!       edendc -radial density dependence exponent 
!               (dens, ane and dust density are a function of r)
!               Rho(r) ~ Ne(r) ~ densdc*(r/rindc)**edendc	                        
!       itdc=1  gas & dust temperatures are constant
!       itdc=2  gas & dust temperatures are a function of r
!               gas:  T(r)=tempdc*(rdc/r)**0.75*(1-(rdc/r)**0.5)**0.25       
!               dust: T(r)=dsttdc*(rdc/r)**0.75*(1-(rdc/r)**0.5)**0.25       
!       itdc=3  gas & dust temperatures as a power law
!               gas:  T(r)=tempdc*(r/rindc)**etmpdc
!               dust: T(r)=dsttdc*(r/rindc)**etmpdc
!       etmpdc  -exponent of the radial temperature dependence
!       dstddc  -dust density at rindc [g/cm^3]
!               you must also set densdc>0. to have an effect
!       dsttdc  -characteristic dust temperature [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbdc -microturbulence [km/s]
!-------envelope around the primary star
!       ienv,emen,ggen,ffen have similar meaning to istar,emstar,qq,ffst
!       ienv=2 envelope is on, has a detached Roche shape
!       ienv=3 envelope is on, has a contact Roche shape
!           (common envelope)
!       emen -mass of the central star [M_sol]
!       qqen -mass ratio (companion/star)
!       ffen<=1 -Roche lobe fill-in factor of the detached envelope. 
!         It is radius of the substellar point of the envelope
!         relative to the radius of the L1. Roche lobe has ffen=1.
!       1<ffen<=2 -Roche lobe fill-out factor of the contact envelope
!               ffen=(C1-C)/(C1-C2)+1
!       hen -vertical limit [R_sol], limits the envelope in 
!         the direction perpendicular to the orbital plane to z<+-hen
!       tempen -constant temperature [K]
!       densen -constant gas density [g/cm^3]
!       aneen -constant electron number density [cm^-3]
!       dstden -constant dust density [g/cm^3]
!               you must also set densen>0. to have an effect
!       dstten -constant dust temperature [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrben -microturbulence [km/s]
!-------nebula (protoplanetary disk/nebula around central object)
!              it is defined in cylindrical coordinates (r,z)
!       inebl not=1  -nebula is off
!       inebl=1 flared protoplanetary or accretion disk
!               vertical scale height is: 
!               H(r)=hcnb*(gamma*k*T_gas(r)/m)**0.5
!             Vertical structure:
!               fdens=dens0*dexp(-erz**2/H**2/2.d0)
!               dens0 is midplane density calculated from surface dens.
!               Density may have wind region.		
!               Gas temperature may have vertical temperature inversion.
!             Radial structure:
!               surface density decreases ~(r/rinnb)**edennb
!               dust dens & electron num. dens are ~ density.
!               Gas and dust temperatures change with radius.
!       aneb    -vertical cutoff of the nebula at particular r in [H]
!               cutoff(r)=+-aneb* H(r)
!       rinnb   -inner radius of the nebula [R_sol]
!       routnb  -outer radius of the nebula [R_sol]
!       emnb    -mass of the object in the nebula center [M_sol]
!       rnb     -radius of the object in the nebula center [R_sol]
!       iinvnb,hinvnb,tinvnb -describe vertical gas temp. inversion
!               (only if itnb=3).
!               There is no inversion in dust temperature.
!       hinvnb  -start of vertical gas temp. inversion in [H]
!               for z(r)>hinvnb*H(r) (except itnb=1 or 2)
!               hinvnb>aneb or itnb=1 or itnb=2 means no inversion
!       tinvnb  -temperature multiplication factor in the inversion
!       iinvnb=0 no inversion
!       iinvnb=1 step fuction inversion
!               gas temp(z,r)=temp0(r)*tinvnb
!       iinvnb=2 linear inversion from temp0 at hinvnb*H(r) to
!               tinvnb*temp0(r) at aneb*H(r)
!               temp0 is midplane gas temperature
!       idennb, hwindnb -describe vertical density profile
!       hwindnb -start of the wind region in vertical scale-heights
!               rho(z)=rho(0)*dexp(-erz**2/H**2/2.d0)
!               but for z>hwindnb*H
!               rho(z)=rho(0)*dexp(-hwindnb**2/2.d0)
!               i.e. rho(z)=rho(hwindnb*H)= const  
!               electron n.d. and dust density are proportional to gas
!               and thus will also have wind region
!               hwindnb>aneb will turn off the wind region
!       idennb  allows to use optional vertical density profile
!       idennb=1 reads file wind_prof with rho=f(z)
!               hwindnb is then ignored
!       hcnb    allows to multiply the classical vertical scale-height
!               by some factor
!       ivelnb, hvelnb, vnb, evelnb -describe radial velocity field
!       ivelnb=1 will add a radial velocity component to Keplerian vel.
!               to simulate disk wind
!       hvelnb  -start of the radial wind in vertical scale-heights
!               will apply for z>=hvelnb*H
!       vnb     -terminal wind velocity [km/s]
!       evelnb  -velocity exponent, radial velocity is:
!               v_rad=vnb*(1.d0-Rnb/er)**evelnb
!               er-is distance from the center=dsqrt(r^2+z^2)
!       ishdnb  =0,1,2,3 (describes shadows for scattering)
!       hshdnb  -start of the scattering region in [H]
!               if z>hshdnb*H then kshade=ishdnb else kshade=0
!       xneb,yneb,zneb -components of the unit aiming vector of 
!               the rotational axis of the Keplerian disc around emnb
!       vxnb, vynb, vznb -net velocity components                 
!               of the center of the nebula [km/s]
!       tempnb  -characteristic gas temperature [K]
!       itnb    -regulates temperature structure
!       itnb=1  nebula gas and dust temperatures are constant
!               gas temp=tempnb, dust temp=dsttnb (no inversion)
!       itnb=2  nebula gas and dust temp. are a function of r only
!               T(r)=tempnb*(Rnb/r)**0.75*(1-(Rnb/r)**0.5)**0.25
!               T(r)=dsttuf*(Rnb/r)**0.75*(1-(Rnb/r)**0.5)**0.25
!               no inversion       
!       itnb=3  disc temperature as a power law (e.g. protopl. discs)
!               gas:  T(r)=tempnb*(r/rinnb)**etmpnb
!               dust: T(r)=dsttnb*(r/rinnb)**etmpnb
!               there may be a gas temperature inversion in z
!       etmpnb  -exponent of radial temperature dependence
!       densnb -gas density at rinnb (at midplane) [g/cm^3]
!       anenb  -electron num. density at rinnb (at midplane) [cm^-3]
!       edennb -radial density dependence exponent of surface density
!               (dens, ane and dust density are a function of r)
!               Ne(r,z) ~ Rho_dust(r,z) ~ Rho_gas(r,z)
!       dstdnb -dust density at rinnb [g/cm^3] (at midplane)
!               you must also set densnb>0. to have an effect
!       dsttnb -characteristic dust temperature [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbnb -microturbulence [km/s]
!-------flow
!         it is identical to the stream but lower priority
!       iflow=0/1  -stream off/on
!       v1fw   -stream velocity at the beginnig of stream [km/s]
!       v2fw   -stream velocity at the end of stream [km/s]
!         velocity is directed from beginning to end 
!       r1fw    -radius of the stream at the beginning [R_sol]
!       r2fw    -radius of the stream at the end [R_sol]
!         notice that although the radius changes the streamlines 
!         are made paralel (contrary to jet) 
!       x1fw,y1fw,z1fw -position of the beginning of the stream [R_sol]
!       x2fw,y2fw,z2fw -position of the end of the stream [R_sol]
!       vxfw, vyfw, vzfw -net velocity [km/s]
!         you can use it also to mimic orbital drag or if the center
!         of rotation is not at the center of coordinates
!       xfw,yfw,zfw -rotational vector of stream
!       pfw -rotational period of stream in days
!       tempfw -temperature [K], constant along the stream 
!       densfw - is density at the beginning [g/cm^3] and scales along
!         the stream to satisfy the continuity equation: 
!         density=densfw*v1fw*r1fw**2/(vfw*rfw**2)*exp(t/rsol*edenfw)
!         where t is distance along the stream
!       anefw - electron number density at the beginning [cm^-3]
!         similar to the density but if ielnd=1 then it is overriden by
!         the calculation from the state quantities
!       edenfw -density dependence exponent to enable the modeling 
!         of additional phenomena
!       dstdfw -dust density [g/cm^3], it changes along the stream like 
!         the gas density
!               you must also set densfw>0. to have an effect
!       dsttfw -dust temperature [K], constant along the stream 
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbfw -microturbulence velocity [km/s]
!-------jet
!       ijet=0   switch off the jet
!       ijet=1  jet has only one -primary cone	
!       ijet=2  jet has two cones: the primary cone and the opposite one
!       ajet -angle halfwidth of the jet cones [deg]
!         streamlines flare according to the opening angle
!       rinjt, routjt -radius boundaries of the jet cones [R_sol]
!       xjt,yjt,zjt -location of the jet origin [R_sol]
!       xjet,yjet,zjet -components of the unit aiming vector 
!               of the primary jet cone
!       ivjt -the velocity field switch
!       ivjt not equal 2 -radial velocity is polynomial 
!                v(r)=vjt*(r/rinjt)**eveljt
!       ivjt=2  radial velocity is: v(r)=vjt*(1-rcjt/r)**eveljt
!       vjt  -radial velocity at the inner edge or terminal velocity
!                depending on ivjt [km/s]
!       eveljt -velocity exponent
!       rcjt -radius of the object in the jet, only if ivjt=2, [R_sol] 
!       vxjt, vyjt, vzjt -net velocity component [km/s] 
!       tempjt -characteristic gas temperature [K]
!               jet gas temperature is a power law 
!               gas:  T(r)=tempjt*(r/rinjt)**etmpjt
!       densjt -gas density [g/cm**3] at rinjt, it scales along the jet
!         to satisfy the continuity equation
!         density=densjt*rinjt**2/routjt**2*v(rinjt)/v(r)
!       anejt - electron number density [cm**-3] at rinjt. It changes
!         along the jet like the gas density but if ielnd=1 then 
!         it is overriden by the calculation from the state quantities
!       dstdjt -dust density [g/cm**3] at rinjt, changes along the jet 
!         like the gas density
!               you must also set densjt>0.0 to have an effect
!       dsttjt -characteristic dust temperature [K]
!               dust temperature is a power law 
!               dust:  T(r)=dsttjt*(r/rinjt)**etmpjt
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbjt -microturbulence [km/s]
!-------ufo 
!         it is identical to DISK (same subroutine) but lower priority
!       iufo=0 switch off the ufo		
!       iufo=1 ufo has the shape of a rotating wedge
!               limited by inner and outer radii (spherical surfaces)
!               tip of the wedge is at the center of the object
!       iufo=2 ufo has the shape of a slab
!               limited by inner and outer radii (spherical surfaces)
!       iufo=3 ufo has the shape of a rotating ellipsoid
!               limited by inner spherical and outer ellipsoidal surface
!       aufo -angular halfwidth of the ufo wedge [deg] 
!              (if iufo=1)
!             -half of the thickness of the ufo slab [R_sol] 
!              (if iufo=2)
!             -semiaxis of the ellipsiod along the rotational axis 
!              [R_sol] (if iufo=3)
!       rinuf -inner radius of the ufo [R_sol]
!       routuf -outer radius of the ufo [R_sol] or
!              -semiaxis of the ellipsoid perpendicular to the rotation 
!               axis, if iufo=3, [R_sol]
!       emuf -mass of the object in the ufo center [M_sol]
!       ruf  -radius of the object in the ufo center [R_sol]
!       xuf,yuf,zuf -location of the disk center in [R_sol]
!       xufo,yufo,zufo -components of the unit aiming vector of 
!               the rotational axis of the Keplerian disc around emuf
!       vxuf, vyuf, vzuf -net velocity components                 
!               of the center of the ufo [km/s]
!       tempuf  -temperature [K]
!       ituf=1  ufo gas and dust temperatures are constant
!       ituf=2  ufo gas and dust temperatures are a function of r
!               T(r)=tempuf*(Ruf/r)**0.75*(1-(Ruf/r)**0.5)**0.25
!               T(r)=dsttuf*(Ruf/r)**0.75*(1-(Ruf/r)**0.5)**0.25       
!       ituf=3  ufo gas and dust temperatures as a power law
!               T(r)=tempuf*(r/rinuf)**etmpuf
!               T(r)=dsttuf*(r/rinuf)**etmpuf
!       etmpuf  -exponent of radial temperature dependence
!       densuf -gas density at rinuf [g/cm^3]
!       aneuf  -electron num. density at rinuf [cm^-3]
!       edenuf -radial density dependence exponent 
!               (dens, ane and dust density are a function of r)
!               Rho(r) ~ Ne(r) ~ densuf*(r/rinuf)**edenuf
!       dstduf -dust density at rinuf [g/cm^3] 
!               you must also set densuf>0. to have an effect
!       dsttuf -characteristic dust temperature [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbuf -microturbulence [km/s]
!-------shell
!       ishell=0   switch off the shell
!       ishell=1   velocity, dens, temp, ane are constant
!       ishell=2   radial velocity is v(r)=vsh*(r/rinsh)**evelsh
!           Ne(r)~Rho(r)=denssh*(rinsh/r)**2*vsh/v(r)), temp=const.
!       ishell=3   radial velocity is v(r)=vsh*(1-rcsh/r)**evelsh
!           Ne(r)~Rho(r)=denssh*(rinsh/r)**2*v(rinsh)/v(r)),temp=const.
!       rinsh, routsh -inner, outer radius of the shell in [R_sol]     
!       vsh -velocity of the uniformly expanding shell [km/s]           
!       evelsh -exponent of velocity dependence
!       rcsh - core/photospheric radius of the star in shell [R_sol]
!       vxsh, vysh, vzsh -net velocity [km/s]
!       tempsh -temperature [K]
!       denssh -gas density at rinsh [g/cm^3]
!       anesh -electron number density at rinsh [cm^-3]
!       dstdsh -dust density at rinsh [g/cm^3], 
!         it changes as the gas density
!               you must also set denssh>0. to have an effect
!       dsttsh -dust temperature [K]
!               dust temperatures must be higher than the condensation
!               temperatures of the species (see dust_opac) 
!               to have an effect
!       vtrbsh -microturbulence [km/s]
!-------background
!       v0 -constant uniformly expanding velocity of background [km/s]
!       temp0 -temperature [K]
!       dens0 -gas density [g/cm^3](note dust density is =0 in the code)
!       ane0 -electron number density [cm^-3]
!-----------------------------------------------------------------------                                                                       
!       If the objects happen to overlap, priority is given by the order
!       of 'if'-blocks in the subroutine smod1 and it is as follows:         
!               star,companion,spot,stream,ring,disc,envelope,nebula,
!               flow,jet,ufo,shell,and background.
!       temp and ane are assumed to have reasonable values all along 
!       the beam. An empty space can be defined as dens<denvac.         
!       Four types of nontransparent objects can be defined as:         
!       dcut1<dens<dcut2 -central star,                                 
!       dcut2<dens<dcut3 -secondary star(=companion), 
!       dcut3<dens<dcutn -3.body (it can be anything)
!       dcutn<dens -any opaque dark matter.                              
!         Note that lunt1, lunt2, lunt3 are in fact associated with 
!       density intervals (<dcut1,dcut2>, <dcut2,dcut3>, <dcut3,dcutn>) 
!       rather then with objects (star,companion,...) and thus can be 
!       used to ascribe the spectrum to any nontransparent object 
!       setting its density within a particular density interval.
!       However, limb darkening is applied to star and companion only
!       and it must be switched off (dlst=dlcp=0.) if you want to use 
!       these density intervals for other objects (without limb dark.).
!       Roche geometry assumes synchronous rotation around z axis with 
!       star in the center and companion at xcp>0,ycp=zcp=0 revolving 
!       towards (0,1,0).
!         Treatment of the scattered light assumes that the medium is 
!       optically thin, more precisely, that there is no significant 
!       obstruction between the source of the light (the two stars) 
!       and scattering medium. Only the radiative transfer along 
!       the line of sight is solved. It means that objects could make 
!       eclipses along the line of sight but would cast no shadows 
!       into other directions (i.e. would be transparent when 
!       considering scattered light from the two stars). 
!       To allow objects to cast shadows we introduced a 3D field 
!       where the user can specify whether a certain space point
!       is in the shadow so that only scattering from the unobscured 
!       source (star) is taken into account.
!         Input variables which are supposed to be components of a unit
!       vector do not need to be normalized.
        dcut1=0.5d15
        dcut2=1.5d15
        dcut3=2.5d15
        dcutn=3.5d15
        denvac=1.d-50
	open(9,file='shellspec.in',status='old')
	read(9,*)
	read(9,*)alam1,alamn,alams,loglam,cutoff
	read(9,*)
	read(9,*)imodel,irotat,ipart,ichemc,ielnd
	read(9,*)
	read(9,*)ithom,irayl,imie,imiepf,ihyd,iopac,iline,eps	
	read(9,*)
	read(9,*)ionu,ior,iot,offset
	read(9,*)
	read(9,*)phase1,phasen,nphase,dinc
	read(9,*)
	read(9,*)dd,vgamma
	read(9,*)
	read(9,*)iext,rv,ebv
	read(9,*)
	read(9,*)
	read(9,*)lunt1,xunt1,yunt1
	read(9,*)
	read(9,*)lunt2,xunt2,yunt2
	read(9,*)
	read(9,*)lunt3,xunt3,yunt3
	read(9,*)
	read(9,*)
	read(9,*)rmdfx1,rmdfx2,rmdfy1,rmdfy2,rmdfz1,rmdfz2
	read(9,*)
	read(9,*)stepfx,stepfy,stepfz,gainfx,gainfy,gainfz
	read(9,*)
	read(9,*)
	read(9,*)rmdfx3,rmdfx4,rmdfy3,rmdfy4,rmdfz3,rmdfz4
	read(9,*)
	read(9,*)stepfxb,stepfyb,stepfzb,gainfxb,gainfyb,gainfzb
	read(9,*)
	read(9,*)
	read(9,*)rmdx1,rmdx2,rmdy1,rmdy2,rmdz1,rmdz2
	read(9,*)
        read(9,*)stepx,stepy,stepz,gainx,gainy,gainz
	read(9,*)
	read(9,*)
	read(9,*)rmdx3,rmdx4,rmdy3,rmdy4,rmdz3,rmdz4
	read(9,*)
        read(9,*)stepxb,stepyb,stepzb,gainxb,gainyb,gainzb
        read(9,*)
        read(9,*)
!
        read(9,*)istar,icomp,ispot,ism,iring,idisc,ienv
	read(9,*)
        read(9,*)inebl,iflow,ijet,iufo,ishell
        read(9,*)
        read(9,*)
	read(9,*)rstar,tstar,emstar
        read(9,*)
        read(9,*)xstar,ystar,zstar,vrotst
        read(9,*)
        read(9,*)idifst,drotst,hst
        read(9,*)
        read(9,*)vxst,vyst,vzst
        read(9,*)
        read(9,*)dlst,dlst2,dgst,ffst
        read(9,*)
        read(9,*)irrst,ialbst,albst,htst,htsta        
	read(9,*)
	read(9,*)ispst,xspst,yspst,zspst,aspst,tspst
	read(9,*)
	read(9,*)
	read(9,*)rcp,tempcp,qq
	read(9,*)
	read(9,*)vrxcp,vrycp,vrzcp,vrotcp
	read(9,*)
	read(9,*)xcp,ycp,zcp
	read(9,*)
	read(9,*)vxcp,vycp,vzcp	
	read(9,*)
	read(9,*)dlcp,dlcp2,dgcp,ffcp
	read(9,*)
	read(9,*)irrcp,ialbcp,albcp,htcp,htcpa	
	read(9,*)
	read(9,*)
	read(9,*)vrxsp,vrysp,vrzsp,vrotsp,rsp
	read(9,*)
	read(9,*)xsp,ysp,zsp,vxsp,vysp,vzsp
	read(9,*)
	read(9,*)tempsp,denssp,anesp,vtrbsp,dstdsp,dsttsp
	read(9,*)
        read(9,*)
        read(9,*)v1sm,v2sm,r1sm,r2sm
        read(9,*)
        read(9,*)x1sm,y1sm,z1sm
        read(9,*)
        read(9,*)x2sm,y2sm,z2sm
        read(9,*)
        read(9,*)vxsm,vysm,vzsm
        read(9,*)
        read(9,*)xsm,ysm,zsm,psm
        read(9,*)
        read(9,*)tempsm,denssm,anesm,vtrbsm,edensm,dstdsm,dsttsm
	read(9,*)
	read(9,*)
	read(9,*)rrg,emrg
	read(9,*)
	read(9,*)b1rg,b2rg
	read(9,*)
	read(9,*)a1rg,a2rg,dr1rg,dr2rg
	read(9,*)
	read(9,*)xrg,yrg,zrg
	read(9,*)
	read(9,*)xpolrg,ypolrg,zpolrg
	read(9,*) 
	read(9,*)vxrg,vyrg,vzrg
	read(9,*) 
	read(9,*)temprg,densrg,anerg,vtrbrg,itrg
	read(9,*)
	read(9,*)edenrg,dstdrg,ede2rg,dst2rg,dsttrg
	read(9,*)	
	read(9,*)
        read(9,*)adisc,rindc,routdc,emdc,rdc
	read(9,*)
        read(9,*)xdc,ydc,zdc
	read(9,*)
        read(9,*)xdisc,ydisc,zdisc
        read(9,*)
        read(9,*)vxdc,vydc,vzdc
        read(9,*)
        read(9,*)tempdc,densdc,anedc,vtrbdc,edendc,itdc,etmpdc
        read(9,*)
        read(9,*)dstddc,dsttdc
	read(9,*)
	read(9,*)
	read(9,*)emen,qqen,aen,ffen,hen
	read(9,*)
	read(9,*)tempen,densen,aneen,vtrben,dstden,dstten        
	read(9,*)	
	read(9,*)
        read(9,*)aneb,rinnb,routnb,emnb,rnb
        read(9,*)
        read(9,*)iinvnb,hinvnb,tinvnb
	read(9,*)
        read(9,*)hwindnb,idennb,hcnb
        read(9,*)
        read(9,*)ivelnb,hvelnb,vnb,evelnb
	read(9,*)
	read(9,*)ishdnb,hshdnb
	read(9,*)
        read(9,*)xneb,yneb,zneb
        read(9,*)
        read(9,*)vxnb,vynb,vznb
        read(9,*)
        read(9,*)tempnb,densnb,anenb,vtrbnb,edennb,itnb,etmpnb
        read(9,*)
        read(9,*)dstdnb,dsttnb
        read(9,*)
        read(9,*)
        read(9,*)v1fw,v2fw,r1fw,r2fw
        read(9,*)
        read(9,*)x1fw,y1fw,z1fw
        read(9,*)
        read(9,*)x2fw,y2fw,z2fw
        read(9,*)
        read(9,*)vxfw,vyfw,vzfw
        read(9,*)
        read(9,*)xfw,yfw,zfw,pfw        
        read(9,*)
        read(9,*)tempfw,densfw,anefw,vtrbfw,edenfw,dstdfw,dsttfw
	read(9,*)
	read(9,*)
        read(9,*)ajet,rinjt,routjt
        read(9,*)
        read(9,*)xjt,yjt,zjt
        read(9,*)
        read(9,*)xjet,yjet,zjet
        read(9,*)
        read(9,*)ivjt,vjt,eveljt,rcjt
        read(9,*)
        read(9,*)vxjt,vyjt,vzjt
        read(9,*)
        read(9,*)tempjt,densjt,anejt,vtrbjt,dstdjt,dsttjt,etmpjt,asymjt
	read(9,*)
	read(9,*)
        read(9,*)aufo,rinuf,routuf,emuf,ruf
	read(9,*)
        read(9,*)xuf,yuf,zuf
	read(9,*)
        read(9,*)xufo,yufo,zufo
        read(9,*)
        read(9,*)vxuf,vyuf,vzuf
        read(9,*)
        read(9,*)tempuf,densuf,aneuf,vtrbuf,edenuf,ituf,etmpuf
        read(9,*)
        read(9,*)dstduf,dsttuf
	read(9,*)
	read(9,*)
	read(9,*)rinsh,routsh,vsh
	read(9,*)
	read(9,*)evelsh,rcsh
        read(9,*)
        read(9,*)vxsh,vysh,vzsh	
	read(9,*)
	read(9,*)tempsh,denssh,anesh,vtrbsh,dstdsh,dsttsh,etmpsh
        read(9,*)
        read(9,*)
        read(9,*)temp0,dens0,ane0,v0
!	only ipart=1 is supported at the moment
	ipart=1
	nfreq=nint((alamn-alam1)/alams)+1
	write(*,*)' nfreq=',nfreq
	nbodf1a=nint((rmdfx2-rmdfx1)/stepfx)+1
	nbodf1b=nint((rmdfx4-rmdfx3)/stepfxb)+1
	nbodf2a=nint((rmdfy2-rmdfy1)/stepfy)+1
	nbodf2b=nint((rmdfy4-rmdfy3)/stepfyb)+1
	nbodf3a=nint((rmdfz2-rmdfz1)/stepfz)+1
	nbodf3b=nint((rmdfz4-rmdfz3)/stepfzb)+1
	nbod1a=nint((rmdx2-rmdx1)/stepx)+1
	nbod1b=nint((rmdx4-rmdx3)/stepxb)+1
	nbod2a=nint((rmdy2-rmdy1)/stepy)+1
	nbod2b=nint((rmdy4-rmdy3)/stepyb)+1
	nbod3a=nint((rmdz2-rmdz1)/stepz)+1
	nbod3b=nint((rmdz4-rmdz3)/stepzb)+1
	if(nbodf1a.lt.3.or.nbodf2a.lt.3.or.nbodf3a.lt.3)then
          write(*,*)' error: nbodfa, stop'
          write(3,*)' 1 error: nbodfa, stop'
	   call exit(1)
        endif
	if(nbod1a.lt.3.or.nbod2a.lt.3.or.nbod3a.lt.3)then
          write(*,*)' error: nboda, stop'
          write(3,*)' 1 error: nboda, stop'
	   call exit(1)
        endif
!       make the number of points odd
	if(nbodf1a/2.eq.(nbodf1a+1)/2)nbodf1a=nbodf1a-1
	if(nbodf1b/2.eq.(nbodf1b+1)/2)nbodf1b=nbodf1b-1
	if(nbodf2a/2.eq.(nbodf2a+1)/2)nbodf2a=nbodf2a-1
	if(nbodf2b/2.eq.(nbodf2b+1)/2)nbodf2b=nbodf2b-1
	if(nbodf3a/2.eq.(nbodf3a+1)/2)nbodf3a=nbodf3a-1
	if(nbodf3b/2.eq.(nbodf3b+1)/2)nbodf3b=nbodf3b-1
	if(nbod1a/2.eq.(nbod1a+1)/2)nbod1a=nbod1a-1
	if(nbod1b/2.eq.(nbod1b+1)/2)nbod1b=nbod1b-1
	if(nbod2a/2.eq.(nbod2a+1)/2)nbod2a=nbod2a-1
	if(nbod2b/2.eq.(nbod2b+1)/2)nbod2b=nbod2b-1
	if(nbod3a/2.eq.(nbod3a+1)/2)nbod3a=nbod3a-1
	if(nbod3b/2.eq.(nbod3b+1)/2)nbod3b=nbod3b-1
	if(nbodf1b.lt.3)nbodf1b=0
	if(nbodf2b.lt.3)nbodf2b=0	
	if(nbodf3b.lt.3)nbodf3b=0
	if(nbod1b.lt.3)nbod1b=0	
	if(nbod2b.lt.3)nbod2b=0	
	if(nbod3b.lt.3)nbod3b=0	
	nbodf1=nbodf1a+nbodf1b
	nbodf2=nbodf2a+nbodf2b
	nbodf3=nbodf3a+nbodf3b
	nbod1=nbod1a+nbod1b
	nbod2=nbod2a+nbod2b
	nbod3=nbod3a+nbod3b
	write(*,*)'nbodf1a nbodf2a nbodf3a  (1st body frozen grid)'
	write(*,'(4i5)')nbodf1a,nbodf2a,nbodf3a
	write(*,*)'nbodf1b nbodf2b nbodf3b  (2nd body frozen grid)'
	write(*,'(4i5)')nbodf1b,nbodf2b,nbodf3b
	write(*,*)'nbodf1  nbodf2  nbodf3   (sum body frozen grid)'
	write(*,'(4i5)')nbodf1,nbodf2,nbodf3
	write(*,*)'nbod1a nbod2a nbod3a     (1st line of sight grid)'
	write(*,'(4i5)')nbod1a,nbod2a,nbod3a
	write(*,*)'nbod1b nbod2b nbod3b     (2nd line of sight grid)'
	write(*,'(4i5)')nbod1b,nbod2b,nbod3b
	write(*,*)'nbod1  nbod2  nbod3      (sum line of sight grid)'
	write(*,'(4i5)')nbod1,nbod2,nbod3
	if(nbodf1.gt.ndimf1.or.nbodf2.gt.ndimf2.or.nbodf3.gt.ndimf3)then
	  write(*,*)' warning: No. grid points may > dimension, check'
	  write(*,*)'ndimf1 ndimf2 ndimf3  (dimension body frozen grid)'
	  write(*,'(4i5)')ndimf1,ndimf2,ndimf3
	  write(3,*)' 1 warning: No. grid points may > dimension, check'
	  write(3,*)'ndimf1 ndimf2 ndimf3  (dimension body frozen grid)'
	  write(3,'(4i5)')ndimf1,ndimf2,ndimf3
!	  nbodf can be < nbodfa+nbodfb that is why only a warning here
!	  call exit(1)
	endif
	if(nbod1.gt.ndim1.or.nbod2.gt.ndim2.or.nbod3.gt.ndim3)then
	  write(*,*)' warning: No. grid points may > dimension, check'
 	  write(*,*)'ndim1 ndim2 ndim3  (dimension line of sight grid)'
	  write(*,'(4i5)')ndim1,ndim2,ndim3
	  write(3,*)' 1 warning: No grid points may > dimension, check'
 	  write(3,*)'ndim1 ndim2 ndim3  (dimension line of sight grid)'
	  write(3,'(4i5)')ndim1,ndim2,ndim3
!	  nbod can be < nboda+nbodb that is why only a warning here
!	  call exit(1)
	endif	
	if(loglam.ne.2)then
	if(nfreq.gt.mfreq)then
	  write(*,*)' error: nfreq>mfreq dimension exceeded, stop'
	  write(*,*)' nfreq, mfreq=',nfreq,mfreq
	  write(3,*)' 1 error: nfreq>mfreq dimension exceeded, stop'
	  write(3,*)' nfreq, mfreq=',nfreq,mfreq
	  call exit(1)
	endif		
	endif		
!	create the set of lambda's	
	if(loglam.eq.0)then
	  do i=1,nfreq
	    alam(i)=alam1+dble(i-1)*alams
          enddo
        elseif(loglam.eq.1)then
          alams=(dlog10(alamn)-dlog10(alam1))/dble(nfreq-1)
          do i=1,nfreq
             alam(i)=dlog10(alam1)+dble(i-1)*alams
             alam(i)=10**(alam(i))
          enddo   
	else
!         including loglam.eq.2
	  write(*,*)' reading file lambda'
	  open(15,file='lambda',status='old')
	  i=1
11	  read(15,*,err=13,end=13)alam(i)
	  i=i+1
	  if(i.gt.mfreq)then
	    write(*,*)' max. nfreq i.e. mfreq reached'
	    goto 13
	  endif		
	  goto 11
13	  nfreq=i-1
	  close(15)
	  write(*,*)' nfreq=',nfreq
        endif
!       check for extinction correction range
        if(iext.eq.1)then
        if(alam(1).lt.1.d3.or.alam(nfreq).gt.1.d5/3.d0)then
	  write(*,*)' error: alam out of extinction corr. range, stop'
	  write(3,*)' 1 error: alam out of extinction corr. range, stop'
          call exit(1)
        endif
        endif        
!	phases        
	if(nphase.eq.0)then
	  write(*,*)' reading file phases'
	  open(15,file='phases',status='old')
	  i=1
15	  read(15,*,err=17,end=17)alpha(i)
	  alpha(i)=(alpha(i)-0.25d0)*2.d0*pi
	  i=i+1
          if (i-1.gt.mphase) then
	    write(*,*)' error: number of phases .gt. mphase = ', mphase
	    write(3,*)' 1 error: number of phases .gt. mphase = ', mphase
            call exit(1)
          endif
	  goto 15
17	  nphase=i-1
	  close(15)
	elseif(nphase.eq.1)then
	  alpha(1)=phase1/180.d0*pi
	elseif(nphase.gt.1)then
	  alphas=(phasen-phase1)/dble(nphase-1)
	  do i=1,nphase
	    alpha(i)=(phase1+dble(i-1)*alphas)/180.d0*pi
	  enddo
	endif
	write(*,*)' nphase=',nphase
!	conversion of time
	psm=psm*24.d0*3600.d0
	pfw=pfw*24.d0*3600.d0
!	conversion of masses	
	emstar=emstar*emsol
	emen=emen*emsol
	emdc=emdc*emsol
	emnb=emnb*emsol
	emuf=emuf*emsol
	emrg=emrg*emsol
!	conversion of length	
	rstar=rstar*rsol
	rcp=rcp*rsol
        xcp=xcp*rsol
        ycp=ycp*rsol
        zcp=zcp*rsol
        aen=aen*rsol
        hen=hen*rsol
	rsp=rsp*rsol
        xsp=xsp*rsol
        ysp=ysp*rsol
        zsp=zsp*rsol
	rinsh=rinsh*rsol
	routsh=routsh*rsol
	rcsh=rcsh*rsol
	rindc=rindc*rsol
	routdc=routdc*rsol
	rdc=rdc*rsol
	xdc=xdc*rsol
	ydc=ydc*rsol
	zdc=zdc*rsol
	rinnb=rinnb*rsol
	routnb=routnb*rsol
	rnb=rnb*rsol
	rinuf=rinuf*rsol
	routuf=routuf*rsol
	ruf=ruf*rsol
	xuf=xuf*rsol
	yuf=yuf*rsol
	zuf=zuf*rsol
	rinjt=rinjt*rsol
	routjt=routjt*rsol
	xjt=xjt*rsol
	yjt=yjt*rsol
	zjt=zjt*rsol
	rcjt=rcjt*rsol
	x1sm=x1sm*rsol
	y1sm=y1sm*rsol
	z1sm=z1sm*rsol
	x2sm=x2sm*rsol
	y2sm=y2sm*rsol
	z2sm=z2sm*rsol
	r1sm=r1sm*rsol
	r2sm=r2sm*rsol
	x1fw=x1fw*rsol
	y1fw=y1fw*rsol
	z1fw=z1fw*rsol
	x2fw=x2fw*rsol
	y2fw=y2fw*rsol
	z2fw=z2fw*rsol
	r1fw=r1fw*rsol
	r2fw=r2fw*rsol
	rrg=rrg*rsol
	a1rg=a1rg*rsol
	a2rg=a2rg*rsol
	dr1rg=dr1rg*rsol
	dr2rg=dr2rg*rsol
        xrg=xrg*rsol
        yrg=yrg*rsol
        zrg=zrg*rsol
        rmdfx1=rmdfx1*rsol
	rmdfx2=rmdfx2*rsol
	rmdfx3=rmdfx3*rsol
	rmdfx4=rmdfx4*rsol
	rmdfy1=rmdfy1*rsol
	rmdfy2=rmdfy2*rsol
	rmdfy3=rmdfy3*rsol
	rmdfy4=rmdfy4*rsol
	rmdfz1=rmdfz1*rsol
	rmdfz2=rmdfz2*rsol
	rmdfz3=rmdfz3*rsol
	rmdfz4=rmdfz4*rsol
	rmdx1=rmdx1*rsol
	rmdx2=rmdx2*rsol
	rmdx3=rmdx3*rsol
	rmdx4=rmdx4*rsol
	rmdy1=rmdy1*rsol
	rmdy2=rmdy2*rsol
	rmdy3=rmdy3*rsol
	rmdy4=rmdy4*rsol
	rmdz1=rmdz1*rsol
	rmdz2=rmdz2*rsol	
	rmdz3=rmdz3*rsol
	rmdz4=rmdz4*rsol
!	conversion of velocities
	vgamma=1.d5*vgamma
	vrotst=1.d5*vrotst
	vxst=1.d5*vxst
        vyst=1.d5*vyst
        vzst=1.d5*vzst
	vxcp=1.d5*vxcp
        vycp=1.d5*vycp
        vzcp=1.d5*vzcp
        vrotcp=1.d5*vrotcp
        vtrben=1.d5*vtrben
	vxsp=1.d5*vxsp
        vysp=1.d5*vysp
        vzsp=1.d5*vzsp
        vrotsp=1.d5*vrotsp
	vtrbsp=1.d5*vtrbsp
        vxrg=1.d5*vxrg        
        vyrg=1.d5*vyrg
        vzrg=1.d5*vzrg
	vtrbrg=1.d5*vtrbrg
	vxdc=1.d5*vxdc
        vydc=1.d5*vydc
        vzdc=1.d5*vzdc
	vtrbdc=1.d5*vtrbdc
	vnb=1.d5*vnb
	vxnb=1.d5*vxnb
        vynb=1.d5*vynb
        vznb=1.d5*vznb
	vtrbnb=1.d5*vtrbnb
	vxuf=1.d5*vxuf
        vyuf=1.d5*vyuf
        vzuf=1.d5*vzuf
	vtrbuf=1.d5*vtrbuf
	vjt=1.d5*vjt
	vxjt=1.d5*vxjt
	vyjt=1.d5*vyjt
	vzjt=1.d5*vzjt
	vtrbjt=1.d5*vtrbjt
	vsh=1.d5*vsh
	vtrbsh=1.d5*vtrbsh
	vxsh=1.d5*vxsh
	vysh=1.d5*vysh
	vzsh=1.d5*vzsh
	v0=1.d5*v0
	v1sm=1.d5*v1sm
	v2sm=1.d5*v2sm
	vtrbsm=1.d5*vtrbsm
	vxsm=1.d5*vxsm
	vysm=1.d5*vysm
	vzsm=1.d5*vzsm
	v1fw=1.d5*v1fw
	v2fw=1.d5*v2fw
	vtrbfw=1.d5*vtrbfw
	vxfw=1.d5*vxfw
	vyfw=1.d5*vyfw
	vzfw=1.d5*vzfw
!	conversion of angles
	aspst=aspst/180.d0*pi
	dinc=dinc/180.d0*pi
	b1rg=b1rg/180.d0*pi
	b2rg=b2rg/180.d0*pi
	if(idisc.eq.1)then
	  adisc=adisc/180.d0*pi
	else
 	  adisc=adisc*rsol
	endif
	if(iufo.eq.1)then
	  aufo=aufo/180.d0*pi
	else
 	  aufo=aufo*rsol
	endif	
	ajet=ajet/180.d0*pi
	dd=dd*pc
	close(9)
!	shift lambda's by the gamma velocity
        do i=1,nfreq
          alam(i)=alam(i)-vgamma/clight*alam(i)
        enddo
!	check for most common errors in input	
        if(istar.gt.1.or.icomp.gt.1)then
          if(xcp.le.0.d0.or.qq.le.0.d0.or.emstar.le.0.d0)then
            write(*,*)' error: inconsistent star/companion input data'
            write(3,*)' 1 error: inconsistent star/companion input data'
	    call exit(1)
          endif
        endif
	if(densdc.gt.dcut1.and.edendc.ne.0.d0)then
          write(*,*)' error: wrong densdc or edendc, stop'
          write(3,*)' 1 error: wrong densdc or edendc, stop'
	  call exit(1)
	endif
	if(densnb.gt.dcut1.and.edennb.ne.0.d0)then
          write(*,*)' error: wrong densnb or edennb, stop'
          write(3,*)' 1 error: wrong densnb or edennb, stop'
	  call exit(1)
	endif
	if(densuf.gt.dcut1.and.edenuf.ne.0.d0)then
          write(*,*)' error: wrong densuf or edenuf, stop'
          write(3,*)' 1 error: wrong densuf or edenuf, stop'
	  call exit(1)
	endif	
	if(ishell.gt.1.and.denssh.gt.dcut1)then
          write(*,*)' error: wrong ishell or denssh, stop'
          write(3,*)' 1 error: wrong ishell or denssh, stop'
	  call exit(1)
	endif
        if((rstar.le.0.d0.or.tstar.le.0.d0)                             &
     &  .and.(istar.gt.0.or.icomp.eq.2))then
          write(*,*)' error: rstar-tstar, stop'
          write(3,*)' 1 error: rstar-tstar, stop'
	  call exit(1)
	endif
        if((rcp.le.0.d0.or.tempcp.le.0.d0)                              &
     &  .and.(istar.gt.0.or.icomp.gt.0))then
          write(*,*)' error: rcp-tempcp, stop'
          write(3,*)' 1 error: rcp-tempcp, stop'
	  call exit(1)
	endif	
	if(rindc.le.rdc.and.itdc.eq.2.and.idisc.gt.0)then
          write(*,*)' error: rdc, rindc, itdc, stop (T-singularity)'
          write(3,*)' 1 error: rdc, rindc, itdc, stop (T-singularity)'
	  call exit(1)
	endif
	if(rinuf.le.ruf.and.ituf.eq.2.and.iufo.gt.0)then
          write(*,*)' error: ruf, rinuf, ituf, stop (T-singularity)'
          write(3,*)' 1 error: ruf, rinuf, ituf, stop (T-singularity)'
	  call exit(1)
	endif	
	if(rmdx2.le.rmdx1.or.rmdy2.le.rmdy1.or.rmdz2.le.rmdz1)then
          write(*,*)' error: rmd(xyz)(21)'
          write(3,*)' 1 error: rmd(xyz)(21)'
	  call exit(1)
	endif	
	if(rmdx4.lt.rmdx3.or.rmdy4.lt.rmdy3.or.rmdz4.lt.rmdz3)then
          write(*,*)' error: rmd(xyz)(43)'
          write(3,*)' 1 error: rmd(xyz)(43)'
	  call exit(1)
	endif	
	if(rmdfx2.le.rmdfx1.or.rmdfy2.le.rmdfy1.or.rmdfz2.le.rmdfz1)then
          write(*,*)' error: rmdf(xyz)(21)'
          write(3,*)' 1 error: rmdf(xyz)(21)'
	  call exit(1)
	endif		
	if(rmdfx4.lt.rmdfx3.or.rmdfy4.lt.rmdfy3.or.rmdfz4.lt.rmdfz3)then
          write(*,*)' error: rmdf(xyz)(43)'
          write(3,*)' 1 error: rmdf(xyz)(43)'
	  call exit(1)
	endif		
	if(vrxsp**2+vrysp**2+vrzsp**2.le.0.d0)then
	  write(*,*)' error: unit vector'
	  write(3,*)' 1 error: unit vector'
	  call exit(1)
	endif		
	if((dlst+dlst2).gt.1.d0)then
	  write(*,*)' error: limb darkening'
	  write(3,*)' 1 error: limb darkening'
	  call exit(1)
	endif		
	if((dlcp+dlcp2).gt.1.d0)then
	  write(*,*)' error: limb darkening'
	  write(3,*)' 1 error: limb darkening'
	  call exit(1)
	endif		
	if(v1sm.le.0.d0.or.v2sm.le.0.d0)then
	  write(*,*)' error: v1sm or v2sm<=0'
	  write(3,*)' 1 error: v1sm or v2sm<=0'
	  call exit(1)
	endif  
	if(imiepf.gt.0.and.imie.lt.1)then
	  write(*,*)' error: Mie scattering'
	  write(3,*)' 1 error: Mie scattering'
	  call exit(1)
	endif
	if(ivjt.eq.2.and.rcjt.ge.rinjt)then
	  write(*,*)' error: jet density singularity'
	  write(3,*)' 1 error: jet density sigularity'
	  call exit(1)
	endif
	if(ishell.eq.3.and.rcsh.ge.rinsh)then
	  write(*,*)' error: shell density singularity'
	  write(3,*)' 1 error: shell density sigularity'
	  call exit(1)
	endif
	if(v1fw*v2fw.lt.0.d0.or.v1sm*v2sm.lt.0.d0)then
	  write(*,*)' error: wrong stream velocities'
	  write(3,*)' 1 error: wrong stream velocities'
	  call exit(1)
	endif
	if(irotat.eq.0)then
	  write(*,*)' warning: irotat=0 does not support shadows'
	endif		
!	reads monochromatic albedos	
	nalb1=0
	nalb2=0
	if(ialbst.eq.1.and.irrst.eq.1)then
	  call albedo(alb1x,alb1y,nalb1,1)
	  write(*,*)' albedo1',nalb1
	endif
	if(ialbcp.eq.1.and.irrcp.eq.1)then
	  call albedo(alb2x,alb2y,nalb2,2)
	  write(*,*)' albedo2',nalb2	  
	endif	
        if(albcp.lt.0.d0.or.albcp.gt.1.d0.or.                           &
     &  htcp.lt.0.d0.or.htcp.gt.1.d0.or.                                &
     &  htcpa.lt.0.d0.or.htcpa.gt.1.d0)then
          write(*,*)' error: alb**, ht**, stop (albedo-heat transport)'
          write(3,*)' 1 error: alb**, ht**, stop (albedo-heat transp.)'
	  call exit(1)
	endif
        if(albst.lt.0.d0.or.albst.gt.1.d0.or.                           &
     &  htst.lt.0.d0.or.htst.gt.1.d0.or.                                &
     &  htsta.lt.0.d0.or.htsta.gt.1.d0)then
          write(*,*)' error: alb**, ht**, stop (albedo-heat transport)'
          write(3,*)' 1 error: alb**, ht**, stop (albedo-heat transp.)'
	  call exit(1)
	endif	
!	reads opacities and phase functions for Mie scattering on dust
	if(imie.gt.0)then
          call mie(opmiex,opmies,opmiea,nmie,pfmiex,pfang,pfmie         &
     &    ,npfmie,imie,imiepf,ndust,dtlow,dthig,drmf)
	endif
!		xsec reads extra gas opacity tables
!       xsecx - freq [Hz]
!       xsecy - cross-section [cm^2]
!	amixr -molecule abundance relative to H nuclei
!		chem reads chemistry table
!	pop -dlog10 of molecule population
	if(iopac.eq.1)then
          call xsec(alam(1),alam(nfreq),cutoff,xsecx,xsecy,txsec        &
     &    ,nxsec,ntxsec,amixr)
          call chem(pop,popt,popd,mtemp,mdens,melm)
	endif	
!		reading the input spectra of nontransparent objects
!	if lunt123>0 and converts them to central intensity
        call untsp(xstar1,star1,xstar2,star2,xstar3,star3               &
     &  ,wstar1,fstar1,tstar,jstar1,wstar2,fstar2,tempcp,jstar2         &
     &  ,xunt1,yunt1,xunt2,yunt2,xunt3,yunt3                            &
     &  ,nstar1,nstar2,nstar3,nspec1,nspec2,tspec1,tspec2               &
     &  ,lunt1,lunt2,lunt3,dlst,dlst2,dlcp,dlcp2,alam,nfreq)
!	reads vertical density profile for nebula
	ndennb=0
	if(idennb.eq.1)then
	  write(*,*)' reading file wind_prof'
	  open(17,file='wind_prof',status='old')
	  read(17,*)
	  read(17,*)
	  i=1
30	  read(17,*,end=40,err=40)denxnb(i),pom,pom,denznb(i)
	  i=i+1
	  goto 30
40	  ndennb=i-1	  
	  write(*,'(a,i6)')' wind_prof:ndennb=',ndennb
	  close(17)	
	endif
!-----------------------------------------------------------------------
!	The model (input) is defined in its body frozen cartesian 
!	coordinates (far,fat,faz). These may rotate relative to your 
!	steady cartesian laboratory frame (ar,at,az) -line of sight 
!	grid. The later has the same center of coordinates but its	
!	'z' axis is always pointed towards the observer so that 
!	the radial velocity of approaching object is positive. 
!	Radiative transfer is then solved in this line of sight grid.
!	Model is first inclined around 'ar' by dinc degrees and 
!	then rotated around 'faz'-its rotational axis.
!	Nbod3 should be large enough so that a huge leaps in opt. depth
!	do not occure because of large velocity and opacity change 
!	between the grid points.
!	ftemp,fdens,fne,fvr,fvt,fvz,fvtrb,fdustd,fdustt,kshade are 
!	the following quantities in body frozen grid:
!	  ftemp -gas temperature
!	  fdens -gas density
!	  fne -electron number density 
!	  fvr -x velocity
!	  fvt -y velocity
!	  fvz -z velocity (radial velocity)
!	  fvtrb -turbulent velocity 
!	  fdustd -density of dust
!         fdustt -dust temperature
!         kshade -integer describing shadows (scattering from 2 stars)
!	atemp,adens,ane,avr,avt,avz,avtrb,adustd,adustt,lshade are 
!	corresponding quantities in the line of sight grid.
!	Body frozen grid is defined in subroutine smod1.
	if (imodel.eq.1) then
         call smod1(nbodf1,nbodf1a,nbodf1b                              &
     &   ,nbodf2,nbodf2a,nbodf2b,nbodf3,nbodf3a,nbodf3b                 &
     &   ,rmdfx1,rmdfx2,rmdfx3,rmdfx4,rmdfy1,rmdfy2,rmdfy3,rmdfy4       &
     &   ,rmdfz1,rmdfz2,rmdfz3,rmdfz4                                   &
     &   ,gainfx,gainfy,gainfz,gainfxb,gainfyb,gainfzb                  &
     &   ,rstar,tstar,emstar,xstar,ystar,zstar,vrotst,idifst,drotst,hst &
     &   ,istar,vxst,vyst,vzst,dgst,ffst,irrst,albst,htst,htsta         &
     &   ,ispst,xspst,yspst,zspst,aspst,tspst                           &
     &   ,icomp,dgcp,ffcp,qq,vrxcp,vrycp,vrzcp,vrotcp,rcp               &
     &   ,xcp,ycp,zcp,vxcp,vycp,vzcp,tempcp,irrcp,albcp,htcp,htcpa      &
     &   ,ienv,emen,qqen,aen,ffen,hen                                   &
     &   ,tempen,densen,aneen,vtrben,dstden,dstten                      &
     &   ,ispot,vrxsp,vrysp,vrzsp,vrotsp,rsp                            &
     &   ,xsp,ysp,zsp,vxsp,vysp,vzsp,tempsp,denssp,anesp,vtrbsp         &
     &   ,dstdsp,dsttsp                                                 &
     &   ,iring,rrg,emrg,b1rg,b2rg,a1rg,a2rg,dr1rg,dr2rg                &
     &   ,xrg,yrg,zrg,xpolrg,ypolrg,zpolrg,vxrg,vyrg,vzrg               &
     &   ,temprg,densrg,anerg,vtrbrg,itrg                               &
     &   ,edenrg,dstdrg,ede2rg,dst2rg,dsttrg                            &
     &   ,idisc,adisc,rindc,routdc,emdc,rdc                             &
     &   ,xdc,ydc,zdc,xdisc,ydisc,zdisc,vxdc,vydc,vzdc                  &
     &   ,tempdc,densdc,anedc,vtrbdc,edendc,itdc,etmpdc,dstddc,dsttdc   &
     &   ,inebl,aneb,rinnb,routnb,emnb,rnb,hcnb,ishdnb,hshdnb           &
     &   ,hinvnb,tinvnb,iinvnb,hwindnb,ndennb,denxnb,denznb             &
     &   ,ivelnb,hvelnb,vnb,evelnb,xneb,yneb,zneb,vxnb,vynb,vznb        &
     &   ,tempnb,densnb,anenb,vtrbnb,edennb,itnb,etmpnb,dstdnb,dsttnb   &
     &   ,x1sm,y1sm,z1sm,x2sm,y2sm,z2sm,v1sm,v2sm,r1sm,r2sm             &
     &   ,ism,vxsm,vysm,vzsm,xsm,ysm,zsm,psm                            &
     &   ,tempsm,denssm,anesm,vtrbsm,edensm,dstdsm,dsttsm               &
     &   ,x1fw,y1fw,z1fw,x2fw,y2fw,z2fw,v1fw,v2fw,r1fw,r2fw             &
     &   ,iflow,vxfw,vyfw,vzfw,xfw,yfw,zfw,pfw                          &
     &   ,tempfw,densfw,anefw,vtrbfw,edenfw,dstdfw,dsttfw               &
     &   ,iufo,aufo,rinuf,routuf,emuf,ruf                               &
     &   ,xuf,yuf,zuf,xufo,yufo,zufo,vxuf,vyuf,vzuf                     &
     &   ,tempuf,densuf,aneuf,vtrbuf,edenuf,ituf,etmpuf,dstduf,dsttuf   &
     &   ,ajet,rinjt,routjt,xjt,yjt,zjt,xjet,yjet,zjet                  &
     &   ,ivjt,vjt,eveljt,rcjt,vxjt,vyjt,vzjt                           &
     &   ,tempjt,densjt,anejt,vtrbjt,ijet,dstdjt,dsttjt,etmpjt,asymjt   &
     &   ,rinsh,routsh,vsh,vxsh,vysh,vzsh                               &
     &   ,tempsh,denssh,anesh,vtrbsh,ishell,evelsh,rcsh,dstdsh,dsttsh   &
     &   ,etmpsh                                                        &
     &   ,v0,temp0,dens0,ane0,dcut1,dcut2,dcut3,dcutn                   &
     &   ,far,fat,faz,ftemp,fdens,fne,fvr,fvt,fvz,fvtrb,fdustd,fdustt   &
     &   ,kshade)
	endif
!		reading the model and bf grid from the file
	if (imodel.eq.2) then
	  write(*,*)' reading file shellspec.mod'
	  open(10,file='shellspec.mod',status='old')
          call smod2(ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3          &
     &    ,far,fat,faz,ftemp,fdens,fne,fvr,fvt,fvz,fvtrb,fdustd,fdustt  &
     &    ,kshade)
	  close(10)
        endif
	write(*,*)'nbodf1 nbodf2 nbodf3  (body frozen grid)'
	write(*,'(3i5)')nbodf1,nbodf2,nbodf3
!		definition of the line of sight grid
        call grid3d(ndim1,ndim2,ndim3                                   &
     &  ,nbod1,nbod1a,nbod1b                                            &
     &  ,nbod2,nbod2a,nbod2b                                            &
     &  ,nbod3,nbod3a,nbod3b                                            &          
     &  ,rmdx1,rmdx2,rmdx3,rmdx4                                        &
     &  ,rmdy1,rmdy2,rmdy3,rmdy4                                        &     
     &  ,rmdz1,rmdz2,rmdz3,rmdz4                                        &          
     &  ,gainx,gainxb,gainy,gainyb,gainz,gainzb,ar,at,az)     
	write(*,*)'nbod1 nbod2 nbod3  (line of sight grid)'
        write(*,'(3i5)')nbod1,nbod2,nbod3
!		reading the atomic data for spectral lines
	if(iline.eq.1)then
          call lindat(iat,iion,wlab,elo,eup,glo,gr0,gs0,gw0,bij         &
     &    ,nline)
        else
          nline=0
        endif
!		reads the atomic data for all elements and ions
!		dyp(*)-name of the element
!		d(1,*)-atom weight
!		d(2,*)-solar abundance relative to hydrogen
!		d(3,*)-number of ions considered
!		xi(*,*)-first 8 ionization potentials
	call state0(dyp,d,xi)
!		extract data for our particular elements of sp.lines 
!		(hmc,nion,abund,zip) and change abundances: 
        call eldat(ichemc,iat,iion,dyp,d,xi,abhyd,abhel,wm              &
     &  ,hmc,nion,abund,zip,nline,ielnd,necod)
!		calculation of the electron number density
        if(ielnd.eq.1) call elnd(d,xi,necod,wm,abhyd,ftemp,fdens,fne    &
     &  ,nbodf1,nbodf2,nbodf3,denvac,dcut1,ane0)
	write(2,130)nbod1,nbod2,nbod3,nfreq
130	format(' nbod1=',i4,' nbod2=',i4,' nbod3=',i4,' nfreq=',i4)
        write(2,'(3(a,i3))')' output for:  ionu=',ionu                  &
     &  ,' ior=',ior,' iot=',iot
!		area(i,j)-area of a projected surface element
	call surf(ndim1,ndim2,nbod1,nbod2,ar,at,area)
!-----------------------------------------------------------------------
!		grand cycle through: angles, x,y, frequencies, and z
!	aint0-incident intesity from behind
	do iang=1,nphase
          write(6,'(2(a,f7.2))')'i=',dinc*180.d0/pi                     &
     &    ,'  alpha=',alpha(iang)*180.d0/pi
          write(2,'(2(a,f7.2))')'i=',dinc*180.d0/pi                     &
     &    ,'  alpha=',alpha(iang)*180.d0/pi
!		cycle through parallel lines of sight (x,y)
	  opdepm=0.d0
	  do 430 ii=1,nbod1
	  do 420 jj=1,nbod2
!	    rotation of the object or ray casting:
!		linear interpolation 
!             (does not support shadows and sets lshade=3)
	    if(irotat.eq.0)then
              call rot1d1(dinc,alpha(iang),dcut1,temp0,ane0,ii,jj       &
     &        ,vxst,vyst,vzst,vxstr,vystr,vzstr                         &
     &        ,xcp,ycp,zcp,xcpr,ycpr,zcpr                               &
     &        ,vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr                         &
     &        ,ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3                &
     &        ,ndim1,ndim2,ndim3,nbod1,nbod2,nbod3                      &
     &        ,ar,at,az,far,fat,faz                                     &
     &        ,atemp,adustt,adustd,adens,ane,avr,avt,avz,avtrb,lshade   &
     &        ,ftemp,fdustt,fdustd,fdens,fne,fvr,fvt,fvz,fvtrb,kshade)
	    else
!		nearest neighbour approximation (supports shadows)
              call rot1d2(dinc,alpha(iang),dcut1,temp0,ane0,ii,jj       &
     &        ,vxst,vyst,vzst,vxstr,vystr,vzstr                         &
     &        ,xcp,ycp,zcp,xcpr,ycpr,zcpr                               &
     &        ,vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr                         &
     &        ,ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3                &
     &        ,ndim1,ndim2,ndim3,nbod1,nbod2,nbod3                      &
     &        ,ar,at,az,far,fat,faz                                     &
     &        ,atemp,adustt,adustd,adens,ane,avr,avt,avz,avtrb,lshade   &
     &        ,ftemp,fdustt,fdustd,fdens,fne,fvr,fvt,fvz,fvtrb,kshade)
	    endif
	    if(ii.eq.ior.and.jj.eq.iot)then
	      iprint=1
	    else
	      iprint=0
	    endif		
!		check for an empty space
	    ivac1=0
	    aint0=0.d0	
	    do kk=1,nbod3
	      if(adens(kk).gt.denvac)then
		ivac1=kk
		goto 150
	      endif
	    enddo
150         ivacn=0
	    do kk=nbod3,1,-1
	      if(adens(kk).gt.denvac)then
		ivacn=kk
		goto 170
	      endif	
            enddo
170	    if(ivac1.eq.0)then
	      do nu=1,nfreq
	        aint(ii,jj,nu)=aint0
	      enddo
	      if(ii.eq.ior.and.jj.eq.iot)then
                write(2,*)' ii=',ii,'   jj=',jj,                        &
     & 	        '   this ray goes through an empty space'
	      endif	
	      goto 420	
	    else
	      do kk=ivac1,ivacn
	        temp(kk-ivac1+1)=atemp(kk)
	        ekonc(kk-ivac1+1)=ane(kk)
	        akonc(kk-ivac1+1)=adens(kk)/wm/hjed
	        dens(kk-ivac1+1)=adens(kk)
	        dustd(kk-ivac1+1)=adustd(kk)
	        dustt(kk-ivac1+1)=adustt(kk)
	        azn(kk-ivac1+1)=az(kk)
	        vr(kk-ivac1+1)=avr(kk)
	        vt(kk-ivac1+1)=avt(kk)
	        vz(kk-ivac1+1)=avz(kk)
	        vtrb(kk-ivac1+1)=avtrb(kk)
		mshade(kk-ivac1+1)=lshade(kk)
	      enddo
	      nivac=ivacn-ivac1+1
	    endif
!
!		you can possibly speed up if you seach for 
!		untransparent object before the above cycle??
!           search for the last untransparent point along the ray
!	    and calculation of the limb darkening factor (iunt,darkl)	
            call bcond1(iunt,darkl,cdelta,vdif,nivac,dens               &
     &      ,dcut1,dcut2,dcut3,dcutn,vr,vt,vz                           &
     &      ,ar(ii),at(jj),azn,dinc,alpha(iang)                         &
     &      ,istar,rstar,vxstr,vystr,vzstr,dlst,dlst2,irrst             &
     &      ,icomp,rcp,xcpr,ycpr,zcpr,vxcpr,vycpr,vzcpr                 &
     &      ,dlcp,dlcp2,irrcp,xcp,qq,iprint)
!
!		hneg - H- hydrogen number density
!		h2 - H2 number density 
!	    	hi -neutral hydrogen number density
!	    	hei -neutral helium number density
!	    	hipf -HI partition function	
!		akonc -number density of all atoms
            call hihei(ipart,nivac,temp,ekonc,akonc                     &
     &      ,hipf,hi,hneg,h2,hei,abhyd,abhel)
!		interpolates amixf from the pop table
            if(iopac.eq.1)then
              call abnd(pop,popt,popd,mtemp,mdens,melm                  &
     &        ,temp,dens,nivac,ndim3,iopac,amixr,amixf)
            endif
	    if(loglam.eq.0)then
!             ophyd calculates at alam1,alamn:     
!	    	ophbf -HI bound-free opacity
!	    	ophff -HI free-free opacity
!	    	ophrs -HI Rayleigh scattering
!		ophn  -H(-) opacity
              call ophyd(nivac,temp,ekonc,akonc,hi,hneg,abhyd           &
     &        ,ophbf1,ophff1,ophrs1,ophn1,hipf,alam1)
              call ophyd(nivac,temp,ekonc,akonc,hi,hneg,abhyd           &
     &        ,ophbf2,ophff2,ophrs2,ophn2,hipf,alamn)
            endif
	    if (iprint.eq.1)then        
              write(2,*)' iunt=',iunt
              write(2,*)' ivac1=',ivac1,'ivacn=',ivacn,' nivac=',nivac
              call tlac1(ndim3,nivac,azn,vr,vt,vz                       &
     &        ,temp,ekonc,akonc,dens,dustd,dustt,mshade                 &
     &        ,hi,hneg,h2,hei,hipf)
	    endif
!
!		calculation of level population and damping for 
!	    the sp. lines (0,nline)
!	    partition functions of the element at different depth
!	    stavs(i,j) :  i-depth point, j-ion
	    if(iline.eq.1)then
	      do ilin=1,nline
	        if(ipart.eq.1)then
                  call pfdwor(ndim3,nivac,mion,nion(ilin),iat(ilin)     &
     &            ,temp,ekonc,stavs)
                else
                  call pfirw(ndim3,nivac,mion,nion(ilin),sta,temp,stavs)
                endif
	        do ll=1,mion-1
	          zipp(ll)=zip(ll,ilin)
	        enddo
!		rr(i,j)- popul of the j-th ion/total element population
!		at i-th depth
	        call zastup(temp,ekonc,zipp,stavs,rr,nivac,nion(ilin))
	        do kk=1,nivac
                  popul(kk,ilin)=bolt(temp(kk),stavs(kk,iion(ilin))     &
     &            ,glo(ilin),elo(ilin),1)                               &
     &            *rr(kk,iion(ilin))*abund(ilin)*abhyd*akonc(kk)
	        enddo
!		calculates damping constants
!	    	aa-combined damping constant
!	    	dop-Doppler halfwidth	
                do kk=1,nivac
                  call utlm(wlab(ilin),eup(ilin),iat(ilin),hmc(ilin)    &
     &            ,iion(ilin),zipp(iion(ilin)),hi(kk),hei(kk)           &
     &            ,ekonc(kk),temp(kk),vtrb(kk),gr(ilin),gr0(ilin)       &
     &            ,gs(ilin),gs0(ilin),gw(ilin),gw0(ilin)                &
     &            ,aa(kk,ilin),dop(kk,ilin))
                enddo
!		prints various quantities	    
	        if (iprint.eq.1)then        
	          write(2,*)' ilin=',ilin	
                  call tlac2(ndim3,nivac,mion                           &
     &            ,nion(ilin),iion(ilin),azn,zipp,eup(ilin)             &
     &            ,temp,stavs,rr,gr(ilin),gs(ilin),gw(ilin))
	        endif
	      enddo
	    endif
!		frequency cycle
	    do 410 nu=1,nfreq
	      if(nu.eq.ionu.and.ii.eq.ior.and.jj.eq.iot)then
	        iprint=1
	      else
	        iprint=0
	      endif		
	      if(loglam.eq.0)then
!		interpolate H continuum opacity to alam(nu)
!		it is 2x faster then calculating continuum opacity
!               directly at alam(nu)
!               It presumes that <alam1,alamn> is short enough.??
                call medop(ndim3,nivac,alam1,alamn,alam(nu)             &
     &          ,ophbf,ophbf1,ophbf2,ophff,ophff1,ophff2                &
     &          ,ophrs,ophrs1,ophrs2,ophn,ophn1,ophn2)
              else
!               or calculate ophbf,ophff,ophrs,ophn at each alam(nu):
                call ophyd(nivac,temp,ekonc,akonc,hi,hneg,abhyd         &
     &          ,ophbf,ophff,ophrs,ophn,hipf,alam(nu))
	      endif	
!		wlab-laboratory lambda of the line center in [A]
!		alam(nu)-solving lambda in [A]
!		freq-frequency of alam(nu)
	      freq=clight/alam(nu)*1.d8
!             	interpolates dust opacity to freq (opmis,opmia)
	      if(imie.gt.0)then
		do idst=1,ndust
		  do ij=1,nmie(idst)
		    pmiex(ij)=opmiex(ij,idst)
		    pmies(ij)=opmies(ij,idst)
		    pmiea(ij)=opmiea(ij,idst)
		  enddo  
	          call intrp(pmiex,pmies,nmie(idst),freq,opmis(idst))
	          call intrp(pmiex,pmiea,nmie(idst),freq,opmia(idst))
	        enddo  
	      else
		do idst=1,ndust
	          opmis(idst)=0.d0
	          opmia(idst)=0.d0
	        enddo  
	      endif  
!             	interpolates dust phase function to freq	      
	      if(imiepf.eq.1)then
	        call intrp2(pfmiex,pfmie,npfmie,freq,pfmief)
	      else
	        do ipfang=1,npfang
		  pfmief(ipfang)=0.d0
	        enddo
	      endif  
!
!		calculation of aintb- boundary condition for intensity
!		behind the last untransparent object along the ray
              call bcond2(iunt,dens,temp,vz,vdif,freq                   &
     &        ,dcut1,dcut2,dcut3,dcutn,darkl,cdelta                     &
     &        ,lunt1,nstar1,xstar1,star1,nspec1,tspec1,alb1x,alb1y,nalb1&
     &        ,tstar,dlst,dlst2,irrst,albst,wstar1,fstar1,jstar1        &
     &        ,lunt2,nstar2,xstar2,star2,nspec2,tspec2,alb2x,alb2y,nalb2&
     &        ,tempcp,dlcp,dlcp2,irrcp,albcp,wstar2,fstar2,jstar2       &
     &        ,lunt3,nstar3,xstar3,star3,aint0,aintb)
              if(iprint.eq.1)then
                write(2,*)' iunt=',iunt,' aintb=',aintb
              endif
!
!	        this is the 1D radiative trasfer routine
              call rte(nivac,alam(nu),aintb,aint(ii,jj,nu),opdep        &
     &        ,iunt,wlab,bij,aa,dop,popul,nline                         &
     &        ,istar,tstar,rstar,wstar1,fstar1,jstar1,lunt1             &
     &        ,vxstr,vystr,vzstr                                        &
     &        ,icomp,tempcp,rcp,wstar2,fstar2,jstar2,lunt2              &
     &        ,vxcpr,vycpr,vzcpr,xcpr,ycpr,zcpr                         &
     &        ,xsecx,xsecy,txsec,nxsec,ntxsec,amixf,amixr,iopac         &
     &        ,ophbf,ophff,ophrs,ophn,ithom,irayl,ihyd                  &
     &        ,opmis,opmia,dtlow,dthig,drmf,ndust,imie                  &
     &        ,imiepf,pfmief,pfang                                      &
     &        ,temp,dens,ekonc,akonc,dustd,dustt,azn,vr,vt,vz,mshade    &
     &        ,ar(ii),at(jj),denvac,iprint,eps)
              if(opdep.gt.opdepm)opdepm=opdep
!             end freq cycle 
410	    continue
!           end x-y cycle
420	  continue
430	  continue
	  write(2,*)' max optical depth for this view angle'
	  write(2,'(a,es10.2)')' is: ',opdepm	  
!
!	  write 2Dimage_xxx files with 2D images for one frequency=ionu
!	  Original alam(nu) were shifted by '-gamma' velocity. 
!	  In the next step alam(nu) will be shifted back by gamma vel.
!	  The 2D images correspond to the curent alam(ionu) in the body
!         frame which corresponds to 
!         alamg=alam(ionu)+vgamma/clight*alam(nu)
!         in the observers frame which corresponds to original alam(nu) 
!          plocha=0.d0
	  image1='2Dimage_'
	  write(image2,'(i3.3)')iang
	  image12=image1 // image2
	  open(100+iang,file=image12,status='unknown')
!         apply reddening/extinction
          if(iext.eq.1)then
            call extccm89(ebv,rv,alam(ionu)*1.d-4,amag,aflx)
	  else
	    amag=0.d0
	    aflx=1.d0
	  endif  
	  do i=1,nbod1
!            write(20+iang,*)
            write(100+iang,*)
	    do 440 j=1,nbod2
!              write(20+iang,455)ar(i),at(j),aint(i,j,ionu)/aflx
              write(100+iang,455)ar(i),at(j),aint(i,j,ionu)/aflx
!	      do 435 k=nbod3,1,-1
!	        if (adens(k).gt.dcut1)then
!		  plocha=area(i,j)+plocha
!		  goto 440
!		endif
!435	      continue
440         continue
          enddo
          close(100+iang)
455       format(2e12.4,e14.5e3)
!	  if(rstar.gt.0.d0)then
!	    ratio=plocha/pi/rstar/rstar
!            write(2,*)' star disk',pi*rstar**2
!            write(2,*)' sum_ij untranarea_ij/star disk',ratio            
!	  endif
!	  integrate flux over the 2D surface
	  do nu=1,nfreq
            flux(nu,iang)=0.d0
	    do i=1,nbod1
	    do j=1,nbod2
              flux(nu,iang)=aint(i,j,nu)*area(i,j)/dd/dd+flux(nu,iang)
	    enddo
	    enddo
!           apply reddening/extinction
            if(iext.eq.1)then
	      call extccm89(ebv,rv,alam(nu)*1.d-4,amag,aflx)
	      flux(nu,iang)=flux(nu,iang)/aflx
	    endif         
	  enddo
	  der=(flux(nfreq,iang)-flux(1,iang))/(alam(nfreq)-alam(1))
	  do nu=1,nfreq
!	    flux in [erg/cm^2/s/Hz], fluxl in [erg/cm^2/s/cm], 
!	    fluxn -normalized flux 
!	    fluxs -shifted fluxn	
!	    alamg -lambda's shifted by gamma velocity [Ang]
	    cont=flux(1,iang)+der*(alam(nu)-alam(1))
	    fluxl=flux(nu,iang)*clight/alam(nu)/alam(nu)*1.d16	
	    fluxn(nu,iang)=flux(nu,iang)/cont
	    fluxs=fluxn(nu,iang)-dble(iang-1)*offset
	    if(nline.ge.1)then
	      refw=wlab(1)	
	    else
              refw=alam1 
	    endif
	    alamg=alam(nu)+vgamma/clight*alam(nu)
	    dvel=(alamg-refw)/refw*clight*1.d-5
!	    shifts lambda & writes shellspectrum
            write(4,510)alamg,dvel,flux(nu,iang),fluxl                  &
     &      ,fluxn(nu,iang),fluxs
!           write(100+iang,'(f10.3,f10.1,5e14.7)')alam(nu),dvel         &
!     &      ,flux(nu,iang),fluxl,fluxn(nu,iang),fluxs
!	    writes ligtcurve
            write(11,530)alpha(iang)/2.d0/pi,dvel                       &
     &      ,-2.5d0*dlog10(fluxl),alamg,flux(nu,iang)
          enddo
	  write(11,*)
	  write(4,*)
	  write(2,*)
510       format(f12.3,1x,f10.1,5es14.5e3)
530 	  format(f7.3,es11.3,f10.4,es12.4,es18.9e3)
!	  end of rotation cycle
	enddo
	write(3,*)' 0 errors found'
	close(2)
	close(3)
	close(4)
	close(11)
!	call exit(0)
!	end program SHELLSPEC
	return
	end subroutine SHELLSPEC
!-----------------------------------------------------------------------
	subroutine exit(icode)
!	icode=0 no error was detected
!	icode><0 some error or warning was detected
!	check the error file or screen for more detailes
	implicit none
        integer, parameter:: i4=4,i8=8
	integer:: icode
	if(icode.ne.0) then
          stop 1
        else
          stop 0
        endif
	end
!-----------------------------------------------------------------------
        subroutine rte(nbod,alam,aintb,aint,opdep                       &
     &  ,iunt,wlab,bij,aa,dop,popul,nline                               &
     &  ,istar,tstar,rstar,xstar1,star1,nstar1,lunt1                    &
     &  ,vxstr,vystr,vzstr                                              &
     &  ,icomp,tempcp,rcp,xstar2,star2,nstar2,lunt2                     &
     &  ,vxcpr,vycpr,vzcpr,xcpr,ycpr,zcpr                               &
     &  ,xsecx,xsecy,txsec,nxsec,ntxsec,amixf,amixr,iopac               &
     &  ,ophbf,ophff,ophrs,ophn,ithom,irayl,ihyd                        &
     &  ,opmis,opmia,dtlow,dthig,drmf,ndust,imie                        &
     &  ,imiepf,pfmief,pfang                                            &
     &  ,temp,dens,ekonc,akonc,dustd,dustt,az,vr,vt,vz,mshade           &
     &  ,ari,ati,denvac,iprint,eps)
!	to solve the radiative tranfer along the line of sight
!	at the wavelength -alam
!	input:
!	  nbod -number of depth points
!	  alam -wavelength [A]
!	  aintb -boundary condition for intensity 
!	  iunt -last untransparent point
!	  wlab -laboratory lambda of the line center [A]
!	  bij - Einstein coefficient B_lu per unit solid angle
!	  aa - combined damping constant
!	  dop - Doppler halfwidth
!	  popul - lower level population of spectral line
!	  nline - number of spectral lines
!	  istar,icomp -on/off switch of the scattered light
!	  tstar,rstar -temperature,radius of the spherical 
!			central star
!	  xstar1,star1,nstar1,lunt1 -spectrum of the central star
!	  vxstr,vystr,vzstr -velocity of the central star
!	  tempcp,rcp -temperature,radius of the spherical 
!			secondary star
!	  xstar2,star2,nstar2,lunt2 -spectrum of the second star
!	  vxcpr,vycpr,vzcpr -velocity of the second star
!	  xcpr,ycpr,zcpr -position of the second star
!	  xsecx,xsecy,txsec,nxsec,ntxsec -extra gas opacity table	
!	  amixf -molecule number density
!	  amixr -constant molecule abundance (mixing ratio)
!               (relative to H nuclei)
!	  iopac -extra gas opacity table switch
!	  ophbf,ophff,ophrs - opacity hydrogen 
!			(bound-free, free-free, Rayleigh scattering)
!         ophn -opacity of H- (bound-fee, free-free)
!	  opmis,opmia - dust opacity (scattering, absorption)
!		    per 1g of dust material for different species
!         dtlow,dthig -temperature interval for different dust species
!	  drmf -relative mas fraction for different dust species
!	  ndust -number of different dust species
!	  imie,imiepf -switch for Mie scattering
!	  pfmief -phase function for dust, function of angles 
!	  pfang -cos(angle) for phase function of dust
!	  ithom,irayl -switch for Thompson & Rayleigh scattering
!         ihyd -switch for hydrogen bf, ff opacity
!	  temp,dens - temperature, density
!	  ekonc,akonc - electron and atom number densities
!	  dustd,dustt - dust density and temperature 
!	  vr,vt,vz - x,y,z velocity along the ray
!         mshade - shadows
!         mshade=1 scattering only from the central star
!         mshade=2 scattering only from the companion star
!         mshade=3 scattering from both stars (default)
!         mshade=0 no scattering
!         az - z coordinate along the ray
!	  ari,ati - x,y coordinates defining the line 
!			of sight along z
!	  denvac - limiting density of vacuum
!	  iprint - signal to print
!	  eps - In LTE eps=1. ( S=eps*B+(1-eps)*J ) 
!	output: 
!	  aint-emerging intensity [erg/cm^2/s/Hz/rad]
!	  opdep -optical depth along the ray after 
!                 	the possible non-transparent object
!
!       temp, ane -are assumed to have reasonable values all along
!       the beam. An empty space can be identified as dens<denvac.
!	Integration starts behind the last untransparent object.
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
        include 'param.inc'
	real(i8):: aa(ndim,mline),dop(ndim,mline),popul(ndim,mline)
        real(i8):: flab(mline),wlab(mline),bij(mline),bij2(mline)
        real(i8):: ophbf(ndim),ophff(ndim),ophrs(ndim),ophn(ndim)
        real(i8):: temp(ndim),dens(ndim),ekonc(ndim),akonc(ndim)
        real(i8):: dustd(ndim),dustt(ndim),amixf(ndim)
        real(i8):: vr(ndim),vt(ndim),vz(ndim),az(ndim)
        integer(i1):: mshade(ndim)
        real(i8):: bnu(ndim),bnud(ndim),ejnu1(ndim),ejnu2(ndim)
        real(i8):: sf2(ndim),opac2(ndim),emis2(ndim),opdept(ndim)
	real(i8):: dtlow(mspecx),dthig(mspecx),drmf(mspecx)
	real(i8):: opmis(mspecx),opmia(mspecx)
        real(i8):: pfang(npfang),pfmief(npfang)
        real(i8):: xstar1(mstarx),star1(mstarx)
        real(i8):: xstar2(mstarx),star2(mstarx)
	real(i8):: xsecx(mstarx),xsecy(mstarx,mspecx)
        real(i8):: txsec(mspecx)
	integer:: nbod,iunt,nline,istar,nstar1,lunt1,icomp,nstar2,lunt2
	integer:: nxsec,ntxsec,iopac,ithom,irayl,ihyd,ndust,imie,imiepf
	integer:: iprint
	real(i8):: alam,aintb,aint,opdep,tstar,rstar,vxstr,vystr,vzstr
	real(i8):: tempcp,rcp,vxcpr,vycpr,vzcpr,xcpr,ycpr,zcpr
	real(i8):: amixr,ari,ati,denvac,eps
	real(i8):: voigt,planck
	integer:: jj1,jj2,ll,idst,id1,iuh
	real(i8):: bol,clight,hjed,plank,pi,EVERG,EVCM,ERGCM,freq,oplin
	real(i8):: flabs,dfreqd,vprof,opmol,opth,opmiej,opmieb
	real(i8):: emhibf,emhiff,emhn,er2,uhcos,dippf1,deruh,dstpf1
	real(i8):: emth1,emra1,emmij1,dippf2,dstpf2,emth2,emra2,emmij2
	real(i8):: emth,emra,emmij,emlinb,emlinj,emmib,emmol
	real(i8):: dzz,opdepd,opdlim,opslim,cf12,cf1,cf2
	real(i8):: aintbt,xsecs
!	CONSTANTS AND UNITS FROM NIST(2002)
	bol=1.3806503d-16
	clight=2.99792458d10
	hjed=1.66053873d-24
	plank=6.62606876d-27
	pi=3.1415926535897931d0
!	1EV=X ERG, X CM^-1; 1ERG= X CM^-1 
	EVERG=1.602176462D-12
	EVCM=8.06554477D3
	ERGCM=5.03411762D15
!		flab -frequency [1/s] of wlab
!		flabs-frequency of the line center in comoving frame
!		alam-solving lambda in [A]
!		freq-frequency of alam
	freq=clight/alam*1.D8
	if(nline.gt.0)then
	  do 10 ll=1,nline
	    flab(ll)=clight/wlab(ll)*1.D8
	    bij2(ll)=bij(ll)*plank*flab(ll)/dsqrt(pi)
10	  continue
	endif
	if (nbod-iunt.lt.2)then
!         aint-is set to boundary condition for intensity
	  aint=aintb
	  opdep=0.d0
          goto 100
        endif  
!		popul-population of the lower level
!		opac-opacity
!		emis-emisivity
!		sf-source function
        if (iprint.eq.1)then
	  write(2,*)' lambda=',alam
          write(2,'(a,a,a)')' id  '                                     &
     &    ,'ophbf     ophff     opth      ophrs     ophn      '         &
     &    ,'opmieb    opmiej    oplin     opmol'
	endif
!		initialize the first value of jj for hunt
	jj1=1
	jj2=1
!		BEGINNING OF THE DEPTH CYCLE ---------------------------
	do 60 id1=iunt+1,nbod
	  if(dens(id1).lt.denvac)then
!	    icheck=0
	    sf2(id1)=0.d0
	    emis2(id1)=0.d0
	    opac2(id1)=0.d0
            if (iprint.eq.1)then
	      write(2,'(i3,a)')id1,' empty space-skipped'
	    endif
	    goto 60
	  endif
!		Line opacity (takes into account the velocity field)	  
	  oplin=0.d0
	  if(nline.gt.0)then
	    do 20 ll=1,nline
	      flabs=vz(id1)/clight*flab(ll)+flab(ll)
	      dfreqd=(freq-flabs)/dop(id1,ll)
!	  	\INT(VOIGT F.)dFREQ 
!		HAS NORMALIZATION= SQRT(PI)*DOP.W.[1/S]
	      vprof=voigt(dfreqd,aa(id1,ll))
              oplin=(1.d0-dexp(-4.7992375d-11*flab(ll)/temp(id1)))      &
     &        *bij2(ll)*popul(id1,ll)*vprof/dop(id1,ll)+oplin
20	    continue
	  endif
!		extra true absorption gas opacity (xsection) from table
!               (takes into account the velocity field)
	  if(iopac.eq.1.and.nxsec.gt.1)then
	    flabs=-vz(id1)/clight*freq+freq
            call int2d(xsecx,txsec,xsecy,nxsec,ntxsec,flabs,temp(id1)   &
     &      ,xsecs)
            opmol=amixf(id1)*xsecs
!            opmol=akonc(id1)*amixr*xsecs
	  else
	    opmol=0.d0
	  endif
!                Hydrogen bf, ff, H- opacity (ignores velocity field)
	  if(ihyd.ne.1)then
	    ophbf(id1)=0.d0
	    ophff(id1)=0.d0
	    ophn(id1)=0.d0
	  endif  
!		Thomson scattering opacity (no velocity dependence)
          if(ithom.eq.1)then
	    opth=ekonc(id1)*6.65d-25
	  else
	    opth=0.d0
	  endif
!               dust opacity (Mie scattering+absorption)
!               (ignores velocity field)
	  opmiej=0.d0
	  opmieb=0.d0
          if(imie.gt.0)then
	    do idst=1,ndust
              if(dustt(id1).gt.dtlow(idst).and.                         &
     &        dustt(id1).le.dthig(idst))then
                opmiej=dustd(id1)*drmf(idst)*opmis(idst)+opmiej
                opmieb=dustd(id1)*drmf(idst)*opmia(idst)+opmieb
              endif
            enddo  
          endif
!		Rayleigh scattering on HI opacity
!               (ignores velocity field)
          if(irayl.ne.1)ophrs(id1)=0.d0
!		total opacity
	  opac2(id1)=oplin+opth+ophbf(id1)+ophff(id1)+ophrs(id1)
	  opac2(id1)=opac2(id1)+ophn(id1)+opmiej+opmieb+opmol
	  if(opac2(id1).le.0.d0)then
	    write(*,*)' error: zero opacity'
	    write(3,*)' 1 error: zero opacity'
	    write(*,*)'id1   dens    dustd   dustt    opac'
	    write(*,30)id1,dens(id1),dustd(id1),dustt(id1),opac2(id1)
30	    format(i4,4es9.1)
	    call exit(1)
	  endif  
	  bnu(id1)=planck(freq,temp(id1),1)
	  bnud(id1)=planck(freq,dustt(id1),1)
!		HI, H-  bound-free and free-fee emissivity
	  emhibf=ophbf(id1)*bnu(id1)
	  emhiff=ophff(id1)*bnu(id1)
	  emhn=ophn(id1)*bnu(id1)
!		scattering from stars
!		ejnu1,2 - mean intensity from the star and companion
!		emth,emra - Thompson and Rayleigh scattering emissivity
!		emmij - Mie scattering emissivity
!               emmib - Mie absorption emissivity
          if(ithom.eq.1.or.irayl.eq.1.or.imie.eq.1)then
	    if(istar.gt.0.and.(mshade(id1).eq.1.or.mshade(id1).eq.3))   &
     &      then	
              call jnu(rstar,tstar,ari,ati,az(id1)                      &
     &        ,vr(id1),vt(id1),vz(id1),freq                             &
     &        ,xstar1,star1,nstar1,lunt1                                &
     &        ,vxstr,vystr,vzstr,jj1,ejnu1(id1),iprint,0.d0,0.d0,0.d0)
	      er2=ari*ari+ati*ati+az(id1)*az(id1)
	      uhcos=az(id1)/dsqrt(er2)	
	      if(er2.gt.(rstar*rstar*9.d0))then
!		dipol phase function scattering
	        dippf1=0.75d0*(1.d0+uhcos**2)
!		dust phase function		
		if(imiepf.eq.1)then
		  call locate(pfang,npfang,uhcos,iuh)
		  deruh=(pfmief(iuh+1)-pfmief(iuh))
		  deruh=deruh/(pfang(iuh+1)-pfang(iuh))
		  dstpf1=deruh*(uhcos-pfang(iuh))+pfmief(iuh)
		else
		  dstpf1=1.d0
		endif
	      else
	        dippf1=1.d0
	        dstpf1=1.d0
	      endif    
	      emth1=dippf1*ejnu1(id1)*opth
	      emra1=dippf1*ejnu1(id1)*ophrs(id1)
	      if(imie.eq.1)then	
	        emmij1=dstpf1*opmiej*ejnu1(id1)
	      else
	        emmij1=0.d0
	      endif  
     	    else
     	      ejnu1(id1)=0.d0
     	      emth1=0.d0
	      emra1=0.d0
	      emmij1=0.d0
     	    endif
	    if(icomp.gt.0.and.(mshade(id1).eq.2.or.mshade(id1).eq.3))   &
     &      then	
              call jnu(rcp,tempcp,ari,ati,az(id1)                       &
     &        ,vr(id1),vt(id1),vz(id1),freq                             &
     &        ,xstar2,star2,nstar2,lunt2                                &
     &        ,vxcpr,vycpr,vzcpr,jj2,ejnu2(id1),iprint,xcpr,ycpr,zcpr)
	      er2=(ari-xcpr)**2+(ati-ycpr)**2+(az(id1)-zcpr)**2
	      uhcos=(az(id1)-zcpr)/dsqrt(er2)
	      if(er2.gt.(rcp*rcp*9.d0))then
!		dipol phase function scattering
	        dippf2=0.75d0*(1.d0+uhcos**2)
!               dust phase function   
		if(imiepf.eq.1)then
                  call locate(pfang,npfang,uhcos,iuh)    
                  deruh=(pfmief(iuh+1)-pfmief(iuh))      
                  deruh=deruh/(pfang(iuh+1)-pfang(iuh))		
                  dstpf2=deruh*(uhcos-pfang(iuh))+pfmief(iuh)
                else
                  dstpf2=1.d0
                endif  
	      else
	        dippf2=1.d0
	        dstpf2=1.d0
	      endif    
	      emth2=dippf2*ejnu2(id1)*opth
	      emra2=dippf2*ejnu2(id1)*ophrs(id1)
	      if(imie.eq.1)then
	        emmij2=dstpf2*opmiej*ejnu2(id1)
	      else
	        emmij2=0.d0
	      endif  
     	    else
     	      ejnu2(id1)=0.d0
     	      emth2=0.d0
	      emra2=0.d0
	      emmij2=0.d0
     	    endif     
	  else
	    ejnu1(id1)=0.d0
	    ejnu2(id1)=0.d0
     	    emth1=0.d0
	    emra1=0.d0
	    emth2=0.d0
	    emra2=0.d0
	    emmij1=0.d0
	    emmij2=0.d0	
	  endif
	  emth=emth1+emth2
	  emra=emra1+emra2
	  emmij=emmij1+emmij2
!		line emissivity
	  emlinb=eps*oplin*bnu(id1)
	  emlinj=(1.d0-eps)*oplin*(ejnu1(id1)+ejnu2(id1))
!		Mie emissivity from dust (thermal)
	  emmib=opmieb*bnud(id1)
!		test for disk	  
	  if(imie.eq.3)emmij=opmiej*bnud(id1)
!		molecule emissivity
	  emmol=opmol*bnu(id1)	  
!		total emissivity
	  emis2(id1)=emlinb+emlinj+emth+emhibf+emhiff+emra+emhn
	  emis2(id1)=emis2(id1)+emmib+emmij+emmol
!		total source function
	  sf2(id1)=emis2(id1)/opac2(id1)
	  if (iprint.eq.1)then
            write(2,'(i3,9es10.2)')id1                                  &
     &      ,ophbf(id1),ophff(id1),opth,ophrs(id1),ophn(id1)            &
     &      ,opmieb,opmiej,oplin,opmol
	  endif
60      continue
!               END OF THE DEPTH CYCLE----------------------------------
!		solver1: for opt thin medium, integrates aint and opdep
!		opdep- optical depth as measured from behind the shell
!		      or from an untransparent object
!	aint-is set to boundary condition for intensity
!	aint=aintb
!        if (iprint.eq.1)then
!	  write(2,'(a)')' id    opdep      aint '  
!	endif
!	opdep=0.d0
!	icheck=0
!        do 65 id1=iunt+1,nbod
!	  if(dens(id1).lt.denvac)then
!	    icheck=0
!	    goto 65
!	  endif
!	  if (icheck.gt.0)then
!	    opac12=(opac1+opac2(id1))*0.5d0
!	    sf12=(sf1+sf2(id1))*0.5d0
!	    dz12=(az(id1)-az(id1-1))
!	    alph=1.d0/opac12/dz12
!	    aint=(aint*(alph-0.5d0)+sf12)/(alph+0.5d0)
!	    opdep=opac12*dz12+opdep
!            if (iprint.eq.1)then
!              write(2,70)id1,opdep,aint
!            endif
!	  endif
!	  opac1=opac2(id1)
!	  sf1=sf2(id1)
!	  icheck=icheck+1
!65	continue
!
!               solver2: formal solution, integrates aint and opdep
!		opdept- optical depth as measured from the observer to
!        behind the shell or up to an untransparent object
        if (iprint.eq.1)then
	  write(2,'(a)')' id   opdep     contrib.f   aint '  
	endif
!	calculate op.depth
	opdept(nbod)=0.d0
	do id1=nbod-1,iunt+1,-1
          dzz=az(id1+1)-az(id1) 
          opdepd=(opac2(id1+1)+opac2(id1))*0.5d0*dzz
          opdept(id1)=opdepd+opdept(id1+1)
        enddo
	aint=0.d0
	do id1=nbod-1,iunt+1,-1
          opdepd=(opdept(id1)-opdept(id1+1))
!         break a huge leap in opt. depth into steps
          opdlim=10.d0
          opslim=0.3d0
          if(opdept(id1+1).lt.opdlim.and.opdepd.gt.opslim)then
	    call split2(opdept(id1),opdept(id1+1),opslim,sf2(id1),      &
     &      sf2(id1+1),cf12)
          else
            cf1=sf2(id1+1)*dexp(-opdept(id1+1))
            cf2=sf2(id1)*dexp(-opdept(id1))
            cf12=(cf1+cf2)*0.5d0*opdepd
          endif  
          aint=cf12+aint
          if (iprint.eq.1)then
            write(2,70)id1,opdept(id1),cf12,aint
          endif
        enddo  
	aintbt=aintb*dexp(-opdept(iunt+1))
        aint=aint+aintbt
        opdep=opdept(iunt+1)
        if(iprint.eq.1)then 
          write(2,*)'iunt+1  opdep   aintbt    aint'
          write(2,70)iunt+1,opdep,aintbt,aint
        endif  
70      format(i3,3es11.3)
!
        if (iprint.eq.1)then
          write(2,'(a,a,a)')' id    opac2     emis2     bnu       bnud' &
     &    ,'      ejnu1     ejnu2     sf2       az        vz'
          do id1=iunt+1,nbod
            if(dens(id1).gt.denvac)then
              write(2,'(i3,10es10.2)')id1,opac2(id1),emis2(id1)         &
     &        ,bnu(id1),bnud(id1),ejnu1(id1),ejnu2(id1),sf2(id1)        &
     &        ,az(id1),vz(id1)
            endif
          enddo  
        endif    
100     continue 
        return   
        end subroutine rte
!-----------------------------------------------------------------------
	FUNCTION VOIGT(X,Y)
!      IMPLICIT REAL*8 (A-H,O-Z)
!      REAL*8 VOIGT
!      REAL*4 X,Y
!	taken from the old Synspec code:
! Hubeny I., Lanz T., Jeffery C.S., 1994, in Newsletter on Analysis
! of Astronomical spectra No.20, ed. C.S. Jeffery (CCP7; St. Andrews:
! St. Andrews Univ.), 30
!   	HUMLICEK VOIGT FUNCTION ALGORITHM  (LOWER PRECISION VERSION)
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: VOIGT
	real(i8):: X,Y
        real(i8):: T(6),C(6),S(6)
	real(i8):: WR,Y1,Y2,R,D,D1,D2,D3,D4
	integer:: i
        data T/.314240376,.947788391,1.59768264,2.27950708,             &
     &  3.02063703,3.8897249/
        data C/1.01172805,-.75197147,1.2557727E-2,                      &
     &  1.00220082E-2,-2.42068135E-4,5.00848061E-7/
        data S/1.393237,.231152406,-.155351466,6.21836624E-3,           &
     &  9.19082986E-5,-6.27525958E-7/
      WR=0.
      Y1=Y+1.5
      Y2=Y1*Y1
      DO 3 I=1,6
        R=X-T(I)
        D=1./(R*R+Y2)
        D1=Y1*D
        D2=R*D
        R=X+T(I)
        D=1./(R*R+Y2)
        D3=Y1*D
        D4=R*D
        WR=WR+C(I)*(D1+D3)-S(I)*(D2-D4)
    3 CONTINUE
      VOIGT=DBLE(WR)
      RETURN
      END FUNCTION VOIGT
!-----------------------------------------------------------------------
	FUNCTION BOLT(TEMP,G1,G2,POT,IRIAD)
!-----------------------------------------------------------------------
!	      TEMP -TEMPERATURE [KELVIN]			       C
!	      G1   -STATISTICAL WEIGHT OF THE LOWER LEVEL-m	       C
!	      G2   -STATISTICAL WEIGHT OF THE UPPER LEVEL-n	       C
!	      POT  -ENERGY DIFFERENCE (UPPER-LOWER LEVEL) [EV], POT>0  C
!	      OUTPUT -POPULATION RATIO N(n)/N(m) 		       C
!---------------------OR IF:            -------------------------------C
!	      TEMP -TEMPERATURE [KELVIN]			       C
!	      G2   -STATISTICAL WEIGHT OF THE LEVEL-n		       C
!	      G1   -PARTITION FUNCTION OF THE ION                      C
!	      POT  -EXCITATION POTENCIAL OF THE n-TH LEVEL [EV]	       C
!		   (WITH RESPECT TO GROUND STATE), POT>0	       C
!	      OUTPUT -RATIO N(n)/N, (n-TH LEVEL POPULATION TO	       C
!		    TO THE ION POPULATION			       C
!-----------------------OR IF------------------------------------------C
!	      IRIAD=1  POT IS IN [1/CM]         		       C
!-----------------------------------------------------------------------
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: BOLT
	real(i8):: TEMP,G1,G2,POT
	integer:: IRIAD
	real(i8):: EVCM
	EVCM=8.06554477D3
	IF(IRIAD.EQ.1)THEN
	  BOLT=G2/(G1*DEXP(2.302585D0*5039.77D0/TEMP*POT/EVCM))
	ELSE
	  BOLT=G2/(G1*DEXP(2.302585D0*5039.77D0/TEMP*POT))
	ENDIF
	RETURN
	END FUNCTION BOLT
!-----------------------------------------------------------------------
        SUBROUTINE UTLM(DL0,EP,IAT,AM,IION,CHI,HI,HEI,EKONC,TEMP,VT     &
     &  ,GR,GR0,GS,GS0,GW,GW0,A,DOP)
!-----------------------------------------------------------------------
!	     DL0   -LINE CENTER [ANGSTROM]			       C
!	     EP    -EXCIT.POTENTIAL OF THE UPPER LEVEL [1/CM]	       C
!	     IAT   -ATOMIC NUMBER 				       C
!	     AM    -MASS NUMBER 				       C
!	     IION  -DEGREE OF IONIZATION, IION=1 FOR NEUTRAL...	       C
!	     CHI   -IONIZ. POTENTIAL OF THE ION [EV] 		       C
!	     HI    -NEUTR. HYDROGEN NUMBER DENSITY [CGS]	       C
!	     HEI   -NEUTR. HELIUM NUMBER DENSITY [CGS]		       C
!	     EKONC -ELECTRON NUMBER DENSITY [CGS]		       C
!	     TEMP   TEMPERATURE [KELVIN]			       C
!	     VT    -TURBULENCE VELOCITY [CGS]    		       C
!	     GR    -RADIATIVE DAMPING				       C
!	     GS    -STARK DAMPING				       C
!	     GW    -VAN DER WAALS DAMPING			       C
!								       C
!	     OUTPUT:                                                   C
!		A   - DAMPING PARAMETER                                C
!		DOP - DOPPLER. HALF-WIDTH [1/S]                        C
!-----------------------------------------------------------------------
	implicit none
        integer, parameter:: i4=4,i8=8
	integer:: IAT,IION
	real(i8):: DL0,EP,AM,CHI,HI,HEI,EKONC,TEMP,VT
	real(i8):: GR,GR0,GS,GS0,GW,GW0,A,DOP
	real(i8):: PI,EVCM,CLIGHT,DLL0,BOL,S,Z,EFF2,R2
	PI=3.1415926535897931d0
	EVCM=8.06554477D3
	CLIGHT=2.99792458D10
	DLL0=1.E-8*DL0
	BOL=1.380622D-16
	S=IAT
	Z=IION
!	     RADIATIVE DAMPING
	IF(GR0.GT.0.D0)THEN
	  GR=DEXP(2.302585D0*GR0)
	ELSE
	  GR=2.4734D-22*clight*clight/DLL0/DLL0
!	     V PRIPADE AK POZNAS SIRKU AUTOIONIZACNEJ HLADINY GRA
!	     V [EV] => GR=1.52E12*GRA
	ENDIF
	EFF2=Z*Z*13.595D0/(CHI-EP/EVCM)
	IF(EFF2.LT.0.D0.OR.EFF2.GT.25.D0)EFF2=25.D0
!	     STARK DAMPING
!	 IF(CHI.LE.EP/EVCM) THEN
!	     AUTOIONIZ. HORNA HLADINA ? =>BERIEME JU AKO EFF2=25
!            AJ VO VAN DER WAALSOVOM PRE S<20 (AT.CISLO)
!	   EFF2=0.D0
!	 ENDIF
	IF(GS0.NE.0.D0)THEN
	  GS=DEXP(2.302585D0*GS0)
	ELSE
	  GS=1.0D-8*EFF2**2.5D0
	ENDIF
!	     VAN DER WAALS DAMPING
	IF (INT(S).LT.21.OR.INT(S).GT.28)THEN
          R2=2.5D0*EFF2*EFF2/Z/Z
	ELSE
	  R2=(45.D0-S)/Z
	ENDIF
	IF(GW0.NE.0.D0)THEN
	  GW=DEXP(2.302585D0*GW0)
	ELSE
	  GW=4.5D-9*R2**0.4D0
	ENDIF
!	     DOPPLER halfwidth in frequency
	DOP=DSQRT(2.D0*BOL*TEMP/(AM*1.672614D-24)+VT*VT)/DLL0
!	     FRAME DAMPING PARAMETER
	A=(GR+GS*EKONC+GW*(HI+0.42D0*HEI)*(TEMP/1.D4)**0.3D0)
        A=A/(4.D0*PI*DOP)
	RETURN
	END
!-----------------------------------------------------------------------
	SUBROUTINE PFIRW(NDIM,NBOD,mion,NION,STA,TEMP,STAVS)
!	partition functions after Irwin
!	input: sta -Irwin's tabulated coeficients for each ion
!		of an atom
!	output: stavs(i,j)-partition function of j-th ion at i-th depth
	implicit none
        integer, parameter:: i4=4,i8=8
!	PARAMETER(mion=9)
	integer:: NDIM,NBOD,mion,NION
        real(i8):: STA(6,mion),STAVS(NDIM,mion),TEMP(NDIM)
	integer:: i,j
	real(i8):: SALA,ARG
        DO 80 I=1,NBOD
	  DO 70 J=1,NION
	    IF (ABS(STA(1,J)-0.D0).LT.1.D-30) GO TO 60
	    SALA=DLOG(TEMP(I))
	    ARG=STA(1,J)+STA(2,J)*SALA+STA(3,J)*SALA*SALA
            ARG=ARG+STA(4,J)*SALA**3+STA(5,J)*SALA**4+STA(6,J)*SALA**5
	    STAVS(I,J)=DEXP(ARG)
	    GO TO 70
60	    STAVS(I,J)=1.D0
70	  CONTINUE
80      CONTINUE
!        WRITE(2,'(A)')' PARTITION FUNCTIONS'
!        WRITE(2,'(A,20I10)')' TEMP',(I,I=1,NION)
!	 DO 140 I=1,NBOD
! 	   WRITE(2,130)TEMP(I),(STAVS(I,J),J=1,NION)
!130       FORMAT(F7.0,20E10.3) 
!140     CONTINUE
	RETURN
	END
!-----------------------------------------------------------------------
        subroutine potlow(temp,ane,eplow)
!       temp - Temperature[K]
!       ane - electron number density [cm^{-3}]
!       lowering of the ionization potential for neutrals, z=1 in [eV]
!       for ions(z):      eplow(z)=z*eplow(1)
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: temp,ane,eplow
	real(i8):: bol,pi,e,pom,debye
        bol=1.3806503d-16
	pi=3.1415926535897931d0
        e=4.80325d-10
        pom=bol*temp/8.d0/pi/ane
!       debye lenght in [CGS]
        debye=dsqrt(pom)/e
        eplow=1.44d-7/debye
        return
        end
!-----------------------------------------------------------------------
	SUBROUTINE PFDWOR(NDIM,NBOD,mion,NION,IAT,TEMP,EKONC,STAVS)
!
!       partition functions from UCLSYN:
! Smith K.C., Dworetsky M.M., 1988, In: Adelman S.J., Lanz T., eds,
! Elemental Abundance Analyses.  Institut d'Astronomie de l'Univ.  de
! Lausanne, Switzerland, p. 32
!	Budaj J., Dworetsky M.M., Smalley B, 2002,
!	Comm. Univ. London Obs. No. 82
!       URL=http://www.ulo.ucl.ac.uk/ulo_comms/82/index.html
!
!	output: stavs(i,j)-partition function of j-th ion at i-th depth
	implicit none
        integer, parameter:: i4=4,i8=8
	integer:: NDIM,NBOD,mion,NION,IAT
	real(i8):: STAVS(NDIM,mion),EKONC(NDIM),TEMP(NDIM)
	real(i8):: U(5),EPLOW,THETA
	integer:: i,j
        DO 80 I=1,NBOD
          CALL POTLOW(TEMP(I),EKONC(I),EPLOW)
          THETA=5039.77D0/TEMP(I)
!		PARTFN IS IN 'pfdwor.inc'
          CALL PARTFN(IAT,THETA,U,EPLOW)
          DO 70 J=1,NION
	    IF(J.LE.5)THEN
	      STAVS(I,J)=U(J)
	    ELSE
	      STAVS(I,J)=1.D0
	    ENDIF
70        CONTINUE
80      CONTINUE
!        WRITE(2,'(A)')' PARTITION FUNCTIONS'
!        WRITE(2,'(A,20I10)')' TEMP',(I,I=1,NION)
!        DO 140 I=1,NBOD
!          WRITE(2,130)TEMP(I),(STAVS(I,J),J=1,NION)
!130       FORMAT(F7.0,20E10.3)
!140     CONTINUE
        RETURN
        END
!-----------------------------------------------------------------------
	SUBROUTINE ZASTUP(TEMP,EKONC,ZIP,STAVS,RR,NBOD,NION)
!	input:  temp,ekonc,nbod,nion
!		zip-ionization potentials
!		stavs-partition functions
!	output:
!       RR(I,J)- POPULATION OF THE J-TH ION/TOTAL ELEMENT POPULATION
!		 AT I-TH DEPTH
!	R1(J)= N(J)/N(1)
!	RION(J)= N(J)/N(J-1)
!	SAH	-CALLING FUNCTION
!	
	implicit none
        integer, parameter:: i4=4,i8=8
	include 'param.inc'
	integer:: NBOD,NION
	real(i8):: ZIP(mion-1),STAVS(NDIM,mion),TEMP(NDIM),EKONC(NDIM)
        real(i8):: RR(NDIM,mion)
	real(i8):: sah
	real(i8):: RION(mion),R1(mion)
	integer:: i,j
	real(i8):: SUM1
        DO 40 I=1,NBOD
	  RION(1)=1.d0
	  R1(1)=1.d0
	  SUM1=1.D0
	  DO 20 J=2,NION
            RION(J)=SAH(TEMP(I),EKONC(I),ZIP(J-1),STAVS(I,J-1)          &
     &      ,STAVS(I,J))
            R1(J)=R1(J-1)*RION(J)
	    SUM1=SUM1+R1(J)
20	  CONTINUE
	  DO 30 J=1,NION
	    RR(I,J)=R1(J)/SUM1
30	  CONTINUE
40	CONTINUE
!         WRITE(2,*)'   N(J)/N  ION POPULATIONS'
!         WRITE(2,'(A,40I11)')' DM',(I,I=1,NION)
!	DO 50 I=1,NBOD
!	  WRITE(2,130)DM(I),(RR(I,J),J=1,NION)
!130	  FORMAT(40(D10.3,1X))
!50	CONTINUE
	RETURN
	END
!-----------------------------------------------------------------------
        FUNCTION SAH(TE,EK,ZP,S1,S2)
        implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: sah
	real(i8):: TE,EK,ZP,S1,S2
!             SAHA EQUATION (RATIO OF THE NUMBER DENSITIES OF
!		TWO ADJACENT IONS)
!             ZP -in [EV], TE -in [K], EK -in [cm^-3]
        SAH=4.8293745D15*S2/S1*DSQRT(TE)**3
        SAH=SAH*DEXP(-1.16045059D4/TE*ZP)/EK
        RETURN   
        END FUNCTION SAH
!-----------------------------------------------------------------------
        function sahan(te,zp,s1,s2)
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: sahan
	real(i8), intent(in):: te,zp,s1,s2 
!             PART OF SAHA EQUATION WITHOUT ELECTRON NUMBER DENSITY
!             ZP -in [EV], TE -in [K]
        SAHAN=4.8293745D15*S2/S1*DSQRT(TE)**3
        SAHAN=SAHAN*DEXP(-1.16045059D4/TE*ZP)
        RETURN   
        end function sahan
!-----------------------------------------------------------------------
        function sah2(te)
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: sah2
	real(i8), intent(in):: te
	real(i8):: bol,everg,bolte	
!             PART OF SAHA EQUATION FOR H2 MOLECULE: N(H2)/N(HI)^2
!       FROM KURUCZ, TE -in [K]
	bol=1.3806503d-16
	everg=1.602176462d-12
	bolte=bol*te/everg
        SAH2=4.477D0/BOLTE-4.6628D1+1.8031D-3*TE-5.0739D-7*TE*TE
        SAH2=SAH2+8.1424D-11*TE**3-5.0501D-15*TE**4-1.5D0*DLOG(TE)
        SAH2=DEXP(SAH2)
        RETURN   
        end function sah2
!-----------------------------------------------------------------------
	FUNCTION PLANCK(FREQ,TEP,IER)
!-----------------------------------------------------------------------
!	CALCULATES INTENSITY OF THE BLACK BODY RADIATION
!		IF   IER<0  FREQ[CM], PLANCK[ERG/S/CM**2/CM]
!		     IER=0  FREQ [ANGST.], PLANCK [ERG/S/CM**2/ANGSTROM]
!		     IER>0  FREQ [1/S], PLANCK [ERG/S/CM**2/HZ]
!-----------------------------------------------------------------------
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: PLANCK
	real(i8):: FREQ,TEP
	integer:: IER
	if(ier.lt.0)then
          PLANCK=1.1904397D-5/FREQ**5/(DEXP(1.438769D0/FREQ/TEP)-1.)
 	elseif(ier.eq.0)then  
	  PLANCK=1.1904397D27/FREQ**5/(DEXP(1.438769D8/FREQ/TEP)-1.)
	else  
	 PLANCK=1.474501D-47*FREQ**3/(DEXP(4.799216D-11*FREQ/TEP)-1.)
	endif
	return 
	END FUNCTION PLANCK
!-----------------------------------------------------------------------
      SUBROUTINE STATE0(Adyp,Ad,Axi)
!		taken from the Synspec code:
! 	Hubeny I., Lanz T., Jeffery C.S., 1994, 
!	in: Newsletter on Analysis of Astronomical spectra No.20, 
!	ed. C.S. Jeffery (CCP7; St. Andrews: St. Andrews Univ.), 30
!
!     Initialization of the basic parameters for the Saha equation
!
	implicit none
        integer, parameter:: i4=4,i8=8
        include 'param.inc'
!	parameter(matom=99)
      	character(4):: ADYP(MATOM)
      	real(i8):: AD(3,MATOM),AXI(8,MATOM)
      	character(4):: DYP(MATOM)
      	real(i8):: D(3,MATOM),XI(8,MATOM)
	integer:: i,j
!
      DATA DYP/' H  ',' He ',' Li ',' Be ',' B  ',' C  ',               &
     &         ' N  ',' O  ',' F  ',' Ne ',' Na ',' Mg ',               &
     &         ' Al ',' Si ',' P  ',' S  ',' Cl ',' Ar ',               &
     &         ' K  ',' Ca ',' Sc ',' Ti ',' V  ',' Cr ',               &
     &         ' Mn ',' Fe ',' Co ',' Ni ',' Cu ',' Zn ',               &
     &         ' Ga ',' Ge ',' As ',' Se ',' Br ',' Kr ',               &
     &         ' Rb ',' Sr ',' Y  ',' Zr ',' Nb ',' Mo ',               &
     &         ' Tc ',' Ru ',' Rh ',' Pd ',' Ag ',' Cd ',               &
     &         ' In ',' Sn ',' Sb ',' Te ',' I  ',' Xe ',               &
     &         ' Cs ',' Ba ',' La ',' Ce ',' Pr ',' Nd ',               &
     &         ' Pm ',' Sm ',' Eu ',' Gd ',' Tb ',' Dy ',               &
     &         ' Ho ',' Er ',' Tm ',' Yb ',' Lu ',' Hf ',               &
     &         ' Ta ',' W  ',' Re ',' Os ',' Ir ',' Pt ',               &
     &         ' Au ',' Hg ',' Tl ',' Pb ',' Bi ',' Po ',               &
     &         ' At ',' Rn ',' Fr ',' Ra ',' Ac ',' Th ',               &
     &         ' Pa ',' U  ',' Np ',' Pu ',' Am ',' Cm ',               &
     &         ' Bk ',' Cf ',' Es '/
!
!    Standard atomic constants for first 99 species 
!      Abundances for the first 30 from Grevesse & Sauval,
!         (1998, Space Sci. Rev. 85, 161)
!
!            Element Atomic  Solar    Std.
!                    weight abundance highest 
!                                     ionization stage 
      DATA D/ 1.008, 1.000D0, 2.,                                       &
     &        4.003, 1.00D-1, 3.,                                       &
     &        6.941, 1.26D-11, 4.,                                      &
     &        9.012, 2.51D-11, 5.,                                      &
     &       10.810, 5.0D-10, 5.,                                       &
     &       12.011, 3.31D-4, 5.,                                       &
     &       14.007, 8.32D-5, 5.,                                       &
     &       16.000, 6.76D-4, 5.,                                       &
     &       18.918, 3.16D-8, 5.,                                       &
     &       20.179, 1.23D-4, 5.,                                       &
     &       22.990, 2.14D-6, 5.,                                       &
     &       24.305, 3.80D-5, 5.,                                       &
     &       26.982, 2.95D-6, 5.,                                       &
     &       28.086, 3.55D-5, 5.,                                       &
     &       30.974, 2.82D-7, 5.,                                       &
     &       32.060, 2.14D-5, 5.,                                       &
     &       35.453, 3.16D-7, 5.,                                       &
     &       39.948, 2.52D-6, 5.,                                       &
     &       39.098, 1.32D-7, 5.,                                       &
     &       40.080, 2.29D-6, 5.,                                       &
     &       44.956, 1.48D-9, 5.,                                       &
     &       47.900, 1.05D-7, 5.,                                       &
     &       50.941, 1.00D-8, 5.,                                       &
     &       51.996, 4.68D-7, 5.,                                       &
     &       54.938, 2.45D-7, 5.,                                       &
     &       55.847, 3.16D-5, 5.,                                       &
     &       58.933, 8.32D-8, 5.,                                       &
     &       58.700, 1.78D-6, 5.,                                       &
     &       63.546, 1.62D-8, 5.,                                       &
     &       65.380, 3.98D-8, 5.,                                       &
     &       69.72 ,   1.34896324e-09  ,  3.,                           &
     &       72.60 ,   4.26579633e-09  ,  3.,                           &
     &       74.92 ,   2.34422821e-10  ,  3.,                           &
     &       78.96 ,   2.23872066e-09  ,  3.,                           &
     &       79.91 ,   4.26579633e-10  ,  3.,                           &
     &       83.80 ,   1.69824373e-09  ,  3.,                           &
     &       85.48 ,   2.51188699e-10  ,  3.,                           &
     &       87.63 ,   8.51138173e-10  ,  3.,                           &
     &       88.91 ,   1.65958702e-10  ,  3.,                           &
     &       91.22 ,   4.07380181e-10  ,  3.,                           &
     &       92.91 ,   2.51188630e-11  ,  3.,                           &
     &       95.95 ,   9.12010923e-11  ,  3.,                           &
     &       99.00 ,   1.00000000e-24  ,  3.,                           &
     &       101.1 ,   6.60693531e-11  ,  3.,                           &
     &       102.9 ,   1.23026887e-11  ,  3.,                           &
     &       106.4 ,   5.01187291e-11  ,  3.,                           &
     &       107.9 ,   1.73780087e-11  ,  3.,                           &
     &       112.4 ,   5.75439927e-11  ,  3.,                           &
     &       114.8 ,   6.60693440e-12  ,  3.,                           &
     &       118.7 ,   1.38038460e-10  ,  3.,                           &
     &       121.8 ,   1.09647810e-11  ,  3.,                           &
     &       127.6 ,   1.73780087e-10  ,  3.,                           &
     &       126.9 ,   3.23593651e-11  ,  3.,                           &
     &       131.3 ,   1.69824373e-10  ,  3.,                           &
     &       132.9 ,   1.31825676e-11  ,  3.,                           &
     &       137.4 ,   1.62181025e-10  ,  3.,                           &
     &       138.9 ,   1.58489337e-11  ,  3.,                           &
     &       140.1 ,   4.07380293e-11  ,  3.,                           &
     &       140.9 ,   6.02559549e-12  ,  3.,                           &
     &       144.3 ,   2.95120943e-11  ,  3.,                           &
     &       147.0 ,   1.00000000e-24  ,  3.,                           &
     &       150.4 ,   9.33254366e-12  ,  3.,                           &
     &       152.0 ,   3.46736869e-12  ,  3.,                           &
     &       157.3 ,   1.17489770e-11  ,  3.,                           &
     &       158.9 ,   2.13796216e-12  ,  3.,                           &
     &       162.5 ,   1.41253747e-11  ,  3.,                           &
     &       164.9 ,   3.16227767e-12  ,  3.,                           &
     &       167.3 ,   8.91250917e-12  ,  3.,                           &
     &       168.9 ,   1.34896287e-12  ,  3.,                           &
     &       173.0 ,   8.91250917e-12  ,  3.,                           &
     &       175.0 ,   1.31825674e-12  ,  3.,                           &
     &       178.5 ,   5.37031822e-12  ,  3.,                           &
     &       181.0 ,   1.34896287e-12  ,  3.,                           &
     &       183.9 ,   4.78630102e-12  ,  3.,                           &
     &       186.3 ,   1.86208719e-12  ,  3.,                           &
     &       190.2 ,   2.39883290e-11  ,  3.,                           &
     &       192.2 ,   2.34422885e-11  ,  3.,                           &
     &       195.1 ,   4.78630036e-11  ,  3.,                           &
     &       197.0 ,   6.76082952e-12  ,  3.,                           &
     &       200.6 ,   1.23026887e-11  ,  3.,                           &
     &       204.4 ,   6.60693440e-12  ,  3.,                           &
     &       207.2 ,   1.12201834e-10  ,  3.,                           &
     &       209.0 ,   5.12861361e-12  ,  3.,                           &
     &       210.0 ,   1.00000000e-24  ,  3.,                           &
     &       211.0 ,   1.00000000e-24  ,  3.,                           &
     &       222.0 ,   1.00000000e-24  ,  3.,                           &
     &       223.0 ,   1.00000000e-24  ,  3.,                           &
     &       226.1 ,   1.00000000e-24  ,  3.,                           &
     &       227.1 ,   1.00000000e-24  ,  3.,                           &
     &       232.0 ,   1.20226443e-12  ,  3.,                           &
     &       231.0 ,   1.00000000e-24  ,  3.,                           &
     &       238.0 ,   3.23593651e-13  ,  3.,                           &
     &       237.0 ,   1.00000000e-24  ,  0.,                           &
     &       244.0 ,   1.00000000e-24  ,  0.,                           &
     &       243.0 ,   1.00000000e-24  ,  0.,                           &
     &       247.0 ,   1.00000000e-24  ,  0.,                           &
     &       247.0 ,   1.00000000e-24  ,  0.,                           &
     &       251.0 ,   1.00000000e-24  ,  0.,                           &
     &       254.0 ,   1.00000000e-24  ,  0./
!
!
!     Ionization potentials for first 99 species:
!
!     Element Ionization potentials (eV) 
!              I     II      III     IV       V     VI     VII    VIII
      DATA XI/                                                          &
     &       13.595,  0.   ,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,&
     &       24.580, 54.400,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,&
     &        5.392, 75.619,122.451,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,&
     &        9.322, 18.206,153.850,217.713,  0.  ,  0.  ,  0.  ,  0.  ,&
     &        8.296, 25.149, 37.920,259.298,340.22,  0.  ,  0.  ,  0.  ,&
     &       11.264, 24.376, 47.864, 64.476,391.99,489.98,  0.  ,  0.  ,&
     &       14.530, 29.593, 47.426, 77.450, 97.86,551.93,667.03,  0.  ,&
     &       13.614, 35.108, 54.886, 77.394,113.87,138.08,739.11,871.39,&
     &       17.418, 34.980, 62.646, 87.140,114.21,157.12,185.14,953.6 ,&
     &       21.559, 41.070, 63.500, 97.020,126.30,157.91,207.21,239.0 ,&
     &        5.138, 47.290, 71.650, 98.880,138.37,172.09,208.44,264.16,&
     &        7.664, 15.030, 80.120,102.290,141.23,186.49,224.9 ,265.96,&
     &        5.984, 18.823, 28.440,119.960,153.77,190.42,241.38,284.53,&
     &        8.151, 16.350, 33.460, 45.140,166.73,205.11,246.41,303.07,&
     &       10.484, 19.720, 30.156, 51.354, 65.01,220.41,263.31,309.26,&
     &       10.357, 23.400, 35.000, 47.290, 72.50, 88.03,280.99,328.8 ,&
     &       12.970, 23.800, 39.900, 53.500, 67.80, 96.7 ,114.27,348.3 ,&
     &       15.755, 27.620, 40.900, 59.790, 75.00, 91.3 ,124.0 ,143.46,&
     &        4.339, 31.810, 46.000, 60.900, 82.6 , 99.7 ,118.0 ,155.0 ,&
     &        6.111, 11.870, 51.210, 67.700, 84.39,109.0 ,128.0 ,147.0 ,&
     &        6.560, 12.890, 24.750, 73.900, 92.0 ,111.1 ,138.0 ,158.7 ,&
     &        6.830, 13.630, 28.140, 43.240, 99.8 ,120.0 ,140.8 ,168.5 ,&
     &        6.740, 14.200, 29.700, 48.000, 65.2 ,128.9 ,151.0 ,173.7 ,&
     &        6.763, 16.490, 30.950, 49.600, 73.0 , 90.6 ,161.1 ,184.7 ,&
     &        7.432, 15.640, 33.690, 53.000, 76.0 , 97.0 ,119.24,196.46,&
     &        7.870, 16.183, 30.652, 54.800, 75.0 , 99.1 ,125.0 ,151.06,&
     &        7.860, 17.060, 33.490, 51.300, 79.5 ,102.0 ,129.0 ,157.0 ,&
     &        7.635, 18.168, 35.170, 54.900, 75.5 ,108.0 ,133.0 ,162.0 ,&
     &        7.726, 20.292, 36.830, 55.200, 79.9 ,103.0 ,139.0 ,166.0 ,&
     &        9.394, 17.964, 39.722, 59.400, 82.6 ,108.0 ,134.0 ,174.0 ,&
     &        6.000,  20.509,   30.700, 99.99,99.99,99.99,99.99,99.99,  &
     &        7.89944,15.93462, 34.058, 45.715,99.99,99.99,99.99,99.99, &
     &        9.7887, 18.5892,  28.351, 99.99,99.99,99.99,99.99,99.99,  &
     &        9.750,21.500, 32.000, 99.99,99.99,99.99,99.99,99.99,      &
     &       11.839,21.600, 35.900, 99.99,99.99,99.99,99.99,99.99,      &
     &       13.995,24.559, 36.900, 99.99,99.99,99.99,99.99,99.99,      &
     &        4.175,27.500, 40.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.692,11.026, 43.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.2171,12.2236, 20.5244,60.607,99.99,99.99,99.99,99.99,   &
     &        6.63390,13.13,23.17,34.418,80.348,99.99,99.99,99.99,      &
     &        6.879,14.319, 25.039, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.099,16.149, 27.149, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.280,15.259, 30.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.364,16.759, 28.460, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.460,18.070, 31.049, 99.99,99.99,99.99,99.99,99.99,      &
     &        8.329,19.419, 32.920, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.574,21.480, 34.819, 99.99,99.99,99.99,99.99,99.99,      &
     &        8.990,16.903, 37.470, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.784,18.860, 28.029, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.342,14.627, 30.490,72.3,99.99,99.99,99.99,99.99,        &
     &        8.639,16.500, 25.299,44.2,55.7,99.99,99.99,99.99,         &
     &        9.0096,18.600, 27.96, 37.4,58.7,99.99,99.99,99.99,        &
     &       10.454,19.090, 32.000, 99.99,99.99,99.99,99.99,99.99,      &
     &       12.12984,20.975,31.05,45.,54.14,99.99,99.99,99.99,         &
     &        3.893,25.100, 35.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.210,10.000, 37.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.580,11.060, 19.169, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.650,10.850, 20.080, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.419,10.550, 23.200, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.490,10.730, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.550,10.899, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.629,11.069, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.680,11.250, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.159,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.849,11.519, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.930,11.670, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.020,11.800, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.099,11.930, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.180,12.050, 23.700, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.250,12.170, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.099,13.899, 19.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.000,14.899, 23.299, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.879,16.200, 24.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        7.86404,17.700, 25.000, 99.99,99.99,99.99,99.99,99.99,    &
     &        7.870,16.600, 26.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        8.500,17.000, 27.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        9.100,20.000, 28.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        8.95868,18.563,33.227, 99.99,99.99,99.99,99.99,99.99,     &
     &        9.220,20.500, 30.000, 99.99,99.99,99.99,99.99,99.99,      &
     &       10.430,18.750, 34.200, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.10829,20.4283,29.852,50.72,99.99,99.99,99.99,99.99,     &
     &        7.416684,15.0325,31.9373,42.33,69.,99.99,99.99,99.99,     &
     &        7.285519,16.679, 25.563,45.32,56.0,88.,99.99,99.99,       &
     &        8.430,19.000, 27.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        9.300,20.000, 29.000, 99.99,99.99,99.99,99.99,99.99,      &
     &       10.745,20.000, 30.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        4.000,22.000, 33.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        5.276,10.144, 34.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.900,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,      &
     &        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99/
!
	do 20 i=1,matom
	  adyp(i)=dyp(i)
	  do 10 j=1,8
	    axi(j,i)=xi(j,i)	
10	  continue
	  do 15 j=1,3
	    ad(j,i)=d(j,i)	
15	  continue
20	continue
        return
        end
!-----------------------------------------------------------------------
        subroutine eldat(ichemc,iat,iion,dyp,d,xi,abhyd,abhel,wm        &
     &  ,hmc,nion,abund,zip,nline,ielnd,necod)
!		extract data for our particular elements of sp. lines
!	and change abundances
!	input: 
!		ichemc-switch to change built in chemical composition
!		iat-atomic number of the sp. line
!		iion-ionization degree of the sp.line
!		dyp,d,xi 
!		nline -number of sp. lines considered
!		ielnd-switch to calculate el.num.d. and read 3.column
!	output: 
!		d -weight,abundance, max ionization stage
!       	abhyd - H abundance N/N(TOT)
!		abhel - He abundance, N(EL)/N(H) 
!		wm -mean molecular weight= sum ai*mi/sum ai
! 		hmc - mass number
!		nion - number of ions of our element considered 
!		abund - abundance, N(EL)/N(H)
!		zip - first 8 ionization potentials [EV]
!		necod - code to include element into the el.num.dens. 
	implicit none
        integer, parameter:: i4=4,i8=8
 	include 'param.inc'
	integer:: ichemc,nline,ielnd
	integer:: nion(mline),iat(mline),iion(mline),necod(matom)
 	character(4):: dyp(matom),elem(mline)
 	real(i8):: d(3,matom),xi(8,matom),zip(mion-1,mline)
	real(i8):: abhyd,abhel,wm
        real(i8):: abund(mline),hmc(mline)
	integer:: i,nichem,ii,necdi,j
	real(i8):: abii
	if(ielnd.eq.1)ichemc=1
	do i=1,matom
	  necod(i)=0
	enddo
	if(ichemc.eq.1)then
	  write(*,*)' reading file abundances'
	  open(7,file='abundances',status='old')
          read(7,*)nichem
	  do 80 i=1,nichem
	    if(ielnd.eq.1)then	
              read(7,*)ii,abii,necdi
              necod(i)=necdi
	    else
              read(7,*)ii,abii
	    endif
	    d(2,ii)=abii
80 	  continue
	  close(7)	
	endif
	abhyd=0.d0
	wm=0.d0
	do 90 i=1,matom
	  abhyd=abhyd+d(2,i)
	  wm=wm+d(1,i)*d(2,i)
90 	continue
	abhyd=1.d0/abhyd
	wm=wm*abhyd
	abhel=d(2,2)
   	write(2,*)' No. of sp. lines considered:',nline
	if(nline.gt.0)then
	  do 110 i=1,nline
	    nion(i)=int(d(3,iat(i))+0.5d0)
	    if(iion(i).gt.nion(i).or.iion(i).gt.mion)then
	      write(*,*)i,'-th sp. line cannot be calculated,'
              write(*,*)' error iion>nion'
              write(3,*)' 1 error iion>nion'
	      call exit(1)
	    endif
	    if(iion(i).gt.5)then
	      write(*,*)i,'-th sp. line, partition f. not available,'
              write(*,*)' warning: p.f. will be set=1'
              write(3,*)' 1 warning: p.f. will be set=1'
	    endif
	    abund(i)=d(2,iat(i))
	    hmc(i)=d(1,iat(i))
	    do 100 j=1,8
	      zip(j,i)=xi(j,iat(i))	
100	    continue
	    elem(i)=dyp(iat(i))
	    write(2,120)elem(i),iat(i),hmc(i),nion(i),abund(i)
110	  continue
	endif
	write(2,130)abhyd,abhel,wm
120     format(' elem=',a4,' iat=',i3,' hmc=',f8.4                      &
     &  ,' nion=',i3,' abund=',e9.3)
130	format(' abhyd=',f5.3,' abhel=',f5.3,' wm=',e9.3)
	return
	end
!-----------------------------------------------------------------------
        subroutine ophyd(nbod,temp,ekonc,akonc,hi,hneg,abhyd            &
     &  ,ophbf,ophff,ophrs,ophn,hipf,wlab)
!	calculates HI bound-free, free-free, Rayleigh opacity
!	input: nbod,temp,ekonc,akonc,hi,hneg,abhyd,hipf,wlab[A]
!	output: ophbf,ophff,ophrs,ophn
	implicit none
        integer, parameter:: i4=4,i8=8
        integer:: nbod
	include 'param.inc'
	real(i8):: temp(ndim),ekonc(ndim),akonc(ndim),hipf(ndim,2)
        real(i8):: hi(ndim),hneg(ndim)
        real(i8):: ophbf(ndim),ophff(ndim),ophrs(ndim),ophn(ndim)
	real(i8):: abhyd,wlab
	real(i8):: hnbf,hnff,gfree,gaunt
	real(i8):: bol,clight,hjed,plank,rydbh,freq0,pot,frray,frjump
	real(i8):: hii,tk,secor,en,chi,en3,gau,gau3,xe,ophnbf,ophnff
	integer:: i,in,j
        bol=1.3806503d-16
        clight=2.99792458d10
        hjed=1.66053873d-24
        plank=6.62606876d-27
	rydbh=109677.585d0
	freq0=clight/wlab*1.D8
	pot=clight*plank*rydbh
	frray=2.463d15
	do 10 i=1,16
	  frjump=clight*rydbh*1.d0/dble(i*i)
	  if(frjump.lt.freq0)goto 20
10	continue
20	in=i
	do 70 i=1,nbod
!		HI bound-free opacity
	  hii=akonc(i)*abhyd-hi(i)
	  tk=bol*temp(i)
	  secor=(1.d0-dexp(-plank*freq0/tk))
	  if(in.lt.16)then
	    ophbf(i)=0.d0
	    do 60 j=1,3
	      en=dble(in+j-1)
	      chi=pot*(1.d0-1.d0/(en*en))
	      en3=en*en*en
	      gau=gaunt(in+j-1,freq0)
	      ophbf(i)=ophbf(i)+gau/en3*dexp(-chi/tk)
60          continue
	    en=dble(in+3)
	    chi=pot*(1.d0-1.d0/(en*en))
            ophbf(i)=tk/pot*0.5d0*(dexp(-chi/tk)-dexp(-pot/tk))         &
     &      +ophbf(i)
!           ophbf(i)=ophbf(i)*1.044d-26*wlab*wlab*wlab*2.d0/hipf(i,1)   &
!     &      *secor*hi(i)
	    ophbf(i)=ophbf(i)*2.8154d29/freq0/freq0/freq0*2.d0
	    ophbf(i)=ophbf(i)/hipf(i,1)*secor*hi(i)
	  else
	    ophbf(i)=0.d0
	  endif
!               HI free-free opacity
	  gau3=gfree(temp(i),freq0)
	  ophff(i)=3.69d8*gau3/(freq0*freq0*freq0*dsqrt(temp(i)))
	  ophff(i)=ophff(i)*ekonc(i)*hii*secor
!               HI Rayleigh scattering on neutral hydrogen
!		Kurucz R.L., 1970, SAO Special Report 309
!		note that scattering has no stimulated emission
	  xe=1.d-16*(min(freq0,frray)/clight)**2
	  ophrs(i)=(5.799d-13+(1.422d-6+2.784d0*xe)*xe)*xe*xe
	  ophrs(i)=ophrs(i)*2.d0*hi(i)/hipf(i,1)
!		H- bound-free opacity
	  ophnbf=hnbf(freq0)*hneg(i)*secor
!		H- free-free opacity
	  ophnff=hnff(freq0,ekonc(i),temp(i))*hi(i)*2.d0/hipf(i,1)
	  ophn(i)=ophnbf+ophnff	  
70	continue
	return
	end
!-----------------------------------------------------------------------
	function hnbf(x)
	implicit none
        integer, parameter:: i4=4,i8=8
        real(i8):: hnbf 
	real(i8):: x
	real(i8):: pom1,pom2,pom3
!	H- bound-free x-section from Kurucz 
!	without the stimulated emission
!	input: x-frequency [Hz]
!	output: cross-section [cm**2] per H- particle
!	works OK
	if(x.gt.2.111d14)then
	  pom1=(-5.519d27+4.808d41/x)
	  pom2=(1.481d13+pom1/x)
	  pom3=(5.358d-3+pom2/x)
	  hnbf=6.801d-20+pom3/x
!	elseif(x.gt.1.8259d14)then
!	this original Kurucz expression has a small bug 
!	it produces negative values at about 16400A 
!	which is near the ionization limit that is why I cut-off
!	before at 1.834d14
	elseif(x.gt.1.834d14)then
	  pom1=-1.251d-1+1.052d13/x
	  hnbf=3.695d-16+pom1/x
	else  
	  hnbf=0.d0
	endif  
	return
	end function hnbf
!-----------------------------------------------------------------------
	function hnff(x,en,t)
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: hnff
	real(i8):: x,en,t
	real(i8):: pom1
!	H- free-free opacity from Kurucz
!	with the stimulated emission
!	input: 
!	x -frequency [Hz], en -electron num. density [cm^-3]
!	t -temperature[K]
!	output: cross-section [cm**2] per HI (in ground state)
!	works OK
	if(t.lt.1200.d0)then
	  hnff=0.d0
	elseif(t.gt.2.d4)then
	  hnff=0.d0
	elseif(x.lt.1.5d12)then
!       artificial constraint at 200mic to reach agreement with Gray 
!	within a factor of 2	
	  hnff=0.d0
	elseif(x.gt.4.d15)then
!	artificial constraint since it is already in Lyman cont.	
	  hnff=0.d0
	else  
	  pom1=4.3748d-10-2.5993d-7/t
	  hnff=(1.3727d-25+pom1/x)/x*en
	endif
	return
	end function hnff
!-----------------------------------------------------------------------
        subroutine hihei(ipart,nbod,temp,ekonc,akonc                    &
     &  ,hipf,hi,hneg,h2,hei,abhyd,abhel)
!	calculates H, He number densities and partition functions
!	input: ipart,nbod,temp,ekonc,akonc
!       	abhyd - H abundance N/N(TOT)
!		abhel - He abundance, N(EL)/N(H) 
!	output: hipf - HI partition function
!		hi,hneg,h2,hei - HI, H-, H2, HeI number densities
	implicit none
        integer, parameter:: i4=4,i8=8
        include 'param.inc'
        real(i8):: hipf(ndim,2),hi(ndim),hneg(ndim),h2(ndim),hei(ndim)
        real(i8):: TEMP(NDIM),EKONC(NDIM),AKONC(NDIM)
	real(i8):: abhyd,abhel
	integer:: ipart,nbod
	real(i8):: sahan,sah2,sah
        real(i8):: hepf(ndim,3),STAHI(6,2),STAHEI(6,3)
	real(i8):: sala,pothm,stavhm,pothi,rt3,rt2,rt1,b,b2,c,d,bot
	integer:: i,j
!		POLYNOMIAL PARTITION FUNCTION COEFFICIENTS FOR 
!		HI,HII   A    HEI,HEII,HEIII 
!		AFTER IRWIN, 1981, ApJS 45,621	
        DATA STAHI/-2.61655891D2,1.63428326D2,-4.06133526D1             &
     &  ,5.03282928D0,-3.10998364D-1,7.66654594D-3,0.0D0,0.0D0,0.0D0    &
     &  ,0.0D0,0.0D0,0.0D0/
        DATA STAHEI/-3.76575219D-1,2.33951687D-1,-5.79755525D-2         &
     &  ,7.16333160D-3,-4.41302573D-4,1.08442997D-5,6.93147179D-1       &
     &  ,9.29636701D-10,-2.30049742D-10,2.83829746D-11,-1.74590774D-12  &
     &  ,4.28355287D-14,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
	if(ipart.eq.1)then
          CALL PFDWOR(NDIM,NBOD,2,2,1,TEMP,EKONC,hipf)
          CALL PFDWOR(NDIM,NBOD,3,3,2,TEMP,EKONC,hepf)
	else
	  do I=1,NBOD
	    SALA=DLOG(TEMP(I))
!		HIPF - H PARTITION FUNCTION
	    do J=1,2
              HIPF(i,j)=DEXP(STAHI(1,J)+STAHI(2,J)*SALA                 &
     &        +STAHI(3,J)*SALA*SALA+STAHI(4,J)*SALA**3                  &
     &        +STAHI(5,J)*SALA**4+STAHI(6,J)*SALA**5)
	    enddo
!               HEPF - HE PARTITION FUNCTION
	    do J=1,3
              HEPF(i,j)=DEXP(STAHEI(1,J)+STAHEI(2,J)*SALA+STAHEI(3,J)   &
     &        *SALA*SALA+STAHEI(4,J)*SALA**3+STAHEI(5,J)*SALA**4        &
     &        +STAHEI(6,J)*SALA**5)
	    enddo
	  enddo
	endif
!		HNEG - H- NUMBER DENSITY
!               HI - HI NUMBER DENSITY
!               HEI - HeI NUMBER DENSITY
!               H2 - H2 NUMBER DENSITY	
        pothm=0.7552d0
        stavhm=1.d0
        pothi=13.598d0
	do i=1,nbod
	  if(temp(i).lt.1.d2)then
!  	    H2
	    hi(i)=0.d0
	    hneg(i)=0.d0
	    h2(i)=akonc(i)*abhyd/2.d0
 	  elseif(temp(i).lt.2.d4)then
!	    HII+HI+H-+H2 	  
            rt3=sahan(temp(i),pothi,hipf(i,1),hipf(i,2))
            rt2=sahan(temp(i),pothm,stavhm,hipf(i,1))          
	    rt1=sah2(temp(i))
	    b=1.d0+ekonc(i)/rt2+rt3/ekonc(i)
	    b=b/2.d0/rt1
	    c=-akonc(i)*abhyd/2.d0/rt1
	    b2=b*b
	    d=b2-4.d0*c
	    if(-4.d0*c/b2.lt.1.d-10)then
!           then avoid substracting two similar numbers -b+b
              hi(i)=-c/b
            else
              hi(i)=(-b+dsqrt(d))/2.d0
            endif
	    hneg(i)=hi(i)*ekonc(i)/rt2
	    h2(i)=rt1*hi(i)*hi(i)
	  else
!	    HII+HI+H-	  	  
            rt2=sahan(temp(i),pothm,stavhm,hipf(i,1))
            rt3=sahan(temp(i),pothi,hipf(i,1),hipf(i,2))
            bot=1.d0+ekonc(i)/rt2+rt3/ekonc(i)
            hi(i)=akonc(i)*abhyd/bot
            hneg(i)=hi(i)*ekonc(i)/rt2
            h2(i)=0.d0
	  endif
          HEI(I)=1.D0/(1.D0+SAH(TEMP(I),EKONC(I),24.587D0,HEPF(i,1)     &
     &    ,HEPF(i,2))*(1.D0+SAH(TEMP(I),EKONC(I),54.416D0,HEPF(i,2)     &
     &    ,HEPF(i,3))))*AKONC(I)*ABHEL*abhyd
	enddo
	return
	end
!-----------------------------------------------------------------------
      SUBROUTINE HUNT(XX,N,X,JLO)
!	taken from Numerical recipes
!	given field XX and a value X subroutine returns JLO such that
!	XX(JLO)<X<XX(JLO+1) 
!	JLO on input is taken as the intial guess for JLO on output
	implicit none
        integer, parameter:: i4=4,i8=8
!?? 	potential problem with dynamical declaration
      	real(i8):: XX(N)
	real(i8):: x
	integer:: n,jlo
	integer:: jhi,inc,jm
	logical:: ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END
!-----------------------------------------------------------------------
        subroutine rot1d1(theta,alpha,dcut1,temp0,ane0,iii,jjj          &
     &  ,vxst,vyst,vzst,vxstr,vystr,vzstr                               &
     &  ,xcp,ycp,zcp,xcpr,ycpr,zcpr                                     &
     &  ,vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr                               &
     &  ,ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3                      &
     &  ,ndim1,ndim2,ndim3,nbod1,nbod2,nbod3                            &
     &  ,ar,at,az,far,fat,faz                                           &
     &  ,atemp,adustt,adustd,adens,ane,avr,avt,avz,avtrb,lshade         &
     &  ,ftemp,fdustt,fdustd,fdens,fne,fvr,fvt,fvz,fvtrb,kshade)
!	linear interpolation in the transparent object
!	nearest neighbour in the nontransparent object     
!	to transform (rotate) the line of sight grid (x,y,z) to the body
!	frozen grid (x'',y'',z'') and interpolate state quantities from 
!	the rotating body frozen coordinates (RBFC) to the line of sight 
!	coordinates (LSC) and transform the velocity field vectors 
!	from RBFC to LSC
!	Note: shadows are not supported/interpolated and lshade=3
!	input: 	
!         theta,alpha,dcut1
!	  temp0,ane0 - out of grid (sphere) reasonable values
!	  iii,jjj -identify x,y the ray
!         vxst,vyst,vzst - net velocity of the primary
!	  xcp,ycp,zcp - position of the secondary
!	  vxcp,vycp,vzcp - velocity of the secondary
!	  ar,at,az -coord. in which the model is necessary (LSC)
!	  far,fat,faz -coord. in which the model is defined (RBFC)
!	  ftemp,fdustt,fdustd,fdens,fne,fvr,fvt,fvz,fvtrb,kshade -model
!	output:
!	  atemp,adustt,adustd,adens,ane,avr,avt,avz,avtrb,lshade
!         (^rotated model^)
!         vxstr,vystr,vzstr- rotated net velocity of the primary
!	  xcpr,ycpr,zcpr -rotated position of the secondary
!	  vxcpr,vycpr,vzcpr - rotated velocity of the secondary 
!	parameter (ndim1=100,ndim2=100,ndim3=200)
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
	real(i8):: ar(ndim1),at(ndim2),az(ndim3)
        real(i8):: far(ndimf1),fat(ndimf2),faz(ndimf3)
        real(i8):: atemp(ndim3),adustt(ndim3),adustd(ndim3)
        real(i8):: adens(ndim3),ane(ndim3)
        real(i8):: avr(ndim3),avt(ndim3)
        real(i8):: avz(ndim3),avtrb(ndim3)
        real(i8):: ftemp(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustt(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustd(ndimf1,ndimf2,ndimf3)
        real(i8):: fdens(ndimf1,ndimf2,ndimf3),fne(ndimf1,ndimf2,ndimf3)
        real(i8):: fvr(ndimf1,ndimf2,ndimf3),fvt(ndimf1,ndimf2,ndimf3)
        real(i8):: fvz(ndimf1,ndimf2,ndimf3),fvtrb(ndimf1,ndimf2,ndimf3)
	integer(i1):: kshade(ndimf1,ndimf2,ndimf3),lshade(ndim3)
	integer:: iii,jjj,ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3
	integer:: ndim1,ndim2,ndim3,nbod1,nbod2,nbod3
	real(i8):: theta,alpha,dcut1,temp0,ane0
	real(i8):: vxst,vyst,vzst,vxstr,vystr,vzstr
	real(i8):: xcp,ycp,zcp,xcpr,ycpr,zcpr
	real(i8):: vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr
	real(i8):: ctheta,stheta,calpha,salpha
	real(i8):: fx,fy,fz,dist1,dist,dx,dy,dz,tt,uu,vv,tt1,uu1,vv1
	integer::ii,jj,kk,i,j,k,ii1,jj1,kk1
	ctheta=dcos(theta)
	stheta=dsin(theta)
	calpha=dcos(alpha)
	salpha=dsin(alpha)
!	initialization of first values to ii,jj,kk (index in RBFC)
	ii=1
	jj=1
	kk=1
	i=iii
	j=jjj
	do 10 k=1,nbod3
!	  calculation of RBFC coordinates corresponding to the grid 
!	  points of the line of sight grid 
          call rotz(ar(i),at(j),az(k),fx,fy,fz                          &
     &    ,ctheta,stheta,calpha,salpha)
!	  	interpolation from body frozen grid to required 
!	  line of sight grid  :
          if(fx.gt.far(1).and.fx.lt.far(nbodf1).and.                    &
     &    fy.gt.fat(1).and.fy.lt.fat(nbodf2).and.                       &
     &    fz.gt.faz(1).and.fz.lt.faz(nbodf3))then
	    call hunt(far,nbodf1,fx,ii)
	    if(ii.lt.1.or.ii.ge.nbodf1)then
	      if(ii.eq.nbodf1.and.fx.eq.far(nbodf1))then
		ii=nbodf1-1
	      else		
	        write(*,*)' compiler/maschine dependent error'
	        write(3,*)' 1 compiler/maschine dependent error'
                write(*,'(a,4i4)')' -out of grid at ii,i,j,k',ii,i,j,k
	      endif
	    endif
	    call hunt(fat,nbodf2,fy,jj)
	    if(jj.lt.1.or.jj.ge.nbodf2)then
	      if(jj.eq.nbodf2.and.fy.eq.fat(nbodf2))then
		jj=nbodf2-1
	      else		
	        write(*,*)' compiler/maschine dependent error'
	        write(3,*)' 1 compiler/maschine dependent error'
                write(*,'(a,4i4)')' -out of grid at jj,i,j,k',jj,i,j,k
	      endif
	    endif
	    call hunt(faz,nbodf3,fz,kk)
	    if(kk.lt.1.or.kk.ge.nbodf3)then
	      if(kk.eq.nbodf3.and.fz.eq.faz(nbodf3))then
		kk=nbodf3-1
	      else		
	        write(*,*)' compiler/maschine dependent error'
	        write(3,*)' 1 compiler/maschine dependent error'
                write(*,'(a,4i4)')' -out of grid at kk,i,j,k',kk,i,j,k
	      endif
	    endif
!		abnormal interpolation-nearest neighbour value
!		in the untransparent object or on its edge.
!		density of an untransparent object dens>dcut1
            if(fdens(ii,jj,kk).gt.dcut1                                 &
     &      .or.fdens(ii,jj,kk+1).gt.dcut1                              &
     &      .or.fdens(ii,jj+1,kk).gt.dcut1                              &
     &      .or.fdens(ii,jj+1,kk+1).gt.dcut1                            &
     &      .or.fdens(ii+1,jj,kk).gt.dcut1                              &
     &      .or.fdens(ii+1,jj,kk+1).gt.dcut1                            &
     &      .or.fdens(ii+1,jj+1,kk).gt.dcut1                            &
     &      .or.fdens(ii+1,jj+1,kk+1).gt.dcut1)then
	      dist1=1.d60
	      do 5 ii1=ii,ii+1
              do 4 jj1=jj,jj+1
              do 3 kk1=kk,kk+1	
		dx=fx-far(ii1)
                dy=fy-fat(jj1)
                dz=fz-faz(kk1)
		dist=dx*dx+dy*dy+dz*dz
		if(dist.lt.dist1)then
	          atemp(k)=ftemp(ii1,jj1,kk1)
	          adustt(k)=fdustt(ii1,jj1,kk1)
	          adustd(k)=fdustd(ii1,jj1,kk1)
	          adens(k)=fdens(ii1,jj1,kk1)
	          ane(k)=fne(ii1,jj1,kk1)
	          avr(k)=fvr(ii1,jj1,kk1)
	          avt(k)=fvt(ii1,jj1,kk1)
	          avz(k)=fvz(ii1,jj1,kk1)
	          avtrb(k)=fvtrb(ii1,jj1,kk1)
		  lshade(k)=3
		  dist1=dist
	        endif
3	      continue
4             continue
5             continue
	      goto 10	
	    endif
!		normal interpolation
	    tt=(fx-far(ii))/(far(ii+1)-far(ii))
	    uu=(fy-fat(jj))/(fat(jj+1)-fat(jj))
	    vv=(fz-faz(kk))/(faz(kk+1)-faz(kk))
	    tt1=1.d0-tt
	    uu1=1.d0-uu
	    vv1=1.d0-vv	
            atemp(k)=tt1*uu1*vv1*ftemp(ii,jj,kk)                        &
     &        +tt*uu1*vv1*ftemp(ii+1,jj,kk)                             &
     &        +tt*uu*vv1*ftemp(ii+1,jj+1,kk)                            &
     &        +tt1*uu*vv1*ftemp(ii,jj+1,kk)                             &
     &        +tt1*uu1*vv*ftemp(ii,jj,kk+1)                             &
     &        +tt*uu1*vv*ftemp(ii+1,jj,kk+1)                            &
     &        +tt*uu*vv*ftemp(ii+1,jj+1,kk+1)                           &
     &        +tt1*uu*vv*ftemp(ii,jj+1,kk+1)
            adustt(k)=tt1*uu1*vv1*fdustt(ii,jj,kk)                      &
     &        +tt*uu1*vv1*fdustt(ii+1,jj,kk)                            &
     &        +tt*uu*vv1*fdustt(ii+1,jj+1,kk)                           &
     &        +tt1*uu*vv1*fdustt(ii,jj+1,kk)                            &
     &        +tt1*uu1*vv*fdustt(ii,jj,kk+1)                            &
     &        +tt*uu1*vv*fdustt(ii+1,jj,kk+1)                           &
     &        +tt*uu*vv*fdustt(ii+1,jj+1,kk+1)                          &
     &        +tt1*uu*vv*fdustt(ii,jj+1,kk+1)
            adustd(k)=tt1*uu1*vv1*fdustd(ii,jj,kk)                      &
     &        +tt*uu1*vv1*fdustd(ii+1,jj,kk)                            &
     &        +tt*uu*vv1*fdustd(ii+1,jj+1,kk)                           &
     &        +tt1*uu*vv1*fdustd(ii,jj+1,kk)                            &
     &        +tt1*uu1*vv*fdustd(ii,jj,kk+1)                            &
     &        +tt*uu1*vv*fdustd(ii+1,jj,kk+1)                           &
     &        +tt*uu*vv*fdustd(ii+1,jj+1,kk+1)                          &
     &        +tt1*uu*vv*fdustd(ii,jj+1,kk+1)
            adens(k)=tt1*uu1*vv1*fdens(ii,jj,kk)                        &
     &        +tt*uu1*vv1*fdens(ii+1,jj,kk)                             &
     &        +tt*uu*vv1*fdens(ii+1,jj+1,kk)                            &
     &        +tt1*uu*vv1*fdens(ii,jj+1,kk)                             &
     &        +tt1*uu1*vv*fdens(ii,jj,kk+1)                             &
     &        +tt*uu1*vv*fdens(ii+1,jj,kk+1)                            &
     &        +tt*uu*vv*fdens(ii+1,jj+1,kk+1)                           &
     &        +tt1*uu*vv*fdens(ii,jj+1,kk+1)
            ane(k)=tt1*uu1*vv1*fne(ii,jj,kk)                            &
     &        +tt*uu1*vv1*fne(ii+1,jj,kk)                               &
     &        +tt*uu*vv1*fne(ii+1,jj+1,kk)                              &
     &        +tt1*uu*vv1*fne(ii,jj+1,kk)                               &
     &        +tt1*uu1*vv*fne(ii,jj,kk+1)                               &
     &        +tt*uu1*vv*fne(ii+1,jj,kk+1)                              &
     &        +tt*uu*vv*fne(ii+1,jj+1,kk+1)                             &
     &        +tt1*uu*vv*fne(ii,jj+1,kk+1)
            avr(k)=tt1*uu1*vv1*fvr(ii,jj,kk)                            &
     &        +tt*uu1*vv1*fvr(ii+1,jj,kk)                               &
     &        +tt*uu*vv1*fvr(ii+1,jj+1,kk)                              &
     &        +tt1*uu*vv1*fvr(ii,jj+1,kk)                               &
     &        +tt1*uu1*vv*fvr(ii,jj,kk+1)                               &
     &        +tt*uu1*vv*fvr(ii+1,jj,kk+1)                              &
     &        +tt*uu*vv*fvr(ii+1,jj+1,kk+1)                             &
     &        +tt1*uu*vv*fvr(ii,jj+1,kk+1)
            avt(k)=tt1*uu1*vv1*fvt(ii,jj,kk)                            &
     &        +tt*uu1*vv1*fvt(ii+1,jj,kk)                               &
     &        +tt*uu*vv1*fvt(ii+1,jj+1,kk)                              &
     &        +tt1*uu*vv1*fvt(ii,jj+1,kk)                               &
     &        +tt1*uu1*vv*fvt(ii,jj,kk+1)                               &
     &        +tt*uu1*vv*fvt(ii+1,jj,kk+1)                              &
     &        +tt*uu*vv*fvt(ii+1,jj+1,kk+1)                             &
     &        +tt1*uu*vv*fvt(ii,jj+1,kk+1)
            avz(k)=tt1*uu1*vv1*fvz(ii,jj,kk)                            &
     &        +tt*uu1*vv1*fvz(ii+1,jj,kk)                               &
     &        +tt*uu*vv1*fvz(ii+1,jj+1,kk)                              &
     &        +tt1*uu*vv1*fvz(ii,jj+1,kk)                               &
     &        +tt1*uu1*vv*fvz(ii,jj,kk+1)                               &
     &        +tt*uu1*vv*fvz(ii+1,jj,kk+1)                              &
     &        +tt*uu*vv*fvz(ii+1,jj+1,kk+1)                             &
     &        +tt1*uu*vv*fvz(ii,jj+1,kk+1)
            avtrb(k)=tt1*uu1*vv1*fvtrb(ii,jj,kk)                        &
     &        +tt*uu1*vv1*fvtrb(ii+1,jj,kk)                             &
     &        +tt*uu*vv1*fvtrb(ii+1,jj+1,kk)                            &
     &        +tt1*uu*vv1*fvtrb(ii,jj+1,kk)                             &
     &        +tt1*uu1*vv*fvtrb(ii,jj,kk+1)                             &
     &        +tt*uu1*vv*fvtrb(ii+1,jj,kk+1)                            &
     &        +tt*uu*vv*fvtrb(ii+1,jj+1,kk+1)                           &
     &        +tt1*uu*vv*fvtrb(ii,jj+1,kk+1)
	    lshade(k)=3
	  else
	    atemp(k)=temp0
	    adustt(k)=temp0
	    adustd(k)=0.d0
	    adens(k)=0.d0
	    ane(k)=ane0
	    avr(k)=0.d0
	    avt(k)=0.d0
	    avz(k)=0.d0
	    avtrb(k)=0.d0
	    lshade(k)=3
	  endif
10	continue
!	back transform of vector quantities
	do k=1,nbod3
          call rotzb(avr(k),avt(k),avz(k)                               &
     &    ,avr(k),avt(k),avz(k),ctheta,stheta,calpha,salpha)
	enddo
!	back transform of the velocity of the center of mass 
!	of the primary explicitely because of the scattered light 
!	treatment
        call rotzb(vxst,vyst,vzst,vxstr,vystr,vzstr                     &
     &  ,ctheta,stheta,calpha,salpha)
!	back transform of the center of the secondary explicitely 
!	because	of the limb darkening treatement
        call rotzb(xcp,ycp,zcp,xcpr,ycpr,zcpr                           &
     &  ,ctheta,stheta,calpha,salpha)
!       back transform of the velocity of the secondary explicitely 
!	because of the scattered light treatment
        call rotzb(vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr                     &
     &  ,ctheta,stheta,calpha,salpha)
        return
        end
!-----------------------------------------------------------------------
        subroutine rot1d2(theta,alpha,dcut1,temp0,ane0,iii,jjj          &
     &  ,vxst,vyst,vzst,vxstr,vystr,vzstr                               &
     &  ,xcp,ycp,zcp,xcpr,ycpr,zcpr                                     &
     &  ,vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr                               &
     &  ,ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3                      &
     &  ,ndim1,ndim2,ndim3,nbod1,nbod2,nbod3                            &
     &  ,ar,at,az,far,fat,faz                                           &
     &  ,atemp,adustt,adustd,adens,ane,avr,avt,avz,avtrb,lshade         &
     &  ,ftemp,fdustt,fdustd,fdens,fne,fvr,fvt,fvz,fvtrb,kshade)
!	nearest neighbour in transparent and nontransparent objects     
!	to transform (rotate) the line of sight grid (x,y,z) to the body
!	frozen grid (x'',y'',z'') and interpolate state quantities from 
!	the rotating body frozen coordinates (RBFC) to the line of sight 
!	coordinates (LSC) and transform the velocity field vectors 
!	from RBFC to LSC
!	input: 	
!         theta,alpha,dcut1
!	  temp0,ane0 - out of grid (sphere) reasonable values
!         iii,jjj -identify x,y the ray
!         vxst,vyst,vzst - net velocity of the primary
!	  xcp,ycp,zcp - position of the secondary
!	  vxcp,vycp,vzcp - velocity of the secondary
!	  ar,at,az -coord. in which the model is necessary (LSC)
!	  far,fat,faz -coord. in which the model is defined (RBFC)
!	  ftemp,fdustt,fdustd,fdens,fne,fvr,fvt,fvz,fvtrb,kshade -model
!	output:
!         atemp,adustt,adustd,adens,ane,avr,avt,avz,avtrb,lshade
!         (^rotated model^)
!         vxstr,vystr,vzstr-rotated net velocity of the primary
!	  xcpr,ycpr,zcpr -rotated position of the secondary
!	  vxcpr,vycpr,vzcpr - rotated velocity of the secondary 
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
	real(i8):: ar(ndim1),at(ndim2),az(ndim3)
        real(i8):: far(ndimf1),fat(ndimf2),faz(ndimf3)
        real(i8):: atemp(ndim3),adustt(ndim3),adustd(ndim3)
        real(i8):: adens(ndim3),ane(ndim3)
        real(i8):: avr(ndim3),avt(ndim3)
        real(i8):: avz(ndim3),avtrb(ndim3)
        real(i8):: ftemp(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustt(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustd(ndimf1,ndimf2,ndimf3)
        real(i8):: fdens(ndimf1,ndimf2,ndimf3),fne(ndimf1,ndimf2,ndimf3)
        real(i8):: fvr(ndimf1,ndimf2,ndimf3),fvt(ndimf1,ndimf2,ndimf3)
        real(i8):: fvz(ndimf1,ndimf2,ndimf3),fvtrb(ndimf1,ndimf2,ndimf3)
	integer(i1):: kshade(ndimf1,ndimf2,ndimf3),lshade(ndim3)
	real(i8):: theta,alpha,dcut1,temp0,ane0
	real(i8):: vxst,vyst,vzst,vxstr,vystr,vzstr
	real(i8):: xcp,ycp,zcp,xcpr,ycpr,zcpr
	real(i8):: vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr
	integer:: iii,jjj,ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3
	integer:: ndim1,ndim2,ndim3,nbod1,nbod2,nbod3
	real(i8):: ctheta,stheta,calpha,salpha
	real(i8):: fx,fy,fz,dist1,dx,dy,dz,dist
	integer:: ii,jj,kk,i,j,k,ii1,jj1,kk1
	ctheta=dcos(theta)
	stheta=dsin(theta)
	calpha=dcos(alpha)
	salpha=dsin(alpha)
!	initialization of first values to ii,jj,kk
	ii=1
	jj=1
	kk=1
	i=iii
	j=jjj
	do k=1,nbod3
!	  calculation of RBFC coordinates corresponding to the grid 
!	  points of the line of sight grid 
          call rotz(ar(i),at(j),az(k),fx,fy,fz                          &
     &    ,ctheta,stheta,calpha,salpha)
!	  	interpolation from body frozen grid to required 
!	  line of sight grid  :
          if(fx.gt.far(1).and.fx.lt.far(nbodf1).and.                    &
     &    fy.gt.fat(1).and.fy.lt.fat(nbodf2).and.                       &
     &    fz.gt.faz(1).and.fz.lt.faz(nbodf3))then
	    call hunt(far,nbodf1,fx,ii)
	    if(ii.lt.1.or.ii.ge.nbodf1)then
	      write(*,*)' compiler/maschine dependent error'
	      write(3,*)' 1 compiler/maschine dependent error'
              write(*,'(a,4i4)')' -out of grid at ii,i,j,k',ii,i,j,k
	    endif
	    call hunt(fat,nbodf2,fy,jj)
	    if(jj.lt.1.or.jj.ge.nbodf2)then
	      write(*,*)' compiler/maschine dependent error'
	      write(3,*)' 1 compiler/maschine dependent error'
              write(*,'(a,4i4)')' -out of grid at jj,i,j,k',jj,i,j,k
	    endif
	    call hunt(faz,nbodf3,fz,kk)
	    if(kk.lt.1.or.kk.ge.nbodf3)then
	      write(*,*)' compiler/maschine dependent error'
	      write(3,*)' 1 compiler/maschine dependent error'
              write(*,'(a,4i4)')' -out of grid at kk,i,j,k',kk,i,j,k
	    endif
!		nearest neighbour value approximation
	    dist1=1.d60
	    do ii1=ii,ii+1
            do jj1=jj,jj+1
            do kk1=kk,kk+1	
	      dx=fx-far(ii1)
              dy=fy-fat(jj1)
              dz=fz-faz(kk1)
	      dist=dx*dx+dy*dy+dz*dz
	      if(dist.lt.dist1)then
	        atemp(k)=ftemp(ii1,jj1,kk1)
	        adustt(k)=fdustt(ii1,jj1,kk1)
	        adustd(k)=fdustd(ii1,jj1,kk1)
	        adens(k)=fdens(ii1,jj1,kk1)
	        ane(k)=fne(ii1,jj1,kk1)
	        avr(k)=fvr(ii1,jj1,kk1)
	        avt(k)=fvt(ii1,jj1,kk1)
	        avz(k)=fvz(ii1,jj1,kk1)
	        avtrb(k)=fvtrb(ii1,jj1,kk1)
		lshade(k)=kshade(ii1,jj1,kk1)
		dist1=dist
	      endif
	    enddo
            enddo
            enddo
	  else
!		out of body frozen grid
	    atemp(k)=temp0
	    adustt(k)=temp0
	    adustd(k)=0.d0
	    adens(k)=0.d0
	    ane(k)=ane0
	    avr(k)=0.d0
	    avt(k)=0.d0
	    avz(k)=0.d0
	    avtrb(k)=0.d0
	    lshade(k)=0	
	  endif
	enddo
!	back transform of vector quantities
	do k=1,nbod3
          call rotzb(avr(k),avt(k),avz(k)                               &
     &    ,avr(k),avt(k),avz(k),ctheta,stheta,calpha,salpha)
	enddo
!	back transform of the velocity of the center of mass 
!	of the primary explicitely because of the scattered light 
!	treatment
        call rotzb(vxst,vyst,vzst,vxstr,vystr,vzstr                     &
     &  ,ctheta,stheta,calpha,salpha)
!	back transform of the center of the secondary explicitely 
!	because	of the limb darkening treatement
        call rotzb(xcp,ycp,zcp,xcpr,ycpr,zcpr                           &
     &  ,ctheta,stheta,calpha,salpha)
!       back transform of the velocity of the secondary explicitely 
!	because of the scattered light treatment
        call rotzb(vxcp,vycp,vzcp,vxcpr,vycpr,vzcpr                     &
     &  ,ctheta,stheta,calpha,salpha)
        return
        end
!-----------------------------------------------------------------------
        subroutine rotzb(xpp,ypp,zpp,x,y,z,ctheta,stheta,calpha,salpha) 
!       transformation from body frozen (x'',y'',z'') to line of sight  
!       coordinates (x,y,z) 
!       i.e. rotation about  z'=z'' and then about x=x'
!       input: xpp,ypp,zpp   
!       output:x,y,z
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: xpp,ypp,zpp,x,y,z,ctheta,stheta,calpha,salpha
	real(i8):: xp,yp,zp
        zp=zpp
        xp=xpp*calpha-ypp*salpha
        yp=ypp*calpha+xpp*salpha
        x=xp
        y=yp*ctheta+zp*stheta
        z=zp*ctheta-yp*stheta
        return
        end   
!-----------------------------------------------------------------------
        subroutine rotz(x,y,z,xpp,ypp,zpp,ctheta,stheta,calpha,salpha)  
!       transformation from line of sight (x,y,z) to body frozen 
!       coordinates (x'',y'',z'')
!       i.e. rotation about x=x' and then about z'=z''
!       input: x,y,z   
!       output:xpp,ypp,zpp
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: x,y,z,xpp,ypp,zpp,ctheta,stheta,calpha,salpha
	real(i8):: xp,yp,zp
        xp=x
        yp=y*ctheta-z*stheta
        zp=z*ctheta+y*stheta
        zpp=zp
        xpp=xp*calpha+yp*salpha
        ypp=yp*calpha-xp*salpha
        return
        end   
!-----------------------------------------------------------------------
	subroutine surf(ndim1,ndim2,nbod1,nbod2,ar,at,area)
!	to calculate area(i,j)-area of a projected surface element
!	input: ar,at
!	output: area
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: ar(ndim1),at(ndim2),area(ndim1,ndim2)
	integer:: ndim1,ndim2,nbod1,nbod2
	real(i8):: pi,pomi12,pominn,dar
	integer:: i,j
	pi=3.1415926535897931d0
	do i=2,nbod1-1
	do j=2,nbod2-1
	  area(i,j)=(ar(i+1)-ar(i-1))*(at(j+1)-at(j-1))*0.25d0
	enddo
	enddo
	pomi12=(ar(2)-ar(1))*0.5d0
	pominn=(ar(nbod1)-ar(nbod1-1))*0.5d0
	do j=2,nbod2-1
	  area(1,j)=pomi12*(at(j+1)-at(j-1))*0.5d0
	  area(nbod1,j)=pominn*(at(j+1)-at(j-1))*0.5d0
	enddo
	do i=2,nbod1-1
          dar=ar(i+1)-ar(i-1)
	  area(i,1)=dar*(at(2)-at(1))*0.25d0
	  area(i,nbod2)=dar*(at(nbod2)-at(nbod2-1))*0.25d0
	enddo
	area(1,1)=pomi12*(at(2)-at(1))*0.5d0
        area(1,nbod2)=pomi12*(at(nbod2)-at(nbod2-1))*0.5d0
        area(nbod1,1)=pominn*(at(2)-at(1))*0.5d0
        area(nbod1,nbod2)=pominn*(at(nbod2)-at(nbod2-1))*0.5d0
	return
	end
!-----------------------------------------------------------------------       
        subroutine tlac1(ndim3,nbod3,az,vr,vt,vz                        &
     &  ,temp,ekonc,akonc,dens,dustd,dustt,mshade                       &
     &  ,hi,hneg,h2,hei,hipf)
!	output of selected quantities for selected ray
!	along the light of sight
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
	integer:: ndim3,nbod3
	real(i8):: az(ndim3),hipf(ndim3,2)
	real(i8):: hi(ndim3),hneg(ndim3),h2(ndim3),hei(ndim3)
        real(i8):: vr(ndim3),vt(ndim3),vz(ndim3)
        real(i8):: temp(ndim3),ekonc(ndim3),akonc(ndim3)
        real(i8):: dens(ndim3),dustd(ndim3),dustt(ndim3)
	integer(i1):: mshade(ndim3)
	integer:: kk
        write(2,'(a,a)')' vel.  vx        vy        vz        dens'     &
     &  ,'      dustd     dustt'
        do kk=1,nbod3
          write(2,20)kk,vr(kk),vt(kk),vz(kk),dens(kk)                   &
     &    ,dustd(kk),dustt(kk)
	enddo
20      format(i4,6es10.2)
        write (2,'(a,a)')'       akonc      ekonc      temp'            &
     &        ,'       az      mshade'
        do kk=1,nbod3
          write(2,40)kk,akonc(kk),ekonc(kk),temp(kk),az(kk),mshade(kk)
	enddo
40      format(i4,4es11.3,i3)
        write(2,'(a,a)')'       AKONC      HI         H-         H2'    &
     &  ,'         HEI        HIPF       HIIPF'
        do kk=1,nbod3
          write(2,60)kk,akonc(kk),hi(kk),hneg(kk),h2(kk),hei(kk)        &
     &    ,hipf(kk,1),hipf(kk,2) 
	enddo
60      format(i4,7es11.3)	
!        write(2,'(a)')'        ophbf1     ophff1     ophrs1'
!        do kk=1,nbod3
!	  write(2,'(i4,3es11.3)')kk,ophbf(kk),ophff(kk),ophrs(kk)
!	enddo
	return
	end
!-----------------------------------------------------------------------
        subroutine tlac2(ndim3,nbod3,mion                               &
     &  ,nion,iion,az,zip,eup                                           &
     &  ,temp,stavs,rr,gr,gs,gw)
!	output of selected quantities for selected ray
!	along the light of sight
	implicit none
        integer, parameter:: i4=4,i8=8
	integer:: ndim3,nbod3,mion,nion,iion
	real(i8):: az(ndim3),temp(ndim3)
        real(i8):: stavs(ndim3,mion),RR(ndim3,mion),zip(mion-1)
	real(i8):: eup,gr,gs,gw
	integer:: kk,jj1
	real(i8):: EVCM
        EVCM=8.06554477D3
        WRITE(2,'(A)')' PARTITION FUNCTIONS'
        WRITE(2,'(A,20I10)')' TEMP',(kk,kk=1,NION)
        do kk=1,nbod3
          WRITE(2,240)TEMP(kk),(STAVS(kk,jj1),jj1=1,NION)
        enddo
240     FORMAT(F7.0,20E10.3)
        WRITE(2,*)'   N(J)/N  ION POPULATIONS'
        WRITE(2,'(A,40I11)')' z',(kk,kk=1,NION)
        do kk=1,nbod3
          WRITE(2,260)az(kk),(RR(kk,jj1),jj1=1,NION)
        enddo
260     FORMAT(40(D10.3,1X))
	write(2,'(a,a)')' i   GR        GS        GW        '
        DO kk=1,nbod3
          IF(ZIP(IION).LT.EUP/EVCM)THEN
	    WRITE(2,310)kk,GR,GS,GW
	  ELSE
	    WRITE(2,320)kk,GR,GS,GW
	  ENDIF
        enddo
310     FORMAT(i3,3E10.3,' AUTO')
320     FORMAT(i3,3E10.3)
	return
	end
!-----------------------------------------------------------------------
         subroutine smod1(nbodf1,nbodf1a,nbodf1b                        &
     &   ,nbodf2,nbodf2a,nbodf2b,nbodf3,nbodf3a,nbodf3b                 &
     &   ,rmdfx1,rmdfx2,rmdfx3,rmdfx4,rmdfy1,rmdfy2,rmdfy3,rmdfy4       &
     &   ,rmdfz1,rmdfz2,rmdfz3,rmdfz4                                   &
     &   ,gainfx,gainfy,gainfz,gainfxb,gainfyb,gainfzb                  &
     &  ,rstar,tstar,emstar,xstar,ystar,zstar,vrotst,idifst,drotst,hst  &
     &  ,istar,vxst,vyst,vzst,dgst,ffst,irrst,albst,htst,htsta          &
     &  ,ispst,xspst,yspst,zspst,aspst,tspst                            &
     &  ,icomp,dgcp,ffcp,qq,vrxcp,vrycp,vrzcp,vrotcp,rcp                &
     &  ,xcp,ycp,zcp,vxcp,vycp,vzcp,tempcp,irrcp,albcp,htcp,htcpa       &
     &  ,ienv,emen,qqen,aen,ffen,hen                                    &
     &  ,tempen,densen,aneen,vtrben,dstden,dstten                       &
     &  ,ispot,vrxsp,vrysp,vrzsp,vrotsp,rsp                             &
     &  ,xsp,ysp,zsp,vxsp,vysp,vzsp,tempsp,denssp,anesp,vtrbsp          &
     &  ,dstdsp,dsttsp                                                  &
     &  ,iring,rrg,emrg,b1rg,b2rg,a1rg,a2rg,dr1rg,dr2rg                 &
     &  ,xrg,yrg,zrg,xpolrg,ypolrg,zpolrg,vxrg,vyrg,vzrg                &
     &  ,temprg,densrg,anerg,vtrbrg,itrg                                &
     &  ,edenrg,dstdrg,ede2rg,dst2rg,dsttrg                             &
     &  ,idisc,adisc,rindc,routdc,emdc,rdc                              &
     &  ,xdc,ydc,zdc,xdisc,ydisc,zdisc,vxdc,vydc,vzdc                   &
     &  ,tempdc,densdc,anedc,vtrbdc,edendc,itdc,etmpdc,dstddc,dsttdc    &
     &  ,inebl,aneb,rinnb,routnb,emnb,rnb,hcnb,ishdnb,hshdnb            &
     &  ,hinvnb,tinvnb,iinvnb,hwindnb,ndennb,denxnb,denznb              &
     &  ,ivelnb,hvelnb,vnb,evelnb,xneb,yneb,zneb,vxnb,vynb,vznb         &
     &  ,tempnb,densnb,anenb,vtrbnb,edennb,itnb,etmpnb,dstdnb,dsttnb    &
     &  ,x1sm,y1sm,z1sm,x2sm,y2sm,z2sm,v1sm,v2sm,r1sm,r2sm              &
     &  ,ism,vxsm,vysm,vzsm,xsm,ysm,zsm,psm                             &
     &  ,tempsm,denssm,anesm,vtrbsm,edensm,dstdsm,dsttsm                &
     &  ,x1fw,y1fw,z1fw,x2fw,y2fw,z2fw,v1fw,v2fw,r1fw,r2fw              &
     &  ,iflow,vxfw,vyfw,vzfw,xfw,yfw,zfw,pfw                           &
     &  ,tempfw,densfw,anefw,vtrbfw,edenfw,dstdfw,dsttfw                &
     &  ,iufo,aufo,rinuf,routuf,emuf,ruf                                &
     &  ,xuf,yuf,zuf,xufo,yufo,zufo,vxuf,vyuf,vzuf                      &
     &  ,tempuf,densuf,aneuf,vtrbuf,edenuf,ituf,etmpuf,dstduf,dsttuf    &
     &  ,ajet,rinjt,routjt,xjt,yjt,zjt,xjet,yjet,zjet                   &
     &  ,ivjt,vjt,eveljt,rcjt,vxjt,vyjt,vzjt                            &
     &  ,tempjt,densjt,anejt,vtrbjt,ijet,dstdjt,dsttjt,etmpjt,asymjt    &
     &  ,rinsh,routsh,vsh,vxsh,vysh,vzsh                                &
     &  ,tempsh,denssh,anesh,vtrbsh,ishell,evelsh,rcsh,dstdsh,dsttsh    &
     &  ,etmpsh                                                         &
     &  ,v0,temp0,dens0,ane0,dcut1,dcut2,dcut3,dcutn                    &
     &  ,far,fat,faz,ftemp,fdens,fne,fvr,fvt,fvz,fvtrb,fdustd,fdustt    &
     &  ,kshade)
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
        include 'param.inc'
	integer:: nbodf1,nbodf1a,nbodf1b
	integer:: nbodf2,nbodf2a,nbodf2b
	integer:: nbodf3,nbodf3a,nbodf3b
	integer:: idifst,istar,irrst,ispst,icomp,irrcp,ienv,ispot,iring
	integer:: itrg,idisc,itdc,inebl,ndennb,iinvnb,itnb,ivelnb,ishdnb
	integer:: ism,iflow,iufo,ituf,ijet,ivjt,ishell
	real(i8):: far(ndimf1),fat(ndimf2),faz(ndimf3)
        real(i8):: ftemp(ndimf1,ndimf2,ndimf3)
        real(i8):: fdens(ndimf1,ndimf2,ndimf3),fne(ndimf1,ndimf2,ndimf3)
        real(i8):: fvr(ndimf1,ndimf2,ndimf3),fvt(ndimf1,ndimf2,ndimf3)
        real(i8):: fvz(ndimf1,ndimf2,ndimf3),fvtrb(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustd(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustt(ndimf1,ndimf2,ndimf3)
        integer(i1):: kshade(ndimf1,ndimf2,ndimf3)
        real(i8):: rosup(ndimf1,ndimf2),rosus(ndimf1,ndimf2)
        real(i8):: rotemp(ndimf1,ndimf2),rotems(ndimf1,ndimf2)
	real(i8):: rosue(ndimf1,ndimf2),roteme(ndimf1,ndimf2)
        real(i8):: trmt(3,3),trmtnb(3,3),trmtuf(3,3)
        real(i8):: trmtst(3,3),trmtrg(3,3)
        real(i8):: denxnb(mstarx),denznb(mstarx)
	real(i8):: rmdfx1,rmdfx2,rmdfx3,rmdfx4
	real(i8):: rmdfy1,rmdfy2,rmdfy3,rmdfy4
        real(i8):: rmdfz1,rmdfz2,rmdfz3,rmdfz4
        real(i8):: gainfx,gainfy,gainfz,gainfxb,gainfyb,gainfzb
        real(i8):: rstar,tstar,emstar,xstar,ystar,zstar,vrotst,drotst
        real(i8):: hst,vxst,vyst,vzst,dgst,ffst,albst,htst,htsta
        real(i8):: xspst,yspst,zspst,aspst,tspst
        real(i8):: dgcp,ffcp,qq,vrxcp,vrycp,vrzcp,vrotcp,rcp
        real(i8):: xcp,ycp,zcp,vxcp,vycp,vzcp,tempcp,albcp,htcp,htcpa
        real(i8):: emen,qqen,aen,ffen,hen
	real(i8):: tempen,densen,aneen,vtrben,dstden,dstten
        real(i8):: vrxsp,vrysp,vrzsp,vrotsp,rsp
        real(i8):: xsp,ysp,zsp,vxsp,vysp,vzsp,tempsp,denssp,anesp,vtrbsp
        real(i8):: dstdsp,dsttsp
	real(i8):: rrg,emrg,b1rg,b2rg,a1rg,a2rg,dr1rg,dr2rg
        real(i8):: xrg,yrg,zrg,xpolrg,ypolrg,zpolrg,vxrg,vyrg,vzrg
        real(i8):: temprg,densrg,anerg,vtrbrg
        real(i8):: edenrg,dstdrg,ede2rg,dst2rg,dsttrg
        real(i8):: adisc,rindc,routdc,emdc,rdc
        real(i8):: xdc,ydc,zdc,xdisc,ydisc,zdisc,vxdc,vydc,vzdc
        real(i8)::tempdc,densdc,anedc,vtrbdc,edendc,etmpdc,dstddc,dsttdc
        real(i8):: aneb,rinnb,routnb,emnb,rnb,hcnb,hshdnb
        real(i8):: hinvnb,tinvnb,hwindnb,hvelnb,vnb,evelnb
	real(i8):: xneb,yneb,zneb,vxnb,vynb,vznb
        real(i8)::tempnb,densnb,anenb,vtrbnb,edennb,etmpnb,dstdnb,dsttnb
        real(i8):: x1sm,y1sm,z1sm,x2sm,y2sm,z2sm,v1sm,v2sm,r1sm,r2sm
        real(i8):: vxsm,vysm,vzsm,xsm,ysm,zsm,psm
	real(i8):: tempsm,denssm,anesm,vtrbsm,edensm,dstdsm,dsttsm
        real(i8):: x1fw,y1fw,z1fw,x2fw,y2fw,z2fw,v1fw,v2fw,r1fw,r2fw
        real(i8):: vxfw,vyfw,vzfw,xfw,yfw,zfw,pfw
        real(i8):: tempfw,densfw,anefw,vtrbfw,edenfw,dstdfw,dsttfw
        real(i8):: aufo,rinuf,routuf,emuf,ruf
        real(i8):: xuf,yuf,zuf,xufo,yufo,zufo,vxuf,vyuf,vzuf
        real(i8)::tempuf,densuf,aneuf,vtrbuf,edenuf,etmpuf,dstduf,dsttuf
        real(i8):: ajet,rinjt,routjt
        real(i8):: xjt,yjt,zjt,xjet,yjet,zjet
        real(i8):: vjt,eveljt,rcjt,vxjt,vyjt,vzjt
        real(i8):: tempjt,densjt,anejt,vtrbjt,dstdjt,dsttjt,etmpjt,asymjt
	real(i8):: rinsh,routsh,vsh,vxsh,vysh,vzsh
        real(i8):: tempsh,denssh,anesh,vtrbsh,evelsh,rcsh,dstdsh,dsttsh
        real(i8):: etmpsh
        real(i8):: v0,temp0,dens0,ane0,dcut1,dcut2,dcut3,dcutn
        real(i8):: clight,pi,grav,stefb,rsol,dunt1,dunt2,a1,a2
        real(i8):: htstb,porbst,omgst,omgcp
        real(i8):: unstar,unspst,omgen,porb
	real(i8):: vxen,vyen,vzen,xen,yen,zen,unsm,dxsm,dysm,dzsm,dsm
        real(i8):: omgsm,unrg,vrg,undisc,unneb,unfw,dxfw,dyfw,dzfw,dfw
        real(i8):: omgfw,unufo,unjet,er2,er,farcp,fatcp,fazcp,ermcp
        real(i8):: unspot,farsp,fatsp,fazsp,ermsp,omega
	real(i8):: arease,areasp,areass
        real(i8):: h00nb,h01nb,si00nb,si01nb,ss01nb
	real(i8):: rfrone,rbacke,rpolee,rsidee,rmeane
	real(i8):: rfronp,rbackp,rpolep,rsidep,rmeanp
	real(i8):: rfrons,rbacks,rpoles,rsides,rmeans
	real(i8):: teffe,teffp,teffs
	integer:: i,j,k,inx
!	Calculates model of a shell (behaviour of state quantities and
!	velocity field in Cartesian coordinates x,y,z=far,fat,faz. 
!	fvr,fvt,fvz are corresponding components of velocity vector.
!	ftemp,fdens,fne are temperature, density
!	and electron number density of gas, respectively.
!	fvtrb-microturbulence, 
!       fdustd - dust density, fdustt - dust temperature.
!       kshade - shadows (i.e. scattering from the two stars)
!       kshade=1 scattering only from the central star
!       kshade=2 scattering only from the companion star
!       kshade=3 scattering from both stars (default)
!       kshade=0 no scattering
!	Reasonable vaules of temp and ekonc must be defined 
!	everywhere but you can use dens<=0.d0 to define empty space
!	or dcut1<dens<dcut2 for central nontransparent objects (star),
!	dcut2<dens<dcut3 for second nontransparent object (companion),
!       dcut3<dens<dcutn for 3. nontransparent object (3.body),
!	dcutn<dens for any dark nontransparent object.
	clight=2.99792458d10
	pi=3.1415926535897931d0
	grav=6.67408d-8
        stefb=5.670400d-5
!        rsol=6.9599d10
!        rsol=6.95508d10
	rsol=6.957d10
!	nbodf1,nbodf2,nbodf3>=1
	dunt1=(dcut1+dcut2)*0.5d0
	dunt2=(dcut2+dcut3)*0.5d0
!		definition of the body frozen grid
        call grid3d(ndimf1,ndimf2,ndimf3                                &
     &  ,nbodf1,nbodf1a,nbodf1b                                         &
     &  ,nbodf2,nbodf2a,nbodf2b                                         &
     &  ,nbodf3,nbodf3a,nbodf3b                                         &          
     &  ,rmdfx1,rmdfx2,rmdfx3,rmdfx4                                    &
     &  ,rmdfy1,rmdfy2,rmdfy3,rmdfy4                                    &     
     &  ,rmdfz1,rmdfz2,rmdfz3,rmdfz4                                    &          
     &  ,gainfx,gainfxb,gainfy,gainfyb,gainfz,gainfzb,far,fat,faz)
!
        if(istar.gt.1.or.icomp.gt.1.or.                                 &
     &	(istar.gt.0.and.icomp.gt.0.and.vxst.gt.clight))then
	  a1=xcp*qq/(1.d0+qq)
	  a2=xcp-a1
	  porbst=2.d0*pi*dsqrt(xcp**3/(grav*emstar*(1.d0+qq)))
!	  a1 is local auxiliary variable	  
!         omgst for istar=1 is temporary and will be overwritten
!	  in the subroutine star
	  omgst=2.d0*pi/porbst
	  omgcp=omgst
	  vxst=0.d0
	  vyst=-omgst*a1
	  vzst=0.d0
   	  vxcp=0.d0
   	  vycp=omgst*a2
   	  vzcp=0.d0
	  ycp=0.d0
	  zcp=0.d0
	else
	  porbst=0.d0  
	endif
	if(istar.eq.1)then
	  unstar=xstar*xstar+ystar*ystar+zstar*zstar
	  unstar=dsqrt(unstar)
	  xstar=xstar/unstar
	  ystar=ystar/unstar
	  zstar=zstar/unstar
	  unspst=xspst*xspst+yspst*yspst+zspst*zspst
	  unspst=dsqrt(unspst)
	  xspst=xspst/unspst
	  yspst=yspst/unspst
	  zspst=zspst/unspst
	  aspst=dcos(aspst)
!	  omgst=vrotst/rstar
!		calculates rotation matrix trmtst
	  call trans(xstar,ystar,zstar,trmtst)
	elseif(istar.eq.2)then
	  xstar=0.d0
	  ystar=0.d0
	  zstar=1.d0
!		rosup>=0 -Roche surface of the primary
          call roche(rosup,rotemp,areasp,teffp                          &
     &    ,rfronp,rbackp,rpolep,rsidep,rmeanp                           &
     &    ,tstar,rstar,dgst,irrst,albst,htst,htsta                      &
     &    ,tempcp,rcp,dgcp,irrcp,albcp,htcp,htcpa                       &
     &    ,far,fat,xcp,qq,ffst,nbodf1,nbodf2,1)
	elseif(istar.eq.3)then
	  xstar=0.d0
	  ystar=0.d0
	  zstar=1.d0
	  icomp=0
!		rosup>=0 -Roche surface of the contact system
          call roche(rosup,rotemp,areasp,teffp                          &
     &    ,rfronp,rbackp,rpolep,rsidep,rmeanp                           &
     &    ,tstar,rstar,dgst,irrst,albst,htst,htsta                      &
     &    ,tempcp,rcp,dgcp,irrcp,albcp,htcp,htcpa                       &
     &    ,far,fat,xcp,qq,ffst,nbodf1,nbodf2,3)     
!          do 34 i=1,nbodf1
!          do 32 j=1,nbodf2
!            write(19,*)far(i),fat(j),rosup(i,j)
!32	   continue
!	   write(19,*)
!34	   continue
	endif
	if(icomp.eq.2)then
	  vrxcp=0.d0
	  vrycp=0.d0
	  vrzcp=1.d0
!  	        rosus>=0 -Roche surface of the secondary
          call roche(rosus,rotems,areass,teffs                          &
     &    ,rfrons,rbacks,rpoles,rsides,rmeans                           &
     &    ,tstar,rstar,dgst,irrst,albst,htst,htsta                      &
     &    ,tempcp,rcp,dgcp,irrcp,albcp,htcp,htcpa                       &
     &    ,far,fat,xcp,qq,ffcp,nbodf1,nbodf2,2)
	endif
	if(ienv.eq.2.or.ienv.eq.3)then
	  a1=aen*qqen/(1.d0+qqen)
	  porb=2.d0*pi*dsqrt(aen**3/(grav*emen*(1.d0+qqen)))
	  omgen=2.d0*pi/porb
	  vxen=0.d0
	  vyen=-omgen*a1
	  vzen=0.d0
	  xen=0.d0
	  yen=0.d0
	  zen=1.d0
	endif
	if(ienv.eq.2)then
!		rosue>=0 -Roche surface of the envelope
          call roche(rosue,roteme,arease,teffe                          &
     &    ,rfrone,rbacke,rpolee,rsidee,rmeane                           &
     &    ,tempen,0.d0,0.d0,0,0.d0,1.d0,1.d0                            &
     &    ,tempen,0.d0,0.d0,0,0.d0,1.d0,1.d0                            &
     &    ,far,fat,aen,qqen,ffen,nbodf1,nbodf2,1)
	elseif(ienv.eq.3)then
!	  icomp=0
!		rosue>=0 -Roche surface of the common envelope
          call roche(rosue,roteme,arease,teffe                          &
     &    ,rfrone,rbacke,rpolee,rsidee,rmeane                           &
     &    ,tempen,0.d0,0.d0,0,0.d0,1.d0,1.d0                            &
     &    ,tempen,0.d0,0.d0,0,0.d0,1.d0,1.d0                            &
     &    ,far,fat,aen,qqen,ffen,nbodf1,nbodf2,3)     
     	endif
	if(ism.eq.1)then
	  unsm=xsm*xsm+ysm*ysm+zsm*zsm
	  unsm=dsqrt(unsm)
	  xsm=xsm/unsm
	  ysm=ysm/unsm
	  zsm=zsm/unsm
	  dxsm=x2sm-x1sm
	  dysm=y2sm-y1sm
	  dzsm=z2sm-z1sm	
	  dsm=dxsm*dxsm+dysm*dysm+dzsm*dzsm
	  dsm=dsqrt(dsm)
	  dxsm=dxsm/dsm
	  dysm=dysm/dsm
	  dzsm=dzsm/dsm	
          omgsm=2.d0*pi/psm
	endif
	if(iring.gt.0)then
	  unrg=xpolrg*xpolrg+ypolrg*ypolrg+zpolrg*zpolrg
	  unrg=dsqrt(unrg)
	  xpolrg=xpolrg/unrg
	  ypolrg=ypolrg/unrg
	  zpolrg=zpolrg/unrg
!		calculates rotation matrix trmtrg
	  call trans(xpolrg,ypolrg,zpolrg,trmtrg)
          vrg=dsqrt(grav*emrg/rrg)/rrg
	endif		
	if(idisc.gt.0)then
	  undisc=xdisc*xdisc+ydisc*ydisc+zdisc*zdisc
	  undisc=dsqrt(undisc)
	  xdisc=xdisc/undisc
	  ydisc=ydisc/undisc
	  zdisc=zdisc/undisc
!		calculates rotation matrix trmt
	  call trans(xdisc,ydisc,zdisc,trmt)
	endif
	if(inebl.eq.1)then
	  unneb=xneb*xneb+yneb*yneb+zneb*zneb
	  unneb=dsqrt(unneb)
	  xneb=xneb/unneb
	  yneb=yneb/unneb
	  zneb=zneb/unneb
!		calculates rotation matrix trmtnb
	  call trans(xneb,yneb,zneb,trmtnb)
	endif
	if(iflow.eq.1)then
	  unfw=xfw*xfw+yfw*yfw+zfw*zfw
	  unfw=dsqrt(unfw)
	  xfw=xfw/unfw
	  yfw=yfw/unfw
	  zfw=zfw/unfw
	  dxfw=x2fw-x1fw
	  dyfw=y2fw-y1fw
	  dzfw=z2fw-z1fw	
	  dfw=dxfw*dxfw+dyfw*dyfw+dzfw*dzfw
	  dfw=dsqrt(dfw)
	  dxfw=dxfw/dfw
	  dyfw=dyfw/dfw
	  dzfw=dzfw/dfw	
          omgfw=2.d0*pi/pfw
	endif
	if(iufo.gt.0)then
	  unufo=xufo*xufo+yufo*yufo+zufo*zufo
	  unufo=dsqrt(unufo)
	  xufo=xufo/unufo
	  yufo=yufo/unufo
	  zufo=zufo/unufo
!		calculates rotation matrix trmtuf
	  call trans(xufo,yufo,zufo,trmtuf)
	endif
        if(ijet.eq.1.or.ijet.eq.2)then	
	  unjet=xjet*xjet+yjet*yjet+zjet*zjet
	  unjet=dsqrt(unjet)
	  xjet=xjet/unjet
	  yjet=yjet/unjet
	  zjet=zjet/unjet  
	endif
	write(2,*)' Porbst[d]  istar  icomp'
	write(2,*)porbst/8.64d4,istar,icomp
        write(2,'(a,a)')'     xstar      ystar      zstar      '        &
     &  ,'vxst       vyst       vzst'
        write(2,'(6e11.3)')xstar,ystar,zstar,vxst,vyst,vzst
        write(2,'(a,a)')'     vrxcp      vrycp      vrzcp      '        &
     &  ,'vxcp       vycp       vzcp'
        write(2,'(6e11.3)')vrxcp,vrycp,vrzcp,vxcp,vycp,vzcp
        write(2,'(a,a)')'       xcp        ycp        zcp      '
        write(2,'(3e11.3)')xcp,ycp,zcp
        if(istar.gt.1)then
          write(2,'(a,a)')' Roche S.  Teffprim '
          write(2,'(2e12.4)')areasp,teffp
          write(2,'(a,a)')'  Rsub     Rback    Rpole    Rside    Rmean' &
     &    ,'    Rsub/Rpole in units of a'
          write(2,'(6f9.5)')rfronp,rbackp,rpolep,rsidep,rmeanp,         &
     &    rfronp/rpolep
        endif  
        if(icomp.eq.2)then
          write(2,'(a,a)')' Roche S. Teffsec'
          write(2,'(2e12.4)')areass,teffs        
          write(2,'(a,a)')'  Rsub     Rback    Rpole    Rside    Rmean' &
     &    ,'    Rsub/Rpole in units of a'
          write(2,'(6f9.5)')rfrons,rbacks,rpoles,rsides,rmeans,         &
     &    rfrons/rpoles
        endif
!        WRITE(2,'(A,A)')'  i  j  k      temp     dustd      dens'      &
!     &  ,'       ane       avr       avt       avz     avtrb'
!	The model with a star, companion, and other objects.
!	The order of objects defines their priority. 
	do 60 i=1,nbodf1
	  do 50 j=1,nbodf2
	    do 40 k=1,nbodf3
	      inx=0
	      er2=far(i)*far(i)+fat(j)*fat(j)+faz(k)*faz(k)
	      er=dsqrt(er2)
!	      default scattering from both stars:
	      kshade(i,j,k)=3	
!	      		inside the star
              call star(istar,rstar,tstar,xstar,ystar,zstar             &
     &        ,vxst,vyst,vzst,rosup(i,j),rotemp(i,j),trmtst             &
     &        ,irrst,htst,htsta,htstb,albst                             &
     &        ,icomp,xcp,ycp,zcp,rcp,tempcp                             &
     &        ,idifst,vrotst,drotst,hst                                 &
     &        ,ispst,xspst,yspst,zspst,aspst,tspst                      &
     &        ,dunt1,ane0,er,er2,porbst                                 &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
              if(inx.eq.1)goto 40
!	      		inside the companion
	      if(icomp.eq.1)then
	        farcp=far(i)-xcp		
	        fatcp=fat(j)-ycp
	        fazcp=faz(k)-zcp
	        ermcp=farcp*farcp+fatcp*fatcp+fazcp*fazcp
	        ermcp=dsqrt(ermcp)	
                if(ermcp.lt.rcp)then
	 	  unspot=vrxcp*vrxcp+vrycp*vrycp+vrzcp*vrzcp
		  unspot=dsqrt(unspot)
		  omgcp=vrotcp/rcp/unspot
	          ftemp(i,j,k)=tempcp
	          fdustt(i,j,k)=tempcp
	          fdustd(i,j,k)=0.d0
	          fdens(i,j,k)=dunt2
	          fne(i,j,k)=ane0
	          fvr(i,j,k)=omgcp*(vrycp*fazcp-vrzcp*fatcp)+vxcp
	          fvt(i,j,k)=omgcp*(vrzcp*farcp-vrxcp*fazcp)+vycp
	          fvz(i,j,k)=omgcp*(vrxcp*fatcp-vrycp*farcp)+vzcp
	          fvtrb(i,j,k)=0.d0
	          goto 40
	        endif	      
	      elseif(icomp.eq.2)then
   	        if(dabs(faz(k)).lt.rosus(i,j))then
   	          farcp=far(i)-xcp                
   	          fatcp=fat(j)-ycp
   	          fazcp=faz(k)-zcp
   	          ftemp(i,j,k)=rotems(i,j)
	          fdustt(i,j,k)=rotems(i,j)
	          fdustd(i,j,k)=0.d0
	          fdens(i,j,k)=dunt2
	          fne(i,j,k)=ane0
	          fvr(i,j,k)=omgcp*(vrycp*fazcp-vrzcp*fatcp)+vxcp
	          fvt(i,j,k)=omgcp*(vrzcp*farcp-vrxcp*fazcp)+vycp
	          fvz(i,j,k)=omgcp*(vrxcp*fatcp-vrycp*farcp)+vzcp
	          fvtrb(i,j,k)=0.d0
	          goto 40
		endif
	      endif	      
!	      		inside the spot
	      call spot(ispot,vrxsp,vrysp,vrzsp,vrotsp,rsp              &
     &        ,xsp,ysp,zsp,vxsp,vysp,vzsp                               &
     &        ,tempsp,denssp,anesp,vtrbsp,dstdsp,dsttsp                 &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
     	      if(inx.eq.1)goto 40
!	      		in the stream
              call stream(ism,dxsm,dysm,dzsm,x1sm,y1sm,z1sm             &
     &        ,r1sm,r2sm,v1sm,v2sm,vxsm,vysm,vzsm                       &
     &        ,omgsm,xsm,ysm,zsm,dsm                                    &
     &        ,vtrbsm,tempsm,dstdsm,dsttsm,denssm,anesm,edensm          &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
	      if(inx.eq.1)goto 40
!	      		in the ring
              call ring(iring,xrg,yrg,zrg                               &
     &        ,a1rg,a2rg,b1rg,b2rg,dr1rg,dr2rg,rrg,vrg,itrg             &
     &        ,vxrg,vyrg,vzrg,xpolrg,ypolrg,zpolrg,trmtrg,dcut1         &
     &        ,vtrbrg,temprg,densrg,anerg                               &
     &        ,dstdrg,dst2rg,edenrg,ede2rg,dsttrg                       &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
     	      if(inx.eq.1)goto 40
!	      		in the disc
              call ufo(idisc,xdisc,ydisc,zdisc                          &
     &        ,adisc,rindc,routdc,emdc,rdc,xdc,ydc,zdc                  &
     &        ,vxdc,vydc,vzdc,trmt                                      &
     &        ,tempdc,densdc,anedc,vtrbdc,edendc,itdc,etmpdc            &
     &        ,dstddc,dsttdc                                            &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
     	      if(inx.eq.1)goto 40
!	      		in the jet 
              call sjet(ijet,ajet,rinjt,routjt                          &
     &        ,xjt,yjt,zjt,xjet,yjet,zjet                               &
     &        ,ivjt,vjt,eveljt,rcjt,vxjt,vyjt,vzjt,vtrbjt               &
     &        ,tempjt,densjt,anejt,dstdjt,dsttjt,etmpjt,asymjt          &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
              if(inx.eq.1)goto 40
!	      		in the nebula
              call nebula(inebl,xneb,yneb,zneb                          &
     &        ,aneb,rinnb,routnb,emnb,rnb,hcnb,ishdnb,hshdnb            &
     &        ,hinvnb,tinvnb,iinvnb,hwindnb,ndennb,denxnb,denznb        &
     &        ,vxnb,vynb,vznb,trmtnb,ivelnb,hvelnb,vnb,evelnb           &
     &        ,tempnb,densnb,anenb,vtrbnb,edennb,itnb,etmpnb            &
     &        ,dstdnb,dsttnb                                            &
     &        ,far(i),fat(j),faz(k),er,inx                              &
     &        ,h00nb,h01nb,si00nb,si01nb,ss01nb                         &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k),kshade(i,j,k)                &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
     	      if(inx.eq.1)goto 40
!	      		in the flow
              call stream(iflow,dxfw,dyfw,dzfw,x1fw,y1fw,z1fw           &
     &        ,r1fw,r2fw,v1fw,v2fw,vxfw,vyfw,vzfw                       &
     &        ,omgfw,xfw,yfw,zfw,dfw                                    &
     &        ,vtrbfw,tempfw,dstdfw,dsttfw,densfw,anefw,edenfw          &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
	      if(inx.eq.1)goto 40
!			inside the envelope
              call envelope(ienv,rosue(i,j),hen,xen,yen,zen             &
     &        ,vxen,vyen,vzen,omgen,inx                                 &
     &        ,tempen,dstden,dstten,densen,aneen,vtrben                 &
     &        ,far(i),fat(j),faz(k)                                     &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
	      if(inx.eq.1)goto 40
!	      		in the ufo
              call ufo(iufo,xufo,yufo,zufo                              &
     &        ,aufo,rinuf,routuf,emuf,ruf,xuf,yuf,zuf                   &
     &        ,vxuf,vyuf,vzuf,trmtuf                                    &
     &        ,tempuf,densuf,aneuf,vtrbuf,edenuf,ituf,etmpuf            &
     &        ,dstduf,dsttuf                                            &
     &        ,far(i),fat(j),faz(k),inx                                 &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
     	      if(inx.eq.1)goto 40
!	      		in the shell
              call shell(ishell,rinsh,routsh,rcsh                       &
     &        ,vsh,vxsh,vysh,vzsh,vtrbsh                                &
     &        ,tempsh,dstdsh,dsttsh,denssh,anesh,evelsh,etmpsh          &
     &        ,far(i),fat(j),faz(k),er,inx                              &
     &        ,ftemp(i,j,k),fdens(i,j,k),fne(i,j,k)                     &
     &        ,fdustd(i,j,k),fdustt(i,j,k)                              &
     &        ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k))
              if(inx.eq.1)goto 40
!	      		in the background	
	      ftemp(i,j,k)=temp0
	      fdens(i,j,k)=dens0
	      fne(i,j,k)=ane0
	      fdustd(i,j,k)=0.d0
	      fdustt(i,j,k)=temp0
	      fvr(i,j,k)=v0*far(i)/er
	      fvt(i,j,k)=v0*fat(j)/er
	      fvz(i,j,k)=v0*faz(k)/er
	      fvtrb(i,j,k)=0.d0	
!	      write(2,70)i,j,k,ftemp(i,j,k),fdustd(i,j,k),fdens(i,j,k)  &
!     &        ,fne(i,j,k),fvr(i,j,k),fvt(i,j,k),fvz(i,j,k)             &
!     &        ,fvtrb(i,j,k)
40	    continue
50 	  continue
60	continue
70 	FORMAT(3i3,8E10.3)
!	to write the model into the file readable by this code (smod2)
!	just uncomment the lines below
!        open(10,file='shellspec.mod',status='unknown')
!        write(10,'(1x,3i5)')nbodf1,nbodf2,nbodf3
!        write(10,'(5e12.4)')(far(i),i=1,nbodf1)
!        write(10,'(5e12.4)')(fat(i),i=1,nbodf2)
!        write(10,'(5e12.4)')(faz(i),i=1,nbodf3)
!        do i=1,nbodf1
!        do j=1,nbodf2
!        do k=1,nbodf3
!          write(10,'(9es11.3,i3)')ftemp(i,j,k)                         &
!     &    ,fdustt(i,j,k),fdustd(i,j,k),fdens(i,j,k),fne(i,j,k)         &
!     &    ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k),kshade(i,j,k)
!        enddo
!        enddo
!        enddo
!        close(10)	
	if(inebl.eq.1)then
	  write(2,110)' inebl=1= pp-disk'
	  write(2,110)' H_scale_in H_scale_out [Rsol] sound_out [km/s]'	  
          write(2,'(4es10.2)')h00nb/rsol,h01nb/rsol,ss01nb/1.d5
	  write(2,110)' Surfdens_in Surfdens_out [g/cm^2]'	  
          write(2,'(4es11.2)')si00nb,si01nb
	endif
110	format(a)	
	return
	end
!-----------------------------------------------------------------------
        subroutine smod2(ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3      &
     &  ,far,fat,faz,ftemp,fdens,fne,fvr,fvt,fvz,fvtrb,fdustd,fdustt    &
     &  ,kshade)
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
	integer:: ndimf1,ndimf2,ndimf3,nbodf1,nbodf2,nbodf3
	real(i8):: far(ndimf1),fat(ndimf2),faz(ndimf3)
        real(i8):: fdustd(ndimf1,ndimf2,ndimf3)
        real(i8):: fdustt(ndimf1,ndimf2,ndimf3)
        real(i8):: ftemp(ndimf1,ndimf2,ndimf3)
        real(i8):: fdens(ndimf1,ndimf2,ndimf3),fne(ndimf1,ndimf2,ndimf3)
        real(i8):: fvr(ndimf1,ndimf2,ndimf3),fvt(ndimf1,ndimf2,ndimf3)
        real(i8):: fvz(ndimf1,ndimf2,ndimf3),fvtrb(ndimf1,ndimf2,ndimf3)
        integer(i1):: kshade(ndimf1,ndimf2,ndimf3)
	integer:: i,j,k
!	reads the input model of the shell from the file
!	nbodf1,nbodf2,nbodf3 -number of x,y,z grid points
!	far,fat,faz -define the x,y,z grid points [cm]
!	ftemp, fdens - gas temperature [K] and density [g/cm^3] 
!	fne -electron number density [cm^-3]
!       fdustd, fdustt -dust density [g/cm^3] and temperature [K]
!	fvr,fvt,fvz,fvtrb -x,y,z components of the velocity field [cm/s]
!	fvtrb -turbulence [cm/s]
!       kshade - shadows 
!       kshade=1 scattering only from the central star
!       kshade=2 scattering only from the companion star
!       kshade=3 scattering from both stars (default)
!       kshade=0 no scattering
	read(10,*)nbodf1,nbodf2,nbodf3
	if(nbodf1.gt.ndimf1.or.nbodf2.gt.ndimf2.or.nbodf3.gt.ndimf3)then
	  write(*,*)' error: space dimension exceeded, stop'
	  write(3,*)' 1 error: space dimension exceeded, stop'
	  goto 100
	endif
	read(10,*)(far(i),i=1,nbodf1)
	read(10,*)(fat(i),i=1,nbodf2)
	read(10,*)(faz(i),i=1,nbodf3)
	do i=1,nbodf1
	do j=1,nbodf2
	do k=1,nbodf3
          read(10,*)ftemp(i,j,k),fdustt(i,j,k),fdustd(i,j,k)            &
     &    ,fdens(i,j,k),fne(i,j,k)                                      &
     &    ,fvr(i,j,k),fvt(i,j,k),fvz(i,j,k),fvtrb(i,j,k),kshade(i,j,k)
        enddo
        enddo
        enddo
100	return
	end
!-----------------------------------------------------------------------	
              subroutine star(istar,rstar,tstar,xstar,ystar,zstar       &
     &        ,vxst,vyst,vzst,rosup,rotemp,trmtst                       &
     &        ,irrst,htst,htsta,htstb,albst                             &
     &        ,icomp,xcp,ycp,zcp,rcp,tempcp                             &
     &        ,idifst,vrotst,drotst,hst                                 &
     &        ,ispst,xspst,yspst,zspst,aspst,tspst                      &
     &        ,dunt1,ane0,er,er2,porbst                                 &
     &        ,far,fat,faz,inst                                         &
     &        ,ftemp,fdens,fne                                          &
     &        ,fdustd,fdustt                                            &
     &        ,fvr,fvt,fvz,fvtrb)
!	central star
        implicit none
        integer, parameter:: i4=4,i8=8
        integer:: istar,irrst,icomp,idifst,ispst,inst
	real(i8):: trmtst(3,3)
	real(i8):: rstar,tstar,xstar,ystar,zstar
	real(i8):: vxst,vyst,vzst,rosup,rotemp
	real(i8):: htst,htsta,htstb,albst
	real(i8):: xcp,ycp,zcp,rcp,tempcp
	real(i8):: vrotst,drotst,hst
	real(i8):: xspst,yspst,zspst,aspst,tspst
	real(i8):: dunt1,ane0,er,er2,porbst
	real(i8):: far,fat,faz
	real(i8):: ftemp,fdens,fne,fdustd,fdustt,fvr,fvt,fvz,fvtrb
	real(i8):: pi,ffx,ffy,ffz,omgste,omgstp,sinlat
	real(i8):: omgst,farcp,fatcp,fazcp,cosir,coslat
	real(i8):: fnlat,pom,t04,firr,tirr4,fxspst,fyspst,fzspst
	real(i8):: cspst
              pi=3.1415926535897931d0
              if(istar.eq.1)then
                if(er.lt.rstar)then
!			new coordinates aligned with the rot.axis		
                  ffx=trmtst(1,1)*far+trmtst(1,2)*fat                   &
     &            +trmtst(1,3)*faz
                  ffy=trmtst(2,1)*far+trmtst(2,2)*fat                   &
     &            +trmtst(2,3)*faz
                  ffz=trmtst(3,1)*far+trmtst(3,2)*fat                   &
     &            +trmtst(3,3)*faz
                  if(idifst.eq.1)then
	            omgste=vrotst/rstar
		    omgstp=drotst*omgste
		    sinlat=ffz**2/(ffx**2+ffy**2+ffz**2)
		    omgst=omgste-(omgste-omgstp)*sinlat
		  elseif(idifst.eq.2)then
		    omgste=vrotst/rstar
		    omgstp=drotst*omgste
		    sinlat=ffz**2/(ffx**2+ffy**2+ffz**2)
		    if(sinlat.lt.hst**2)then
		      omgst=omgste
		    else
		      omgst=omgstp
		    endif  
		  else
		    omgst=vrotst/rstar
		  endif  
	          ftemp=tstar
		  fdustt=tstar
                  if(er2.le.0.d0)goto 35
!                 correct ftemp for irradiation effect                  
	          if(irrst.eq.1.and.icomp.gt.0)then
	            farcp=far-xcp
	            fatcp=fat-ycp
	            fazcp=faz-zcp
		    cosir=far*xcp+fat*ycp+faz*zcp
		    cosir=cosir/dsqrt(xcp**2+ycp**2+zcp**2)/er
     		    coslat=(ffx**2+ffy**2)/(ffx**2+ffy**2+ffz**2)
		    coslat=dsqrt(coslat)
    		    htstb=4.d0*(1.d0-htsta)/pi
    		    fnlat=htsta+htstb*coslat
    		    pom=rcp**2/(xcp**2+ycp**2+zcp**2)	
                    t04=0.25d0*htst*(1.d0-albst)*pom*tempcp**4
		    if(cosir.gt.0.d0)then
!		      dayside
		      pom=rcp**2/(farcp**2+fatcp**2+fazcp**2)	
	              firr=pom*tempcp**4*cosir
	              tirr4=(1.d0-htst)*(1.d0-albst)*firr
	              ftemp=(tirr4+t04*fnlat+tstar**4)**0.25d0
	            else  
!		      nightside	            
                      ftemp=(t04*fnlat+tstar**4)**0.25d0
		    endif
	          endif
!                 add a spot
                  if(ispst.eq.1)then
                    fxspst=trmtst(1,1)*xspst+trmtst(1,2)*yspst          &
     &              +trmtst(1,3)*zspst
                    fyspst=trmtst(2,1)*xspst+trmtst(2,2)*yspst          &
     &              +trmtst(2,3)*zspst
                    fzspst=trmtst(3,1)*xspst+trmtst(3,2)*yspst          &
     &              +trmtst(3,3)*zspst
                    cspst=fxspst*ffx+fyspst*ffy+fzspst*ffz
                    cspst=cspst/er
!                   inside the spot                    
                    if(cspst.gt.aspst)then
                      ftemp=ftemp*tspst
                    endif
                  endif
35	          fdustd=0.d0
	          fdens=dunt1
	          fne=ane0
	          fvr=omgst*(ystar*faz-zstar*fat)+vxst
	          fvt=omgst*(zstar*far-xstar*faz)+vyst
	          fvz=omgst*(xstar*fat-ystar*far)+vzst
	          fvtrb=0.d0
	          inst=1
                endif  
              elseif(istar.eq.2.or.istar.eq.3)then
                if(dabs(faz).lt.rosup)then
                  omgst=2.d0*pi/porbst
	          ftemp=rotemp
	          fdustt=rotemp
	          fdustd=0.d0
	          fdens=dunt1
	          fne=ane0
	          fvr=omgst*(ystar*faz-zstar*fat)+vxst
	          fvt=omgst*(zstar*far-xstar*faz)+vyst
	          fvz=omgst*(xstar*fat-ystar*far)+vzst
	          fvtrb=0.d0
	          inst=1
                endif  
	      endif	      
	      return
	      end
!-----------------------------------------------------------------------
	subroutine spot(ispot,vrxsp,vrysp,vrzsp,vrotsp,rsp              &
     &  ,xsp,ysp,zsp,vxsp,vysp,vzsp                                     &
     &  ,tempsp,denssp,anesp,vtrbsp,dstdsp,dsttsp                       &
     &  ,far,fat,faz,insp                                               &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt                                                  &
     &  ,fvr,fvt,fvz,fvtrb)
!	spot=sphere
	implicit none
        integer, parameter:: i4=4,i8=8
	integer:: ispot,insp
	real(i8):: vrxsp,vrysp,vrzsp,vrotsp,rsp
	real(i8):: xsp,ysp,zsp,vxsp,vysp,vzsp
	real(i8):: tempsp,denssp,anesp,vtrbsp,dstdsp,dsttsp
	real(i8):: far,fat,faz
	real(i8):: ftemp,fdens,fne,fdustd,fdustt,fvr,fvt,fvz,fvtrb
	real(i8):: farsp,fatsp,fazsp,ermsp,unspot,omega
	if(ispot.eq.1)then
	  farsp=far-xsp		
	  fatsp=fat-ysp
	  fazsp=faz-zsp
	  ermsp=farsp*farsp+fatsp*fatsp+fazsp*fazsp
	  ermsp=dsqrt(ermsp)	
          if(ermsp.lt.rsp)then
	    unspot=vrxsp*vrxsp+vrysp*vrysp+vrzsp*vrzsp
	    unspot=dsqrt(unspot)
	    omega=vrotsp/rsp/unspot
	    ftemp=tempsp
	    fdens=denssp
	    fne=anesp
	    fdustd=dstdsp
	    fdustt=dsttsp
	    fvr=omega*(vrysp*fazsp-vrzsp*fatsp)+vxsp
	    fvt=omega*(vrzsp*farsp-vrxsp*fazsp)+vysp
	    fvz=omega*(vrxsp*fatsp-vrysp*farsp)+vzsp
	    fvtrb=vtrbsp
	    insp=1
	  endif	      
	endif	  
	return
	end    
!-----------------------------------------------------------------------
        subroutine envelope(ienv,rosue,hen,xen,yen,zen                  &
     &  ,vxen,vyen,vzen,omgen,inen                                      &
     &  ,tempen,dstden,dstten,densen,aneen,vtrben                       &
     &  ,far,fat,faz                                                    &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt                                                  &
     &  ,fvr,fvt,fvz,fvtrb)
!	envelope	      	     
	implicit none
        integer, parameter:: i4=4,i8=8
!     	include 'param.inc'
	integer:: ienv,inen
	real(i8):: rosue,hen,xen,yen,zen,vxen,vyen,vzen,omgen
	real(i8):: tempen,dstden,dstten,densen,aneen,vtrben
	real(i8):: far,fat,faz
	real(i8):: ftemp,fdens,fne,fdustd,fdustt,fvr,fvt,fvz,fvtrb
        if(ienv.eq.2.or.ienv.eq.3)then
          if(dabs(faz).lt.rosue.and.dabs(faz).lt.hen)then
	    ftemp=tempen
	    fdustt=dstten
	    fdustd=dstden
	    fdens=densen
	    fne=aneen
	    fvr=omgen*(yen*faz-zen*fat)+vxen
	    fvt=omgen*(zen*far-xen*faz)+vyen
	    fvz=omgen*(xen*fat-yen*far)+vzen
	    fvtrb=vtrben
	    inen=1
          endif  
	endif	
	return
	end      
!-----------------------------------------------------------------------
        subroutine stream(ism,dxsm,dysm,dzsm,x1sm,y1sm,z1sm             &
     &  ,r1sm,r2sm,v1sm,v2sm,vxsm,vysm,vzsm                             &
     &  ,omgsm,xsm,ysm,zsm,dsm                                          &
     &  ,vtrbsm,tempsm,dstdsm,dsttsm,denssm,anesm,edensm                &
     &  ,far,fat,faz,insm                                               &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt                                                  &
     &  ,fvr,fvt,fvz,fvtrb)
!	stream
!	dxsm,dysm,dzsm,dsm -are unit vectors and length of the stream
!	notice that although stream radius changes the streamlines 
!       do not flare but are parallel in this approximation
	implicit none
        integer, parameter:: i4=4,i8=8
	integer:: ism,insm
	real(i8):: dxsm,dysm,dzsm,x1sm,y1sm,z1sm
	real(i8):: r1sm,r2sm,v1sm,v2sm,vxsm,vysm,vzsm
	real(i8):: omgsm,xsm,ysm,zsm,dsm
	real(i8):: vtrbsm,tempsm,dstdsm,dsttsm,denssm,anesm,edensm
	real(i8):: far,fat,faz
	real(i8):: ftemp,fdens,fne,fdustd,fdustt,fvr,fvt,fvz,fvtrb
	real(i8):: rsol,tsm,xsm3,ysm3,zsm3,dxsm3,dysm3,dzsm3,er2sm,rsm
	real(i8):: vsm,fact
!	rsol=6.9599d10
!	rsol=6.95508d10
	rsol=6.957d10
	if(ism.eq.1)then
!         tsm is projection of the point on the stream vector
	  tsm=dxsm*(far-x1sm)+dysm*(fat-y1sm)+dzsm*(faz-z1sm)
!	  point1/2 are at the beginning/end of stream	  
!	  point 3 is at the distance t on the axis
	  xsm3=dxsm*tsm+x1sm
	  ysm3=dysm*tsm+y1sm
	  zsm3=dzsm*tsm+z1sm
	  dxsm3=xsm3-far
	  dysm3=ysm3-fat
	  dzsm3=zsm3-faz
	  er2sm=dxsm3*dxsm3+dysm3*dysm3+dzsm3*dzsm3
	  rsm=(r2sm-r1sm)*tsm/dsm+r1sm
	  if(tsm.gt.0.d0.and.tsm.lt.dsm.and.er2sm.lt.rsm*rsm)then
!	    vsm \sim dsqrt(tsm)+v1sm     or
!           vsm \sim (2/r-1/a)          might be more physical
	    vsm=(v2sm-v1sm)*tsm/dsm+v1sm
!	    edensm determines density drops on 1 Rsol
	    fact=v1sm/vsm*r1sm**2/rsm**2*dexp(tsm/rsol*edensm)
	    ftemp=tempsm
	    fdustt=dsttsm
	    fdustd=dstdsm*fact
	    fdens=denssm*fact
	    fne=anesm*fact
!		internal velocity + rotational drag + net velocity
	    fvr=vsm*dxsm+omgsm*(ysm*faz-zsm*fat)+vxsm
	    fvt=vsm*dysm+omgsm*(zsm*far-xsm*faz)+vysm
	    fvz=vsm*dzsm+omgsm*(xsm*fat-ysm*far)+vzsm
	    fvtrb=vtrbsm
            insm=1
	  endif
        endif
	return
	end
!-----------------------------------------------------------------------
        subroutine ring(iring,xrg,yrg,zrg                               &
     &  ,a1rg,a2rg,b1rg,b2rg,dr1rg,dr2rg,rrg,vrg,itrg                   &
     &  ,vxrg,vyrg,vzrg,xpolrg,ypolrg,zpolrg,trmtrg,dcut1               &
     &  ,vtrbrg,temprg,densrg,anerg                                     &
     &  ,dstdrg,dst2rg,edenrg,ede2rg,dsttrg                             &
     &  ,far,fat,faz,inrg                                               &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt                                                  &
     &  ,fvr,fvt,fvz,fvtrb)
!	ring
!	b1rg, b2rg -beginning, end of arc
!	a1rg, a2rg -vertical half width at beginning, end
!	dr1rg, dr2rg -horizontal half thickness the ring at beg., end
!	trmtrg is transformation matrix
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: trmtrg(3,3)
	integer:: iring,itrg,inrg
	real(i8):: xrg,yrg,zrg,a1rg,a2rg,b1rg,b2rg,dr1rg,dr2rg,rrg,vrg
	real(i8):: vxrg,vyrg,vzrg,xpolrg,ypolrg,zpolrg,dcut1
	real(i8):: vtrbrg,temprg,densrg,anerg,dstdrg,dst2rg
	real(i8):: edenrg,ede2rg,dsttrg,far,fat,faz
	real(i8):: ftemp,fdens,fne,fdustd,fdustt,fvr,fvt,fvz,fvtrb
 	real(i8):: pi,grav,dxrg,dyrg,dzrg,ffx,ffy,ffz,ffxyrg,angle
	real(i8):: ffzrg,trg,arg,drrg,trgpi,pomrg
	pi=3.1415926535897931d0
	grav=6.67408d-8
	if(iring.gt.0)then	
	  dxrg=far-xrg
	  dyrg=fat-yrg
	  dzrg=faz-zrg
!	  rotation of the coordinates to align with the ring
	  ffx=trmtrg(1,1)*dxrg+trmtrg(1,2)*dyrg+trmtrg(1,3)*dzrg
	  ffy=trmtrg(2,1)*dxrg+trmtrg(2,2)*dyrg+trmtrg(2,3)*dzrg
	  ffz=trmtrg(3,1)*dxrg+trmtrg(3,2)*dyrg+trmtrg(3,3)*dzrg
	  ffxyrg=dsqrt(ffx*ffx+ffy*ffy)
!         you may try angle=dacos(ffy)
          angle=dacos(ffx/ffxyrg)
          if(ffy.lt.0.d0)angle=2.d0*pi-angle
          if((angle.ge.b1rg.and.angle.le.b2rg).or.                      &
     &    (angle.le.b1rg.and.angle.ge.b2rg))then 
!	    within the arc
	    ffzrg=dsqrt(ffz*ffz)
	    trg=angle-b1rg
	    arg=(a2rg-a1rg)/(b2rg-b1rg)*trg+a1rg
	    if(ffzrg.lt.arg)then
!             within the vertical limit	    
	      drrg=(dr2rg-dr1rg)/(b2rg-b1rg)*trg+dr1rg
	      if(drrg.lt.rrg.and.drrg.gt.0.d0)then
	        if(ffxyrg.lt.rrg+drrg.and.ffxyrg.gt.rrg-drrg)then
!		  within the radial limit	        
	          ftemp=temprg
	          fdustt=dsttrg
		  if(densrg.lt.dcut1)then
		    trg=dsqrt(trg*trg)
		    if(itrg.eq.1)then
	             trgpi=trg/pi+1.d0
	             pomrg=dstdrg*trgpi**edenrg
	             pomrg=pomrg+dst2rg*trgpi**ede2rg
		     fdustd=pomrg*a1rg/arg*dr1rg/drrg
                     fdens=densrg*trgpi**edenrg
                     fdens=fdens*a1rg/arg*dr1rg/drrg
	             fne=anerg*trgpi**edenrg
	             fne=fne*a1rg/arg*dr1rg/drrg
	            else
	             trgpi=trg/pi 
	             pomrg=dstdrg*dexp(trgpi*edenrg)
	             pomrg=pomrg+dst2rg*dexp(trgpi*ede2rg)
		     fdustd=pomrg*a1rg/arg*dr1rg/drrg
                     fdens=densrg*dexp(trgpi*edenrg)
                     fdens=fdens*a1rg/arg*dr1rg/drrg
	             fne=anerg*dexp(trgpi*edenrg)
	             fne=fne*a1rg/arg*dr1rg/drrg
		    endif
	          else
	            fdustd=dstdrg
	            fdens=densrg
	            fne=anerg
	          endif  
	          fvr=vrg*(ypolrg*dzrg-zpolrg*dyrg)+vxrg
	          fvt=vrg*(zpolrg*dxrg-xpolrg*dzrg)+vyrg
	          fvz=vrg*(xpolrg*dyrg-ypolrg*dxrg)+vzrg
	          fvtrb=vtrbrg
	          inrg=1
!                 within the radial limit	          
	        endif
	      endif
!             within the vertical limit
	    endif      
!           within the arc	    
	  endif
	endif	
	return
	end
!-----------------------------------------------------------------------
        subroutine ufo(iufo,xufo,yufo,zufo                              &
     &  ,aufo,rinuf,routuf,emuf,ruf,xuf,yuf,zuf                         &
     &  ,vxuf,vyuf,vzuf,trmtuf                                          &
     &  ,tempuf,densuf,aneuf,vtrbuf,edenuf,ituf,etmpuf                  &
     &  ,dstduf,dsttuf                                                  &
     &  ,far,fat,faz,inuf                                               &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt                                                  &
     &  ,fvr,fvt,fvz,fvtrb)
!	disk=ufo
!	trmtuf is transformation matrix
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: trmtuf(3,3)
	integer:: iufo,ituf,inuf
	real(i8):: xufo,yufo,zufo,aufo,rinuf,routuf,emuf,ruf,xuf,yuf,zuf
	real(i8):: vxuf,vyuf,vzuf,tempuf,densuf,aneuf,vtrbuf
	real(i8):: edenuf,etmpuf,dstduf,dsttuf,far,fat,faz
	real(i8):: ftemp,fdens,fne,fdustd,fdustt,fvr,fvt,fvz,fvtrb
	real(i8):: pi,grav,dxuf,dyuf,dzuf,druf,angle,vuf,ratio,pom
	real(i8):: ffx,ffy,ffz,er2uf
	pi=3.1415926535897931d0
	grav=6.67408d-8
!	place these 4 lines inside the if blocks?? to speed up
        dxuf=far-xuf
	dyuf=fat-yuf
	dzuf=faz-zuf
	druf=dsqrt(dxuf*dxuf+dyuf*dyuf+dzuf*dzuf)
	if(iufo.eq.1)then
!         flared disk	      
!	  ffx=trmtuf(1,1)*dxuf+trmtuf(1,2)*dyuf+trmtuf(1,3)*dzuf
!	  ffy=trmtuf(2,1)*dxuf+trmtuf(2,2)*dyuf+trmtuf(2,3)*dzuf
	  ffz=trmtuf(3,1)*dxuf+trmtuf(3,2)*dyuf+trmtuf(3,3)*dzuf
	  angle=dabs(ffz)/druf
	  angle=dasin(angle)
          if(dabs(angle).lt.aufo                                        &
     &    .and.druf.gt.rinuf.and.druf.lt.routuf)then
	    vuf=dsqrt(grav*emuf/druf)/druf
	    if(ituf.eq.1)then
	      ftemp=tempuf
	      fdustt=dsttuf
	    elseif(ituf.eq.2)then
	      ratio=ruf/druf
	      pom=dsqrt(dsqrt(1.d0-dsqrt(ratio)))
	      ftemp=tempuf*ratio**0.75d0*pom
	      fdustt=dsttuf*ratio**0.75d0*pom
	    else
	      ratio=druf/rinuf
	      ftemp=tempuf*ratio**etmpuf
	      fdustt=dsttuf*ratio**etmpuf  	      
	    endif
	    ratio=druf/rinuf
	    fdens=densuf*ratio**edenuf
	    fne=aneuf*ratio**edenuf
	    fdustd=dstduf*ratio**edenuf
	    fvr=vuf*(yufo*dzuf-zufo*dyuf)+vxuf
	    fvt=vuf*(zufo*dxuf-xufo*dzuf)+vyuf
	    fvz=vuf*(xufo*dyuf-yufo*dxuf)+vzuf
	    fvtrb=vtrbuf
	    inuf=1
	  endif
	elseif(iufo.eq.2)then	
!         slab	      
!	  ffx=trmtuf(1,1)*dxuf+trmtuf(1,2)*dyuf+trmtuf(1,3)*dzuf
!	  ffy=trmtuf(2,1)*dxuf+trmtuf(2,2)*dyuf+trmtuf(2,3)*dzuf
	  ffz=trmtuf(3,1)*dxuf+trmtuf(3,2)*dyuf+trmtuf(3,3)*dzuf
          if(dabs(ffz).lt.aufo                                          &
     &    .and.druf.gt.rinuf.and.druf.lt.routuf)then
	    vuf=dsqrt(grav*emuf/druf)/druf
	    if(ituf.eq.1)then
	      ftemp=tempuf
	      fdustt=dsttuf
	    elseif(ituf.eq.2)then
	      ratio=ruf/druf
	      pom=dsqrt(dsqrt(1.d0-dsqrt(ratio)))
	      ftemp=tempuf*ratio**0.75d0*pom
	      fdustt=dsttuf*ratio**0.75d0*pom
	    else
	      ratio=druf/rinuf
	      ftemp=tempuf*ratio**etmpuf
	      fdustt=dsttuf*ratio**etmpuf  
	    endif
	    ratio=druf/rinuf
	    fdens=densuf*ratio**edenuf
	    fne=aneuf*ratio**edenuf
	    fdustd=dstduf*ratio**edenuf
	    fvr=vuf*(yufo*dzuf-zufo*dyuf)+vxuf
	    fvt=vuf*(zufo*dxuf-xufo*dzuf)+vyuf
	    fvz=vuf*(xufo*dyuf-yufo*dxuf)+vzuf
	    fvtrb=vtrbuf
	    inuf=1
	  endif
	elseif(iufo.eq.3)then	
!	  rotation ellipsoid	      
!	  rotation of the coordinates to align with the disc
	  ffx=trmtuf(1,1)*dxuf+trmtuf(1,2)*dyuf+trmtuf(1,3)*dzuf
	  ffy=trmtuf(2,1)*dxuf+trmtuf(2,2)*dyuf+trmtuf(2,3)*dzuf
	  ffz=trmtuf(3,1)*dxuf+trmtuf(3,2)*dyuf+trmtuf(3,3)*dzuf
	  er2uf=ffx*ffx/routuf/routuf+ffy*ffy/routuf/routuf
          er2uf=er2uf+ffz*ffz/aufo/aufo
          if(er2uf.lt.1.d0.and.druf.gt.rinuf)then
	    vuf=dsqrt(grav*emuf/druf)/druf
	    if(ituf.eq.1)then
	      ftemp=tempuf
	      fdustt=dsttuf
	    elseif(ituf.eq.2)then
	      ratio=ruf/druf
	      pom=dsqrt(dsqrt(1.d0-dsqrt(ratio)))
	      ftemp=tempuf*ratio**0.75d0*pom
	      fdustt=dsttuf*ratio**0.75d0*pom
	    else
	      ratio=druf/rinuf
	      ftemp=tempuf*ratio**etmpuf
	      fdustt=dsttuf*ratio**etmpuf
	    endif
	    ratio=druf/rinuf
	    fdens=densuf*ratio**edenuf
	    fne=aneuf*ratio**edenuf
	    fdustd=dstduf*ratio**edenuf
	    fvr=vuf*(yufo*dzuf-zufo*dyuf)+vxuf
	    fvt=vuf*(zufo*dxuf-xufo*dzuf)+vyuf
	    fvz=vuf*(xufo*dyuf-yufo*dxuf)+vzuf
	    fvtrb=vtrbuf
	    inuf=1
	  endif
	endif	                
	return
	end
!-----------------------------------------------------------------------??-
        subroutine nebula(inebl,xneb,yneb,zneb                          &
     &  ,aneb,rinnb,routnb,emnb,rnb,hcnb,ishdnb,hshdnb                  &
     &  ,hinvnb,tinvnb,iinvnb,hwindnb,ndennb,denxnb,denznb              &
     &  ,vxnb,vynb,vznb,trmtnb,ivelnb,hvelnb,vnb,evelnb                 &
     &  ,tempnb,densnb,anenb,vtrbnb,edennb,itnb,etmpnb                  &
     &  ,dstdnb,dsttnb                                                  &
     &  ,far,fat,faz,er,innb                                            &
     &  ,hsca00,hsca01,sig00,sig01,soun01                               &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt,kshade                                           &
     &  ,fvr,fvt,fvz,fvtrb)
!       protoplanetary or accretion disk
!       input: emnb, rnb -mass and radius of the central object
!	output: hsca00,hsca01 -vertical presure scale height
!         sig00,sig01 -surface density at rin, rout
!         soun01 -speed of sound
!         innb -indicator that state quantities were updated
!	trmtnb is transformation matrix
	implicit none
        integer, parameter:: i1=1,i4=4,i8=8
	include 'param.inc'
	real(i8):: trmtnb(3,3),denxnb(mstarx),denznb(mstarx)
	real(i8):: xneb,yneb,zneb,aneb,rinnb,routnb,emnb,rnb
	real(i8):: hinvnb,tinvnb,hwindnb,vxnb,vynb,vznb,hcnb,hshdnb
	real(i8):: hvelnb,vnb,evelnb
	real(i8):: tempnb,densnb,anenb,vtrbnb,edennb,etmpnb
	real(i8):: dstdnb,dsttnb,far,fat,faz,er
	real(i8):: hsca00,hsca01,sig00,sig01,soun01,ftemp,fdens,fne
	real(i8):: fdustd,fdustt,fvr,fvt,fvz,fvtrb
	integer(i1):: kshade
	integer:: inebl,ndennb,iinvnb,ivelnb,ishdnb,itnb,innb
	real(i8):: pi,grav,bol,hjed,ffx,ffy,ffz,erxy,erz,temp0,ratio
	real(i8):: pom,gamma,em,sound,vuf0,hscale,ftem00,ftem01,ftemp0
	real(i8):: soun00,vuf00,vuf01,vuf,sca01
	real(i8):: sig0,dens0,hh,densw0,vrx,vry,vrz,vrnb
	pi=3.1415926535897931d0
	grav=6.67408d-8
	bol=1.3806503d-16
	hjed=1.66053873d-24
	if(inebl.eq.1)then	
!	  protoplanetary disk
!	  rotation of the coordinates to align with the disc		
	  ffx=trmtnb(1,1)*far+trmtnb(1,2)*fat+trmtnb(1,3)*faz
	  ffy=trmtnb(2,1)*far+trmtnb(2,2)*fat+trmtnb(2,3)*faz
	  ffz=trmtnb(3,1)*far+trmtnb(3,2)*fat+trmtnb(3,3)*faz
          erxy=dsqrt(ffx*ffx+ffy*ffy)
!          erz=dsqrt(ffz*ffz)
          erz=dabs(ffz)
!	  define mid plane temperature temp0 depending on itnb          
          if(itnb.eq.1)then
	    temp0=tempnb
	  elseif(itnb.eq.2)then
	    ratio=rnb/erxy
	    pom=dsqrt(dsqrt(1.d0-dsqrt(ratio)))
	    temp0=tempnb*ratio**0.75d0*pom
	  else
	    ratio=erxy/rinnb
	    temp0=tempnb*ratio**etmpnb
	  endif
!         monoatommic=5/3, diatomic=7/5, many degrees of freedom=1
!         m(H2)=2, m(He)=4
!         gamma=5.d0/3.d0
!         gamma=7.d0/5.d0
          gamma=1.d0
          em=(45.d0*2.d0+10.d0*4.d0)/(45.d0+10.d0)*hjed
	  sound=dsqrt(gamma*bol*temp0/em)        
!	  velocity in cylindrical coord.
          vuf0=dsqrt(grav*emnb/erxy)
!         scale height H(r)          
          hscale=sound/vuf0*erxy*hcnb
!??	  z=hscale*2 corresponds to 10^-1 drop in density          
!	  z=hscale*3 corresponds to 10^-2 drop in density          
!	  z=hscale*4 corresponds  3x10^-4 drop in density          
!         z=hscale*7 corresponds to 10^-10 drop in density
!         erz.lt.hscale*6.d0 may be required to reach wind (1d-8rho0)
!	  soun* -sound speed, hsca* -disk scale hight
!	  vuf* -Keplerian velocity
!         *00 -quantity at rinnb, *01 -quantity at routnb
!	  you may slightly speed up if you take the calculations of
!         *00, *01 quantities out of this subroutine
!	  inside the nebula:
          if(erxy.gt.rinnb.and.erxy.lt.routnb.and.erz.lt.hscale*aneb)   &
     &    then
	    if(itnb.eq.1)then
!		constant temperatures	    
	      ftem00=tempnb
	      ftem01=tempnb	
	      ftemp=tempnb
	      fdustt=dsttnb
	    elseif(itnb.eq.2)then
!		radially dependent temperatures (Pringle 1981)	    
	      ratio=rnb/rinnb
	      pom=dsqrt(dsqrt(1.d0-dsqrt(ratio)))
	      ftem00=tempnb*ratio**0.75d0*pom
	      ratio=rnb/routnb
	      pom=dsqrt(dsqrt(1.d0-dsqrt(ratio)))
	      ftem01=tempnb*ratio**0.75d0*pom
	      ratio=rnb/erxy
	      pom=dsqrt(dsqrt(1.d0-dsqrt(ratio)))
	      ftemp=tempnb*ratio**0.75d0*pom
	      fdustt=dsttnb*ratio**0.75d0*pom
	    else
!	      radial power law & vertical temperature dependence	    
	      ftem00=tempnb
	      ratio=routnb/rinnb
	      ftem01=tempnb*ratio**etmpnb
	      ratio=erxy/rinnb
!	      mid-plane temperature	      
	      ftemp0=tempnb*ratio**etmpnb
!             gas temperature inversion	      
	      if(erz.gt.hscale*hinvnb)then	
!	        ftemp=2.d0*ftemp0/hscale*(erz-hscale*5.5d0)+ftemp0
		if(iinvnb.eq.1)then
	          ftemp=tinvnb*ftemp0
		elseif(iinvnb.eq.2)then
   		  ftemp=(tinvnb-1.d0)/(aneb-hinvnb)*ftemp0/hscale
		  ftemp=ftemp*(erz-hinvnb*hscale)+ftemp0
	        else
		  ftemp=ftemp0
		endif
	      else
	        ftemp=ftemp0
	      endif  
	      fdustt=dsttnb*ratio**etmpnb
	    endif
	    soun00=dsqrt(gamma*bol*ftem00/em)
	    soun01=dsqrt(gamma*bol*ftem01/em)
	    vuf00=dsqrt(grav*emnb/rinnb)
	    vuf01=dsqrt(grav*emnb/routnb)
	    hsca00=soun00/vuf00*rinnb*hcnb
	    hsca01=soun01/vuf01*routnb*hcnb
	    sig00=densnb*hsca00*dsqrt(2.d0*pi)
	    sig01=sig00*(routnb/rinnb)**edennb
	    ratio=(erxy/rinnb)**edennb
            sig0=sig00*ratio
            dens0=sig0/hscale/dsqrt(2.d0*pi)
!	    wind region in densities
	    if(ndennb.eq.0)then
!             gaussian-like wind profile	    
              if(erz.lt.hscale*hwindnb)then
	        fdens=dens0*dexp(-erz**2/hscale**2/2.d0)
	      else
	        fdens=dens0*dexp(-hwindnb**2/2.d0)
	      endif  
	    else
!             precalculated wind profile	    
	      hh=hsca01/5.374d10
	      densw0=dens0*hh
	      call intrp(denxnb,denznb,ndennb,erz,fdens)
	      fdens=fdens/denznb(1)*densw0
	    endif  
	    fne=anenb/densnb*fdens
	    fdustd=dstdnb/densnb*fdens
	    if(erz.lt.hscale*hshdnb)then
	      kshade=0
	    else
	      kshade=ishdnb
	    endif
!??	    write(19,*)
!	    write(19,'(6es10.2)')erxy,erz,rinnb,hscale,hsca00
!	    write(19,'(6es10.2)')dens0,sig0,sig00
!	    write(19,'(6es10.2)')fdens,fne,fdustd
!		velocities
!	    v = omega x r = omega x erxy 
!           but v is not the same as for the DISK since 
!           omega=sqrt(GM/erxy**3)=vuf
!	    radial vel. component
            if(ivelnb.eq.1)then
              if(erz.ge.hscale*hvelnb)then
                vrnb=vnb*(1.d0-rnb/er)**evelnb
                vrx=vrnb*ffx/er
                vry=vrnb*ffy/er
                vrz=vrnb*ffz/er
              else
                vrx=0.d0
                vry=0.d0
                vrz=0.d0
              endif  
            else
              vrx=0.d0
              vry=0.d0
              vrz=0.d0
            endif	    
	    vuf=vuf0/erxy
	    fvr=vuf*(yneb*faz-zneb*fat)+vrx+vxnb
	    fvt=vuf*(zneb*far-xneb*faz)+vry+vynb
	    fvz=vuf*(xneb*fat-yneb*far)+vrz+vznb
	    fvtrb=vtrbnb
	    innb=1
!	    inside the nebula
	  endif
!         inside inebl	  
	endif	                
	return
	end
!-----------------------------------------------------------------------
        subroutine sjet(ijet,ajet,rinjt,routjt                          &
     &  ,xjt,yjt,zjt,xjet,yjet,zjet                                     &
     &  ,ivjt,vjt,eveljt,rcjt,vxjt,vyjt,vzjt,vtrbjt                     &
     &  ,tempjt,densjt,anejt,dstdjt,dsttjt,etmpjt,asymjt                &
     &  ,far,fat,faz,injt                                               &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt                                                  &
     &  ,fvr,fvt,fvz,fvtrb)
!	jet with one or 2 cones
!	streamlines diverge
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: ajet,rinjt,routjt,xjt,yjt,zjt,xjet,yjet,zjet
	real(i8):: dxjt,dyjt,dzjt,drjt
	real(i8):: vjt,eveljt,rcjt,vxjt,vyjt,vzjt,vtrbjt
	real(i8):: tempjt,densjt,anejt,etmpjt,asymjt
	real(i8):: dstdjt,dsttjt,far,fat,faz,ftemp,fdens,fne
	real(i8):: fdustd,fdustt,fvr,fvt,fvz,fvtrb
	integer:: ijet,injt,ivjt
	real(i8):: pi,angle,rr2,rr,vel,vel1,velr
	pi=3.1415926535897931d0
	if(ijet.eq.1)then
!         a single jet cone	
	  dxjt=far-xjt
	  dyjt=fat-yjt
	  dzjt=faz-zjt
	  drjt=dsqrt(dxjt*dxjt+dyjt*dyjt+dzjt*dzjt)
	  angle=(xjet*dxjt+yjet*dyjt+zjet*dzjt)/drjt
	  angle=dacos(angle)
          if(dabs(angle-0.d0).lt.ajet                                   &
     &    .and.drjt.ge.rinjt.and.drjt.lt.routjt)then
            rr2=(rinjt/drjt)**2
	    rr=drjt/rinjt
	    ftemp=tempjt*rr**etmpjt
	    fdustt=dsttjt*rr**etmpjt
	    if(ivjt.eq.2)then
	      vel=vjt*(1.d0-rcjt/drjt)**eveljt
	      vel1=vjt*(1.d0-rcjt/rinjt)**eveljt
	    else
	      vel=vjt*rr**eveljt
	      vel1=vjt
	    endif
	    velr=vel1/vel
	    fdustd=dstdjt*rr2*velr
	    fdens=densjt*rr2*velr
	    fne=anejt*rr2*velr
	    fvr=vel*dxjt/drjt+vxjt
	    fvt=vel*dyjt/drjt+vyjt
	    fvz=vel*dzjt/drjt+vzjt
	    fvtrb=vtrbjt
	    injt=1
	  endif
	elseif(ijet.eq.2)then
!         two jet cones		
	  dxjt=far-xjt
	  dyjt=fat-yjt
	  dzjt=faz-zjt
	  drjt=dsqrt(dxjt*dxjt+dyjt*dyjt+dzjt*dzjt)
	  angle=(xjet*dxjt+yjet*dyjt+zjet*dzjt)/drjt
	  angle=dacos(angle)	
          if((dabs(angle-0.d0).lt.ajet.or.dabs(angle-pi).lt.ajet)       &
     &    .and.drjt.ge.rinjt.and.drjt.lt.routjt)then
            rr2=(rinjt/drjt)**2
            rr=drjt/rinjt
	    ftemp=tempjt*rr**etmpjt
	    fdustt=dsttjt*rr**etmpjt
	    if(ivjt.eq.2)then
	      vel=vjt*(1.d0-rcjt/drjt)**eveljt
	      vel1=vjt*(1.d0-rcjt/rinjt)**eveljt
	    else
	      vel=vjt*rr**eveljt
	      vel1=vjt
	    endif
	    velr=vel1/vel
	    fdustd=dstdjt*rr2*velr
	    fdens=densjt*rr2*velr
            if (dabs(angle-0.d0).lt.ajet) then
              fdens=fdens*(1.0d0+asymjt)
            else
              fdens=fdens*(1.0d0-asymjt)
            endif
	    fne=anejt*rr2*velr
	    fvr=vel*dxjt/drjt+vxjt
	    fvt=vel*dyjt/drjt+vyjt
	    fvz=vel*dzjt/drjt+vzjt
	    fvtrb=vtrbjt
	    injt=1
	  endif
	endif	
	return
	end
!-----------------------------------------------------------------------
        subroutine shell(ishell,rinsh,routsh,rcsh                       &
     &  ,vsh,vxsh,vysh,vzsh,vtrbsh                                      &
     &  ,tempsh,dstdsh,dsttsh,denssh,anesh,evelsh,etmpsh                &
     &  ,far,fat,faz,er,insh                                            &
     &  ,ftemp,fdens,fne                                                &
     &  ,fdustd,fdustt                                                  &
     &  ,fvr,fvt,fvz,fvtrb)
!	shell
	implicit none
        integer, parameter:: i4=4,i8=8
	real(i8):: rinsh,routsh,rcsh,vsh,vxsh,vysh,vzsh,vtrbsh
	real(i8):: tempsh,dstdsh,dsttsh,denssh,anesh,evelsh,etmpsh
	real(i8):: far,fat,faz,er,ftemp,fdens,fne,fdustd,fdustt
	real(i8):: fvr,fvt,fvz,fvtrb
	integer:: ishell,insh
	real(i8):: vshr,rr2vv,vshin
        if(ishell.eq.1.and.er.ge.rinsh.and.er.lt.routsh)then
          ftemp=tempsh
          fdustt=dsttsh
          fdustd=dstdsh
          fdens=denssh
          fne=anesh   
          fvr=vsh*far/er+vxsh
          fvt=vsh*fat/er+vysh
          fvz=vsh*faz/er+vzsh
          fvtrb=vtrbsh
          insh=1
        elseif(ishell.eq.2.and.er.ge.rinsh.and.er.lt.routsh)then
	  vshr=vsh*(er/rinsh)**(evelsh)
	  vshin=vsh
	  rr2vv=rinsh*rinsh/(er*er)*vshin/vshr
	  fvr=vshr*far/er+vxsh
	  fvt=vshr*fat/er+vysh
	  fvz=vshr*faz/er+vzsh
	  fvtrb=vtrbsh
	  ftemp=tempsh
	  fdustt=dsttsh
	  fdustd=dstdsh*rr2vv
	  fdens=denssh*rr2vv
	  fne=anesh*rr2vv
	  insh=1
	elseif(ishell.eq.3.and.er.ge.rinsh.and.er.lt.routsh)then
	  vshr=vsh*(1.d0-rcsh/er)**evelsh
	  vshin=vsh*(1.d0-rcsh/rinsh)**evelsh
	  rr2vv=rinsh*rinsh/(er*er)*vshin/vshr
	  fvr=vshr*far/er+vxsh
	  fvt=vshr*fat/er+vysh
	  fvz=vshr*faz/er+vzsh
	  fvtrb=vtrbsh
	  ftemp=tempsh*(er/rinsh)**etmpsh
	  fdustt=dsttsh
	  fdustd=dstdsh*rr2vv
	  fdens=denssh*rr2vv
	  fne=anesh*rr2vv
	  insh=1
	endif	  
	return
	end
!-----------------------------------------------------------------------
        subroutine grid3d(ndim1,ndim2,ndim3                             &
     &  ,nbod1,nbod1a,nbod1b                                            &
     &  ,nbod2,nbod2a,nbod2b                                            &
     &  ,nbod3,nbod3a,nbod3b                                            &          
     &  ,rmdx1,rmdx2,rmdx3,rmdx4                                        &
     &  ,rmdy1,rmdy2,rmdy3,rmdy4                                        &     
     &  ,rmdz1,rmdz2,rmdz3,rmdz4                                        &          
     &  ,gainx,gainxb,gainy,gainyb,gainz,gainzb,ar,at,az)
!       definition of the 3D grid
!	it allows to merge two grids 
!	input: interval & number of points & gain
!              code adjusts steps accordingly
!	output: ar,at,az,nbod1,nbod2,nbod3
        implicit none
        integer, parameter:: i4=4,i8=8
        integer:: ndim1,ndim2,ndim3
        real(i8):: ar(ndim1),at(ndim2),az(ndim3)
        integer:: nbod1,nbod1a,nbod1b
        integer:: nbod2,nbod2a,nbod2b
        integer:: nbod3,nbod3a,nbod3b
        real(i8):: rmdx1,rmdx2,rmdx3,rmdx4
        real(i8):: rmdy1,rmdy2,rmdy3,rmdy4
        real(i8):: rmdz1,rmdz2,rmdz3,rmdz4
        real(i8):: gainx,gainy,gainz,gainxb,gainyb,gainzb  
        call grid1d(ndim1,nbod1,nbod1a,nbod1b                           &
     &  ,rmdx1,rmdx2,rmdx3,rmdx4,gainx,gainxb,ar)
        call grid1d(ndim2,nbod2,nbod2a,nbod2b                           &
     &  ,rmdy1,rmdy2,rmdy3,rmdy4,gainy,gainyb,at)
        call grid1d(ndim3,nbod3,nbod3a,nbod3b                           &
     &  ,rmdz1,rmdz2,rmdz3,rmdz4,gainz,gainzb,az)
        return
        end
!-----------------------------------------------------------------------
        subroutine grid1d(ndim1,nbod1,nbod1a,nbod1b                     &
     &  ,rmdx1,rmdx2,rmdx3,rmdx4,gainx,gainxb,ar)
!       definition of the x coordinate grid
!	it allows to merge two grids 
!	gain?=1 grid is equidistant
!       gain?>1. step increases geometrically from the center 
!          of the interval to the left and to the right
!	input: interval & number of points & gain
!              code adjusts steps accordingly
!	output: ar,nbod1
	implicit none
        integer, parameter:: i4=4,i8=8
        integer:: ndim1,nbod1,nbod1a,nbod1b
        real(i8):: ar(ndim1),ar1(ndim1),ar2(ndim1)
	real(i8):: rmdx1,rmdx2,rmdx3,rmdx4
	real(i8):: gainx,gainxb
	real(i8):: rstepa,rstepb,st1
	integer:: i,k1,k2,n2
!	geom. sequence assumes that nbod1a,nbod1b are odd
!		x
!	1st grid
   	if(gainx.eq.1.d0)then
!	  equidistant step   	
          rstepa=(rmdx2-rmdx1)/dble(nbod1a-1)
          do i=1,nbod1a
            ar1(i)=rmdx1+rstepa*dble(i-1)
          enddo
	else
!	  geometrical sequence of steps
          rstepa=(rmdx2-rmdx1)/dble(nbod1a-1)
   	  n2=nbod1a/2
   	  st1=0.5d0*(rmdx2-rmdx1)*(1.d0-gainx)/(1.d0-gainx**n2)
   	  ar1(n2+1)=0.5d0*(rmdx1+rmdx2)   	  
   	  do i=1,n2  
   	    ar1(n2+1+i)=st1*(1.d0-gainx**i)/(1.d0-gainx)+ar1(n2+1)
   	    ar1(n2+1-i)=-st1*(1.d0-gainx**i)/(1.d0-gainx)+ar1(n2+1)
 	  enddo
        endif
!	2nd grid
        if(nbod1b.gt.2)then
          if(gainxb.eq.1.d0)then
!	    equidistant step          
            rstepb=(rmdx4-rmdx3)/dble(nbod1b-1)
            do i=1,nbod1b
              ar2(i)=rmdx3+rstepb*dble(i-1)
            enddo           
	  else
!           geometrical sequence of steps
            rstepb=(rmdx4-rmdx3)/dble(nbod1b-1)
   	    n2=nbod1b/2
   	    st1=0.5d0*(rmdx4-rmdx3)*(1.d0-gainxb)/(1.d0-gainxb**n2)
   	    ar2(n2+1)=0.5d0*(rmdx3+rmdx4)   	  
   	    do i=1,n2  
   	      ar2(n2+1+i)=st1*(1.d0-gainxb**i)/(1.d0-gainxb)
   	      ar2(n2+1-i)=-st1*(1.d0-gainxb**i)/(1.d0-gainxb)
 	    enddo
 	  endif   	  
 	else
 	  rstepb=1.d300  
 	endif
!       combine the two grids
        if(rstepa.le.rstepb)then
!         first grid is finer          
          if(nbod1b.gt.2)then  
            k1=0
            do i=1,nbod1b
              if(ar2(i).lt.ar1(1))then
                k1=k1+1
              endif
            enddo    
!           place 2nd grid points first              
            if(k1.gt.0)then
              do i=1,k1
                ar(i)=ar2(i)
              enddo
            endif
!           add 1st grid points
            do i=1+k1,nbod1a+k1
              ar(i)=ar1(i-k1)
            enddo
            k2=0
            do i=1,nbod1b
              if(ar2(i).gt.ar1(nbod1a))then
                k2=k2+1
              endif
            enddo    
!           add 2nd grid points              
            if(k2.gt.0)then
              do i=1,k2
                ar(nbod1a+k1+i)=ar2(nbod1b-k2+i)
              enddo
            endif  
            nbod1=nbod1a+k1+k2
          else  
            do i=1,nbod1a
              ar(i)=ar1(i)
            enddo
            nbod1=nbod1a
          endif  
        else
!	  second grid is finer
          if(nbod1b.gt.2)then  
            k1=0
            do i=1,nbod1a
              if(ar1(i).lt.ar2(1))then
                k1=k1+1
              endif
            enddo    
!           place 1st grid points first              
            if(k1.gt.0)then
              do i=1,k1
                ar(i)=ar1(i)
              enddo
            endif
!           add 2nd grid points
            do i=1+k1,nbod1b+k1
              ar(i)=ar2(i-k1)
            enddo
            k2=0
            do i=1,nbod1a
              if(ar1(i).gt.ar2(nbod1b))then
                k2=k2+1
              endif
            enddo    
!           add 1st grid points              
            if(k2.gt.0)then
              do i=1,k2
                ar(nbod1b+k1+i)=ar1(nbod1a-k2+i)
              enddo
            endif  
            nbod1=nbod1b+k1+k2
          else  
            do i=1,nbod1a
              ar(i)=ar1(i)
            enddo
            nbod1=nbod1a
          endif  
        endif  
	return
	end	

!-----------------------------------------------------------------------
      	FUNCTION GAUNT(I,FR)
!	taken from the Syspec code:
! Hubeny I., Lanz T., Jeffery C.S., 1994, in Newsletter on Analysis
! of Astronomical spectra No.20, ed. C.S. Jeffery (CCP7; St. Andrews:
! St. Andrews Univ.), 30
!
!     Hydrogenic bound-free Gaunt factor for the principal quantum
!     number I and frequency FR
!
	implicit none 
        integer, parameter:: i4=4,i8=8
	real(i8):: FR
	integer:: I
	real(i8):: GAUNT
	real(i8):: X
      X=FR/2.99793E14
      GAUNT=1.
      IF(I.GT.10) GO TO 16
      GO TO (1,2,3,4,5,6,7,8,9,10),I
    1 GAUNT=1.2302628+X*(-2.9094219E-3+X*(7.3993579E-6-8.7356966E-9*X)) &
     & +(12.803223/X-5.5759888)/X
      GO TO 16
    2 GAUNT=1.1595421+X*(-2.0735860E-3+2.7033384E-6*X)+(-1.2709045+     &
     & (-2.0244141/X+2.1325684)/X)/X
      GO TO 16
    3 GAUNT=1.1450949+X*(-1.9366592E-3+2.3572356E-6*X)+(-0.55936432+    &
     & (-0.23387146/X+0.52471924)/X)/X
      GO TO 16
    4 GAUNT=1.1306695+X*(-1.3482273E-3+X*(-4.6949424E-6+2.3548636E-8*X))&
     & +(-0.31190730+(0.19683564-5.4418565E-2/X)/X)/X
      GO TO 16
    5 GAUNT=1.1190904+X*(-1.0401085E-3+X*(-6.9943488E-6+2.8496742E-8*X))&
     & +(-0.16051018+(5.5545091E-2-8.9182854E-3/X)/X)/X
      GO TO 16
    6 GAUNT=1.1168376+X*(-8.9466573E-4+X*(-8.8393133E-6+3.4696768E-8*X))&
     & +(-0.13075417+(4.1921183E-2-5.5303574E-3/X)/X)/X
      GO TO 16
    7 GAUNT=1.1128632+X*(-7.4833260E-4+X*(-1.0244504E-5+3.8595771E-8*X))&
     & +(-9.5441161E-2+(2.3350812E-2-2.2752881E-3/X)/X)/X
      GO TO 16
    8 GAUNT=1.1093137+X*(-6.2619148E-4+X*(-1.1342068E-5+4.1477731E-8*X))&
     & +(-7.1010560E-2+(1.3298411E-2 -9.7200274E-4/X)/X)/X
      GO TO 16
    9 GAUNT=1.1078717+X*(-5.4837392E-4+X*(-1.2157943E-5+4.3796716E-8*X))&
     & +(-5.6046560E-2+(8.5139736E-3-4.9576163E-4/X)/X)/X
      GO TO 16
   10 GAUNT=1.1052734+X*(-4.4341570E-4+X*(-1.3235905E-5+4.7003140E-8*X))&
     & +(-4.7326370E-2+(6.1516856E-3-2.9467046E-4/X)/X)/X
   16 RETURN
	END FUNCTION GAUNT
!-----------------------------------------------------------------------
 	FUNCTION GFREE(T,FR)
!	taken from the Syspec code:
! Hubeny I., Lanz T., Jeffery C.S., 1994, in Newsletter on Analysis
! of Astronomical spectra No.20, ed. C.S. Jeffery (CCP7; St. Andrews:
! St. Andrews Univ.), 30
!
!	Hydrogenic free-free Gaunt factor, for temperature T and
!	frequency FR
!
	implicit none 
        integer, parameter:: i4=4,i8=8
	real(i8):: T,FR
	real(i8):: GFREE
	real(i8):: THET,X,C1,C2,C3,C4
 	THET=5040.4/T
 	IF(THET.LT.4.E-2) THET=4.E-2
 	X=FR/2.99793E14
 	IF(X.GT.1) GO TO 10
 	IF(X.LT.0.2) X=0.2 
 	GFREE=(1.0823+2.98E-2/THET)+(6.7E-3+1.12E-2/THET)/X
 	RETURN
10 	C1=(3.9999187E-3-7.8622889E-5/THET)/THET+1.070192
      	C2=(6.4628601E-2-6.1953813E-4/THET)/THET+2.6061249E-1
      	C3=(1.3983474E-5/THET+3.7542343E-2)/THET+5.7917786E-1
      	C4=3.4169006E-1+1.1852264E-2/THET
      	GFREE=((C4/X-C3)/X+C2)/X+C1
      	RETURN
      	END FUNCTION GFREE
!-----------------------------------------------------------------------
	SUBROUTINE LOCATE(XX,N,X,J)
!	taken from Numerical recipes, 
!	search an ordered table by bisection
!	given array xx and value x returns j such that 
!	x is between xx(j),xx(j+1)
      	implicit none 
        integer, parameter:: i4=4,i8=8
!      	include 'param.inc'
!      	DIMENSION XX(mstarx)
!?? potential problem with dynamic declaration
	real(i8):: XX(N),x
	integer:: n,j
	integer:: jl,ju,jm
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GE.XX(1)).EQV.(X.GE.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
        GO TO 10
      ENDIF
      if(x.eq.xx(1))then
        j=1
      elseif(x.eq.xx(n))then
       j=n-1
      else
        J=JL
      endif  
      RETURN
      END
!-----------------------------------------------------------------------
      	SUBROUTINE LOCATE2(XX,N,X,J,k)
!	taken from Numerical recipes, 
!	search an ordered table by bisection
	implicit none 
        integer, parameter:: i4=4,i8=8
      	include 'param.inc'
      	real(i8):: xx(mstarx,mspecx),x
	integer:: n,j,k
	integer:: jl,ju,jm
!      	DIMENSION XX(N,k)
      	JL=0
      	JU=N+1
10    	IF(JU-JL.GT.1)THEN
          JM=(JU+JL)/2
          IF((XX(N,k).GT.XX(1,k)).EQV.(X.GT.XX(JM,k)))THEN
            JL=JM
          ELSE
            JU=JM
          ENDIF
          GO TO 10
        ENDIF
        J=JL
        RETURN
        END      
!-----------------------------------------------------------------------
        subroutine bcond1(iunt,darkl,cdelta,vdif,nbod,dens              &
     &  ,dcut1,dcut2,dcut3,dcutn,vr,vt,vz                               &
     &  ,ari,ati,azn,dinc,alpha                                         &
     &  ,istar,rstar,vxstr,vystr,vzstr,dlst,dlst2,irrst                 &
     &  ,icomp,rcp,xcpr,ycpr,zcpr,vxcpr,vycpr,vzcpr                     &
     &  ,dlcp,dlcp2,irrcp,xcp,qq,iprint)
!       Identifies the last untransparent point along the ray and
!	calculates quantities which are not frequency dependent:
!	limb darkening factor if applicable, cos of the irr. angle,
!	radial velocities.
!	The limb darkening applies only to 
!	'star' and 'companion' and you must set dlst=dlcp=0.d0 if you
!	want to use their density intervals for other objects
!	(without limb darkening). 
!       input:  nbod -number of points
!       	dens - density
!       	dcut1,2,3,n -density interval boundaries
!		vr,vt,vz -x,y,z velocity along the ray
!		ari,ati -x,y coordinates of the ray (line of sight)
!		azn - set of z coordinates
!		dinc - inclination
!		alpha - phase angle
!		istar - shape of the primary (1=sphere, 
!			2=Roche detached, 3=Roche contact)
!		rstar -radius of the 'star'
!               vxstr,vystr,vzstr -velocity of the center of 'star'
!		dlst -limb darkening coefficients of the 'star'
!		icomp - shape of the secondary (1=sphere,
!			2=Roche detached)
!		rcp -radius of the 'companion'
!		xcpr,ycpr,zcpr -position of the center of 'companion'
!			in the rotated line of sight coordinates
!		vxcpr,vycpr,vzcpr -velocity of the center of 'companion'
!		dlcp -limb darkening coeficients of the 'companion'
!		xcp - scaling factor (separation of stars)
!		qq - mass ratio
!		iprint - printing switch (if =1)
!       output: 
!	  iunt -number of the last untransparent point
!	  darkl -limb darkening factor
!	  cdelta -cos(irradiation angle)*R**2/separation**2
!	  vdif -differential radial velocity between the reflecting 
!		surface and the source of light
!              
	implicit none   
        integer, parameter:: i4=4,i8=8
        include 'param.inc'
        real(i8):: dens(ndim),azn(ndim),vr(ndim),vt(ndim),vz(ndim)
	integer:: iunt,nbod,istar,irrst,icomp,irrcp,iprint
	real(i8):: darkl,cdelta,vdif,dcut1,dcut2,dcut3,dcutn,ari,ati
	real(i8):: dinc,alpha,rstar,vxstr,vystr,vzstr,dlst,dlst2
	real(i8):: rcp,xcpr,ycpr,zcpr,vxcpr,vycpr,vzcpr,dlcp,dlcp2
	real(i8):: xcp,qq
	real(i8):: grmin,cdinc,sdinc,calpha,salpha,fx,fy,fz
	real(i8):: xlos,ylos,zlos,potx,poty,potz,gravt,ctheta,rvec,sth2
	real(i8):: rr2,vdif1,vdif2,vdif3,dari,dati,rvecx,rvecy,rvecz
	integer:: i
        iunt=0
        darkl=1.d0
        cdelta=0.d0
        vdif=0.d0
!       icut=0
!       grmin -optional small value??
	grmin=1.d-20
        do 10 i=nbod,1,-1
          if(dens(i).gt.dcut1.and.dens(i).le.dcut2)then
            iunt=i
	    if(istar.gt.1)then
              cdinc=dcos(dinc)
              sdinc=dsin(dinc)
              calpha=dcos(alpha)
              salpha=dsin(alpha)
!		coordinates of the last untransparent point 
!		in the body frozen frame
              call rotz(ari,ati,azn(iunt),fx,fy,fz                      &
     &        ,cdinc,sdinc,calpha,salpha)
	      fx=fx/xcp
	      fy=fy/xcp
	      fz=fz/xcp	
!		line of sight coordinates in the body frozen frame
              call rotz(0.d0,0.d0,1.d0,xlos,ylos,zlos                   &
     &        ,cdinc,sdinc,calpha,salpha)
!		gravity vector in the last untransparent point
	      call gravn(fx,fy,fz,qq,potx,poty,potz,gravt)
	      if(gravt.lt.grmin)then
		ctheta=0.d0
	      else
 	        ctheta=-(xlos*potx+ylos*poty+zlos*potz)/gravt
		if(ctheta.lt.0.d0)ctheta=0.d0
	      endif
! 		reflection off the roche surface
	      if(irrst.eq.1.and.istar.eq.2)then
		rvec=(fx-1.d0)**2+fy**2+fz**2
		if(gravt.lt.grmin)then
		  cdelta=1.d0
		else
		  cdelta=((fx-1.d0)*potx+fy*poty+fz*potz)
		  cdelta=cdelta/dsqrt(rvec)/gravt
		endif
	        cdelta=cdelta*(rcp/xcp)**2/rvec
	      else
	        cdelta=0.d0
	      endif	
	      vdif=0.d0
	    else	
              sth2=(ari*ari+ati*ati)/(rstar*rstar)
              if(sth2.lt.1.d0)then
                ctheta=dsqrt(1.d0-sth2)
              else
                ctheta=0.d0
              endif
!		reflection off the sphere
	      if(irrst.eq.1.and.istar.eq.1)then
	        cdelta=ari*xcpr+ati*ycpr+azn(iunt)*zcpr
	        cdelta=cdelta/dsqrt(ari**2+ati**2+azn(iunt)**2)
	        cdelta=cdelta/dsqrt(xcpr**2+ycpr**2+zcpr**2)
	        rr2=(ari-xcpr)**2+(ati-ycpr)**2+(azn(iunt)-zcpr)**2
	        cdelta=cdelta*rcp**2/rr2
     	        vdif1=(ari-xcpr)*(vr(iunt)-vxcpr)
     	        vdif2=(ati-ycpr)*(vt(iunt)-vycpr)
     	        vdif3=(azn(iunt)-zcpr)*(vz(iunt)-vzcpr)
     	        vdif=-(vdif1+vdif2+vdif3)/dsqrt(rr2)	        
	      else
	        cdelta=0.d0
		vdif=0.d0
	      endif
	    endif
            darkl=1.d0-dlst*(1.d0-ctheta)-dlst2*(1.d0-ctheta)**2
	    goto 20
!           icut=1
          elseif(dens(i).gt.dcut2.and.dens(i).le.dcut3)then
            iunt=i
	    if(icomp.gt.1)then
              cdinc=dcos(dinc)
              sdinc=dsin(dinc)
              calpha=dcos(alpha)
              salpha=dsin(alpha)
!		coordinates of the last untransparent point 
!		in the body frozen frame
              call rotz(ari,ati,azn(iunt),fx,fy,fz                      &
     &        ,cdinc,sdinc,calpha,salpha)
	      fx=fx/xcp
	      fy=fy/xcp
	      fz=fz/xcp	
!		line of sight coordinates in the body frozen frame
              call rotz(0.d0,0.d0,1.d0,xlos,ylos,zlos                   &
     &        ,cdinc,sdinc,calpha,salpha)
!		gravity vector in the last untransparent point
              call gravn(fx,fy,fz,qq,potx,poty,potz,gravt)
	      if(gravt.lt.grmin)then
		ctheta=0.d0
	      else
 	        ctheta=-(xlos*potx+ylos*poty+zlos*potz)/gravt
		if(ctheta.lt.0.d0)ctheta=0.d0
	      endif
! 		reflection off the roche surface
	      if(irrcp.eq.1.and.icomp.eq.2)then
		rvec=fx**2+fy**2+fz**2
		if(gravt.lt.grmin)then
		  cdelta=1.d0
		else
		  cdelta=(fx*potx+fy*poty+fz*potz)/dsqrt(rvec)/gravt
		endif
	        cdelta=cdelta*(rstar/xcp)**2/rvec
	      else
	        cdelta=0.d0
	      endif	
	      vdif=0.d0
	    else	
	      dari=ari-xcpr
	      dati=ati-ycpr
	      sth2=(dari*dari+dati*dati)/(rcp*rcp)
	      if(sth2.lt.1.d0)then
	        ctheta=dsqrt(1.d0-sth2)
	      else
	        ctheta=0.d0
	      endif
!		reflection off the sphere
	      if(irrcp.eq.1.and.icomp.eq.1)then
	        rvecx=ari-xcpr
	        rvecy=ati-ycpr
	        rvecz=azn(iunt)-zcpr
	        cdelta=-rvecx*xcpr-rvecy*ycpr-rvecz*zcpr
	        cdelta=cdelta/dsqrt(rvecx**2+rvecy**2+rvecz**2)
	        cdelta=cdelta/dsqrt(xcpr**2+ycpr**2+zcpr**2)
	        rr2=ari**2+ati**2+azn(iunt)**2	
	        cdelta=cdelta*rstar**2/rr2
     	        vdif1=ari*(vr(iunt)-vxstr)
     	        vdif2=ati*(vt(iunt)-vystr)
     	        vdif3=azn(iunt)*(vz(iunt)-vzstr)
     	        vdif=-(vdif1+vdif2+vdif3)/dsqrt(rr2)	      	        
	      else
	        cdelta=0.d0
	        vdif=0.d0
	      endif
	    endif
	    darkl=1.d0-dlcp*(1.d0-ctheta)-dlcp2*(1.d0-ctheta)**2
	    goto 20
!           icut=2
          elseif(dens(i).gt.dcut3.and.dens(i).le.dcutn)then
            iunt=i
	    darkl=1.d0
	    goto 20
!           icut=3
          elseif(dens(i).gt.dcutn)then
            iunt=i
	    darkl=1.d0
	    goto 20
!           icut=4
          endif
10      continue
20	return
	end
!-----------------------------------------------------------------------
        subroutine bcond2(iunt,dens,temp,vz,vdif,freq                   &
     &  ,dcut1,dcut2,dcut3,dcutn,darkl,cdelta                           &
     &  ,lunt1,nstar1,xstar1,star1,nspec1,tspec1,alb1x,alb1y,nalb1      &
     &  ,tstar,dlst,dlst2,irrst,albst,wstar1,fstar1,jstar1              &
     &  ,lunt2,nstar2,xstar2,star2,nspec2,tspec2,alb2x,alb2y,nalb2      &
     &  ,tempcp,dlcp,dlcp2,irrcp,albcp,wstar2,fstar2,jstar2             &
     &  ,lunt3,nstar3,xstar3,star3,aint0,aintb)
!	Calculation of the boundary condition of intensity behind 
!	the last untransparent object along the line of sight
!	input:  
!         iunt-last untransparent point
!	  dens,temp -density, temperature
!	  vz - radial velocity field
!         vdif -differential radial velocity between the reflecting 
!               surface and the source of light
!	  freq - frequency at which we are about to solve RTE
!	  dcut1,dcut2,dcut3,dcutn -boundaries of density intervals
!	  darkl-limb darkening factor of the particular ray
!         cdelta -cos(irradiation angle)*R**2/separation**2
!	  lunt1,nstar1,xstar1,star1,nspec1,tspec1 -refere to 
!		'star' i.e. <dcut1,dcut2>
!	  tstar -polar temperature of the 'star'
!	  dlst - limb darkening coefficients of the `star`
!	  irrst -irradiation switch
!	  albst -Bond albedo of the `star`
!	  alb1x,alb1y,nalb1 -monochromatic albedo of the 'star'
!	  wstar1,fstar1,jstar1 -spectrum of the star with temp=tstar
!	  lunt2,nstar2,xstar2,star2,nspec2,tspec2 -refere to 
! 		'companion' i.e. <dcut2,dcut3>
!	  tempcp -polar temperature of the 'companion'
!         dlcp - limb darkening coefficients of the 'companion'
!	  irrcp -irradiation switch
!         albcp -Bond albedo of the `companion`
!	  alb2x,alb2y,nalb2 -monochromatic albedo of the 'companion'
!         wstar2,fstar2,jstar2 -spectrum of `companion` with temp=tempcp
!	  lunt3,nstar3,xstar3,star3 -refere to <dcut3,dcutn>
!	  aint0 -incident intensity from behind the shell
!	output:
!	  aintb- boundary condition for intensity, intensity heading 
!	    towards you along the line of sight after the untransparent 
!	    object
	implicit none
        integer, parameter:: i4=4,i8=8
        include 'param.inc'
	real(i8):: temp(ndim),dens(ndim),vz(ndim)
        real(i8):: xstar1(mstarx,mspecx),star1(mstarx,mspecx)
        real(i8):: xstar2(mstarx,mspecx),star2(mstarx,mspecx)
        real(i8):: xstar3(mstarx),star3(mstarx)
        integer:: nstar1(mspecx),nstar2(mspecx)
	real(i8):: tspec1(mspecx),yy(mspecx),tspec2(mspecx)
        real(i8):: wstar1(mstarx),fstar1(mstarx)
        real(i8):: wstar2(mstarx),fstar2(mstarx)
        real(i8):: alb1x(mstarx),alb1y(mstarx)
        real(i8):: alb2x(mstarx),alb2y(mstarx)
	real(i8):: vdif,freq,dcut1,dcut2,dcut3,dcutn,darkl,cdelta
	real(i8):: tstar,dlst,dlst2,albst,tempcp,dlcp,dlcp2,albcp
	real(i8):: aint0,aintb
	integer:: iunt,lunt1,nspec1,nalb1,irrst,jstar1,lunt2,nspec2
	integer:: nalb2,irrcp,jstar2,lunt3,nstar3
	real(i8):: clight,pi,freqs,freqr,alamf,der,bunt,albnu,buntr
	real(i8):: alamr
	real(i8):: planck
	integer:: j,jj
	clight=2.99792458d10
	pi=3.1415926535897931d0
	aintb=aint0
!       Four types of untransparent objects can be recognized here:
!       dcut1<dens<dcut2 -central star,
!       dcut2<dens<dcut3 -secondary star(=companion) 
!       dcut3<dens<dcutn -3.body (it can be anything)
!       dcutn<dens -any opaque dark matter
!	Note that lunt1, lunt2, lunt3 are in fact associated with 
!	particular density intervals rather then with objects. 
!       They can be used to ascribe the spectrum to any untransparent 
!       object setting its density within a particular density interval.
!	Be cautious to turn off the limb darkening in that case.
	if(iunt.gt.0)then
  	  freqs=freq*(1.d0-vz(iunt)/clight)
	  freqr=freq*(1.d0-(vdif+vz(iunt))/clight)
!	  behind the primary
	  if(dens(iunt).gt.dcut1.and.dens(iunt).le.dcut2)then
	    if(lunt1.gt.0)then
	      alamf=clight/freqs*1.d8
	      do j=1,nspec1
	        call locate2(xstar1,nstar1(j),alamf,jj,j)
	        if(jj.lt.1)then
	          yy(j)=star1(1,j)
	        elseif(jj.ge.nstar1(j))then
	          yy(j)=star1(nstar1(j),j)
	        else
	          der=star1(jj+1,j)-star1(jj,j)
                  der=der/(xstar1(jj+1,j)-xstar1(jj,j))
	          yy(j)=(alamf-xstar1(jj,j))*der+star1(jj,j)
	        endif
              enddo
              call locate(tspec1,nspec1,temp(iunt),jj)
	      if(jj.lt.1)then
                bunt=yy(1)*darkl
              elseif(jj.ge.nspec1)then 
                bunt=yy(nspec1)*darkl
              else
                der=(yy(jj+1)-yy(jj))/(tspec1(jj+1)-tspec1(jj))
                bunt=(temp(iunt)-tspec1(jj))*der+yy(jj)
                bunt=bunt*darkl                
              endif
            else
              bunt=planck(freqs,temp(iunt),1)*darkl
              bunt=bunt/(1.d0-dlst/3.d0-dlst2/6.d0)
	    endif
!		reflection off the primary
	    if(irrst.eq.1.and.albst.gt.0.d0.and.cdelta.gt.0.d0)then
	      if(nalb1.gt.1)then
                alamf=clight/freqs*1.d8
!                write(*,*)'before locate alb1x',(alb1x(i),i=1,nalb1)
	        call locate(alb1x,nalb1,alamf,jj)
	        if(jj.lt.1)then
	          albnu=alb1y(1)
	        elseif(jj.ge.nalb1)then
	          albnu=alb1y(nalb1)
	        else
	          der=alb1y(jj+1)-alb1y(jj)
	          der=der/(alb1x(jj+1)-alb1x(jj))
	          albnu=(alamf-alb1x(jj))*der+alb1y(jj)
	        endif
	      else
	        albnu=albst
	      endif
	    endif
            if(lunt2.gt.0)then
              if(irrst.eq.1.and.albst.gt.0.d0.and.cdelta.gt.0.d0)then
	        alamr=clight/freqr*1.d8              
                call locate(wstar2,jstar2,alamr,jj)
                if(jj.lt.1)then
                  buntr=fstar2(1)
                elseif(jj.ge.jstar2)then
                  buntr=fstar2(jstar2)
                else
                  der=fstar2(jj+1)-fstar2(jj)
                  der=der/(wstar2(jj+1)-wstar2(jj))
                  buntr=(alamr-wstar2(jj))*der+fstar2(jj)
                endif
                buntr=albnu*cdelta*buntr*(1.d0-dlcp/3.d0-dlcp2/6.d0)
              else
                buntr=0.d0
              endif  
            else
              if(irrst.eq.1.and.albst.gt.0.d0.and.cdelta.gt.0.d0)then
                buntr=albnu*cdelta*planck(freqr,tempcp,1)
              else
                buntr=0.d0
              endif  
            endif
            bunt=bunt+buntr
!         behind the secondary
	  elseif(dens(iunt).gt.dcut2.and.dens(iunt).le.dcut3)then
	    if(lunt2.gt.0)then
	      alamf=clight/freqs*1.d8
! 	      darkg=planck(freqs,temp(iunt),1)/planck(freqs,tempcp,1)
              do j=1,nspec2
  	        call locate2(xstar2,nstar2(j),alamf,jj,j)
	        if(jj.lt.1)then
	          yy(j)=star2(1,j)
	        elseif(jj.ge.nstar2(j))then
	          yy(j)=star2(nstar2(j),j)
	        else
	          der=star2(jj+1,j)-star2(jj,j)
                  der=der/(xstar2(jj+1,j)-xstar2(jj,j))
	          yy(j)=(alamf-xstar2(jj,j))*der+star2(jj,j)
	        endif
              enddo
              call locate(tspec2,nspec2,temp(iunt),jj)
              if(jj.lt.1)then
!             if(temp(iunt).lt.tspec2(1))then              
                bunt=yy(1)*darkl
              elseif(jj.ge.nspec2)then
!             elseif(temp(iunt).ge.tspec2(nspec2))then              
                bunt=yy(nspec2)*darkl
              else
                der=(yy(jj+1)-yy(jj))/(tspec2(jj+1)-tspec2(jj))
                bunt=(temp(iunt)-tspec2(jj))*der+yy(jj)
                bunt=bunt*darkl
              endif              
	    else
              bunt=planck(freqs,temp(iunt),1)*darkl
              bunt=bunt/(1.d0-dlcp/3.d0-dlcp2/6.d0)
	    endif
!		reflection off the secoondary
	    if(irrcp.eq.1.and.albcp.gt.0.d0.and.cdelta.gt.0.d0)then
	      if(nalb2.gt.1)then
                alamf=clight/freqs*1.d8
	        call locate(alb2x,nalb2,alamf,jj)
	        if(jj.lt.1)then
	          albnu=alb2y(1)
	        elseif(jj.ge.nalb2)then
	          albnu=alb2y(nalb2)
	        else
	          der=alb2y(jj+1)-alb2y(jj)
	          der=der/(alb2x(jj+1)-alb2x(jj))
	          albnu=(alamf-alb2x(jj))*der+alb2y(jj)
	        endif
	      else
	        albnu=albcp
	      endif
	    endif
            if(lunt1.gt.0)then
              if(irrcp.eq.1.and.albcp.gt.0.d0.and.cdelta.gt.0.d0)then
	        alamr=clight/freqr*1.d8              
                call locate(wstar1,jstar1,alamr,jj)
                if(jj.lt.1)then
                  buntr=fstar1(1)
                elseif(jj.ge.jstar1)then
                  buntr=fstar1(jstar1)
                else
                  der=fstar1(jj+1)-fstar1(jj)
                  der=der/(wstar1(jj+1)-wstar1(jj))
                  buntr=(alamr-wstar1(jj))*der+fstar1(jj)
                endif
                buntr=albnu*cdelta*buntr*(1.d0-dlst/3.d0-dlst2/6.d0)
              else
                buntr=0.d0
              endif  
            else
              if(irrcp.eq.1.and.albcp.gt.0.d0.and.cdelta.gt.0.d0)then
                buntr=albnu*cdelta*planck(freqr,tstar,1)
              else
                buntr=0.d0
              endif  
            endif
            bunt=bunt+buntr
!         behind the third body
	  elseif(dens(iunt).gt.dcut3.and.dens(iunt).le.dcutn)then
	    if(lunt3.gt.0)then
	      alamf=clight/freqs*1.d8
	      call locate(xstar3,nstar3,alamf,jj)
	      if(jj.lt.1)then
	        bunt=star3(1)
	      elseif(jj.ge.nstar3)then
	        bunt=star3(nstar3)
	      else
	        der=(star3(jj+1)-star3(jj))/(xstar3(jj+1)-xstar3(jj))
	        bunt=(alamf-xstar3(jj))*der+star3(jj)
	      endif
	    else
              bunt=planck(freqs,temp(iunt),1)
	    endif
	  else
            bunt=0.d0
	  endif	
          aintb=bunt
	endif
	return
	end
!-----------------------------------------------------------------------
        subroutine jnu(rstar,tstar,ari,ati,az                           &
     &  ,vr,vt,vz,freq,xstar1,star1,nstar1,lunt1                        &
     &  ,vxstr,vystr,vzstr,jj,ejnu,iprint,xcpr,ycpr,zcpr)
!	calculation of the mean intensity
!	from a spherical object in the opt. thin medium
!	input:
!	  rstar,tstar -radius, temperature of the object
!         lunt1 -black body switch
!	  ari,ati,az -coordinates of the point
!	  vr,vt,vz -velocity at the point
!	  freq -frequency
!	  xstar1,star1,nstar1 -spectrum of the object
!	  vxstr,vystr,vzstr -velocity of the object (source)
!	  iprint -dummy
!	  xcpr,ycpr,zcpr - coordinates of the object (source)
!	output:
!		ejnu = mean intensity = I_nu*omega/4/pi
!		jj- input for subseqent call of hunt
	implicit none
        integer, parameter:: i4=4,i8=8
        include 'param.inc'
	real(i8):: xstar1(mstarx),star1(mstarx)
	real(i8):: rstar,tstar,ari,ati,az,vr,vt,vz,freq
	real(i8):: vxstr,vystr,vzstr,ejnu,xcpr,ycpr,zcpr 
	integer:: nstar1,lunt1,jj,iprint
	real(i8):: clight,dx,dy,dz,er2,ratio2,omega4,rdotv,dvel,freqs
	real(i8):: alamf,bstar,der
	real(i8):: planck
	clight=2.99792458d10
	dx=ari-xcpr
	dy=ati-ycpr
	dz=az-zcpr
	er2=dx**2+dy**2+dz**2
	ratio2=rstar*rstar/er2
	if(ratio2.lt.1.d0)then	
	  omega4=0.5d0*(1.d0-dsqrt(1.d0-ratio2))
	else
          omega4=0.d0
	endif
	rdotv=(dx*(vr-vxstr)+dy*(vt-vystr)+dz*(vz-vzstr))/dsqrt(er2)
	dvel=-rdotv+vz
!		freqs -frequency of alam in comoving frame
	  freqs=-dvel/clight*freq+freq
!		in ejnu the rotation of the star is ignored so far
	  if(lunt1.gt.0)then
	    alamf=clight/freqs*1.d8	
            call hunt(xstar1,nstar1,alamf,jj)
	    if(jj.lt.1)then
	      bstar=star1(1)
	    elseif(jj.ge.nstar1)then
	      bstar=star1(nstar1)
	    else
	      der=(star1(jj+1)-star1(jj))/(xstar1(jj+1)-xstar1(jj))
	      bstar=(alamf-xstar1(jj))*der+star1(jj)
	    endif
	  else
	    bstar=planck(freqs,tstar,1)
	  endif
          ejnu=bstar*omega4
	  return
	  end
!-----------------------------------------------------------------------
        subroutine roche(rosu,rotem,areas,teff                          &
     &  ,rfront,rback,rpole,rside,rmean                                 &
     &  ,tstar,rstar,dgst,irrst,albst,htst,htsta                        &
     &  ,tempcp,rcp,dgcp,irrcp,albcp,htcp,htcpa                         &
     &  ,x,y,xcomp,q,ff,nbod1,nbod2,ipsc)
!	Subroutine calculates the Roche surface (RS), normalized 
!       gravity on RS of the primary, secondary or contact system,
!	irradiation angle, and temperature on RS including gravity 
!	darkening and an approximate irradiation effect with
!	reflection, heating, and heat redistribution for detached 
!	components.
!	Code works in normalized scale where center of 
!	mass of the primary is at x=0, center of mass of the secondary 
!	at x=1 but the input/output is scaled/unscaled. 
!	CGS units are employed.
!		input: 
!	x,y -grid points in orbital plane
!	nbod1,nbod2 -number of grid points in x and y
!	xcomp -is the scale=x coordinate of the secondary 
!	q -mass ratio (q<1 -more massive star is in the center,
!		q>1 -less massive star is in the center)
!	ff -Roche lobe filling factor of the primary (if ipsc=1). It is
!	   the distance of the inner substellar point of the primary 
!	   (between the stars) from the center of the primary relative 
!	   to the distance to L1, 
!	   ff<=1, the Roche lobe is reproduced if ff=1 
!	ff -Roche lobe filling factor of the secondary (if ipsc=2).
!	   It is the radius of the secondary at the substellar point 
!	   relative to 1-L1, 
!	   ff<=1, the Roche lobe is reproduced if ff=1 
!	ff -Roche lobe fill-out factor of the contact system (if ipsc=3)
!		ff=(P1-P)/(P1-P2)+1   where P,P1,P2 are potentials
!	ipsc=1 -calculates RS of a detached central star (primary)
!	ipsc=2 -calculates RS of a detached secondary
!	ipsc=3 -calculates RS of a contact system
!       tstar -effective temperature of the central star in [K]      
!          in the absence of irradiation effect
!          if ipsc=1 it is the temperature at the rotation pole 
!          if ipsc=3 it is the temperature at the rotation pole of
!          the more massive star
!       rstar -radius of the central star used for the irradiation 
!		of other object
!       dgst -gravity darkening coefficient (beta) of the central star
!               (0.25 for radiative, 0.08 for convective atmospheres) 
!       irrst=0 -irradiation (reflection effect) is off
!           (albstp,htst,htsta have no meaning in this case)
!       irrst=1 -irradiation of the object from the companion is on.
!       albst  -Bond albedo <0,1> of the central star
!       htst   -heat redistribution parameter in case of the irradiation, 
!           <0,1>, 0-nothing is redistributed over the day side and 
!          nothing goes to the night, 1-all the heat is evenly 
!          redistributed over the day and night sides 
!       htsta  -zonal temperature redistribution parameter i.e.
!	   a degree of the homegenity of the heat transport  
!           <0,1>. 1-homegeneous, 0-cosine dependence (zonal)
!           T**4=T0**4(htsta+(1-htsta)*cos_latitude*constant)
!       tempcp -effective temperature at the rotation pole of 
!	   the secondary in [K] in the absence of irradiation effect
!       rcp -radius of the secondary used for the irradiation
!               of other object
!	irrcp,albcp,htcp,htcpa -the same as in the case of primary
!
!		output: 
!	rosu(x,y) -Roche surface
!	rotem(x,y) -surface temperature
!       areas -area of the RS
!	teff -effective temperature of the object including irradiation
!		and gravity darkening
! 	rfront,rback,rpole,rside,rmean -radius relative to the semimajor
!		axis at the substellar point, antistellar point, 
!		rotation pole, side point, and radius of a sphere with
!		 the same volume
!	rogrv(x,y) -normalized gravity on the RS relative to the gravity
!		on the rotation pole (of more massive star in case of
!		contact system)
!	rocos(x,y) -cosine of the angle between gravity vector
!		(-1*normal to the surface) and direction towards
!		towards the irradiating star i.e. cosine of the zenit 
!		distance of the irradiating star as seen from 
!		the irradiated surface). It is used for 
!		an approximate reflection effect (only for ipsc=1,2). 
!		rocos=0 on the night side.
!
!       Author/contact: Jan Budaj, http://www.ta3.sk/~budaj
!                       budaj@ta3.sk
!
	implicit none
	integer, parameter:: i4=4,i8=8
 	include 'param.inc'
!	parameter(ndimf1=6001,ndimf2=6001)
	real(i8):: x(ndimf1),y(ndimf2),xn(ndimf1),yn(ndimf2)
        real(i8):: rosuy(ndimf1)
	real(i8):: rosu(ndimf1,ndimf2),rogrv(ndimf1,ndimf2)
        real(i8):: rocos(ndimf1,ndimf2),rotem(ndimf1,ndimf2)
        real(i8):: area(ndimf1,ndimf2),areaz(ndimf1,ndimf2)
        real(i8):: arear(ndimf1,ndimf2)
	real(i8):: pot,dxpot,dypot,dzpot,dxxpot
	external pot,dxpot,dypot,dzpot,dxxpot
	real(i8):: raphx,raphy,raphz,radlob
	real(i8):: areas,teff,rfront,rback,rpole,rside,rmean
	real(i8):: tstar,rstar,dgst,albst,htst,htsta
	real(i8):: tempcp,rcp,dgcp,albcp,htcp,htcpa
	real(i8):: xcomp,q,ff
	real(i8):: pi,stefb,prc,prc2,el1,el2,el3
	real(i8):: pinn,pout,rpot,ystart,xpole,ypole,zpole
	real(i8):: potel1,potel2,potel3,pom,sout,sinn,rvec,vol,sum1,sum2
	real(i8):: fxy2,coslat,surf2,rr2,firr,surf1,t04,t04new
	real(i8):: htb,htbnew,fnlat,tirr4,potx,poty,potz,gravec,gravnp
	integer:: irrst,irrcp,nbod1,nbod2,ipsc,i,j
	pi=3.1415926535897931d0
	stefb=5.670400d-5
        if((ipsc.eq.1.or.ipsc.eq.2).and.(ff.gt.1.d0.or.ff.lt.1.d-2))then
	  write(*,130)
	  write(3,131)
	  call exit(1)
!	  goto 120
	endif
	if(ipsc.eq.3.and.(ff.le.1.d0.or.ff.gt.2.d0))then
	  write(*,130)
	  write(3,131)
	  call exit(1)
!	  goto 120
	endif
	if(.not.(ipsc.eq.1.or.ipsc.eq.2.or.ipsc.eq.3))then
	  write(*,130)
	  write(3,131)
	  call exit(1)
!	  goto 120
	endif
        if(htst.lt.0.d0.or.htst.gt.1.d0.or.                             &
     &    htsta.lt.0.d0.or.htsta.gt.1.d0.or.                            &
     &    htcp.lt.0.d0.or.htcp.gt.1.d0.or.                              &
     &    htcpa.lt.0.d0.or.htcpa.gt.1.d0)then
	  write(*,130)
	  write(3,131)
	  call exit(1)
!	  goto 120
	endif
!	precision required
	prc=1.d-7
	prc2=1.d-9
!	normalization of coordinates
        do 10 i=1,nbod1
          xn(i)=x(i)/xcomp
10      continue
        do 20 i=1,nbod2
          yn(i)=y(i)/xcomp
20      continue
!		el1- the x coordinate of the L1 point 
!	it may be quicker to start with 0.5 then with radlob function
	el1=raphx(dxpot,dxxpot,0.d0,q,0.5d0,prc2)
!       	el2,el3 -the L2,L3 points (they are reversed if q>1)
        el2=raphx(dxpot,dxxpot,0.d0,q,2.d0-el1,prc2)
        el3=raphx(dxpot,dxxpot,0.d0,q,-el1,prc2)
!
!			primary, ff<=1
	if(ipsc.eq.1.and.ff.gt.0.d0.and.ff.le.1.d0)then
	  pinn=ff*el1
	  pout=-pinn
	  rpot=pot(pinn,0.d0,0.d0,q)
!	  pinn-the inner substellar point of the surface of the primary
!	  pout-the outer point of the primary on the averted hemisphere
	  if(pout.gt.el3)then
   	    pout=pout
	  else
	    pout=el3*0.5d0
	  endif
	  pout=raphx(pot,dxpot,rpot,q,pout,prc2)
!	  rosuy(x)-the shape of the primary Roche surface for x=x(i),z=0
!	  ystart=el1*0.2d0
!	  expression below might help for small ff
	  ystart=pinn*0.5d0
	  do 50 i=1,nbod1
	    if(xn(i).gt.pout.and.xn(i).lt.pinn)then
	      rosuy(i)=dabs(raphy(pot,dypot,rpot,q,xn(i),ystart,prc))
	      ystart=rosuy(i)
	    else
	      rosuy(i)=0.d0
	    endif
50	  continue
!	  coordinates of the rotation pole
	  xpole=0.d0
	  ypole=0.d0
!	  zpole=raphz(pot,dzpot,rpot,q,xpole,ypole,radlob(1.d0/q),prc2)
	  zpole=raphz(pot,dzpot,rpot,q,xpole,ypole,dabs(pout),prc2)
	  rfront=pinn
	  rback=-pout
	  rpole=zpole
	  rside=dabs(raphy(pot,dypot,rpot,q,0.d0,pout,prc2))
	endif
!
!			contact, 1<ff<=2
	if(ipsc.eq.3.and.ff.gt.1.d0.and.ff.le.2.d0)then
	  potel1=pot(el1,0.d0,0.d0,q)	
	  potel2=pot(el2,0.d0,0.d0,q)
	  potel3=pot(el3,0.d0,0.d0,q)
	  if(potel2.lt.potel3)then
	    pom=potel2
	    potel2=potel3
	    potel3=pom
	  endif
	  rpot=potel1-(ff-1.d0)*(potel1-potel2)
	  pout=raphx(pot,dxpot,rpot,q,el3*0.5d0,prc2)
	  pinn=raphx(pot,dxpot,rpot,q,0.5d0+el2*0.5d0,prc2)
!	  pinn-the right hand side edge of the contact Roche surface
!	  pout-the left hand side edge of the contact Roche surface
!	  rosuy(x)-the shape of the contact Roche surface for x=x(i),z=0
	  ystart=radlob(1.d0/q)*0.2d0
	  do 70 i=1,nbod1
	    if(xn(i).gt.pout.and.xn(i).lt.pinn)then
	      rosuy(i)=dabs(raphy(pot,dypot,rpot,q,xn(i),ystart,prc))
	      ystart=rosuy(i)
	    else
	      rosuy(i)=0.d0
	    endif
70	  continue
!	  coordinates of the rotation pole of the more massive star
	  if(q.lt.1.d0)then
	  xpole=0.d0
	  ypole=0.d0
	  zpole=raphz(pot,dzpot,rpot,q,xpole,ypole,radlob(1.d0/q),prc2)
          rside=dabs(raphy(pot,dypot,rpot,q,xpole,radlob(1.d0/q),prc2))
	  else
	  xpole=1.d0
	  ypole=0.d0
	  zpole=raphz(pot,dzpot,rpot,q,xpole,ypole,radlob(q),prc2)
	  rside=dabs(raphy(pot,dypot,rpot,q,xpole,radlob(q),prc2))
	  endif
	  rfront=pinn
	  rback=-pout
	  rpole=zpole
	endif
!
!			secondary, ff<=1
	if(ipsc.eq.2.and.ff.gt.0.d0.and.ff.le.1.d0)then
	  sout=ff*(1.d0-el1)
	  sinn=1.d0-sout
	  rpot=pot(sinn,0.d0,0.d0,q)
!	  sinn-inner substellar point of the surface of the secondary
!	  sout-outer point of the secondary on the averted hemisphere
	  if((1.d0+sout).lt.el2)then
   	    sout=1.d0+sout
	  else
	    sout=1.d0+(el2-1.d0)*0.5d0
	  endif
	  sout=raphx(pot,dxpot,rpot,q,sout,prc2)
!	  rosuy(x)-the shape of the Roche surface for x=x(i),z=0
	  ystart=(sout-1.d0)*0.2d0
	  do 90 i=1,nbod1
	    if(xn(i).gt.sinn.and.xn(i).lt.sout)then
	      rosuy(i)=dabs(raphy(pot,dypot,rpot,q,xn(i),ystart,prc))
	      ystart=rosuy(i)
	    else
	      rosuy(i)=0.d0
	    endif
90	  continue
!	  coordinates of the rotation pole
	  xpole=1.d0
	  ypole=0.d0
	  zpole=radlob(q)*ff
	  zpole=raphz(pot,dzpot,rpot,q,xpole,ypole,zpole,prc2)
	  rfront=1.d0-sinn
	  rback=sout-1.d0
	  rpole=zpole
	  rside=dabs(raphy(pot,dypot,rpot,q,1.d0,zpole,prc2))
	endif
!	check mainly this
  	do 95 i=1,nbod1
          if(rosuy(i).gt.0.d0.and.                                      &
     &    dypot(xn(i),rosuy(i),0.d0,q).ge.0.d0)then
            write(*,*)' error: Newton-Raphson slipped in subr.: Roche,',&
     &      ' try to modify slightly: nbod,q,ff'
            write(3,*)' 1 error: Newton-Raphson slipped in subr. Roche'
            call exit(1)
!           goto 120
	  endif
95	continue
!
!	gravnp-normalized gravity at the rotation pole
	call gravn(xpole,ypole,zpole,q,potx,poty,potz,gravnp)
!
!		rosu(x,y)-the shape of the Roche surface
!		rogrv(x,y) -normalized gravity relative to the pole
!		rocos(i,j) -cosine of the irradiation angle
	do  i=1,nbod1
          do  j=1,nbod2
            if(dabs(yn(j)).lt.rosuy(i))then
              rosu(i,j)=raphz(pot,dzpot,rpot,q,xn(i),yn(j),rosuy(i),prc)
	      call gravn(xn(i),yn(j),rosu(i,j),q,potx,poty,potz,gravec)
	      rogrv(i,j)=gravec/gravnp
	      areaz(i,j)=dabs(potz/gravec)
	      if(ipsc.eq.1)then
	        rvec=dsqrt((xn(i)-1.d0)**2+yn(j)**2+rosu(i,j)**2)
                rocos(i,j)=((xn(i)-1.d0)*potx+yn(j)*poty+rosu(i,j)*potz)
                rocos(i,j)=rocos(i,j)/rvec/gravec
       	        if(rogrv(i,j).lt.1.d-4)then
	          rogrv(i,j)=1.d-4
	          rocos(i,j)=0.84d0
	        endif
		rotem(i,j)=tstar*rogrv(i,j)**dgst
	      elseif(ipsc.eq.2)then
	        rvec=dsqrt(xn(i)*xn(i)+yn(j)*yn(j)+rosu(i,j)*rosu(i,j))
	        rocos(i,j)=(xn(i)*potx+yn(j)*poty+rosu(i,j)*potz)
	        rocos(i,j)=rocos(i,j)/rvec/gravec
!		take care of the unimportant case of singularity
     	        if(rogrv(i,j).lt.1.d-4)then
	          rogrv(i,j)=1.d-4
	          rocos(i,j)=0.84d0
	        endif
		rotem(i,j)=tempcp*rogrv(i,j)**dgcp
              else
	        if(rogrv(i,j).lt.1.d-4)rogrv(i,j)=1.d-4
                rocos(i,j)=0.d0
		rotem(i,j)=tstar*rogrv(i,j)**dgst
	      endif
    	      if(rocos(i,j).lt.0.d0)rocos(i,j)=0.d0
            else
	      rosu(i,j)=0.d0
	      rogrv(i,j)=0.d0
	      rocos(i,j)=0.d0
	      rotem(i,j)=0.d0
	      areaz(i,j)=0.d0
            endif
!		unscale the result
	    rosu(i,j)=xcomp*rosu(i,j)
	  enddo
	enddo
!	
!		irradiation and heat redistribution
	call surf(ndim1,ndim2,nbod1,nbod2,x,y,area)
!	area-projection of the surface element to xy
!	arear-area of the surface element
!	vol-volume under the Roche surface
        areas=0.d0
	vol=0.d0
        sum1=0.d0
        sum2=0.d0        
        do i=1,nbod1  
          do j=1,nbod2
            if(dabs(yn(j)).lt.rosuy(i))then  
              arear(i,j)=area(i,j)/areaz(i,j)
              if(irrst.eq.1.and.ipsc.eq.1)then
                fxy2=x(i)**2+y(j)**2
                coslat=fxy2/(fxy2+rosu(i,j)**2)
                coslat=dsqrt(coslat)
                surf2=coslat*arear(i,j)
                if(rocos(i,j).gt.0.d0)then
!                 dayside
                  rr2=rcp**2/((x(i)-xcomp)**2+y(j)**2+rosu(i,j)**2)
                  firr=stefb*rr2*tempcp**4*rocos(i,j)
		  surf1=firr*arear(i,j)
                endif
              endif
              if(irrcp.eq.1.and.ipsc.eq.2)then
                fxy2=(x(i)-xcomp)**2+y(j)**2
                coslat=fxy2/(fxy2+rosu(i,j)**2)
                coslat=dsqrt(coslat)
                surf2=coslat*arear(i,j)
                if(rocos(i,j).gt.0.d0)then
!                 day-side
                  rr2=rstar*rstar/(x(i)**2+y(j)**2+rosu(i,j)**2)
                  firr=stefb*rr2*tstar**4*rocos(i,j)
                  surf1=firr*arear(i,j)
                endif    
              endif              
              areas=areas+arear(i,j)
              vol=vol+area(i,j)*rosu(i,j)
	      sum1=sum1+surf1
	      sum2=sum2+surf2
            endif
          enddo
        enddo
        areas=areas*2.d0
        rmean=(vol*2.d0/(4.d0/3.d0*pi))**(1.d0/3.d0)/xcomp
        sum1=sum1*2.d0
        sum2=sum2*2.d0
        if(irrst.eq.1.and.ipsc.eq.1)then
          rr2=rcp**2/xcomp**2
          t04=0.25d0*htst*(1.d0-albst)*rr2*tempcp**4
          t04new=htst*(1.d0-albst)*sum1/areas/stefb
!            t04=htst*(1.d0-albst)*sum1/areas/stefb          
          htb=4.d0*(1.d0-htsta)/pi
          htbnew=(1.d0-htsta)*areas/sum2
!            htb=(1.d0-htsta)*areas/sum2          
        endif  
        if(irrcp.eq.1.and.ipsc.eq.2)then
          rr2=rstar*rstar/xcomp**2  
          t04=0.25d0*htcp*(1.d0-albcp)*rr2*tstar**4            
          t04new=htcp*(1.d0-albcp)*sum1/areas/stefb
!            t04=htcp*(1.d0-albcp)*sum1/areas/stefb          
          htb=4.d0*(1.d0-htcpa)/pi
          htbnew=(1.d0-htcpa)*areas/sum2
!            htb=(1.d0-htcpa)*areas/sum2          
        endif  
!	sum1,t04new may not be precise enough
!	(polar or cylindrical coordinates might be better).
!	However, if you want to use them simply uncomment 
!	the four lines above and replace
!	t04 by t04new and htb by htbnew.
!	write(*,*)'T0=',t04new**0.25d0,t04**0.25d0
!	write(*,*)'b=',htbnew,htb
!	write(*,*)'sum1,2=',sum1,sum2
	teff=0.d0
	do i=1,nbod1
          do j=1,nbod2
            if(dabs(yn(j)).lt.rosuy(i))then
              if(irrst.eq.1.and.ipsc.eq.1)then
                fxy2=x(i)**2+y(j)**2
                coslat=fxy2/(fxy2+rosu(i,j)**2)
                coslat=dsqrt(coslat)
                fnlat=htsta+htb*coslat
                if(rocos(i,j).gt.0.d0)then
!                 dayside
                  rr2=rcp**2/((x(i)-xcomp)**2+y(j)**2+rosu(i,j)**2)
                  firr=stefb*rr2*tempcp**4*rocos(i,j)
                  tirr4=(1.d0-htst)*(1.d0-albst)*firr/stefb
                  rotem(i,j)=(tirr4+t04*fnlat+rotem(i,j)**4)
		  rotem(i,j)=dsqrt(dsqrt(rotem(i,j)))
                else 
!                 nightside             
                  rotem(i,j)=(t04*fnlat+rotem(i,j)**4)
                  rotem(i,j)=dsqrt(dsqrt(rotem(i,j)))
                endif
              endif
              if(irrcp.eq.1.and.ipsc.eq.2)then
                fxy2=(x(i)-xcomp)**2+y(j)**2
                coslat=fxy2/(fxy2+rosu(i,j)**2)
                coslat=dsqrt(coslat)
                fnlat=htcpa+htb*coslat  
                if(rocos(i,j).gt.0.d0)then
!                 day-side
                  rr2=rstar*rstar/(x(i)**2+y(j)**2+rosu(i,j)**2)
                  firr=stefb*rr2*tstar**4*rocos(i,j)
                  tirr4=(1.d0-htcp)*(1.d0-albcp)*firr/stefb
                  rotem(i,j)=(tirr4+t04*fnlat+rotem(i,j)**4)
                  rotem(i,j)=dsqrt(dsqrt(rotem(i,j)))
                else
!                 night-side
                  rotem(i,j)=(t04*fnlat+rotem(i,j)**4)
		  rotem(i,j)=dsqrt(dsqrt(rotem(i,j)))
                endif    
              endif
              teff=teff+arear(i,j)*rotem(i,j)**4
            endif
	  enddo
	enddo
	teff=(teff*2.d0/areas)**0.25d0
!
!     	open(4,file='roche.out',status='unknown')
!	do i=1,nbod1
!         do j=1,nbod2
!           write(4,180)x(i),y(j),rosu(i,j),rogrv(i,j),rocos(i,j)       &
!     &      ,rotem(i,j)
!	  enddo
!	  write(4,*)
!	enddo
!180	format(6e14.5)
!	close(4)	
!
120	return
130	format(' error: inconsistent input in subroutine: roche')
131	format(' 1 error: inconsistent input in subroutine: roche')
	end
!-----------------------------------------------------------------------
        function raphx(fun,dfun,fun0,q,xstart,xacc)
!	Newton-Raphson Method in 1D (x coordinate)
!	to solve the equation fun(x)=fun0
!	input:  fun-external function of x and parameters: y,z,q
!		dfun-external function, partial derivative df/dx
!		fun0-value of the function 
!		q-mass ratio
!	        xstart-first estimate
!		xacc-accuracy
!	output:  raphx-root for which fun=fun0
        implicit none
	integer, parameter:: i4=4,i8=8
        real(i8):: raphx
        real(i8):: fun,dfun
	real(i8):: fun0,q,xstart,xacc,dx
	integer:: i
	raphx=xstart
	do i=1,22
	  dx=(fun(raphx,0.d0,0.d0,q)-fun0)/dfun(raphx,0.d0,0.d0,q)
	  raphx=raphx-dx
	  if(dabs(dx).lt.xacc)return
 	enddo
	return
	end
!-----------------------------------------------------------------------
        function raphy(fun,dfun,fun0,q,x0,ystart,yacc)
!	Newton-Raphson Method in 1D (y coordinate)
!	to solve the equation fun(y)=fun0
!	input:  fun-external function of y and parameters: x,z,q
!		dfun-external function, partial derivative df/dy
!		fun0-value of the function 
!		q-mass ratio
!	        ystart-first estimate
!		yacc-accuracy
!	output:  raphy-root for which fun=fun0
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: raphy
	real(i8):: fun,dfun
	real(i8):: fun0,q,x0,ystart,yacc,dy
	integer:: i
	raphy=ystart
	do i=1,22
	  dy=(fun(x0,raphy,0.d0,q)-fun0)/dfun(x0,raphy,0.d0,q)
	  raphy=raphy-dy
	  if(dabs(dy).lt.yacc)return
 	enddo
	return
	end function raphy
!-----------------------------------------------------------------------
        function raphz(fun,dfun,fun0,q,x0,y0,zstart,zacc)
!	Newton-Raphson Method in 1D (z coordinate)
!	to solve the equation fun(z)=fun0
!	input:  fun-external function of z and parameters: x,y,q
!		dfun-external function, partial derivative df/dz
!		fun0-value of the function 
!		q-mass ratio
!	        zstart-first estimate
!		zacc-accuracy
!	output:  raphz-root for which fun=fun0
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: raphz
	real(i8):: fun,dfun
	real(i8):: fun0,q,x0,y0,zstart,zacc,dz
	integer:: i		
	raphz=zstart
	do i=1,22
	  dz=(fun(x0,y0,raphz,q)-fun0)/dfun(x0,y0,raphz,q)
	  raphz=raphz-dz
	  if(dabs(dz).lt.zacc)return
 	enddo
	return
	end function raphz
!-----------------------------------------------------------------------
	function pot(x,y,z,q)
!	input: x,y,z,q
!	output: normalized Roche potential 
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: pot
	real(i8):: x,y,z,q,q1,x1,r1,r2
	q1=1.d0+q
	x1=x-1.d0
	r1=dsqrt(x*x+y*y+z*z)
	r2=dsqrt(x1*x1+y*y+z*z)
	pot=2.d0/q1/r1+2.d0*q/q1/r2+(x-q/q1)*(x-q/q1)+y*y
	return
	end function pot
!-----------------------------------------------------------------------
	function dxpot(x,y,z,q)
!	input: x,y,z,q
!	output: partial derivation of the Roche potential after x
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: dxpot
	real(i8):: x,y,z,q,q1,x1,r1,r2
	q1=1.d0+q
	x1=x-1.d0
	r1=dsqrt(x*x+y*y+z*z)
	r2=dsqrt(x1*x1+y*y+z*z)
	dxpot=-x/q1/r1**3-q*x1/q1/r2**3+x-q/q1
	dxpot=2.d0*dxpot
	return
	end function dxpot
!-----------------------------------------------------------------------
	function dypot(x,y,z,q)
!	input: x,y,z,q
!	output: partial derivation of the Roche potential after y
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: dypot
	real(i8):: x,y,z,q,q1,x1,r1,r2
	q1=1.d0+q
	x1=x-1.d0
	r1=dsqrt(x*x+y*y+z*z)
	r2=dsqrt(x1*x1+y*y+z*z)
	dypot=-2.d0/q1/r1**3-2.d0*q/q1/r2**3+2.d0
	dypot=dypot*y
	return
	end function dypot
!-----------------------------------------------------------------------
	function dzpot(x,y,z,q)
!	input: x,y,z,q
!	output: partial derivation of the Roche potential after z
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: dzpot
	real(i8):: x,y,z,q,q1,x1,r1,r2
	q1=1.d0+q
	x1=x-1.d0
	r1=dsqrt(x*x+y*y+z*z)
	r2=dsqrt(x1*x1+y*y+z*z)
	dzpot=-2.d0/q1/r1**3-2.d0*q/q1/r2**3
	dzpot=dzpot*z
	return
	end function dzpot
!-----------------------------------------------------------------------
	function dxxpot(x,y,z,q)
!	input: x,y,z,q
!	output: second partial derivation of the Roche potential after x
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: dxxpot
	real(i8):: x,y,z,q,q1,x1,r1,r2
	q1=1.d0+q
	x1=x-1.d0
	r1=dsqrt(x*x+y*y+z*z)
	r2=dsqrt(x1*x1+y*y+z*z)
	dxxpot=6.d0*x*x/q1/r1**5+6.d0*q*x1*x1/q1/r2**5-2.d0/q1/r1**3
        dxxpot=dxxpot-2.d0*q/q1/r2**3+2.d0
	return
	end function dxxpot
!-----------------------------------------------------------------------
	function radlob(q)
!	returns the volumee radius of the Roche lobe, 0<q<infty
!	q<1 -radius of the less massive star
!	q>1 -radius of the more massive star
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: radlob
	real(i8):: q,q13
	q13=q**(1.d0/3.d0)
	radlob=0.49d0*q13*q13/(0.69d0*q13*q13+dlog(1.d0+q13))
	return
	end function radlob
!-----------------------------------------------------------------------
        subroutine gravn(x,y,z,q,potx,poty,potz,pot)
!       calculates the normalized gravity vector
!       input: x,y,z,q
!       output: potx,poty,potz,pot
        implicit none
	integer, parameter:: i4=4,i8=8
	interface
	  function dxpot(x,y,z,q)
	    implicit none
            integer, parameter:: i4=4,i8=8
            real(i8):: dxpot
            real(i8):: x,y,z,q,q1,x1,r1,r2
	  endfunction dxpot
	  function dypot(x,y,z,q)
            implicit none
            integer, parameter:: i4=4,i8=8
            real(i8):: dypot
            real(i8):: x,y,z,q,q1,x1,r1,r2
          endfunction dypot
	  function dzpot(x,y,z,q)
            implicit none
            integer, parameter:: i4=4,i8=8
            real(i8):: dzpot
            real(i8):: x,y,z,q,q1,x1,r1,r2
          endfunction dzpot
	endinterface
	real(i8):: x,y,z,q,potx,poty,potz,pot
        potx=dxpot(x,y,z,q)
        poty=dypot(x,y,z,q)
        potz=dzpot(x,y,z,q)
        pot=dsqrt(potx*potx+poty*poty+potz*potz)
        return
        end
!-----------------------------------------------------------------------
        subroutine untsp(xstar1,star1,xstar2,star2,xstar3,star3         &
     &  ,wstar1,fstar1,tstar1,jstar1,wstar2,fstar2,tstar2,jstar2        &
     &  ,xunt1,yunt1,xunt2,yunt2,xunt3,yunt3                            &
     &  ,nstar1,nstar2,nstar3,nspec1,nspec2,tspec1,tspec2               &
     &  ,lunt1,lunt2,lunt3,dlst,dlst2,dlcp,dlcp2,alam,nfreq)
!	reading the input spectra of nontransparent objects
!	input: lunt1,lunt2,lunt3,xunt1,yunt1,xunt2,yunt2,xunt3,yunt3,
!		dlst,dlst2,dlcp,dlcp2,tstar1,tstar2,alam,nfreq
!	output:xstar1,star1,xstar2,star2,xstar3,star3,
!		nstar1,nstar2,nstar3,nspec1,nspec2,tspec1,tspec2,
!		wstar1,fstar1,jstar1,wstar2,fstar2,jstar2
        implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	character(50):: in1(mspecx),in2(mspecx)
	real(i8):: xstar1(mstarx,mspecx),star1(mstarx,mspecx)
        real(i8):: wstar1(mstarx),fstar1(mstarx)
        real(i8):: xstar2(mstarx,mspecx),star2(mstarx,mspecx)
        real(i8):: wstar2(mstarx),fstar2(mstarx) 
        real(i8):: xstar3(mstarx),star3(mstarx)
        integer:: nstar1(mspecx),nstar2(mspecx)
	real(i8):: tspec1(mspecx),tspec2(mspecx)     
        real(i8):: fjj(mstarx),wjj1(mstarx),fjj1(mstarx),fjj1n(mstarx)
        real(i8):: alam(mfreq)
	real(i8):: pi,clight,xunt1,yunt1,xunt2,yunt2,xunt3,yunt3
	real(i8):: dlst,dlst2,dlcp,dlcp2,tstar1,tstar2,dummy,der
	integer:: lunt1,lunt2,lunt3,nspec1,nspec2,nstar3,nfreq
	integer:: i,j,idummy,jstar1,jstar2,jj
	pi=3.1415926535897931d0
	clight=2.99792458d10
!	our units will be:	
! 	star1 -intrinsic nonrotated spectrum of the central star
!		intensity for angle=0 in [erg/cm^2/s/Hz/sterad]
!	xstar1(i) -lambda in [A]
	if(lunt1.gt.0)then
	  write(*,*)' reading file starspec1'
	  open(12,file='starspec1',STATUS='OLD')
	  read(12,*)
	  read(12,*)nspec1
	  read(12,*)
	  do i=1,nspec1
	    read(12,*)tspec1(i)
	  enddo
	  read(12,*)
	  do i=1,nspec1
	    read(12,'(a)')in1(i)
	  enddo
	  close(12)
  	  if(nspec1.gt.1)then
	    do i=1,nspec1-1
	      if(tspec1(i+1).le.tspec1(i))then 
	        write(*,*)' error: tspec1 in starspec1, stop'
	        write(3,*)' 1 error: tspec1 in starspec1, stop'
	        call exit(1)
	      endif  
            enddo
	  endif
	  do j=1,nspec1
	    open(12,file=in1(j),status='old')
	    i=1
            if(lunt1.eq.3)then
              read(12,*,end=20)
              read(12,*,end=20)
            endif
10	    if(lunt1.eq.1)then
	      read(12,*,end=20)xstar1(i,j),star1(i,j)
	      if(i.gt.1)then
!		fix synspec bug	      
	        if(xstar1(i,j).eq.xstar1(i-1,j))goto 10
	        if(.not.(xstar1(i,j).gt.xstar1(i-1,j)))then
                  write(*,*)' input error in starspec1 on line ',i,j
                  write(3,*)' 1 input error in starspec1 on line ',i,j
                endif
	        if(i+1.gt.mstarx)then
	          write(*,*)' error in starspec1: '
                  write(*,*)'dimension mstarx almost exceeded'
	          write(3,*)' 1 error in starspec1: '
                  write(3,*)'dimension mstarx almost exceeded'
                endif
	      endif
!		convert input of Hlambda from synspec to our units
	      xstar1(i,j)=xstar1(i,j)*xunt1
	      star1(i,j)=star1(i,j)*yunt1/(1.d0-dlst/3.d0-dlst2/6.d0)
	      star1(i,j)=4.d-8*star1(i,j)*xstar1(i,j)**2/clight
	    endif 
 	    if(lunt1.eq.2)then
 	      read(12,*,end=20)xstar1(i,j),star1(i,j)
	      if(i.gt.1)then
	        if(.not.(xstar1(i,j).gt.xstar1(i-1,j)))then
                  write(*,*)' input error in starspec1 on line ',i,j
                  write(3,*)' 1 input error in starspec1 on line ',i,j
                endif
	        if(i+1.gt.mstarx)then
	          write(*,*)' error in starspec1: '
	          write(*,*)'dimension mstarx almost exceeded'
	          write(3,*)' 1 error in starspec1: '
	          write(3,*)'dimension mstarx almost exceeded'
	        endif  
	      endif
	      xstar1(i,j)=xstar1(i,j)*xunt1
	      star1(i,j)=star1(i,j)*yunt1
	    endif
 	    if(lunt1.eq.3)then
 	      read(12,*,end=20)idummy,xstar1(i,j),dummy,star1(i,j)
	      if(i.gt.1)then
	        if(.not.(xstar1(i,j).gt.xstar1(i-1,j)))then
                  write(*,*)' input error in starspec1 on line ',i,j
                  write(3,*)' 1 input error in starspec1 on line ',i,j
                endif
	        if(i+1.gt.mstarx)then
	          write(*,*)' error in starspec1: '
                  write(*,*)'dimension mstarx almost exceeded'
	          write(3,*)' 1 error in starspec1: '
                  write(3,*)'dimension mstarx almost exceeded'
                endif  
	      endif
!		convert input of nu, Fnu from coolTlusty.21 to our units
	      xstar1(i,j)=xunt1*clight/xstar1(i,j)*1.d8
	      star1(i,j)=star1(i,j)*yunt1
	      star1(i,j)=star1(i,j)/pi/(1.d0-dlst/3.d0-dlst2/6.d0)
	    endif
	    i=i+1
	    goto 10
20	    nstar1(j)=i-1
	    write(6,*)' nstar1=',j,nstar1(j)
	    if(alam(1).lt.xstar1(1,j).or.                               &
	    alam(nfreq).gt.xstar1(nstar1(j),j))then
	      write(*,*)' error in starspec1: '
              write(*,*)'short range of lambdas'
	      write(3,*)' 1 error in starspec1: '
              write(3,*)'short range of lambdas'
	    endif	    
	    close(12)
	  enddo  
!		an interpolation to tstar1 to get wstar1,fstar1 
!		these will be used e.g. for scattering or reflection
	  if(tstar1.ge.tspec1(nspec1))then
	    do i=1,nstar1(nspec1)
	      wstar1(i)=xstar1(i,nspec1)
	      fstar1(i)=star1(i,nspec1)
	      jstar1=nstar1(nspec1)
	    enddo  
	  elseif(tstar1.le.tspec1(1))then
	    do i=1,nstar1(1)
	      wstar1(i)=xstar1(i,1)
	      fstar1(i)=star1(i,1)
	      jstar1=nstar1(1)
            enddo
	  else
	    call locate(tspec1,nspec1,tstar1,jj)
	    jstar1=nstar1(jj)
	    do i=1,nstar1(jj)
	      wstar1(i)=xstar1(i,jj)
	      fjj(i)=star1(i,jj)
	    enddo
	    do i=1,nstar1(jj+1)
	      wjj1(i)=xstar1(i,jj+1) 
              fjj1(i)=star1(i,jj+1) 
            enddo
	    call interp(wjj1,fjj1,wstar1,fjj1n,nstar1(jj+1),nstar1(jj))
	    do i=1,nstar1(jj)
	      der=(fjj1n(i)-fjj(i))/(tspec1(jj+1)-tspec1(jj))
	      fstar1(i)=der*(tstar1-tspec1(jj))+fjj(i)
	    enddo
	  endif
	endif
! 	star2 -intrinsic nonrotated spectrum of the secondary star
!		intensity for angle=0 in erg/cm^2/s/Hz/sterad
!	xstar2(i) -lambda in [A]
	if(lunt2.gt.0)then
	  write(*,*)' reading file starspec2'
	  open(13,file='starspec2',status='old')
	  read(13,*)
	  read(13,*)nspec2
	  read(13,*)
	  do i=1,nspec2
	    read(13,*)tspec2(i)
	  enddo
	  read(13,*)
	  do i=1,nspec2
	    read(13,'(a)')in2(i)
	  enddo
	  close(13)
  	  if(nspec2.gt.1)then
	    do i=1,nspec2-1
	      if(tspec2(i+1).le.tspec2(i))then 
	        write(*,*)' error: tspec2 in starspec2, stop'
	        write(3,*)' 1 error: tspec2 in starspec2, stop'
	        call exit(1)
	      endif  
            enddo
	  endif
          do j=1,nspec2
            open(13,file=in2(j),status='old')
  	    i=1
            if(lunt2.eq.3)then
              read(13,*,end=60)
              read(13,*,end=60)
            endif
50	    if(lunt2.eq.1)then
	      read(13,*,end=60)xstar2(i,j),star2(i,j)
	      if(i.gt.1)then
!		fix synspec bug
		if(xstar2(i,j).eq.xstar2(i-1,j))goto 50
	        if(.not.(xstar2(i,j).gt.xstar2(i-1,j)))then
                  write(*,*)' input error in starspec2 on line ',i,j
                  write(3,*)' 1 input error in starspec2 on line ',i,j
                endif  
	        if(i+1.gt.mstarx)then
	          write(*,*)' error in starspec2: '
                  write(*,*)'dimension mstarx almost exceeded'
	          write(3,*)' 1 error in starspec2: '
                  write(3,*)'dimension mstarx almost exceeded'
                endif  
	      endif
!		convert input of Hlambda from synspec to our units
  	      xstar2(i,j)=xstar2(i,j)*xunt2
	      star2(i,j)=star2(i,j)*yunt2/(1.d0-dlcp/3.d0-dlcp2/6.d0)
	      star2(i,j)=4.d-8*star2(i,j)*xstar2(i,j)**2/clight
	    endif
	    if(lunt2.eq.2)then
	      read(13,*,end=60)xstar2(i,j),star2(i,j)
	      if(i.gt.1)then
	        if(.not.(xstar2(i,j).gt.xstar2(i-1,j)))then
                  write(*,*)' input error in starspec2 on line ',i,j
                  write(3,*)' 1 input error in starspec2 on line ',i,j
                endif  
	        if(i+1.gt.mstarx)then
	          write(*,*)' error in starspec2: '
                  write(*,*)'dimension mstarx almost exceeded'
	          write(3,*)' 1 error in starspec2: '
                  write(3,*)'dimension mstarx almost exceeded'
                endif  
	      endif
	      xstar2(i,j)=xstar2(i,j)*xunt2
	      star2(i,j)=star2(i,j)*yunt2
	    endif
	    if(lunt2.eq.3)then
	      read(13,*,end=60)idummy,xstar2(i,j),dummy,star2(i,j)
	      if(i.gt.1)then
	        if(.not.(xstar2(i,j).gt.xstar2(i-1,j)))then
                  write(*,*)' input error in starspec2 on line ',i,j
                  write(3,*)' 1 input error in starspec2 on line ',i,j
                endif  
	        if(i+1.gt.mstarx)then
	          write(*,*)' error in starspec2: '
                  write(*,*)'dimension mstarx almost exceeded'
	          write(3,*)' 1 error in starspec2: '
                  write(3,*)'dimension mstarx almost exceeded'
                endif  
	      endif
!               convert input of nu, Fnu from coolTlusty.21 to our units	    
	      xstar2(i,j)=xunt2*clight/xstar2(i,j)*1.d8
	      star2(i,j)=star2(i,j)*yunt2
              star2(i,j)=star2(i,j)/pi/(1.d0-dlcp/3.d0-dlcp2/6.d0)
	    endif
	    i=i+1
	    goto 50
60	    nstar2(j)=i-1
	    write(6,*)' nstar2=',j,nstar2(j)
	    if(alam(1).lt.xstar2(1,j).or.                               &
	    alam(nfreq).gt.xstar2(nstar2(j),j))then
	      write(*,*)' error in starspec2: '
              write(*,*)'short range of lambdas'
	      write(3,*)' 1 error in starspec2: '
              write(3,*)'short range of lambdas'
	    endif	       
	    close(13)
          enddo
!		an interpolation to tstar2
	  if(tstar2.ge.tspec2(nspec2))then
	    do i=1,nstar2(nspec2)
	      wstar2(i)=xstar2(i,nspec2)
	      fstar2(i)=star2(i,nspec2)
	      jstar2=nstar2(nspec2)
	    enddo  
	  elseif(tstar2.le.tspec2(1))then
	    do i=1,nstar2(1)
	      wstar2(i)=xstar2(i,1)
	      fstar2(i)=star2(i,1)
	      jstar2=nstar2(1)
            enddo
	  else
	    call locate(tspec2,nspec2,tstar2,jj)
	    jstar2=nstar2(jj)
	    do i=1,nstar2(jj)
	      wstar2(i)=xstar2(i,jj)
	      fjj(i)=star2(i,jj)
	    enddo
	    do i=1,nstar2(jj+1)
	      wjj1(i)=xstar2(i,jj+1) 
              fjj1(i)=star2(i,jj+1) 
            enddo
	    call interp(wjj1,fjj1,wstar2,fjj1n,nstar2(jj+1),nstar2(jj))
	    do i=1,nstar2(jj)
	      der=(fjj1n(i)-fjj(i))/(tspec2(jj+1)-tspec2(jj))
	      fstar2(i)=der*(tstar2-tspec2(jj))+fjj(i)
	    enddo
	  endif          
	endif
! 	star3 -intrinsic nonrotated 3. intensity spectrum 
!		in erg/cm^2/s/Hz/sterad
!	xstar3(i) -lambda in [A]
	if(lunt3.gt.0)then
	  write(*,*)' reading file starspec3'
	  OPEN(14,FILE='starspec3',STATUS='OLD')
	  i=1
          if(lunt3.eq.3)then
            read(14,*,end=100)
            read(14,*,end=100)
          endif
90	  if(lunt3.eq.1)then
	    read(14,*,end=100)xstar3(i),star3(i)
	    if(i.gt.1)then
!	      fix synspec bug		    	    
	      if(xstar3(i).eq.xstar3(i-1))goto 90
	      if(.not.(xstar3(i).gt.xstar3(i-1)))then
                write(*,*)' input error in starspec3 on line ',i
                write(3,*)' 1 input error in starspec3 on line ',i
              endif  
	      if(i+1.gt.mstarx)then
	        write(*,*)' error in starspec3: '
                write(*,*)'dimension mstarx almost exceeded'
	        write(3,*)' 1 error in starspec3: '
                write(3,*)'dimension mstarx almost exceeded'
              endif  
	    endif
!		convert input of Hlambda from synspec to our units
!		with no limb darkening
            xstar3(i)=xstar3(i)*xunt3
	    star3(i)=star3(i)*yunt3*4.d-8
	    star3(i)=star3(i)*xstar3(i)*xstar3(i)/clight
	  endif
	  if(lunt3.eq.2)then
	    read(14,*,end=100)xstar3(i),star3(i)
	    if(i.gt.1)then
	      if(.not.(xstar3(i).gt.xstar3(i-1)))then
                write(*,*)' input error in starspec3 on line ',i
                write(3,*)' 1 input error in starspec3 on line ',i
              endif  
	      if(i+1.gt.mstarx)then
	        write(*,*)' error in starspec3: '
                write(*,*)'dimension mstarx almost exceeded'
	        write(3,*)' 1 error in starspec3: '
                write(3,*)'dimension mstarx almost exceeded'
              endif  
	    endif
	    xstar3(i)=xstar3(i)*xunt3
	    star3(i)=star3(i)*yunt3
	  endif
	  if(lunt3.eq.3)then
	    read(14,*,end=100)idummy,xstar3(i),dummy,star3(i)
	    if(i.gt.1)then
	      if(.not.(xstar3(i).gt.xstar3(i-1)))then
                write(*,*)' input error in starspec3 on line ',i
                write(3,*)' 1 input error in starspec3 on line ',i
              endif  
	      if(i+1.gt.mstarx)then
	        write(*,*)' error in starspec3: '
                write(*,*)'dimension mstarx almost exceeded'
	        write(3,*)' 1 error in starspec3: '
                write(3,*)'dimension mstarx almost exceeded'
              endif  
	    endif
!               convert input of nu, Fnu from coolTlusty.21 to our units
	    xstar3(i)=xunt3*clight/xstar3(i)*1.d8
	    star3(i)=star3(i)*yunt3/pi
	  endif
	  i=i+1
	  goto 90
100	  nstar3=i-1
	  write(6,*)' nstar3=',nstar3
          if(alam(1).lt.xstar3(1).or.                                   &
	  alam(nfreq).gt.xstar3(nstar3))then
	      write(*,*)' error in starspec3: '
              write(*,*)'short range of lambdas'
	      write(3,*)' 1 error in starspec3: '
              write(3,*)'short range of lambdas'
	  endif	       
	  close(14)
	endif
	return
	end
!-----------------------------------------------------------------------
        subroutine lindat(iat,iion,wlab,elo,eup,glo,gr0,gs0,gw0,bij     &
     &  ,nline)
! 		Line data as in Kurucz linelist or in the SYNSPEC code
!       dll -wavelength [nm]
!       cod -element.ion cod, e.g. 26.02. It is interpreted as:
!         26=atomic number=iron, 02=2xtimes ionized i.e. FeIII line
!       gf -log_10 (gf) 
!       elo,eup -energy of the lower and upper level in [1/cm]
!       qlo -quantum number -J of the lower level[=>stat.weight=2*J+1] 
!       qup -quantum number -J of the upper level[=>stat.weight=2*J+1] 
!       gr0,gs0,gw0-radiative, Stark, Van der Waals damping constants
! 	iat - element atomic number
! 	wlab - lambda [A] 
!	freq0 -frequency [1/s]
! 	stio-ionization degree of the ion [1.-neutral,...] 
! 	bij-transition probability [1/s]=Einstein coef. per unit solid
!		angle
        implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	integer:: iat(mline),iion(mline)
	real(i8):: wlab(mline),elo(mline),eup(mline),glo(mline)
        real(i8):: gr0(mline),gs0(mline),gw0(mline),bij(mline)
	real(i8):: dll,cod,gf,qlo,qup,stio
	real(i8):: eloi,eupi,gr0i,gs0i,gw0i
	integer:: nline,i
        write(*,*)' reading file line.dat'
	open(8,file='line.dat',status='old')
	i=1
10      read(8,*,err=20,end=20)dll,cod,gf,eloi,qlo,eupi,qup             &
     &  ,gr0i,gs0i,gw0i
	if(i.gt.mline)then
	  write(*,*)' error: No. of lines>mline in line.dat. Stop.'
	  write(3,*)' 1 error: No. of lines>mline in line.dat. Stop.'
	  call exit(1)
	endif
	elo(i)=eloi
	eup(i)=eupi
	gr0(i)=gr0i
	gs0(i)=gs0i
	gw0(i)=gw0i
	wlab(i)=dll*10.d0
	iat(i)=int(cod)
	stio=(cod-dint(cod))*100.d0+1.d0
	iion(i)=int(stio+0.5d0)
	glo(i)=2.d0*qlo+1.d0
	bij(i)=dexp(2.302585d0*gf)/glo(i)*wlab(i)/7.484d-7
	i=i+1
	goto 10
20	nline=i-1
	close(8)
	write(*,*)' No. of sp. lines read from line.dat=',nline
	return
	end
!-----------------------------------------------------------------------
        subroutine medop(ndim3,nivac,alam1,alamn,alam                   &
     &  ,ophbf,ophbf1,ophbf2,ophff,ophff1,ophff2                        &
     &  ,ophrs,ophrs1,ophrs2,ophn,ophn1,ophn2)
!	calculates the continuum opacities at the particular wavelength
!	via a simple linear interpolation in between the boundaries
!	input:  ndim3,nivac,alam1,alamn
!		alam -actual value of lambda[A]
!		ophbf1,ophbf2,ophff1,ophff2,ophrs1,ophrs2,ophn1,ophn2
!	output: ophbf,ophff,ophrs,ophn
 	implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: ophbf(ndim3),ophbf1(ndim3),ophbf2(ndim3)   
        real(i8):: ophff(ndim3),ophff1(ndim3),ophff2(ndim3)
        real(i8):: ophrs(ndim3),ophrs1(ndim3),ophrs2(ndim3)
        real(i8):: ophn(ndim3),ophn1(ndim3),ophn2(ndim3)
	real(i8):: alam1,alamn,alam,falam
	integer:: ndim3,nivac,kk
	falam=(alam-alam1)/(alamn-alam1)
	do kk=1,nivac
	  ophbf(kk)=(ophbf2(kk)-ophbf1(kk))*falam+ophbf1(kk)
	  ophff(kk)=(ophff2(kk)-ophff1(kk))*falam+ophff1(kk)
	  ophrs(kk)=(ophrs2(kk)-ophrs1(kk))*falam+ophrs1(kk)
	  ophn(kk)=(ophn2(kk)-ophn1(kk))*falam+ophn1(kk)
        enddo
	return
	end
!-----------------------------------------------------------------------
	subroutine fend(temp,ekonc,entot,zip,stavs,nion,endf,endfp)
!	subroutine calculates the electron number density function
!	and its first parcial derivative (after electron num.dens.) 
!	for one chemical element (valid also for H without H-, H2)
!	input:  temp-temperature
!		ekonc-electron number density
!		entot-number density of particular chemical element
!		zip-ionization potentials
!		stavs-partition functions
!		nion-the highest ion of the element considered
!	output:
!		endf-electron number density function
!		endfp-  \parc endf / \parc Ne
!	notation:
!       rr(j)- population of the J-th ion/total element population
!	rt(j)= Ne*N(j)/N(j-1)
!	rt1(j)= Ne**(j-1)*N(j)/N(1)
!	sahan	-calling function 
!		(Saha equation without el.num. dens. term)
!	
 	implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	interface
	  function sahan(te,zp,s1,s2)                 
            implicit none
            integer, parameter:: i4=4,i8=8
            real(i8):: sahan
            real(i8), intent(in):: te,zp,s1,s2
	  end function sahan
	endinterface
	real(i8):: zip(mion-1),stavs(mion),rt(mion),rt1(mion)
	real(i8):: temp,ekonc,entot,endf,endfp,a,b,rr,rrp,pom
	integer:: nion,j,nj1
!	if(nion.lt.2)then
!	  write(*,*)' error in subroutine: fend'
!	  write(3,*)' 1 error in subroutine: fend'
!	  call exit(1)
!	endif
	rt1(1)=1.d0
	a=ekonc**(nion)
	b=dble(nion)*ekonc**(nion-1)
	do 10 j=2,nion
	  rt(j)=sahan(temp,zip(j-1),stavs(j-1),stavs(j))
          rt1(j)=rt1(j-1)*rt(j)
	  nj1=nion-j+1
	  a=a+ekonc**nj1*rt1(j)
	  b=b+dble(nj1)*ekonc**(nion-j)*rt1(j)
10	continue
	endf=0.d0
	endfp=0.d0
	do 20 j=2,nion
	  nj1=nion-j+1
	  rr=ekonc**nj1*rt1(j)/a
	  rrp=(dble(nj1)*ekonc**(nion-j)*rt1(j)-rr*b)/a
	  pom=dble(j-1)*entot
	  endf=endf+pom*rr
	  endfp=endfp+pom*rrp
20	continue
	return
	end
!-----------------------------------------------------------------------
	subroutine fendh(temp,ekonc,entot,zip,stavs,endfh,endfhp)
!	subroutine calculates the electron number density function
!	and its first parcial derivative (after electron num.dens.) 
!	for Hydrogen including H-
!	input:  temp -temperature
!		ekonc -electron number density
!		entot -number density of particular chemical element
!		zip -ionization potentials
!		stavs -partition functions
!	output:
!		endfh -electron number density function
!		endfhp -  \parc endfh / \parc Ne
!	notation:
!       rr(j) - population of the J-th ion/total element population
!	rt(j)= Ne*N(j)/N(j-1)
!	sahan	-calling function 
!		(Saha equation without el.num. dens. term)
!	
 	implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	interface
	  function sahan(te,zp,s1,s2)                 
            implicit none
            integer, parameter:: i4=4,i8=8
            real(i8):: sahan
            real(i8), intent(in):: te,zp,s1,s2
	  end function sahan
	endinterface
	real(i8):: zip(mion-1),stavs(mion)
	real(i8):: temp,ekonc,entot,endfh,endfhp
	real(i8):: pothm,stavhm,rt2,rt3,bot,rr1,rr3,rr1p,rr3p
        endfh=0.d0
        endfhp=0.d0
!	H- ionization potential	and partition function
	pothm=0.7552d0
	stavhm=1.d0
	rt2=sahan(temp,pothm,stavhm,stavs(1))
	rt3=sahan(temp,zip(1),stavs(1),stavs(2))
	bot=ekonc*ekonc+ekonc*rt2+rt3*rt2
!       n(H-)=rr1*n(H), n(HI)=rr2*n(H), n(HII)=rr3*n(H)	
	rr1=ekonc*ekonc/bot
!	rr2=ekonc*rt2/bot
	rr3=rt3*rt2/bot
        endfh=(rr3-rr1)*entot
	rr1p=(ekonc**2*rt2+2.d0*ekonc*rt3*rt2)/bot**2
	rr3p=-rt3*rt2*(2.d0*ekonc+rt2)/bot**2
        endfhp=(rr3p-rr1p)*entot
	return
	end
!-----------------------------------------------------------------------
	subroutine fendh2 (temp,ekonc,entot,zip,stavs,endfh,endfhp)
!	subroutine calculates the electron number density function
!	and its first parcial derivative (after electron num.dens.) 
!	for Hydrogen including HII,HI,H-,H2
!	input:  temp-temperature
!		ekonc-electron number density
!		entot-number density of particular chemical element
!		zip-ionization potentials
!		stavs-partition functions
!	output:
!		endfh-electron number density function
!		endfhp-  \parc endfh / \parc Ne
!       rr(j)- population of the J-th ion/total element population
!	rt(j)= Ne*N(j)/N(j-1)
!	sahan	-calling function 
!		(Saha equation without el.num. dens. term)
!	
 	implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	interface
	  function sah2(te)        
            implicit none
            integer, parameter:: i4=4,i8=8
            real(i8):: sah2
            real(i8), intent(in):: te
	  end function sah2
	  function sahan(te,zp,s1,s2)                 
            implicit none
            integer, parameter:: i4=4,i8=8
            real(i8):: sahan
            real(i8), intent(in):: te,zp,s1,s2
	  end function sahan
	endinterface
	real(i8):: zip(mion-1),stavs(mion)
	real(i8):: temp,ekonc,entot,endfh,endfhp,pothm,stavhm
	real(i8):: rt1,rt2,rt3,rt3ek,b,c,b2,d,hi,hib,bne,hine,bot,bot2
	real(i8):: rr1,rr3,rr1p,rr3p,ekrt2
        endfh=0.d0
        endfhp=0.d0
!	H- ionization potential	and partition function
	pothm=0.7552d0
	stavhm=1.d0
	if(temp.lt.1.d2)then
!	  rt1 will exceed 1.d300 i.e. infti for T<70K	
	  return
	endif
	if(temp.lt.2.d4)then
!         HII+HI+H-+H2	
!	  for densities< 1.d-20 there are problems for 1d4<T<2d4 K
!	  so you may use if(temp.lt.1.d4)then
	  rt1=sah2(temp)
	  rt2=sahan(temp,pothm,stavhm,stavs(1))
	  rt3=sahan(temp,zip(1),stavs(1),stavs(2))
	  rt3ek=rt3/ekonc
!	  write(*,'(a,3es15.5)')'r123',rt1,rt2,rt3
	  ekrt2=ekonc/rt2
	  b=0.5d0/rt1*(1.d0+rt3ek+ekrt2)
	  c=-0.5d0*entot/rt1
	  b2=b*b
	  d=dsqrt(b2-4.d0*c)
!	  write(*,*)'-4.d0*c/b2',-4.d0/b*c/b
	  if(-4.d0/b*c/b.lt.1.d-10)then
!           then avoid substracting two similar numbers -b+b
	    hi=-c/b
	    endfh=(rt3ek-ekrt2)*hi
	    hib=(c/b)/b
	    bne=0.5d0/rt1*(1.d0/rt2-rt3ek/ekonc)
	    hine=hib*bne
	    endfhp=(-rt3ek/ekonc-1.d0/rt2)*hi+(rt3ek-ekrt2)*hine
	  else
	    hi=0.5d0*(-b+d)
            endfh=(rt3ek-ekrt2)*hi
!	    hib=parc HI / parc b,  bne= parc b / parc ne
!	    hine= parc HI / parc ne 	
	    hib=0.5d0*(-1.d0+b/d)
	    bne=0.5d0/rt1*(1.d0/rt2-rt3ek/ekonc)
	    hine=hib*bne
            endfhp=(-rt3ek/ekonc-1.d0/rt2)*hi+(rt3ek-ekrt2)*hine
          endif
        else
!         HII+HI+H-
!         n(H-)=rr1*n(H), n(HI)=rr2*n(H), n(HII)=rr3*n(H)
 	  rt2=sahan(temp,pothm,stavhm,stavs(1))
	  rt3=sahan(temp,zip(1),stavs(1),stavs(2))
	  bot=ekonc*ekonc+ekonc*rt2+rt3*rt2
	  rr1=ekonc*ekonc/bot
!	  rr2=ekonc*rt2/bot
	  rr3=rt3*rt2/bot
          endfh=(rr3-rr1)*entot
	  bot2=bot**2
	  rr1p=(ekonc**2*rt2+2.d0*ekonc*rt3*rt2)/bot2
	  rr3p=-rt3*rt2*(2.d0*ekonc+rt2)/bot2
          endfhp=(rr3p-rr1p)*entot
        endif  
	return
	end
!-----------------------------------------------------------------------
        subroutine elnd(d,xi,necod,wm,abhyd,ftemp,fdens,fne             &
     &  ,nbod1,nbod2,nbod3,denvac,dcut1,ane0)
!       calculates electron number density from the temperature, density
!	and abundances using Newton-Raphson iteration method 
!	input:  
!		d(2,*) -abundance with respect to H
!		xi -first 8 ionization potentials
!		necod=1 -code for taking the element into account
!		wm -mean molecular weight in [hjed]
!		abhyd -H abundance relative to total elem.num.dens.
!		ftemp,fdens -temperature,density
!		nbod1,nbod2,nbod3 -number of xyz grid points
!		denvac,dcut1 -density limits for vacuum and opaque obj.
!		ane0 -fake electron num.dens beyond the limits
!	output:
!		fne-electron number density
!	notation:
!		xacc-accuracy
        implicit none
	integer, parameter:: i4=4,i8=8
        include 'param.inc'     
        real(i8):: zip(mion-1),stavs(mion),u(5)
        real(i8):: d(3,matom),xi(8,matom)
	integer:: necod(matom)
        real(i8):: ftemp(ndimf1,ndimf2,ndimf3)
        real(i8):: fdens(ndimf1,ndimf2,ndimf3)
        real(i8):: fne(ndimf1,ndimf2,ndimf3)
	real(i8):: wm,abhyd,denvac,dcut1,ane0,hjed,wmhj,xacc
	real(i8):: temp,ekonc,pom,theta,eplow,endf,endfp,entot
	real(i8):: endfi,endfpi,dx
	integer:: nbod1,nbod2,nbod3,niter,l,m,n,i,j,iat,nion
	hjed=1.66053873d-24
	wmhj=hjed*wm
	xacc=1.d-3
	niter=35
	do 70 l=1,nbod1
	do 60 m=1,nbod2
	do 50 n=1,nbod3
	  if(fdens(l,m,n).lt.denvac.or.fdens(l,m,n).gt.dcut1)then
	    fne(l,m,n)=ane0
	    goto 50
	  endif
	  temp=ftemp(l,m,n)
	  ekonc=fdens(l,m,n)/wmhj
	  pom=ekonc*abhyd
 	  theta=5039.77d0/temp
!	  	iteration cycle
	  do i=1,niter
 	    call potlow(temp,ekonc,eplow)
	    endf=0.d0
	    endfp=0.d0
	    do iat=1,matom
	      if(necod(iat).eq.1)then
	      nion=int(d(3,iat)+0.5d0)
              do j=1,8
  	        zip(j)=xi(j,iat)
	      enddo
!	        PARTFN is in 'pfdwor.inc'
!	      mind, for temperatures over a few 1.d5 deg. partition f.
!	      might do strange things 
 	      call partfn(iat,theta,u,eplow) 
 	      do j=1,nion
 	        if(j.le.5)then
 	          stavs(j)=u(j)
 	        else
 	          stavs(j)=1.d0
 	        endif
	      enddo
	      entot=pom*d(2,iat)
	      if(iat.eq.1)then
!	        call fendh(temp,ekonc,entot,zip,stavs,endfi,endfpi)
	        call fendh2(temp,ekonc,entot,zip,stavs,endfi,endfpi)
	      else
                call fend(temp,ekonc,entot,zip,stavs,nion,endfi,endfpi)
              endif
	      endf=endf+endfi
	      endfp=endfp+endfpi
	      endif
            enddo
!	  	heart of iteration
	    dx=(endf-ekonc)/(endfp-1.d0)
	    if((ekonc-dx).lt.1.d-3*ekonc)then
	      ekonc=1.d-3*ekonc
	    else
	      ekonc=ekonc-dx
	    endif
	    if(dabs(dx/ekonc).lt.xacc)goto 45
	    if(i.eq.niter)then
	      write(*,*)' error, slow/no convergence in sub.: elnd'
	      write(3,*)' 1 error, slow/no convergence in sub.: elnd'
	      call exit(1)
	    endif
          enddo
45	  fne(l,m,n)=ekonc
50 	continue
60 	continue
70 	continue
	return
	end
!-----------------------------------------------------------------------
	subroutine trans(x,y,z,trmt)
!	Assumes coordinate system e1,e2,e3.
!	Defines new coordinate system ep1,ep2,ep3 such that ep(3)
!	is along the chosen direction ep3=(x,y,z).
!	Optionaly defines ep2 axes such that ep2 is perpendicular 
!	to e3,ep3 i.e. ep2=e3xep3.
!	Calculates the transformation matrix trmt between 
!	the old and new rotated coordinates such that the new 
!	coordinates (Vp) are related to old coordinates (V):
!	Vp=trmt \times V.  
!	input: x,y,z   (chosen direction)
! 	output:trmt
	implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: trmt(3,3),e(3,3),ep(3,3)
	real(i8):: x,y,z,pom1,pom2,pom3
	integer:: i,j
	e(1,1)=1.d0
	e(1,2)=0.d0
	e(1,3)=0.d0
	e(2,1)=0.d0
	e(2,2)=1.d0
	e(2,3)=0.d0
	e(3,1)=0.d0
	e(3,2)=0.d0
	e(3,3)=1.d0
	pom1=x*x+y*y
	pom2=dsqrt(pom1)
	pom3=dsqrt(x*x*z*z+y*y*z*z+pom1*pom1)
	if(pom2.lt.1.d-30)then
	  ep(1,1)=1.d0
	  ep(1,2)=0.d0
	  ep(1,3)=0.d0
          ep(2,1)=0.d0
          ep(2,2)=1.d0
          ep(2,3)=0.d0
 	else
	  ep(1,1)=x*z/pom3
	  ep(1,2)=y*z/pom3
	  ep(1,3)=-pom1/pom3
	  ep(2,1)=-y/pom2
	  ep(2,2)=x/pom2
	  ep(2,3)=0.d0
	endif
	ep(3,1)=x
	ep(3,2)=y
	ep(3,3)=z
	do i=1,3
	  do j=1,3
	    trmt(i,j)=ep(i,1)*e(j,1)+ep(i,2)*e(j,2)+ep(i,3)*e(j,3)
 	  enddo
	enddo
	return
	end
!-----------------------------------------------------------------------
        subroutine albedo(albx,alby,nalb,ialb)
!	subroutine reads the albedo of the star and companion from 
!	the table
!	albx(i) -wavelength in Angstrom
!	alby(i) -monochromatic albedo
        implicit none
	integer, parameter:: i4=4,i8=8
        include 'param.inc'
        real(i8):: albx(mstarx),alby(mstarx)
	integer:: nalb,ialb,i
        if(ialb.eq.1)then
          write(*,*)'reading file albedo1'
          open(12,file='albedo1',status='old')
          i=1
10        read(12,*,err=20,end=20)albx(i),alby(i)
          i=i+1
          goto 10
20	  nalb=i-1          
          close(12)
        endif
        if(ialb.eq.2)then
          write(*,*)'reading file albedo2'
          open(13,file='albedo2',status='old')
          i=1
30        read(13,*,err=40,end=40)albx(i),alby(i)
          i=i+1
          goto 30
40	  nalb=i-1          
          close(13)
        endif        
        return   
        end          
!-----------------------------------------------------------------------
        subroutine mie(opmiex,opmies,opmiea,nmie,pfmiex,pfang,pfmie     &
     &  ,npfmie,imie,imiepf,ndust,dtlow,dthig,drmf)
!       subroutine reads the Mie data for the dust from the tables:
!       dust_opac:
!       ndust -number of dust species (or separate input files)
!       dtlow, dthig -temperature range of each species [K]
!       drmf -species relative mas fraction within the dust <0,1>
!       file -file names of the tables with: 
!       opmiex(i) -frequency in Hz
!       opmies(i) -scattering opacity per gram of dust material
!       opmiea(i) -absorption opacity per gram of dust material
!       mie_phase:
!	pfmiex(i) -frequency in Hz 
!	pfang(j)  -angle in deg. (cos of the angle on output)
!	pfmie(i,j) -phase function at frequency -i, angle -j,
!		normalised over all space angles to 4pi
!
!       nmie -number of frequencies in dust_opac files
!       npfmie -number of frequencies in mie_phase
!       npfang -number of angles for the Mie phase function 
!	imie -option of threating dust opacity & emissivity
!       imiepf -extra option for treating angle dependent scattering
        implicit none
	integer, parameter:: i4=4,i8=8
        include 'param.inc'
        character(50):: file(mspecx)
        real(i8):: dtlow(mspecx),dthig(mspecx),drmf(mspecx)
        integer:: nmie(mspecx)
	real(i8):: opmiex(mstarx,mspecx)
	real(i8):: opmies(mstarx,mspecx),opmiea(mstarx,mspecx) 
	real(i8):: pfmiex(mstarx),pfang(npfang),pfmie(mstarx,npfang)
	real(i8):: pi,clight,a,b,c,sum,sump
	integer:: npfmie,imie,imiepf,ndust,idtab,i,j
	pi=3.1415926535897931d0
	clight=2.99792458d10
	write(*,*)'reading file dust_opac'
	open(12,file='dust_opac',status='old')
	read(12,*)
	read(12,*)ndust,idtab
!	idtab -table format (1-Budaj, 2-Semenov)	
	read(12,*)
	do i=1,ndust
	  read(12,*)dtlow(i),dthig(i),drmf(i)
	enddo 
	read(12,*)
	do i=1,ndust
	  read(12,'(a)')file(i)
        enddo 
        close(12)
	if(idtab.eq.1)then
!	Budaj format, opacities [cm**2/g] per gram of dust
	  do j=1,ndust
            open(12,file=file(j),status='old')
            i=1
10    	    read(12,*,err=20,end=20)a,b,opmies(i,j),opmiea(i,j)
	    opmiex(i,j)=clight/b*1.d4
            i=i+1
            goto 10
20          nmie(j)=i-1 
            close(12)
          enddo
        else
!	Semenov format, opacities [cm**2/g] per gram of dust+gas        
!	multiply drmf by 100. (gas/dust ratio) in dust_opac
	  do j=1,ndust
            open(12,file=file(j),status='old')
            i=1
23    	    read(12,*,err=25,end=25)a,b,c,opmiea(i,j),opmies(i,j)
	    opmiex(i,j)=clight/a*1.d4
            i=i+1
            goto 23
25          nmie(j)=i-1 
            close(12)
          enddo
        endif  
        write(*,'(a,i3)')' dust species=',ndust
        do i=1,ndust
          write(*,'(a,i3,i5)')' freq. points=',i,nmie(i)
        enddo  
!        	phase functions
	if(imiepf.eq.1)then
	  write(*,*)'reading file mie_phase'
          open(13,file='mie_phase',status='old')
          i=1
30        read(13,*,err=40,end=40)pfmiex(i)
	  do j=1,npfang
	    read(13,*)pfang(j),pfmie(i,j)
	    pfang(j)=pi/180.d0*pfang(j)
	  enddo
	  sum=0.d0
	  do j=2,npfang
	    sump=pfmie(i,j-1)*dsin(pfang(j-1))
	    sump=sump+pfmie(i,j)*dsin(pfang(j))
	    sump=sump*(pfang(j)-pfang(j-1))
	    sum=sum+sump
	  enddo
	  sum=sum/2.d0
	  do j=1,npfang
	    pfang(j)=dcos(pfang(j))
	    pfmie(i,j)=2.d0*pfmie(i,j)/sum
	  enddo
          i=i+1
          goto 30
40        npfmie=i-1 
          close(13) 
        endif  
	return   
        end 
!-----------------------------------------------------------------------
	subroutine interp(x,y,xnew,ynew,ndat,nnew)
!	interpolates in y(x) and returns ynew=y(xnew) 	
!	assumes that x is monotonically increasing
!	input fields: x,y,xnew
!	output field: ynew
!	??-possible problem with dynamical declaration
        implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: x(ndat),y(ndat),xnew(nnew),ynew(nnew)
	real(i8):: der
	integer:: ndat,nnew,jj,i
        jj=1
	do i=1,nnew
	  if(xnew(i).ge.x(ndat))then
	    ynew(i)=y(ndat)
	  elseif(xnew(i).le.x(1))then
	    ynew(i)=y(1)
	  else
	    call hunt(x,ndat,xnew(i),jj)
	    der=(y(jj+1)-y(jj))/(x(jj+1)-x(jj))
	    ynew(i)=der*(xnew(i)-x(jj))+y(jj)
	  endif  
	enddo
        return       
        end
!-----------------------------------------------------------------------
	subroutine intrp(x,y,ndat,xvalue,yvalue)
!	general subroutine for linear interpolation	
!	x may be increasing or decreasing
!	input: x, y, xvalue
!	output: yvalue
	implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	real(i8):: x(mstarx),y(mstarx)
	real(i8):: xvalue,yvalue,der
	integer:: ndat,jj
	if(x(1).lt.x(ndat).and.xvalue.le.x(1))then
	  yvalue=y(1)
	elseif(x(1).lt.x(ndat).and.xvalue.ge.x(ndat))then
          yvalue=y(ndat)
	elseif(x(1).gt.x(ndat).and.xvalue.ge.x(1))then
	  yvalue=y(1)
	elseif(x(1).gt.x(ndat).and.xvalue.le.x(ndat))then
	  yvalue=y(ndat)
	else
          call locate(x,ndat,xvalue,jj)
          der=(y(jj+1)-y(jj))/(x(jj+1)-x(jj))
          yvalue=der*(xvalue-x(jj))+y(jj)
	endif                         
	return     
        end 
!-----------------------------------------------------------------------
	subroutine intrp2(x,pfmie,ndat,xvalue,pfmief)
!	special linear interpolation
!	x may be increasing or decreasing
!	input: x, pfmie, xvalue
!	output: pfmief
	implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	real(i8):: x(mstarx),pfmie(mstarx,npfang),pfmief(npfang)
	real(i8):: xvalue,der
	integer:: ndat,i,jj
	if(x(1).le.x(ndat).and.xvalue.le.x(1))then
	  do i=1,npfang
	    pfmief(i)=pfmie(1,i)
	  enddo  
	elseif(x(1).le.x(ndat).and.xvalue.ge.x(ndat))then
	  do i=1,npfang
            pfmief(i)=pfmie(ndat,i)
          enddo  
	elseif(x(1).ge.x(ndat).and.xvalue.ge.x(1))then
	  do i=1,npfang
             pfmief(i)=pfmie(1,i)
          enddo  
	elseif(x(1).ge.x(ndat).and.xvalue.le.x(ndat))then
          do i=1,npfang
            pfmief(i)=pfmie(ndat,i) 
          enddo         
	else
	  do i=1,npfang
            call locate(x,ndat,xvalue,jj)
            der=(pfmie(jj+1,i)-pfmie(jj,i))/(x(jj+1)-x(jj))
            pfmief(i)=der*(xvalue-x(jj))+pfmie(jj,i)
          enddo  
	endif                         
	return     
        end         
!-----------------------------------------------------------------------
        subroutine xsec(alam1,alamn,cutoff,xsecx,xsecy,txsec,nxsec,ny   &
     &  ,amixr)
!	reads extra cross-section from EXOMOL for gas molecules
!	x -wavenumber=1/lamda is increasing
!       y -cross-section for several temperatures
!       xsecx - frequency [Hz]
!       xsecy -cross-section [cm^2] for several temperatures
!       txsec -temperatures [K]
!	ny -muber of temperatures
!       nxsec -number of frequency points
!	amixr -abundance (mixing ratio) of the molecules 
!		relative to the H nuclei number density
!	input: 'gas_opac', alam1, alamn, cutoff [Ang]
!	output: xsecx,xsecy,txsec,nxsec,ny,amixr
	implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	real(i8):: xsecx(mstarx),xsecy(mstarx,mspecx)
	real(i8):: txsec(mspecx)
	real(i8):: alam1,alamn,cutoff
	real(i8):: clight,xsec1,xsecn,amixr,x
	integer:: nx,ny,nxsec,i,j
	clight=2.99792458d10
	xsec1=1.d8/(alam1-cutoff)
	xsecn=1.d8/(alamn+cutoff)
	write(*,*)'reading file gas_opac'
	open(16,file='gas_opac',status='old')
	read(16,*)nx,ny,amixr
	if(ny.lt.2.or.ny.gt.mspecx)then
	  write(*,*)'1 error in gas_opac: ny.lt.2.or.ny.gt.mspecx'
	  write(3,*)'1 error in gas_opac: ny.lt.2.or.ny.gt.mspecx'
	  call exit(1)
	endif
	if(xsecn.gt.xsec1)then
	  write(*,*)' 1 error in gas_opac: xsecn.gt.xsec1'
	  write(3,*)' 1 error in gas_opac: xsecn.gt.xsec1'
	  call exit(1)
        endif
	read(16,*)(txsec(i),i=1,ny)
!	write(*,*)(txsec(i),i=1,ny)
	nxsec=0
	i=1
	j=1
10	read(16,*,end=30,err=30)x
	if(j.eq.1.and.x.gt.xsec1)then
	  write(*,*)' 1 error in gas_opac: x.gt.xsec1'
	  write(3,*)' 1 error in gas_opac: x.gt.xsec1'
	  call exit(1)
	endif  
	j=j+1
	if(x.ge.xsecn)then
          if(i.eq.1)then
	    if(j.gt.2)then
	      backspace(16)
              backspace(16)
            else
              backspace(16)
            endif
          endif  
20        read(16,*,end=30,err=30)xsecx(i),(xsecy(i,j),j=1,ny)
!		write(*,*)i,xsecx(i)
	  if(xsecx(i).gt.xsec1)then
            nxsec=i
	    goto 40
	  endif  
	  i=i+1
          goto 20
	endif
	goto 10
30      nxsec=i-1
40	write(6,*)'nxsec=',nxsec
	close(16)
	if(nxsec.gt.mstarx)then
	  write(*,*)' 1 error in gas_opac: nxsec.gt.mstarx'
	  write(3,*)' 1 error in gas_opac: nxsec.gt.mstarx'
	  call exit(1)
	endif
	do i=1,nxsec
	  xsecx(i)=xsecx(i)*clight
	enddo
	return     
        end         
!-----------------------------------------------------------------------
	subroutine int2d(x,y,z,nx,ny,xvalue,yvalue,zvalue)
!	subroutine for interpolation in 2D EXOMOL data	
!	x may be increasing or decreasing
!       y may be increasing or decreasing
!	input: x, y, z, xvalue, yvalue
!	output: zvalue
	implicit none
	integer, parameter:: i4=4,i8=8
	include 'param.inc'
	real(i8):: x(mstarx),y(mspecx),z(mstarx,mspecx)
	real(i8):: xvalue,yvalue,zvalue,dx,dz,zz1,zz2,a,b
	integer:: nx,ny,i,j
	if(x(1).lt.x(nx).and.xvalue.le.x(1))then
	  zvalue=0.d0
	elseif(x(1).lt.x(nx).and.xvalue.ge.x(nx))then
          zvalue=0.d0
	elseif(x(1).gt.x(nx).and.xvalue.ge.x(1))then
	  zvalue=0.d0
	elseif(x(1).gt.x(nx).and.xvalue.le.x(nx))then
	  zvalue=0.d0
	elseif(y(1).lt.y(ny).and.yvalue.le.y(1))then
	  zvalue=0.d0
	elseif(y(1).lt.y(ny).and.yvalue.ge.y(ny))then
          zvalue=0.d0
!		extend the table beyond the edge
!          call locate(x,nx,xvalue,i)
!          dz=(z(i+1,ny)-z(i,ny))
!          dx=(xvalue-x(i))/(x(i+1)-x(i))
!          zvalue=z(i,ny)+dz*dx
	elseif(y(1).gt.y(ny).and.yvalue.ge.y(1))then
	  zvalue=0.d0
	elseif(y(1).gt.y(ny).and.yvalue.le.x(ny))then
	  zvalue=0.d0
	else
          call locate(x,nx,xvalue,i)
          call locate(y,ny,yvalue,j)
!		linear interp., 30% faster
!          dzdx=(z(i+1,j)-z(i,j))/(x(i+1)-x(i))
!          dzdy=(z(i,j+1)-z(i,j))/(y(j+1)-y(j))
!          dx=xvalue-x(i)
!          dy=yvalue-y(j)
!          zvalue=z(i,j)+dzdx*dx+dzdy*dy
!		linear in lambda
	  dz=(z(i+1,j)-z(i,j))
	  dx=(xvalue-x(i))/(x(i+1)-x(i))
	  zz1=z(i,j)+dz*dx
	  dz=(z(i+1,j+1)-z(i,j+1))
	  zz2=z(i,j+1)+dz*dx 
!		exponential in temp
	  b=dlog(zz1/zz2)/(1.d0/y(j+1)-1.d0/y(j))
	  a=zz1*dexp(b/y(j))
	  zvalue=a*dexp(-b/yvalue)
	endif          
	return     
        end 
!-----------------------------------------------------------------------
 	subroutine chem(pop,popt,popd,mtemp,mdens,melm)
!	reads chemistry table with equilibrium molecular populations
!	pop,popt,popd -log10 of population, temperature, density
	implicit none
	integer, parameter:: i4=4,i8=8
!	integer, parameter:: mtemp=57,mdens=73,melm=11
	integer:: mtemp,mdens,melm
	integer:: i,j,k,l
	integer:: index(melm)
	character(18):: a18(melm)
	real(i8):: popt(mtemp),popd(mdens),pop(mtemp,mdens,melm)
	write(*,*)' reading file chem_eq_tab'
	open(16,file='chem_eq_tab',status='unknown')
        do i=1,mtemp
          read(16,*)
          read(16,30)popt(i),(popd(j),j=1,mdens)
30        format(6x,18x,f18.5,73f14.2)
          do k=1,melm
            read(16,40)index(k),a18(k),(pop(i,j,k),j=1,mdens)
          enddo  
        enddo
40      format(i6,18x,a18,73f14.5)
	close(16)
 	return
 	end
!-----------------------------------------------------------------------        
        subroutine abnd(pop,popt,popd,mtemp,mdens,melm                  &
     &  ,temp,dens,ndepth,mdepth,iopac,amixr,amixf)
!	CO molecule population from the table
!	input: amixr -molecule abundance relative to H nuclei
!         pop(mtemp,mdens,melm) -log10num.dens(log10temp,log10dens)
!	  popt -log10 temperature, popd -log10 density
!	output: amixf -number density of molecule imol in cm^-3
	implicit none
	integer, parameter:: i4=4,i8=8
	integer:: mtemp,mdens,melm,ndepth,mdepth,i,j,iopac,imol
	real(i8):: t,d,p
	real(i8):: hjed,wm,enh,amixr
	real(i8):: popt(mtemp),popd(mdens),pop(mtemp,mdens,melm)
	real(i8):: dens(mdepth),temp(mdepth)
	real(i8):: amixf(mdepth),popmol(mtemp,mdens)
!	index of HI,H2,CO,H2O	
!	ih=1 ih2=5 ico=6 iw=8
	imol=6
        hjed=1.66053873d-24
!       mean atomic weight for solar composition gas
        wm=1.3d0
!       calculate C,O abundance relative to H
!       abc=2.7d-4 abo=4.9d-4
	do i=1,mtemp
	do j=1,mdens
	  popmol(i,j)=pop(i,j,imol)
	enddo
	enddo
	do i=1,ndepth
	  if(temp(i).lt.255.d0)then
	    enh=0.9d0*dens(i)/wm/hjed
	    amixf(i)=amixr*enh
	  elseif(temp(i).gt.6300.d0)then
	    amixf(i)=0.d0
	  elseif(dens(i).lt.1.d-19.or.dens(i).gt.1.d-1)then
	    enh=0.9d0*dens(i)/wm/hjed
	    amixf(i)=amixr*enh
	  else
	    t=dlog10(temp(i))
	    d=dlog10(dens(i))
	    call int2b(popt,popd,popmol,mtemp,mdens,mtemp,mdens,t,d,p)
	    amixf(i)=10**p
	  endif  
	enddo
        return
        end
!-----------------------------------------------------------------------
	subroutine int2b(x,y,z,nx,ny,mx,my,xvalue,yvalue,zvalue)
!	subroutine for interpolation in 2D
!	x may be increasing or decreasing
!       y may be increasing or decreasing
!	input: x, y, z, xvalue, yvalue
!	output: zvalue
	implicit none
	integer, parameter:: i4=4,i8=8
	integer:: nx,ny,mx,my,i,j
	real(i8):: x(mx),y(my),z(mx,my)
	real(i8):: xvalue,yvalue,zvalue
	real(i8):: dzdx,dzdy,dx,dy
	if(x(1).lt.x(nx).and.xvalue.le.x(1))then
	  zvalue=0.d0
	elseif(x(1).lt.x(nx).and.xvalue.ge.x(nx))then
          zvalue=0.d0
	elseif(x(1).gt.x(nx).and.xvalue.ge.x(1))then
	  zvalue=0.d0
	elseif(x(1).gt.x(nx).and.xvalue.le.x(nx))then
	  zvalue=0.d0
	elseif(y(1).lt.y(ny).and.yvalue.le.y(1))then
	  zvalue=0.d0
	elseif(y(1).lt.y(ny).and.yvalue.ge.y(ny))then
          zvalue=0.d0
!		extend the table beyond the edge
!          call locate(x,nx,xvalue,i)
!          dz=(z(i+1,ny)-z(i,ny))
!          dx=(xvalue-x(i))/(x(i+1)-x(i))
!          zvalue=z(i,ny)+dz*dx
	elseif(y(1).gt.y(ny).and.yvalue.ge.y(1))then
	  zvalue=0.d0
	elseif(y(1).gt.y(ny).and.yvalue.le.x(ny))then
	  zvalue=0.d0
	else
          call locate(x,nx,xvalue,i)
          call locate(y,ny,yvalue,j)
!		linear interp.
          dzdx=(z(i+1,j)-z(i,j))/(x(i+1)-x(i))
          dzdy=(z(i,j+1)-z(i,j))/(y(j+1)-y(j))
          dx=xvalue-x(i)
          dy=yvalue-y(j)
          zvalue=z(i,j)+dzdx*dx+dzdy*dy
	endif          
	return     
        end 
!-----------------------------------------------------------------------
	subroutine split1(x1,x2,dmax,s1,s2,cf12)
!	splits a single huge step into modp smaller steps
!	input: 
!		x1,x2 -op.depths
!		dmax -max.step in op.depth
!		s1,s2 -source functions
!	output: cf12 -contrib function
	implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: x1,x2,dmax,s1,s2
!       modp -maximum optical dept points to break a single huge step
	integer, parameter:: modp=2000
	real(i8):: si(modp),xi(modp),cfi(modp)
	integer:: nodp,i
	real(i8):: d,step,der,cf12
	d=x1-x2
	nodp=d/dmax+1
	if(nodp+1.gt.modp)then
          write(*,*)' error: nodp+1>modp in split1',nodp,modp
          write(*,*)'steep/high density gradient on the edge'
          write(3,*)' 1 error: nodp+1>modp in split1',nodp,modp
          write(3,*)'steep/high density gradient on the edge'
          call exit(1)
        endif 
        step=d/dble(nodp)
	der=(s1-s2)/d	
	do i=1,nodp+1
	  si(i)=der*step*dble(i-1)+s2
	  xi(i)=step*dble(i-1)+x2
	  cfi(i)=si(i)*dexp(-xi(i))
	enddo
	cf12=0.d0
	do i=1,nodp
	  cf12=cf12+cfi(i)+cfi(i+1)
	enddo	
	cf12=cf12/2.d0*step
	return     
        end 
!-----------------------------------------------------------------------
	subroutine split2(x1,x2,dmax,s1,s2,cf12)
!	splits a single huge step into n smaller steps
!	input: 
!		x1,x2 -op.depths
!		dmax -1st step in op.depth
!		s1,s2 -source functions
!	output: cf12 -contrib function
	implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: x1,x2,dmax,s1,s2,cf12
!       m -maximum number of optical dept points to break into
	integer, parameter:: m=200
	real(i8):: si(m),xi(m),dxi(m),sxi(m),cfi(m)
	integer:: n,i
	real(i8):: d,der,r,dmax0,r1,rn,ri
	d=x1-x2
	r=1.2d0
	r1=1.d0-r
	n=dlog(1.d0-d*r1/dmax)/dlog(r)+1
	rn=1.d0-r**n
	if((n+1).gt.m)then
          write(*,*)' error: n+1>m in split2',n,m
          write(*,*)'steep/high density gradient on the edge'
          write(3,*)' 1 error: n+1>m in split2',n,m
          write(3,*)'steep/high density gradient on the edge'
          call exit(1)
        endif 
!	modified 1.step
	dmax0=d*r1/rn
	der=(s1-s2)/d
	xi(1)=x2
	dxi(1)=0.d0
	sxi(1)=0.d0
	si(1)=s2
	cfi(1)=si(1)*dexp(-xi(1))	
	ri=1.d0
	do i=2,n+1
	  dxi(i)=dmax0*ri
	  sxi(i)=sxi(i-1)+dxi(i)
	  xi(i)=sxi(i)+x2
	  si(i)=der*sxi(i)+s2
	  cfi(i)=si(i)*dexp(-xi(i))
	  ri=ri*r
	enddo
	cf12=0.d0
	do i=2,n+1
	  cf12=cf12+(cfi(i-1)+cfi(i))*dxi(i)
	enddo	
	cf12=cf12/2.d0
	return     
        end 
!-----------------------------------------------------------------------
	subroutine extccm89(ebv,rv,w,amag,aflx)
!	Extincion curve according to Cardelli, Clayton & Mathis 1989
!       covers: 0.1<=lambda[mic]<=10/3 
!	input: ebv=E(B-V), rv=A(V)/E(B-V), w=wavelength[mic]
!	output: amag=A(lambda)[mag], aflx=linear dimming factor
!       To deredden your observations flux_obs or mag_obs use:
!               flux=flux_obs*aflx
!               mag=mag_obs-amag
	implicit none
	integer, parameter:: i4=4,i8=8
	real(i8):: ebv,rv,w,amag,aflx
	real(i8):: x,av,aa,a,b,y,fa,fb
	av=rv*ebv
!	aa=A(lambda)/A(V)
	x=1.d0/w
!	infrared
	if(x.ge.0.3d0.and.x.lt.1.1d0)then
          a=0.574d0*x**1.61d0
          b=-0.527d0*x**1.61d0
!       optical          
        elseif(x.ge.1.1d0.and.x.lt.3.3d0)then
          y=x-1.82d0
          a=1.d0+0.17699d0*y-0.50447d0*y**2-0.02427d0*y**3
          a=a+0.72085d0*y**4+0.01979d0*y**5-0.77530d0*y**6
          a=a+0.32999d0*y**7
          b=1.41338d0*y+2.28305d0*y**2+1.07233d0*y**3-5.38434d0*y**4
          b=b-0.62251d0*y**5+5.30260d0*y**6-2.09002d0*y**7
!       UV
        elseif(x.ge.3.3d0.and.x.lt.8.d0)then
          if(x.lt.8.d0.and.x.ge.5.9d0)then
            y=x-5.9d0
            fa=-0.04473d0*y**2-0.009779d0*y**3
            fb=+0.2130d0*y**2+0.1207d0*y**3
          else
            fa=0.d0
            fb=0.d0
          endif  
          a=+1.752d0-0.316d0*x-0.104d0/((x-4.67d0)**2+0.341d0)+fa          
          b=-3.090d0+1.825d0*x+1.206d0/((x-4.62d0)**2+0.263d0)+fb
!       far UV
        elseif(x.ge.8.d0.and.x.le.10.d0)then
          y=x-8.d0
          a=-1.073d0-0.628d0*y+0.137d0*y**2-0.070d0*y**3
          b=13.670d0+4.257d0*y-0.420d0*y**2+0.374d0*y**3
!       out of range
        else
          a=0.d0
          b=0.d0      
        endif  
        aa=a+b/rv
        amag=aa*av
        aflx=10.d0**(amag*0.4d0)
	return
	end
!-----------------------------------------------------------------------
	include 'pfdwor_inc.f90'
