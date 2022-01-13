!	include file with the maximum field dimensions
!	most memory consuming parameters are: ndim,ndim1-3
!       this was used on a PC with 4GB RAM
!	(301,301,301) grid can readily be used on a PC with 8GB RAM
!
!       ndimf1,2,3 -number of grid points in 3D (body frozen grid)
!       ndim1,2,3 -number of grid points in 3D (line of sight grid)
!	mion, matom -number of ions and atoms considered
!	mfreq -number of frequences
!	mline -number of spectral lines
!	mphase -number of phases
!	mstarx -number of frequences for stellar spectra, albedo, 
!		dust, and gas opacity tables
!	mspecx -number of spectra in the temperature grid of
!               stellar spectra, gas opacities
!	npfang -number of angles for the Mie scattering phase function
!
        integer, parameter:: ndimf1=401,ndimf2=401,ndimf3=401
        integer, parameter:: ndim1=401,ndim2=401,ndim3=401,ndim=ndim3
        integer, parameter:: mion=9,matom=99,mfreq=1001,mline=10
	integer, parameter:: mphase=200
        integer, parameter:: mstarx=500000,mspecx=36,npfang=65
