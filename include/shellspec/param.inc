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
        integer, parameter:: ndimf1=501,ndimf2=501,ndimf3=501
        integer, parameter:: ndim1=501,ndim2=501,ndim3=501,ndim=ndim3
        integer, parameter:: mion=9,matom=99,mfreq=1001,mline=10
	integer, parameter:: mphase=251
        integer, parameter:: mstarx=500000,mspecx=36,npfang=65
