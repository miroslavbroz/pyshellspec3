c limcof_test.f
c Quadrilinear interpolation of limb-darkening coefficients.
c Miroslav Broz (miroslav.broz@email.cz), Jun 6th 2016

      program limcof_test

      implicit none
      include 'limcof.inc'
c constants
      real*8 Ang
      parameter(Ang=1.d-10)
c parameters
      real*8 Teff, logg, Z
      real*8 lambda1, lambda2, dlambda
c internal
      integer i, j, n, m
      real*8 lambda, u
      real*8 lambda_(LAMBDAMAX), u_(LAMBDAMAX)
      character*255 str
c functions
      real*8 interp

c read input parameters
      if (iargc().lt.6) then
        write(*,*) "Usage: limcof_test Teff logg Z lambda1 lambda2",
     :    " dlambda [Angstroms]"
        stop 1
      endif
      call getarg(1,str)
      read(str,*) Teff
      call getarg(2,str)
      read(str,*) logg
      call getarg(3,str)
      read(str,*) Z
      call getarg(4,str)
      read(str,*) lambda1
      call getarg(5,str)
      read(str,*) lambda2
      call getarg(6,str)
      read(str,*) dlambda

      lambda1 = lambda1*Ang
      lambda2 = lambda2*Ang
      dlambda = dlambda*Ang

c write input
      write(*,*)
      write(*,*) "# Teff = ", Teff
      write(*,*) "# logg = ", logg
      write(*,*) "# Z = ", Z
      write(*,*) "# lambda1 = ", lambda1, " m"
      write(*,*) "# lambda2 = ", lambda2, " m"
      write(*,*) "# dlambda = ", dlambda, " m"
      write(*,*)

      call limcof_read("limcof.dat")
      call limcof_interp(Teff, logg, Z, lambda_, u_, m)

      write(*,*)
      write(*,*) "# lambda [m] & u [] linear limb-darkening coef. & id"
      do j = 1, m
        write(*,*) lambda_(j), u_(j), j
      enddo
      write(*,*)

      i = 0
      j = 2
      lambda = lambda1
      do while (lambda.le.lambda2)
        i = i+1
        do while ((lambda_(j).lt.lambda).and.(j.lt.m))
          j = j+1
        enddo
        if (lambda.lt.lambda_(j-1)) then
          write(*,*) "# Warning: Extrapolation lambda = ", lambda, " m"
          u = u_(j-1)
        else if (lambda.gt.lambda_(j)) then
          write(*,*) "# Warning: Extrapolation lambda = ", lambda, " m"
          u = u_(j)
        else
          u = interp(lambda_(j-1), lambda_(j), u_(j-1), u_(j), lambda)
        endif
        write(*,*) lambda, u, i, j
        lambda = lambda + dlambda
      enddo

      stop 0
      end


