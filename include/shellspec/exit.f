
c-----------------------------------------------------------------------
        subroutine exit(exitcode)
c
c the simplest exit (but one can add something in the future)
c
c input:
c   exitcode ... dtto
c
c   0 ... success
c   1 ... error reading input parameters
c   2 ... convergence problems during RTE integration
c   3 ... interpolation or extrapolation problems
c   4 ... array dimensions too small
c
        implicit none
        integer exitcode
        if (exitcode.eq.4) then
          stop 4
        else if (exitcode.eq.3) then
          stop 3
        else if (exitcode.eq.2) then
          stop 2
        else if (exitcode.eq.1) then
          stop 1
        else
          stop 0
        endif
        end
