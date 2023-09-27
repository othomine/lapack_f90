c     cabs1sub.f
c
c     The program is a fortran wrapper for cabs1.
c
      subroutine scabs1sub(z, cabs1)
c
      complex z
      real cabs1
c
      cabs1=ABS(REAL(Z)) + ABS(AIMAG(Z))
      return
      end
