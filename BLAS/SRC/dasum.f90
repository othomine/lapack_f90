!> \brief \b DASUM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DASUM takes the sum of the absolute values.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Olivier Thomine
!
!> \ingroup asum
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER INCX,N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   DOUBLE PRECISION DTEMP
   INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC DABS,MOD
!     ..
   DASUM = 0.0d0
   DTEMP = 0.0d0
   IF (N <= 0 .OR. INCX <= 0) RETURN
   IF (INCX == 1) THEN
!
!
!        clean-up loop
!
      DASUM = sum(DABS(DX(1:N)))
! stop 'has to be changed in original file'
      RETURN
   ELSE
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO I = 1,NINCX,INCX
         DTEMP = DTEMP + DABS(DX(I))
      END DO
      DASUM = DTEMP
      RETURN
   END IF
!
!     End of DASUM
!
END
