!> \brief \b SASUM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SASUM(N,SX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SASUM takes the sum of the absolute values.
!>    uses unrolled loops for increment equal to one.
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
!> \param[in] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
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
   REAL FUNCTION SASUM(N,SX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER INCX,N
!     ..
!     .. Array Arguments ..
   REAL SX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   REAL STEMP
   INTEGER I,NINCX
!     ..
   SASUM = 0.0e0
   STEMP = 0.0e0
   IF (N <= 0 .OR. INCX <= 0) RETURN
   IF (INCX == 1) THEN
!        code for increment equal to 1
!
!
!        clean-up loop
!
      STEMP = sum(ABS(SX(1:N)))
   ELSE
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO I = 1,NINCX,INCX
         STEMP = STEMP + ABS(SX(I))
      END DO
   END IF
   SASUM = STEMP
   RETURN
!
!     End of SASUM
!
END
