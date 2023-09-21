!> \brief \b ZDSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZDSCAL scales a vector by a constant.
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
!> \param[in] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!>           On entry, DA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in,out] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
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
!> \ingroup scal
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION DA
   INTEGER INCX,N
!     ..
!     .. Array Arguments ..
   COMPLEX*16 ZX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   INTEGER NINCX
!     ..
   IF (N <= 0 .OR. INCX <= 0 .OR. DA == 1.0D+0) RETURN
   IF (INCX == 1) THEN
!
!        code for increment equal to 1
!
      ZX(1:N) = DCMPLX(DA*DBLE(ZX(1:N)),DA*DIMAG(ZX(1:N)))
   ELSE
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      ZX(1:NINCX:INCX) = DCMPLX(DA*DBLE(ZX(1:NINCX:INCX)),DA*DIMAG(ZX(1:NINCX:INCX)))
   END IF
   RETURN
!
!     End of ZDSCAL
!
END
