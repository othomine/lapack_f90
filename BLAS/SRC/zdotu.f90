!> \brief \b ZDOTU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*),ZY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDOTU forms the dot product of two complex vectors
!>      ZDOTU = X^T * Y
!>
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
!> \param[in] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
!> \endverbatim
!>
!> \param[in] ZY
!> \verbatim
!>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of ZY
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
!> \ingroup dot
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
   COMPLEX*16 ZX(*),ZY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   COMPLEX*16 ZTEMP
   INTEGER I,IX,IY
!     ..
   ZTEMP = (0.0d0,0.0d0)
   ZDOTU = (0.0d0,0.0d0)
   IF (N <= 0) RETURN
   IF (INCX == 1 .AND. INCY == 1) THEN
!
!        code for both increments equal to 1
!
      ZTEMP = sum(ZX(1:N)*ZY(1:N))
   ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX < 0) IX = (-N+1)*INCX + 1
      IF (INCY < 0) IY = (-N+1)*INCY + 1
      DO I = 1,N
         ZTEMP = ZTEMP + ZX(IX)*ZY(IY)
         IX = IX + INCX
         IY = IY + INCY
      END DO
   END IF
   ZDOTU = ZTEMP
   RETURN
!
!     End of ZDOTU
!
END
