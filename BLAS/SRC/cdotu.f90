!> \brief \b CDOTU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       COMPLEX CX(*),CY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CDOTU forms the dot product of two complex vectors
!>      CDOTU = X^T * Y
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
!> \param[in] CX
!> \verbatim
!>          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of CX
!> \endverbatim
!>
!> \param[in] CY
!> \verbatim
!>          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of CY
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
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!>
!>     converted to F90 and optimized 2023, olivier thomine
!> \endverbatim
!>
!  =====================================================================
   COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER, INTENT(IN) :: INCX,INCY,N
!     ..
!     .. Array Arguments ..
   COMPLEX, INTENT(IN) :: CX(*),CY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   COMPLEX CTEMP
   INTEGER I,IX,IY
!     ..
   CDOTU = (0.0,0.0)
   IF (N <= 0) RETURN
   IF (INCX == 1 .AND. INCY == 1) THEN
!
!        code for both increments equal to 1
!
      CDOTU = sum(CX(1:N)*CY(1:N))
      RETURN
   ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      CTEMP = (0.0,0.0)
      IX = 1
      IY = 1
      IF (INCX < 0) IX = (-N+1)*INCX + 1
      IF (INCY < 0) IY = (-N+1)*INCY + 1
      DO I = 1,N
         CTEMP = CTEMP + CX(IX)*CY(IY)
         IX = IX + INCX
         IY = IY + INCY
      END DO
      CDOTU = CTEMP
      RETURN
   END IF
!
!     End of CDOTU
!
END
