!> \brief \b ICAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ICAMAX(N,CX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX CX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup iamax
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
   INTEGER FUNCTION ICAMAX(N,CX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER INCX,N
!     ..
!     .. Array Arguments ..
   COMPLEX CX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   INTEGER I,IX
!     ..
!     .. External Functions ..
   REAL SCABS1
   EXTERNAL SCABS1
!     ..
   ICAMAX = 0
   IF (N < 1 .OR. INCX <= 0) RETURN
   ICAMAX = 1
   IF (N == 1) RETURN
   IF (INCX == 1) THEN
!
!        code for increment equal to 1
!
      ICAMAX = maxloc(ABS(REAL(CX(1:N))) + ABS(AIMAG(CX(1:N))),1)
   ELSE
!
!        code for increment not equal to 1
!
      IX = 1
      SMAX = SCABS1(CX(1))
      IX = IX + INCX
      DO I = 2,N
         IF (SCABS1(CX(IX)) > SMAX) THEN
            ICAMAX = I
            SMAX = SCABS1(CX(IX))
         END IF
         IX = IX + INCX
      END DO
   END IF
   RETURN
!
!     End of ICAMAX
!
END
