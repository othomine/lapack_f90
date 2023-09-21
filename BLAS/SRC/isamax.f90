!> \brief \b ISAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ISAMAX(N,SX,INCX)
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
!>    ISAMAX finds the index of the first element having maximum absolute value.
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
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   INTEGER FUNCTION ISAMAX(N,SX,INCX)
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
   REAL SMAX
   INTEGER I,IX
!     ..
   ISAMAX = 0
   IF (N < 1 .OR. INCX <= 0) RETURN
   ISAMAX = 1
   IF (N == 1) RETURN
   IF (INCX == 1) THEN
!
!        code for increment equal to 1
!
      ISAMAX = maxloc(ABS(SX(1:N)),1)
   ELSE
!
!        code for increment not equal to 1
!
      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO I = 2,N
         IF (ABS(SX(IX)) > SMAX) THEN
            ISAMAX = I
            SMAX = ABS(SX(IX))
         END IF
         IX = IX + INCX
      END DO
   END IF
   RETURN
!
!     End of ISAMAX
!
END
