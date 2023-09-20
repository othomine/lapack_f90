!> \brief \b SROTM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SPARAM(5),SX(*),SY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
!>
!>    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
!>    (SX**T)
!>
!>    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX  >=  0, ELSE
!>    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
!>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!>
!>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
!>
!>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
!>    H=(          )    (          )    (          )    (          )
!>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
!>    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
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
!> \param[in,out] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
!> \endverbatim
!>
!> \param[in,out] SY
!> \verbatim
!>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of SY
!> \endverbatim
!>
!> \param[in] SPARAM
!> \verbatim
!>          SPARAM is REAL array, dimension (5)
!>     SPARAM(1)=SFLAG
!>     SPARAM(2)=SH11
!>     SPARAM(3)=SH21
!>     SPARAM(4)=SH12
!>     SPARAM(5)=SH22
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
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!
!> \ingroup rotm
!
!  =====================================================================
   SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
   REAL SPARAM(5),SX(*),SY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   REAL SFLAG,SH11,SH12,SH21,SH22,W,Z
   INTEGER I,KX,KY,NSTEPS
!     ..
!
   SFLAG = SPARAM(1)
   IF (N <= 0 .OR. (SFLAG+2.E0 == 0.E0)) RETURN
   IF (INCX == INCY.AND.INCX > 0) THEN
!
      NSTEPS = N*INCX
      IF (SFLAG < 0.E0) THEN
         SH11 = SPARAM(2)
         SH12 = SPARAM(4)
         SH21 = SPARAM(3)
         SH22 = SPARAM(5)
         DO I = 1,NSTEPS,INCX
            W = SX(I)
            Z = SY(I)
            SX(I) = W*SH11 + Z*SH12
            SY(I) = W*SH21 + Z*SH22
         END DO
      ELSE IF (SFLAG == 0.E0) THEN
         SH12 = SPARAM(4)
         SH21 = SPARAM(3)
         DO I = 1,NSTEPS,INCX
            W = SX(I)
            Z = SY(I)
            SX(I) = W + Z*SH12
            SY(I) = W*SH21 + Z
         END DO
      ELSE
         SH11 = SPARAM(2)
         SH22 = SPARAM(5)
         DO I = 1,NSTEPS,INCX
            W = SX(I)
            Z = SY(I)
            SX(I) = W*SH11 + Z
            SY(I) = -W + SH22*Z
         END DO
      END IF
   ELSE
      KX = 1
      KY = 1
      IF (INCX < 0) KX = 1 + (1-N)*INCX
      IF (INCY < 0) KY = 1 + (1-N)*INCY
!
      IF (SFLAG < 0.E0) THEN
         SH11 = SPARAM(2)
         SH12 = SPARAM(4)
         SH21 = SPARAM(3)
         SH22 = SPARAM(5)
         DO I = 1,N
            W = SX(KX)
            Z = SY(KY)
            SX(KX) = W*SH11 + Z*SH12
            SY(KY) = W*SH21 + Z*SH22
            KX = KX + INCX
            KY = KY + INCY
         END DO
      ELSE IF (SFLAG == 0.E0) THEN
         SH12 = SPARAM(4)
         SH21 = SPARAM(3)
         DO I = 1,N
            W = SX(KX)
            Z = SY(KY)
            SX(KX) = W + Z*SH12
            SY(KY) = W*SH21 + Z
            KX = KX + INCX
            KY = KY + INCY
         END DO
      ELSE
          SH11 = SPARAM(2)
          SH22 = SPARAM(5)
          DO I = 1,N
             W = SX(KX)
             Z = SY(KY)
             SX(KX) = W*SH11 + Z
             SY(KY) = -W + SH22*Z
             KX = KX + INCX
             KY = KY + INCY
         END DO
      END IF
   END IF
   RETURN
!
!     End of SROTM
!
END
