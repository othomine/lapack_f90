!> \brief \b DROTMG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DD1,DD2,DX1,DY1
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DPARAM(5)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
!>    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*>    DY2)**T.
!>    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!>
!>    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
!>
!>      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
!>    H=(          )    (          )    (          )    (          )
!>      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
!>    LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
!>    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
!>    VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
!>
!>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
!>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
!>    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] DD1
!> \verbatim
!>          DD1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] DD2
!> \verbatim
!>          DD2 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] DX1
!> \verbatim
!>          DX1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] DY1
!> \verbatim
!>          DY1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] DPARAM
!> \verbatim
!>          DPARAM is DOUBLE PRECISION array, dimension (5)
!>     DPARAM(1)=DFLAG
!>     DPARAM(2)=DH11
!>     DPARAM(3)=DH21
!>     DPARAM(4)=DH12
!>     DPARAM(5)=DH22
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
!> \ingroup rotmg
!
!  =====================================================================
   SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION DD1,DD2,DX1,DY1
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION DPARAM(5)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP, &
                    DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC DABS
!     ..
!     .. Data statements ..
!
   DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
!     ..

   IF (DD1 < 0.D0) THEN
!        GO ZERO-H-D-AND-DX1..
      DFLAG = -1.D0
      DH11 = 0.D0
      DH12 = 0.D0
      DH21 = 0.D0
      DH22 = 0.D0
!
      DD1 = 0.D0
      DD2 = 0.D0
      DX1 = 0.D0
   ELSE
!        CASE-DD1-NONNEGATIVE
      DP2 = DD2*DY1
      IF (DP2 == 0.D0) THEN
         DFLAG = -2.D0
         DPARAM(1) = DFLAG
         RETURN
      END IF
!        REGULAR-CASE..
      DP1 = DD1*DX1
      DQ2 = DP2*DY1
      DQ1 = DP1*DX1
!
      IF (DABS(DQ1) > DABS(DQ2)) THEN
         DH21 = -DY1/DX1
         DH12 = DP2/DP1
!
         DU = 1.D0 - DH12*DH21
!
        IF (DU > 0.D0) THEN
          DFLAG = 0.D0
          DD1 = DD1/DU
          DD2 = DD2/DU
          DX1 = DX1*DU
        ELSE
!            This code path if here for safety. We do not expect this
!            condition to ever hold except in edge cases with rounding
!            errors. See DOI: 10.1145/355841.355847
          DFLAG = -1.D0
          DH11 = 0.D0
          DH12 = 0.D0
          DH21 = 0.D0
          DH22 = 0.D0
!
          DD1 = 0.D0
          DD2 = 0.D0
          DX1 = 0.D0
        END IF
      ELSE

         IF (DQ2 < 0.D0) THEN
!              GO ZERO-H-D-AND-DX1..
            DFLAG = -1.D0
            DH11 = 0.D0
            DH12 = 0.D0
            DH21 = 0.D0
            DH22 = 0.D0
!
            DD1 = 0.D0
            DD2 = 0.D0
            DX1 = 0.D0
         ELSE
            DFLAG = 1.D0
            DH11 = DP1/DP2
            DH22 = DX1/DY1
            DU = 1.D0 + DH11*DH22
            DTEMP = DD2/DU
            DD2 = DD1/DU
            DD1 = DTEMP
            DX1 = DY1*DU
         END IF
      END IF

!     PROCEDURE..SCALE-CHECK
      IF (DD1 /= 0.D0) THEN
         DO WHILE ((DD1 <= RGAMSQ) .OR. (DD1 >= GAMSQ))
            IF (DFLAG == 0.D0) THEN
               DH11 = 1.D0
               DH22 = 1.D0
               DFLAG = -1.D0
            ELSE
               DH21 = -1.D0
               DH12 = 1.D0
               DFLAG = -1.D0
            END IF
            IF (DD1 <= RGAMSQ) THEN
               DD1 = DD1*GAM**2
               DX1 = DX1/GAM
               DH11 = DH11/GAM
               DH12 = DH12/GAM
            ELSE
               DD1 = DD1/GAM**2
               DX1 = DX1*GAM
               DH11 = DH11*GAM
               DH12 = DH12*GAM
            END IF
         ENDDO
      END IF

      IF (DD2 /= 0.D0) THEN
         DO WHILE ( (DABS(DD2) <= RGAMSQ) .OR. (DABS(DD2) >= GAMSQ) )
            IF (DFLAG == 0.D0) THEN
               DH11 = 1.D0
               DH22 = 1.D0
               DFLAG = -1.D0
            ELSE
               DH21 = -1.D0
               DH12 = 1.D0
               DFLAG = -1.D0
            END IF
            IF (DABS(DD2) <= RGAMSQ) THEN
               DD2 = DD2*GAM**2
               DH21 = DH21/GAM
               DH22 = DH22/GAM
            ELSE
               DD2 = DD2/GAM**2
               DH21 = DH21*GAM
               DH22 = DH22*GAM
            END IF
         END DO
      END IF

   END IF

   IF (DFLAG < 0.D0) THEN
      DPARAM(2) = DH11
      DPARAM(3) = DH21
      DPARAM(4) = DH12
      DPARAM(5) = DH22
   ELSE IF (DFLAG == 0.D0) THEN
      DPARAM(3) = DH21
      DPARAM(4) = DH12
   ELSE
      DPARAM(2) = DH11
      DPARAM(5) = DH22
   END IF

   DPARAM(1) = DFLAG
   RETURN
!
!     End of DROTMG
!
END
