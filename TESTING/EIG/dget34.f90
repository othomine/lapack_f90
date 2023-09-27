!> \brief \b DGET34
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET34( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX
!       DOUBLE PRECISION   RMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            NINFO( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET34 tests DLAEXC, a routine for swapping adjacent blocks (either
!> 1 by 1 or 2 by 2) on the diagonal of a matrix in real Schur form.
!> Thus, DLAEXC computes an orthogonal matrix Q such that
!>
!>     Q' * [ A B ] * Q  = [ C1 B1 ]
!>          [ 0 C ]        [ 0  A1 ]
!>
!> where C1 is similar to C and A1 is similar to A.  Both A and C are
!> assumed to be in standard form (equal diagonal entries and
!> offdiagonal with differing signs) and A1 and C1 are returned with the
!> same properties.
!>
!> The test code verifies these last assertions, as well as that
!> the residual in the above equation is small.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is DOUBLE PRECISION
!>          Value of the largest test ratio.
!> \endverbatim
!>
!> \param[out] LMAX
!> \verbatim
!>          LMAX is INTEGER
!>          Example number where largest test ratio achieved.
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (2)
!>          NINFO(J) is the number of examples where INFO=J occurred.
!> \endverbatim
!>
!> \param[out] KNT
!> \verbatim
!>          KNT is INTEGER
!>          Total number of examples tested.
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
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DGET34( RMAX, LMAX, NINFO, KNT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX
   DOUBLE PRECISION   RMAX
!     ..
!     .. Array Arguments ..
   INTEGER            NINFO( 2 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, HALF, ONE
   PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
   DOUBLE PRECISION   TWO, THREE
   PARAMETER          ( TWO = 2.0D0, THREE = 3.0D0 )
   INTEGER            LWORK
   PARAMETER          ( LWORK = 32 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IA, IA11, IA12, IA21, IA22, IAM, IB, IC, &
                      IC11, IC12, IC21, IC22, ICM, INFO, J
   DOUBLE PRECISION   BIGNUM, EPS, RES, SMLNUM, TNRM
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   Q( 4, 4 ), RESULT( 2 ), T( 4, 4 ), T1( 4, 4 ), &
                      VAL( 9 ), VM( 2 ), WORK( LWORK )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DHST01, DLAEXC
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = DLAMCH( 'P' )
   SMLNUM = DLAMCH( 'S' ) / EPS
   BIGNUM = ONE / SMLNUM
!
!     Set up test case parameters
!
   VAL( 1 ) = ZERO
   VAL( 2 ) = SQRT( SMLNUM )
   VAL( 3 ) = ONE
   VAL( 4 ) = TWO
   VAL( 5 ) = SQRT( BIGNUM )
   VAL( 6 ) = -SQRT( SMLNUM )
   VAL( 7 ) = -ONE
   VAL( 8 ) = -TWO
   VAL( 9 ) = -SQRT( BIGNUM )
   VM( 1 ) = ONE
   VM( 2 ) = ONE + TWO*EPS
   CALL DCOPY( 16, VAL( 4 ), 0, T( 1, 1 ), 1 )
!
   NINFO( 1 ) = 0
   NINFO( 2 ) = 0
   KNT = 0
   LMAX = 0
   RMAX = ZERO
!
!     Begin test loop
!
   DO IA = 1, 9
      DO IAM = 1, 2
         DO IB = 1, 9
            DO IC = 1, 9
               T( 1, 1 ) = VAL( IA )*VM( IAM )
               T( 2, 2 ) = VAL( IC )
               T( 1, 2 ) = VAL( IB )
               T( 2, 1 ) = ZERO
               TNRM = MAX( ABS( T( 1, 1 ) ), ABS( T( 2, 2 ) ), &
                      ABS( T( 1, 2 ) ) )
               CALL DCOPY( 16, T, 1, T1, 1 )
               CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
               CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
               CALL DLAEXC( .TRUE., 2, T, 4, Q, 4, 1, 1, 1, WORK, &
                            INFO )
               IF( INFO /= 0 ) &
                  NINFO( INFO ) = NINFO( INFO ) + 1
               CALL DHST01( 2, 1, 2, T1, 4, T, 4, Q, 4, WORK, LWORK, &
                            RESULT )
               RES = RESULT( 1 ) + RESULT( 2 )
               IF( INFO /= 0 ) &
                  RES = RES + ONE / EPS
               IF( T( 1, 1 ) /= T1( 2, 2 ) ) &
                  RES = RES + ONE / EPS
               IF( T( 2, 2 ) /= T1( 1, 1 ) ) &
                  RES = RES + ONE / EPS
               IF( T( 2, 1 ) /= ZERO ) &
                  RES = RES + ONE / EPS
               KNT = KNT + 1
               IF( RES > RMAX ) THEN
                  LMAX = KNT
                  RMAX = RES
               END IF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
!
   DO IA = 1, 5
      DO IAM = 1, 2
         DO IB = 1, 5
            DO IC11 = 1, 5
               DO IC12 = 2, 5
                  DO IC21 = 2, 4
                     DO IC22 = -1, 1, 2
                        T( 1, 1 ) = VAL( IA )*VM( IAM )
                        T( 1, 2 ) = VAL( IB )
                        T( 1, 3 ) = -TWO*VAL( IB )
                        T( 2, 1 ) = ZERO
                        T( 2, 2 ) = VAL( IC11 )
                        T( 2, 3 ) = VAL( IC12 )
                        T( 3, 1 ) = ZERO
                        T( 3, 2 ) = -VAL( IC21 )
                        T( 3, 3 ) = VAL( IC11 )*DBLE( IC22 )
                        TNRM = MAX( ABS( T( 1, 1 ) ), &
                               ABS( T( 1, 2 ) ), ABS( T( 1, 3 ) ), &
                               ABS( T( 2, 2 ) ), ABS( T( 2, 3 ) ), &
                               ABS( T( 3, 2 ) ), ABS( T( 3, 3 ) ) )
                        CALL DCOPY( 16, T, 1, T1, 1 )
                        CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
                        CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
                        CALL DLAEXC( .TRUE., 3, T, 4, Q, 4, 1, 1, 2, &
                                     WORK, INFO )
                        IF( INFO /= 0 ) &
                           NINFO( INFO ) = NINFO( INFO ) + 1
                        CALL DHST01( 3, 1, 3, T1, 4, T, 4, Q, 4, &
                                     WORK, LWORK, RESULT )
                        RES = RESULT( 1 ) + RESULT( 2 )
                        IF( INFO == 0 ) THEN
                           IF( T1( 1, 1 ) /= T( 3, 3 ) ) &
                              RES = RES + ONE / EPS
                           IF( T( 3, 1 ) /= ZERO ) &
                              RES = RES + ONE / EPS
                           IF( T( 3, 2 ) /= ZERO ) &
                              RES = RES + ONE / EPS
                           IF( T( 2, 1 ) /= 0 .AND. &
                               ( T( 1, 1 ) /= T( 2, &
                               2 ) .OR. SIGN( ONE, T( 1, &
                               2 ) ) == SIGN( ONE, T( 2, 1 ) ) ) ) &
                               RES = RES + ONE / EPS
                        END IF
                        KNT = KNT + 1
                        IF( RES > RMAX ) THEN
                           LMAX = KNT
                           RMAX = RES
                        END IF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
!
   DO IA11 = 1, 5
      DO IA12 = 2, 5
         DO IA21 = 2, 4
            DO IA22 = -1, 1, 2
               DO ICM = 1, 2
                  DO IB = 1, 5
                     DO IC = 1, 5
                        T( 1, 1 ) = VAL( IA11 )
                        T( 1, 2 ) = VAL( IA12 )
                        T( 1, 3 ) = -TWO*VAL( IB )
                        T( 2, 1 ) = -VAL( IA21 )
                        T( 2, 2 ) = VAL( IA11 )*DBLE( IA22 )
                        T( 2, 3 ) = VAL( IB )
                        T( 3, 1 ) = ZERO
                        T( 3, 2 ) = ZERO
                        T( 3, 3 ) = VAL( IC )*VM( ICM )
                        TNRM = MAX( ABS( T( 1, 1 ) ), &
                               ABS( T( 1, 2 ) ), ABS( T( 1, 3 ) ), &
                               ABS( T( 2, 2 ) ), ABS( T( 2, 3 ) ), &
                               ABS( T( 3, 2 ) ), ABS( T( 3, 3 ) ) )
                        CALL DCOPY( 16, T, 1, T1, 1 )
                        CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
                        CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
                        CALL DLAEXC( .TRUE., 3, T, 4, Q, 4, 1, 2, 1, &
                                     WORK, INFO )
                        IF( INFO /= 0 ) &
                           NINFO( INFO ) = NINFO( INFO ) + 1
                        CALL DHST01( 3, 1, 3, T1, 4, T, 4, Q, 4, &
                                     WORK, LWORK, RESULT )
                        RES = RESULT( 1 ) + RESULT( 2 )
                        IF( INFO == 0 ) THEN
                           IF( T1( 3, 3 ) /= T( 1, 1 ) ) &
                              RES = RES + ONE / EPS
                           IF( T( 2, 1 ) /= ZERO ) &
                              RES = RES + ONE / EPS
                           IF( T( 3, 1 ) /= ZERO ) &
                              RES = RES + ONE / EPS
                           IF( T( 3, 2 ) /= 0 .AND. &
                               ( T( 2, 2 ) /= T( 3, &
                               3 ) .OR. SIGN( ONE, T( 2, &
                               3 ) ) == SIGN( ONE, T( 3, 2 ) ) ) ) &
                               RES = RES + ONE / EPS
                        END IF
                        KNT = KNT + 1
                        IF( RES > RMAX ) THEN
                           LMAX = KNT
                           RMAX = RES
                        END IF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
   DO IA11 = 1, 5
      DO IA12 = 2, 5
         DO IA21 = 2, 4
            DO IA22 = -1, 1, 2
               DO IB = 1, 5
                  DO IC11 = 3, 4
                     DO IC12 = 3, 4
                        DO IC21 = 3, 4
                           DO IC22 = -1, 1, 2
                              DO ICM = 5, 7
                                 IAM = 1
                                 T( 1, 1 ) = VAL( IA11 )*VM( IAM )
                                 T( 1, 2 ) = VAL( IA12 )*VM( IAM )
                                 T( 1, 3 ) = -TWO*VAL( IB )
                                 T( 1, 4 ) = HALF*VAL( IB )
                                 T( 2, 1 ) = -T( 1, 2 )*VAL( IA21 )
                                 T( 2, 2 ) = VAL( IA11 )* &
                                             DBLE( IA22 )*VM( IAM )
                                 T( 2, 3 ) = VAL( IB )
                                 T( 2, 4 ) = THREE*VAL( IB )
                                 T( 3, 1 ) = ZERO
                                 T( 3, 2 ) = ZERO
                                 T( 3, 3 ) = VAL( IC11 )* &
                                             ABS( VAL( ICM ) )
                                 T( 3, 4 ) = VAL( IC12 )* &
                                             ABS( VAL( ICM ) )
                                 T( 4, 1 ) = ZERO
                                 T( 4, 2 ) = ZERO
                                 T( 4, 3 ) = -T( 3, 4 )*VAL( IC21 )* &
                                             ABS( VAL( ICM ) )
                                 T( 4, 4 ) = VAL( IC11 )* &
                                             DBLE( IC22 )* &
                                             ABS( VAL( ICM ) )
                                 TNRM = ZERO
                                 DO I = 1, 4
                                    DO J = 1, 4
                                       TNRM = MAX( TNRM, &
                                              ABS( T( I, J ) ) )
                                       ENDDO
                                    ENDDO
                                 CALL DCOPY( 16, T, 1, T1, 1 )
                                 CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
                                 CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
                                 CALL DLAEXC( .TRUE., 4, T, 4, Q, 4, &
                                              1, 2, 2, WORK, INFO )
                                 IF( INFO /= 0 ) &
                                    NINFO( INFO ) = NINFO( INFO ) + 1
                                 CALL DHST01( 4, 1, 4, T1, 4, T, 4, &
                                              Q, 4, WORK, LWORK, &
                                              RESULT )
                                 RES = RESULT( 1 ) + RESULT( 2 )
                                 IF( INFO == 0 ) THEN
                                    IF( T( 3, 1 ) /= ZERO ) &
                                       RES = RES + ONE / EPS
                                    IF( T( 4, 1 ) /= ZERO ) &
                                       RES = RES + ONE / EPS
                                    IF( T( 3, 2 ) /= ZERO ) &
                                       RES = RES + ONE / EPS
                                    IF( T( 4, 2 ) /= ZERO ) &
                                       RES = RES + ONE / EPS
                                    IF( T( 2, 1 ) /= 0 .AND. &
                                        ( T( 1, 1 ) /= T( 2, &
                                        2 ) .OR. SIGN( ONE, T( 1, &
                                        2 ) ) == SIGN( ONE, T( 2, &
                                        1 ) ) ) )RES = RES + &
                                        ONE / EPS
                                    IF( T( 4, 3 ) /= 0 .AND. &
                                        ( T( 3, 3 ) /= T( 4, &
                                        4 ) .OR. SIGN( ONE, T( 3, &
                                        4 ) ) == SIGN( ONE, T( 4, &
                                        3 ) ) ) )RES = RES + &
                                        ONE / EPS
                                 END IF
                                 KNT = KNT + 1
                                 IF( RES > RMAX ) THEN
                                    LMAX = KNT
                                    RMAX = RES
                                 END IF
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
   RETURN
!
!     End of DGET34
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
