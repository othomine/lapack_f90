!> \brief \b SGET33
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET33( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NINFO
!       REAL               RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET33 tests SLANV2, a routine for putting 2 by 2 blocks into
!> standard form.  In other words, it computes a two by two rotation
!> [[C,S];[-S,C]] where in
!>
!>    [ C S ][T(1,1) T(1,2)][ C -S ] = [ T11 T12 ]
!>    [-S C ][T(2,1) T(2,2)][ S  C ]   [ T21 T22 ]
!>
!> either
!>    1) T21=0 (real eigenvalues), or
!>    2) T11=T22 and T21*T12<0 (complex conjugate eigenvalues).
!> We also  verify that the residual is small.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is REAL
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
!>          NINFO is INTEGER
!>          Number of examples returned with INFO  /=  0.
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SGET33( RMAX, LMAX, NINFO, KNT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX, NINFO
   REAL               RMAX
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
   REAL               TWO, FOUR
   PARAMETER          ( TWO = 2.0E0, FOUR = 4.0E0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I1, I2, I3, I4, IM1, IM2, IM3, IM4, J1, J2, J3
   REAL               BIGNUM, CS, EPS, RES, SMLNUM, SN, SUM, TNRM, &
                      WI1, WI2, WR1, WR2
!     ..
!     .. Local Arrays ..
   REAL               Q( 2, 2 ), T( 2, 2 ), T1( 2, 2 ), T2( 2, 2 ), &
                      VAL( 4 ), VM( 3 )
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SLANV2
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, SIGN
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = SLAMCH( 'P' )
   SMLNUM = SLAMCH( 'S' ) / EPS
   BIGNUM = ONE / SMLNUM
!
!     Set up test case parameters
!
   VAL( 1 ) = ONE
   VAL( 2 ) = ONE + TWO*EPS
   VAL( 3 ) = TWO
   VAL( 4 ) = TWO - FOUR*EPS
   VM( 1 ) = SMLNUM
   VM( 2 ) = ONE
   VM( 3 ) = BIGNUM
!
   KNT = 0
   NINFO = 0
   LMAX = 0
   RMAX = ZERO
!
!     Begin test loop
!
   DO I1 = 1, 4
      DO I2 = 1, 4
         DO I3 = 1, 4
            DO I4 = 1, 4
               DO IM1 = 1, 3
                  DO IM2 = 1, 3
                     DO IM3 = 1, 3
                        DO IM4 = 1, 3
                           T( 1, 1 ) = VAL( I1 )*VM( IM1 )
                           T( 1, 2 ) = VAL( I2 )*VM( IM2 )
                           T( 2, 1 ) = -VAL( I3 )*VM( IM3 )
                           T( 2, 2 ) = VAL( I4 )*VM( IM4 )
                           TNRM = MAX( ABS( T( 1, 1 ) ), &
                                  ABS( T( 1, 2 ) ), ABS( T( 2, 1 ) ), &
                                  ABS( T( 2, 2 ) ) )
                           T1( 1, 1 ) = T( 1, 1 )
                           T1( 1, 2 ) = T( 1, 2 )
                           T1( 2, 1 ) = T( 2, 1 )
                           T1( 2, 2 ) = T( 2, 2 )
                           Q( 1, 1 ) = ONE
                           Q( 1, 2 ) = ZERO
                           Q( 2, 1 ) = ZERO
                           Q( 2, 2 ) = ONE
!
                           CALL SLANV2( T( 1, 1 ), T( 1, 2 ), &
                                        T( 2, 1 ), T( 2, 2 ), WR1, &
                                        WI1, WR2, WI2, CS, SN )
                           DO J1 = 1, 2
                              RES = Q( J1, 1 )*CS + Q( J1, 2 )*SN
                              Q( J1, 2 ) = -Q( J1, 1 )*SN + &
                                           Q( J1, 2 )*CS
                              Q( J1, 1 ) = RES
                           ENDDO
!
                           RES = ZERO
                           RES = RES + ABS( Q( 1, 1 )**2+ &
                                 Q( 1, 2 )**2-ONE ) / EPS
                           RES = RES + ABS( Q( 2, 2 )**2+ &
                                 Q( 2, 1 )**2-ONE ) / EPS
                           RES = RES + ABS( Q( 1, 1 )*Q( 2, 1 )+ &
                                 Q( 1, 2 )*Q( 2, 2 ) ) / EPS
                           DO J1 = 1, 2
                              DO J2 = 1, 2
                                 T2( J1, J2 ) = ZERO
                                 DO J3 = 1, 2
                                    T2( J1, J2 ) = T2( J1, J2 ) + &
                                                   T1( J1, J3 )* &
                                                   Q( J3, J2 )
                                 ENDDO
                              ENDDO
                           ENDDO
                           DO J1 = 1, 2
                              DO J2 = 1, 2
                                 SUM = T( J1, J2 )
                                 DO J3 = 1, 2
                                    SUM = SUM - Q( J3, J1 )* &
                                          T2( J3, J2 )
                                 ENDDO
                                 RES = RES + ABS( SUM ) / EPS / TNRM
                              ENDDO
                           ENDDO
                           IF( T( 2, 1 ) /= ZERO .AND. &
                               ( T( 1, 1 ) /= T( 2, &
                               2 ) .OR. SIGN( ONE, T( 1, &
                               2 ) )*SIGN( ONE, T( 2, &
                               1 ) ) > ZERO ) )RES = RES + ONE / EPS
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
!
   RETURN
!
!     End of SGET33
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
