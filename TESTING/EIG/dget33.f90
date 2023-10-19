!> \brief \b DGET33
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET33( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NINFO
!       DOUBLE PRECISION   RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET33 tests DLANV2, a routine for putting 2 by 2 blocks into
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DGET33( RMAX, LMAX, NINFO, KNT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX, NINFO
   DOUBLE PRECISION   RMAX
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I1, I2, I3, I4, IM1, IM2, IM3, IM4, J1, J2, J3
   DOUBLE PRECISION   BIGNUM, CS, EPS, RES, SMLNUM, SN, SUM, TNRM, &
                      WI1, WI2, WR1, WR2
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   Q( 2, 2 ), T( 2, 2 ), T1( 2, 2 ), T2( 2, 2 ), &
                      VAL( 4 ), VM( 3 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLANV2
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = DLAMCH( 'P' )
   SMLNUM = DLAMCH( 'S' ) / EPS
   BIGNUM = 1.0D0 / SMLNUM
!
!     Set up test case parameters
!
   VAL( 1 ) = 1.0D0
   VAL( 2 ) = 1.0D0 + 2.0D0*EPS
   VAL( 3 ) = 2.0D0
   VAL( 4 ) = 2.0D0 - 4.0D0*EPS
   VM( 1 ) = SMLNUM
   VM( 2 ) = 1.0D0
   VM( 3 ) = BIGNUM
!
   KNT = 0
   NINFO = 0
   LMAX = 0
   RMAX = 0.0D0
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
                           T1( 1:2, 1:2 ) = T( 1:2, 1:2 )
                           Q( 1, 1 ) = 1.0D0
                           Q( 1, 2 ) = 0.0D0
                           Q( 2, 1 ) = 0.0D0
                           Q( 2, 2 ) = 1.0D0
!
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL DLANV2( T( 1, 1 ), T( 1, 2 ), &
                                        T( 2, 1 ), T( 2, 2 ), WR1, &
                                        WI1, WR2, WI2, CS, SN )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : DLANV2 : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                           DO J1 = 1, 2
                              RES = Q( J1, 1 )*CS + Q( J1, 2 )*SN
                              Q( J1, 2 ) = -Q( J1, 1 )*SN + &
                                           Q( J1, 2 )*CS
                              Q( J1, 1 ) = RES
                           ENDDO
!
                           RES = 0.0D0
                           RES = RES + ABS( Q( 1, 1 )**2+ &
                                 Q( 1, 2 )**2-1.0D0 ) / EPS
                           RES = RES + ABS( Q( 2, 2 )**2+ &
                                 Q( 2, 1 )**2-1.0D0 ) / EPS
                           RES = RES + ABS( Q( 1, 1 )*Q( 2, 1 )+ &
                                 Q( 1, 2 )*Q( 2, 2 ) ) / EPS
                           DO J1 = 1, 2
                              DO J2 = 1, 2
                                 T2( J1, J2 ) = 0.0D0
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
                                    SUM = SUM - Q( J3, J1 )* T2( J3, J2 )
                                 ENDDO
                                 RES = RES + ABS( SUM ) / EPS / TNRM
                              ENDDO
                           ENDDO
                           IF( T( 2, 1 ) /= 0.0D0 .AND. ( T( 1, 1 ) /= T( 2, 2 ) .OR. SIGN( 1.0D0, T( 1, &
                               2 ) )*SIGN( 1.0D0, T( 2, 1 ) ) > 0.0D0 ) )RES = RES + 1.0D0 / EPS
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
!     End of DGET33
!
END




