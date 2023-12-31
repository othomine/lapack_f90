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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
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
   INTEGER            LWORK
   PARAMETER          ( LWORK = 32 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IA, IA11, IA12, IA21, IA22, IAM, IB, IC, &
                      IC11, IC12, IC21, IC22, ICM, INFO, J
   DOUBLE PRECISION   BIGNUM, EPS, RES, SMLNUM, TNRM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
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
   VAL( 1 ) = 0.0D0
   VAL( 2 ) = SQRT( SMLNUM )
   VAL( 3 ) = 1.0D0
   VAL( 4 ) = 2.0D0
   VAL( 5 ) = SQRT( BIGNUM )
   VAL( 6 ) = -SQRT( SMLNUM )
   VAL( 7 ) = -1.0D0
   VAL( 8 ) = -2.0D0
   VAL( 9 ) = -SQRT( BIGNUM )
   VM( 1 ) = 1.0D0
   VM( 2 ) = 1.0D0 + 2.0D0*EPS
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DCOPY( 16, VAL( 4 ), 0, T( 1, 1 ), 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   NINFO( 1 ) = 0
   NINFO( 2 ) = 0
   KNT = 0
   LMAX = 0
   RMAX = 0.0D0
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
               T( 2, 1 ) = 0.0D0
               TNRM = MAX( ABS( T( 1, 1 ) ), ABS( T( 2, 2 ) ), &
                      ABS( T( 1, 2 ) ) )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( 16, T, 1, T1, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLAEXC( .TRUE., 2, T, 4, Q, 4, 1, 1, 1, WORK, &
                            INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLAEXC : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( INFO /= 0 ) &
                  NINFO( INFO ) = NINFO( INFO ) + 1
               CALL DHST01( 2, 1, 2, T1, 4, T, 4, Q, 4, WORK, LWORK, &
                            RESULT )
               RES = RESULT( 1 ) + RESULT( 2 )
               IF( INFO /= 0 ) &
                  RES = RES + 1.0D0 / EPS
               IF( T( 1, 1 ) /= T1( 2, 2 ) ) &
                  RES = RES + 1.0D0 / EPS
               IF( T( 2, 2 ) /= T1( 1, 1 ) ) &
                  RES = RES + 1.0D0 / EPS
               IF( T( 2, 1 ) /= 0.0D0 ) &
                  RES = RES + 1.0D0 / EPS
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
                        T( 1, 3 ) = -2.0D0*VAL( IB )
                        T( 2, 1 ) = 0.0D0
                        T( 2, 2 ) = VAL( IC11 )
                        T( 2, 3 ) = VAL( IC12 )
                        T( 3, 1 ) = 0.0D0
                        T( 3, 2 ) = -VAL( IC21 )
                        T( 3, 3 ) = VAL( IC11 )*DBLE( IC22 )
                        TNRM = MAX( ABS( T( 1, 1 ) ), &
                               ABS( T( 1, 2 ) ), ABS( T( 1, 3 ) ), &
                               ABS( T( 2, 2 ) ), ABS( T( 2, 3 ) ), &
                               ABS( T( 3, 2 ) ), ABS( T( 3, 3 ) ) )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DCOPY( 16, T, 1, T1, 1 )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DLAEXC( .TRUE., 3, T, 4, Q, 4, 1, 1, 2, &
                                     WORK, INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DLAEXC : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        IF( INFO /= 0 ) &
                           NINFO( INFO ) = NINFO( INFO ) + 1
                        CALL DHST01( 3, 1, 3, T1, 4, T, 4, Q, 4, &
                                     WORK, LWORK, RESULT )
                        RES = RESULT( 1 ) + RESULT( 2 )
                        IF( INFO == 0 ) THEN
                           IF( T1( 1, 1 ) /= T( 3, 3 ) ) &
                              RES = RES + 1.0D0 / EPS
                           IF( T( 3, 1 ) /= 0.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           IF( T( 3, 2 ) /= 0.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           IF( T( 2, 1 ) /= 0 .AND. &
                               ( T( 1, 1 ) /= T( 2, &
                               2 ) .OR. SIGN( 1.0D0, T( 1, &
                               2 ) ) == SIGN( 1.0D0, T( 2, 1 ) ) ) ) &
                               RES = RES + 1.0D0 / EPS
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
                        T( 1, 3 ) = -2.0D0*VAL( IB )
                        T( 2, 1 ) = -VAL( IA21 )
                        T( 2, 2 ) = VAL( IA11 )*DBLE( IA22 )
                        T( 2, 3 ) = VAL( IB )
                        T( 3, 1 ) = 0.0D0
                        T( 3, 2 ) = 0.0D0
                        T( 3, 3 ) = VAL( IC )*VM( ICM )
                        TNRM = MAX( ABS( T( 1, 1 ) ), &
                               ABS( T( 1, 2 ) ), ABS( T( 1, 3 ) ), &
                               ABS( T( 2, 2 ) ), ABS( T( 2, 3 ) ), &
                               ABS( T( 3, 2 ) ), ABS( T( 3, 3 ) ) )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DCOPY( 16, T, 1, T1, 1 )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DLAEXC( .TRUE., 3, T, 4, Q, 4, 1, 2, 1, &
                                     WORK, INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DLAEXC : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        IF( INFO /= 0 ) &
                           NINFO( INFO ) = NINFO( INFO ) + 1
                        CALL DHST01( 3, 1, 3, T1, 4, T, 4, Q, 4, &
                                     WORK, LWORK, RESULT )
                        RES = RESULT( 1 ) + RESULT( 2 )
                        IF( INFO == 0 ) THEN
                           IF( T1( 3, 3 ) /= T( 1, 1 ) ) &
                              RES = RES + 1.0D0 / EPS
                           IF( T( 2, 1 ) /= 0.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           IF( T( 3, 1 ) /= 0.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           IF( T( 3, 2 ) /= 0 .AND. &
                               ( T( 2, 2 ) /= T( 3, &
                               3 ) .OR. SIGN( 1.0D0, T( 2, &
                               3 ) ) == SIGN( 1.0D0, T( 3, 2 ) ) ) ) &
                               RES = RES + 1.0D0 / EPS
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
                                 T( 1, 3 ) = -2.0D0*VAL( IB )
                                 T( 1, 4 ) = 0.5D0*VAL( IB )
                                 T( 2, 1 ) = -T( 1, 2 )*VAL( IA21 )
                                 T( 2, 2 ) = VAL( IA11 )* &
                                             DBLE( IA22 )*VM( IAM )
                                 T( 2, 3 ) = VAL( IB )
                                 T( 2, 4 ) = 3.0D0*VAL( IB )
                                 T( 3, 1 ) = 0.0D0
                                 T( 3, 2 ) = 0.0D0
                                 T( 3, 3 ) = VAL( IC11 )* &
                                             ABS( VAL( ICM ) )
                                 T( 3, 4 ) = VAL( IC12 )* &
                                             ABS( VAL( ICM ) )
                                 T( 4, 1 ) = 0.0D0
                                 T( 4, 2 ) = 0.0D0
                                 T( 4, 3 ) = -T( 3, 4 )*VAL( IC21 )* &
                                             ABS( VAL( ICM ) )
                                 T( 4, 4 ) = VAL( IC11 )* &
                                             DBLE( IC22 )* &
                                             ABS( VAL( ICM ) )
                                 TNRM = 0.0D0
                                 DO I = 1, 4
                                    DO J = 1, 4
                                       TNRM = MAX( TNRM, &
                                              ABS( T( I, J ) ) )
                                       ENDDO
                                    ENDDO
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                                 CALL DCOPY( 16, T, 1, T1, 1 )
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S2_time)
                                 open(file='results.out', unit=10, position = 'append')
                                 write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                                       real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                                 close(10)
#endif
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                                 CALL DCOPY( 16, VAL( 1 ), 0, Q, 1 )
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S2_time)
                                 open(file='results.out', unit=10, position = 'append')
                                 write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                                       real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                                 close(10)
#endif
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                                 CALL DCOPY( 4, VAL( 3 ), 0, Q, 5 )
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S2_time)
                                 open(file='results.out', unit=10, position = 'append')
                                 write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                                       real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                                 close(10)
#endif
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                                 CALL DLAEXC( .TRUE., 4, T, 4, Q, 4, &
                                              1, 2, 2, WORK, INFO )
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S2_time)
                                 open(file='results.out', unit=10, position = 'append')
                                 write(10,'(A,F16.10,A)') 'Total time : DLAEXC : ',&
                                       real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                                 close(10)
#endif
                                 IF( INFO /= 0 ) &
                                    NINFO( INFO ) = NINFO( INFO ) + 1
                                 CALL DHST01( 4, 1, 4, T1, 4, T, 4, &
                                              Q, 4, WORK, LWORK, &
                                              RESULT )
                                 RES = RESULT( 1 ) + RESULT( 2 )
                                 IF( INFO == 0 ) THEN
                                    IF( T( 3, 1 ) /= 0.0D0 ) &
                                       RES = RES + 1.0D0 / EPS
                                    IF( T( 4, 1 ) /= 0.0D0 ) &
                                       RES = RES + 1.0D0 / EPS
                                    IF( T( 3, 2 ) /= 0.0D0 ) &
                                       RES = RES + 1.0D0 / EPS
                                    IF( T( 4, 2 ) /= 0.0D0 ) &
                                       RES = RES + 1.0D0 / EPS
                                    IF( T( 2, 1 ) /= 0 .AND. &
                                        ( T( 1, 1 ) /= T( 2, &
                                        2 ) .OR. SIGN( 1.0D0, T( 1, &
                                        2 ) ) == SIGN( 1.0D0, T( 2, &
                                        1 ) ) ) )RES = RES + &
                                        1.0D0 / EPS
                                    IF( T( 4, 3 ) /= 0 .AND. &
                                        ( T( 3, 3 ) /= T( 4, &
                                        4 ) .OR. SIGN( 1.0D0, T( 3, &
                                        4 ) ) == SIGN( 1.0D0, T( 4, &
                                        3 ) ) ) )RES = RES + &
                                        1.0D0 / EPS
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




