!> \brief \b CSYL01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYL01( THRESH, NFAIL, RMAX, NINFO, KNT )
!
!     .. Scalar Arguments ..
!     INTEGER            KNT
!     REAL               THRESH
!     ..
!     .. Array Arguments ..
!     INTEGER            NFAIL( 3 ), NINFO( 2 )
!     REAL               RMAX( 2 )
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYL01 tests CTRSYL and CTRSYL3, routines for solving the Sylvester matrix
!> equation
!>
!>    op(A)*X + ISGN*X*op(B) = scale*C,
!>
!> where op(A) and op(B) are both upper triangular form, op() represents an
!> optional conjugate transpose, and ISGN can be -1 or +1. Scale is an output
!> less than or equal to 1, chosen to avoid overflow in X.
!>
!> The test code verifies that the following residual does not exceed
!> the provided threshold:
!>
!>    norm(op(A)*X + ISGN*X*op(B) - scale*C) /
!>        (EPS*max(norm(A),norm(B))*norm(X))
!>
!> This routine complements CGET35 by testing with larger,
!> random matrices, of which some require rescaling of X to avoid overflow.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          A test will count as "failed" if the residual, computed as
!>          described above, exceeds THRESH.
!> \endverbatim
!>
!> \param[out] NFAIL
!> \verbatim
!>          NFAIL is INTEGER array, dimension (3)
!>          NFAIL(1) = No. of times residual CTRSYL exceeds threshold THRESH
!>          NFAIL(2) = No. of times residual CTRSYL3 exceeds threshold THRESH
!>          NFAIL(3) = No. of times CTRSYL3 and CTRSYL deviate
!> \endverbatim
!>
!> \param[out] RMAX
!> \verbatim
!>          RMAX is DOUBLE PRECISION array, dimension (2)
!>          RMAX(1) = Value of the largest test ratio of CTRSYL
!>          RMAX(2) = Value of the largest test ratio of CTRSYL3
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (2)
!>          NINFO(1) = No. of times CTRSYL where INFO is nonzero
!>          NINFO(2) = No. of times CTRSYL3 where INFO is nonzero
!> \endverbatim
!>
!> \param[out] KNT
!> \verbatim
!>          KNT is INTEGER
!>          Total number of examples tested.
!> \endverbatim

!
!  -- LAPACK test routine --
   SUBROUTINE CSYL01( THRESH, NFAIL, RMAX, NINFO, KNT )
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            NFAIL( 3 ), NINFO( 2 )
   REAL               RMAX( 2 )
!     ..
!
!  =====================================================================
!     ..
!     .. Parameters ..
   INTEGER            MAXM, MAXN, LDSWORK
   PARAMETER          ( MAXM = 101, MAXN = 138, LDSWORK = 18 )
!     ..
!     .. Local Scalars ..
   CHARACTER          TRANA, TRANB
   INTEGER            I, INFO, IINFO, ISGN, ITRANA, ITRANB, J, KLA, &
                      KUA, KLB, KUB, M, N
   REAL               ANRM, BNRM, BIGNUM, EPS, RES, RES1, &
                      SCALE, SCALE3, SMLNUM, TNRM, XNRM
   COMPLEX            RMUL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   COMPLEX            DUML( MAXM ), DUMR( MAXN ), &
                      D( MAX( MAXM, MAXN ) )
   REAL               DUM( MAXN ), VM( 2 )
   INTEGER            ISEED( 4 ), IWORK( MAXM + MAXN + 2 )
!     ..
!     .. Allocatable Arrays ..
   INTEGER            AllocateStatus
   COMPLEX, DIMENSION(:,:), ALLOCATABLE :: A, B, C, CC, X
   REAL,    DIMENSION(:,:), ALLOCATABLE :: SWORK
!     ..
!     .. External Functions ..
   LOGICAL            SISNAN
   REAL               SLAMCH, CLANGE
   EXTERNAL           SISNAN, SLAMCH, CLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLATMR, CLACPY, CGEMM, CTRSYL, CTRSYL3
!     ..
!     .. Allocate memory dynamically ..
   ALLOCATE ( A( MAXM, MAXM ), STAT = AllocateStatus )
   IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
   ALLOCATE ( B( MAXN, MAXN ), STAT = AllocateStatus )
   IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
   ALLOCATE ( C( MAXM, MAXN ), STAT = AllocateStatus )
   IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
   ALLOCATE ( CC( MAXM, MAXN ), STAT = AllocateStatus )
   IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
   ALLOCATE ( X( MAXM, MAXN ), STAT = AllocateStatus )
   IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
   ALLOCATE ( SWORK( LDSWORK, 54 ), STAT = AllocateStatus )
   IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = SLAMCH( 'P' )
   SMLNUM = SLAMCH( 'S' ) / EPS
   BIGNUM = 1.0E+0 / SMLNUM
!
!     Expect INFO = 0
   VM( 1 ) = 1.0E+0
!     Expect INFO = 1
   VM( 2 ) = 0.5E+0
!
!     Begin test loop
!
   NINFO( 1:2 ) = 0
   NFAIL( 1:3 ) = 0
   RMAX( 1:2 ) = 0.0E+0
   KNT = 0
   ISEED( 1:4 ) = 1
   SCALE = 1.0E+0
   SCALE3 = 1.0E+0
   DO J = 1, 2
      DO ISGN = -1, 1, 2
!           Reset seed (overwritten by LATMR)
         ISEED( 1:4 ) = 1
         DO M = 32, MAXM, 23
            KLA = 0
            KUA = M - 1
            CALL CLATMR( M, M, 'S', ISEED, 'N', D, &
                         6, 1.0E+0, (1.0E+0,0.0E+0), 'T', 'N', &
                         DUML, 1, 1.0E+0, DUMR, 1, 1.0E+0, &
                         'N', IWORK, KLA, KUA, 0.0E+0, &
                         1.0E+0, 'NO', A, MAXM, IWORK, &
                         IINFO )
            FORALL (I = 1:M) A( I, I ) = A( I, I ) * VM( J )
            ANRM = CLANGE( 'M', M, M, A, MAXM, DUM )
            DO N = 51, MAXN, 29
               KLB = 0
               KUB = N - 1
               CALL CLATMR( N, N, 'S', ISEED, 'N', D, &
                            6, 1.0E+0, (1.0E+0,0.0E+0), 'T', 'N', &
                            DUML, 1, 1.0E+0, DUMR, 1, 1.0E+0, &
                            'N', IWORK, KLB, KUB, 0.0E+0, &
                            1.0E+0, 'NO', B, MAXN, IWORK, &
                            IINFO )
               FORALL (I = 1:N) B( I, I ) = B( I, I ) * VM ( J )
               BNRM = CLANGE( 'M', N, N, B, MAXN, DUM )
               TNRM = MAX( ANRM, BNRM )
               CALL CLATMR( M, N, 'S', ISEED, 'N', D, &
                            6, 1.0E+0, (1.0E+0,0.0E+0), 'T', 'N', &
                            DUML, 1, 1.0E+0, DUMR, 1, 1.0E+0, &
                            'N', IWORK, M, N, 0.0E+0, 1.0E+0, &
                            'NO', C, MAXM, IWORK, IINFO )
               DO ITRANA = 1, 2
                  IF( ITRANA == 1 ) TRANA = 'N'
                  IF( ITRANA == 2 ) TRANA = 'C'
                  DO ITRANB = 1, 2
                     IF( ITRANB == 1 ) TRANB = 'N'
                     IF( ITRANB == 2 ) TRANB = 'C'
                     KNT = KNT + 1
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'All', M, N, C, MAXM, X, MAXM)
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'All', M, N, C, MAXM, CC, MAXM)
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CTRSYL( TRANA, TRANB, ISGN, M, N, &
                                  A, MAXM, B, MAXN, X, MAXM, &
                                  SCALE, IINFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CTRSYL : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     IF( IINFO /= 0 ) NINFO( 1 ) = NINFO( 1 ) + 1
                     XNRM = CLANGE( 'M', M, N, X, MAXM, DUM )
                     RMUL = (1.0E+0,0.0E+0)
                     IF( XNRM > 1.0E+0 .AND. TNRM > 1.0E+0 ) THEN
                        IF( XNRM > BIGNUM / TNRM ) THEN
                           RMUL = (1.0E+0,0.0E+0) / MAX( XNRM, TNRM )
                        END IF
                     END IF
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CGEMM( TRANA, 'N', M, N, M, RMUL, &
                                 A, MAXM, X, MAXM, -SCALE*RMUL, CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CGEMM( 'N', TRANB, M, N, N, &
                                 REAL( ISGN )*RMUL, X, MAXM, B, &
                                 MAXN, (1.0E+0,0.0E+0), CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     RES1 = CLANGE( 'M', M, N, CC, MAXM, DUM )
                     RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, &
                           ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )
                     IF( RES > THRESH ) NFAIL( 1 ) = NFAIL( 1 ) + 1
                     IF( RES > RMAX( 1 ) ) RMAX( 1 ) = RES
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'All', M, N, C, MAXM, X, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'All', M, N, C, MAXM, CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CTRSYL3( TRANA, TRANB, ISGN, M, N, &
                                   A, MAXM, B, MAXN, X, MAXM, &
                                   SCALE3, SWORK, LDSWORK, INFO)
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CTRSYL3 : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     IF( INFO /= 0 ) NINFO( 2 ) = NINFO( 2 ) + 1
                     XNRM = CLANGE( 'M', M, N, X, MAXM, DUM )
                     RMUL = (1.0E+0,0.0E+0)
                     IF( XNRM > 1.0E+0 .AND. TNRM > 1.0E+0 ) THEN
                        IF( XNRM > BIGNUM / TNRM ) THEN
                           RMUL = (1.0E+0,0.0E+0) / MAX( XNRM, TNRM )
                        END IF
                     END IF
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CGEMM( TRANA, 'N', M, N, M, RMUL, &
                                 A, MAXM, X, MAXM, -SCALE3*RMUL, &
                                 CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CGEMM( 'N', TRANB, M, N, N, &
                                 REAL( ISGN )*RMUL, X, MAXM, B, &
                                 MAXN, (1.0E+0,0.0E+0), CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     RES1 = CLANGE( 'M', M, N, CC, MAXM, DUM )
                     RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, &
                                ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )
!                       Verify that TRSYL3 only flushes if TRSYL flushes (but
!                       there may be cases where TRSYL3 avoid flushing).
                     IF( SCALE3 == 0.0E+0 .AND. SCALE > 0.0E+0 .OR. &
                         IINFO /= INFO ) THEN
                        NFAIL( 3 ) = NFAIL( 3 ) + 1
                     END IF
                     IF( RES > THRESH .OR. SISNAN( RES ) ) &
                        NFAIL( 2 ) = NFAIL( 2 ) + 1
                     IF( RES > RMAX( 2 ) ) &
                        RMAX( 2 ) = RES
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
!
   DEALLOCATE (A, STAT = AllocateStatus)
   DEALLOCATE (B, STAT = AllocateStatus)
   DEALLOCATE (C, STAT = AllocateStatus)
   DEALLOCATE (CC, STAT = AllocateStatus)
   DEALLOCATE (X, STAT = AllocateStatus)
   DEALLOCATE (SWORK, STAT = AllocateStatus)
!
   RETURN
!
!     End of CSYL01
!
END


