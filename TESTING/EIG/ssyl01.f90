!> \brief \b SSYL01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYL01( THRESH, NFAIL, RMAX, NINFO, KNT )
!
!      .. Scalar Arguments ..
!      INTEGER            KNT
!      REAL               THRESH
!      ..
!      .. Array Arguments ..
!      INTEGER            NFAIL( 3 ), NINFO( 2 )
!      REAL               RMAX( 2 )
!      ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYL01 tests STRSYL and STRSYL3, routines for solving the Sylvester matrix
!> equation
!>
!>    op(A)*X + ISGN*X*op(B) = scale*C,
!>
!> A and B are assumed to be in Schur canonical form, op() represents an
!> optional transpose, and ISGN can be -1 or +1.  Scale is an output
!> less than or equal to 1, chosen to avoid overflow in X.
!>
!> The test code verifies that the following residual does not exceed
!> the provided threshold:
!>
!>    norm(op(A)*X + ISGN*X*op(B) - scale*C) /
!>        (EPS*max(norm(A),norm(B))*norm(X))
!>
!> This routine complements SGET35 by testing with larger,
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
!>          NFAIL(1) = No. of times residual STRSYL exceeds threshold THRESH
!>          NFAIL(2) = No. of times residual STRSYL3 exceeds threshold THRESH
!>          NFAIL(3) = No. of times STRSYL3 and STRSYL deviate
!> \endverbatim
!>
!> \param[out] RMAX
!> \verbatim
!>          RMAX is REAL, dimension (2)
!>          RMAX(1) = Value of the largest test ratio of STRSYL
!>          RMAX(2) = Value of the largest test ratio of STRSYL3
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (2)
!>          NINFO(1) = No. of times STRSYL returns an expected INFO
!>          NINFO(2) = No. of times STRSYL3 returns an expected INFO
!> \endverbatim
!>
!> \param[out] KNT
!> \verbatim
!>          KNT is INTEGER
!>          Total number of examples tested.
!> \endverbatim

!
!  -- LAPACK test routine --
   SUBROUTINE SSYL01( THRESH, NFAIL, RMAX, NINFO, KNT )
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
                      KUA, KLB, KUB, LIWORK, M, N
   REAL               ANRM, BNRM, BIGNUM, EPS, RES, RES1, RMUL, &
                      SCALE, SCALE3, SMLNUM, TNRM, XNRM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   REAL               DUML( MAXM ), DUMR( MAXN ), &
                      D( MAX( MAXM, MAXN ) ), DUM( MAXN ), &
                      VM( 2 )
   INTEGER            ISEED( 4 ), IWORK( MAXM + MAXN + 2 )
!     ..
!     .. Allocatable Arrays ..
   INTEGER            AllocateStatus
   REAL, DIMENSION(:,:), ALLOCATABLE :: A, B, C, CC, X, &
                      SWORK
!     ..
!     .. External Functions ..
   LOGICAL            SISNAN
   REAL               SLAMCH, SLANGE
   EXTERNAL           SISNAN, SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SLATMR, SLACPY, SGEMM, STRSYL, STRSYL3
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
   VM( 1 ) = 1.0E+0
   VM( 2 ) = 0.05E+0
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
   LIWORK = MAXM + MAXN + 2
   DO J = 1, 2
      DO ISGN = -1, 1, 2
!           Reset seed (overwritten by LATMR)
         ISEED( 1:4 ) = 1
         DO M = 32, MAXM, 71
            KLA = 0
            KUA = M - 1
            CALL SLATMR( M, M, 'S', ISEED, 'N', D, &
                         6, 1.0E+0, 1.0E+0, 'T', 'N', &
                         DUML, 1, 1.0E+0, DUMR, 1, 1.0E+0, &
                         'N', IWORK, KLA, KUA, 0.0E+0, &
                         1.0E+0, 'NO', A, MAXM, IWORK, IINFO )
            FORALL (I = 1:M) A( I, I ) = A( I, I ) * VM( J )
            ANRM = SLANGE( 'M', M, M, A, MAXM, DUM )
            DO N = 51, MAXN, 47
               KLB = 0
               KUB = N - 1
               CALL SLATMR( N, N, 'S', ISEED, 'N', D, &
                            6, 1.0E+0, 1.0E+0, 'T', 'N', &
                            DUML, 1, 1.0E+0, DUMR, 1, 1.0E+0, &
                            'N', IWORK, KLB, KUB, 0.0E+0, &
                            1.0E+0, 'NO', B, MAXN, IWORK, IINFO )
               BNRM = SLANGE( 'M', N, N, B, MAXN, DUM )
               TNRM = MAX( ANRM, BNRM )
               CALL SLATMR( M, N, 'S', ISEED, 'N', D, &
                            6, 1.0E+0, 1.0E+0, 'T', 'N', &
                            DUML, 1, 1.0E+0, DUMR, 1, 1.0E+0, &
                            'N', IWORK, M, N, 0.0E+0, 1.0E+0, &
                            'NO', C, MAXM, IWORK, IINFO )
               DO ITRANA = 1, 2
                  IF( ITRANA == 1 ) THEN
                     TRANA = 'N'
                  END IF
                  IF( ITRANA == 2 ) THEN
                     TRANA = 'T'
                  END IF
                  DO ITRANB = 1, 2
                     IF( ITRANB == 1 ) THEN
                        TRANB = 'N'
                     END IF
                     IF( ITRANB == 2 ) THEN
                        TRANB = 'T'
                     END IF
                     KNT = KNT + 1
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'All', M, N, C, MAXM, X, MAXM)
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'All', M, N, C, MAXM, CC, MAXM)
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL STRSYL( TRANA, TRANB, ISGN, M, N, &
                                  A, MAXM, B, MAXN, X, MAXM, &
                                  SCALE, IINFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : STRSYL : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     IF( IINFO /= 0 ) &
                        NINFO( 1 ) = NINFO( 1 ) + 1
                     XNRM = SLANGE( 'M', M, N, X, MAXM, DUM )
                     RMUL = 1.0E+0
                     IF( XNRM > 1.0E+0 .AND. TNRM > 1.0E+0 ) THEN
                        IF( XNRM > BIGNUM / TNRM ) THEN
                           RMUL = 1.0E+0 / MAX( XNRM, TNRM )
                        END IF
                     END IF
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SGEMM( TRANA, 'N', M, N, M, RMUL, &
                                 A, MAXM, X, MAXM, -SCALE*RMUL, &
                                 C, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SGEMM( 'N', TRANB, M, N, N, &
                                  REAL( ISGN )*RMUL, X, MAXM, B, &
                                  MAXN, 1.0E+0, C, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     RES1 = SLANGE( 'M', M, N, C, MAXM, DUM )
                     RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, &
                                 ( ( RMUL*TNRM )*EPS )*XNRM )
                     IF( RES > THRESH ) &
                        NFAIL( 1 ) = NFAIL( 1 ) + 1
                     IF( RES > RMAX( 1 ) ) &
                        RMAX( 1 ) = RES
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'All', M, N, C, MAXM, X, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'All', M, N, C, MAXM, CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL STRSYL3( TRANA, TRANB, ISGN, M, N, &
                                   A, MAXM, B, MAXN, X, MAXM, &
                                   SCALE3, IWORK, LIWORK, &
                                   SWORK, LDSWORK, INFO)
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : STRSYL3 : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     IF( INFO /= 0 ) &
                        NINFO( 2 ) = NINFO( 2 ) + 1
                     XNRM = SLANGE( 'M', M, N, X, MAXM, DUM )
                     RMUL = 1.0E+0
                     IF( XNRM > 1.0E+0 .AND. TNRM > 1.0E+0 ) THEN
                        IF( XNRM > BIGNUM / TNRM ) THEN
                           RMUL = 1.0E+0 / MAX( XNRM, TNRM )
                        END IF
                     END IF
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SGEMM( TRANA, 'N', M, N, M, RMUL, &
                                 A, MAXM, X, MAXM, -SCALE3*RMUL, &
                                 CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SGEMM( 'N', TRANB, M, N, N, &
                                 REAL( ISGN )*RMUL, X, MAXM, B, &
                                 MAXN, 1.0E+0, CC, MAXM )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     RES1 = SLANGE( 'M', M, N, CC, MAXM, DUM )
                     RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, &
                                ( ( RMUL*TNRM )*EPS )*XNRM )
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
!     End of SSYL01
!
END



