!> \brief \b SCHKQ3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKQ3( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL,
!                          THRESH, A, COPYA, S, TAU, WORK, IWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NM, NN, NNB, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NVAL( * ),
!      $                   NXVAL( * )
!       REAL               A( * ), COPYA( * ), S( * ),
!      $                   TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKQ3 tests SGEQP3.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          The matrix types to be used for testing.  Matrices of type j
!>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
!>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
!> \endverbatim
!>
!> \param[in] NM
!> \verbatim
!>          NM is INTEGER
!>          The number of values of M contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB and NX contained in the
!>          vectors NBVAL and NXVAL.  The blocking parameters are used
!>          in pairs (NB,NX).
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NNB)
!>          The values of the blocksize NB.
!> \endverbatim
!>
!> \param[in] NXVAL
!> \verbatim
!>          NXVAL is INTEGER array, dimension (NNB)
!>          The values of the crossover point NX.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (MMAX*NMAX)
!>          where MMAX is the maximum value of M in MVAL and NMAX is the
!>          maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] COPYA
!> \verbatim
!>          COPYA is REAL array, dimension (MMAX*NMAX)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension
!>                      (min(MMAX,NMAX))
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (MMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                      (MMAX*NMAX + 4*NMAX + MMAX)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*NMAX)
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
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
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SCHKQ3( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, &
                      THRESH, A, COPYA, S, TAU, WORK, IWORK, &
                      NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NM, NN, NNB, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NVAL( * ), &
                      NXVAL( * )
   REAL               A( * ), COPYA( * ), S( * ), &
                      TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 6 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 3 )
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
!     ..
!     .. Local Scalars ..
   CHARACTER*3        PATH
   INTEGER            I, IHIGH, ILOW, IM, IMODE, IN, INB, INFO, &
                      ISTEP, K, LDA, LW, LWORK, M, MNMIN, MODE, N, &
                      NB, NERRS, NFAIL, NRUN, NX
   REAL               EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SQPT01, SQRT11, SQRT12
   EXTERNAL           SLAMCH, SQPT01, SQRT11, SQRT12
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAHD, ALASUM, ICOPY, SGEQP3, SLACPY, SLAORD, &
                      SLASET, SLATMS, XLAENV
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, IOUNIT
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Single precision'
   PATH( 2: 3 ) = 'Q3'
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
   EPS = SLAMCH( 'Epsilon' )
   INFOT = 0
!
   DO IM = 1, NM
!
!        Do for each value of M in MVAL.
!
      M = MVAL( IM )
      LDA = MAX( 1, M )
!
      DO IN = 1, NN
!
!           Do for each value of N in NVAL.
!
         N = NVAL( IN )
         MNMIN = MIN( M, N )
         LWORK = MAX( 1, M*MAX( M, N )+4*MNMIN+MAX( M, N ), &
                      M*N + 2*MNMIN + 4*N )
!
         DO IMODE = 1, NTYPES
            IF( .NOT.DOTYPE( IMODE ) ) &
               GO TO 70
!
!              Do for each type of matrix
!                 1:  zero matrix
!                 2:  one small singular value
!                 3:  geometric distribution of singular values
!                 4:  first n/2 columns fixed
!                 5:  last n/2 columns fixed
!                 6:  every second column fixed
!
            MODE = IMODE
            IF( IMODE > 3 ) &
               MODE = 1
!
!              Generate test matrix of size m by n using
!              singular value distribution indicated by `mode'.
!
            DO I = 1, N
               IWORK( I ) = 0
            ENDDO
            IF( IMODE == 1 ) THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLASET( 'Full', M, N, ZERO, ZERO, COPYA, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               DO I = 1, MNMIN
                  S( I ) = ZERO
               ENDDO
            ELSE
               CALL SLATMS( M, N, 'Uniform', ISEED, 'Nonsymm', S, &
                            MODE, ONE / EPS, ONE, M, N, 'No packing', &
                            COPYA, LDA, WORK, INFO )
               IF( IMODE >= 4 ) THEN
                  IF( IMODE == 4 ) THEN
                     ILOW = 1
                     ISTEP = 1
                     IHIGH = MAX( 1, N / 2 )
                  ELSE IF( IMODE == 5 ) THEN
                     ILOW = MAX( 1, N / 2 )
                     ISTEP = 1
                     IHIGH = N
                  ELSE IF( IMODE == 6 ) THEN
                     ILOW = 1
                     ISTEP = 2
                     IHIGH = N
                  END IF
                  DO I = ILOW, IHIGH, ISTEP
                     IWORK( I ) = 1
                  ENDDO
               END IF
               CALL SLAORD( 'Decreasing', MNMIN, S, 1 )
            END IF
!
            DO INB = 1, NNB
!
!                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
!
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
               NX = NXVAL( INB )
               CALL XLAENV( 3, NX )
!
!                 Get a working copy of COPYA into A and a copy of
!                 vector IWORK.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( 'All', M, N, COPYA, LDA, A, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               CALL ICOPY( N, IWORK( 1 ), 1, IWORK( N+1 ), 1 )
!
!                 Compute the QR factorization with pivoting of A
!
               LW = MAX( 1, 2*N+NB*( N+1 ) )
!
!                 Compute the QP3 factorization of A
!
               SRNAMT = 'SGEQP3'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGEQP3( M, N, A, LDA, IWORK( N+1 ), TAU, WORK, &
                            LW, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGEQP3 : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Compute norm(svd(a) - svd(r))
!
               RESULT( 1 ) = SQRT12( M, N, A, LDA, S, WORK, &
                             LWORK )
!
!                 Compute norm( A*P - Q*R )
!
               RESULT( 2 ) = SQPT01( M, N, MNMIN, COPYA, A, LDA, TAU, &
                             IWORK( N+1 ), WORK, LWORK )
!
!                 Compute Q'*Q
!
               RESULT( 3 ) = SQRT11( M, MNMIN, A, LDA, TAU, WORK, &
                             LWORK )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 1, NTESTS
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9999 )'SGEQP3', M, N, NB, &
                        IMODE, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
               ENDDO
               NRUN = NRUN + NTESTS
!
            ENDDO
70       CONTINUE
         ENDDO
      ENDDO
   ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A, ' M =', I5, ', N =', I5, ', NB =', I4, ', type ', &
         I2, ', test ', I2, ', ratio =', G12.5 )
!
!     End of SCHKQ3
!
END
                                                                                                                                                                                                                                                                                                            




