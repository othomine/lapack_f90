!> \brief \b SCHKPO
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKPO( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                          THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
!                          XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
!       REAL               A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKPO tests SPOTRF, -TRI, -TRS, -RFS, and -CON
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
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix dimension N.
!> \endverbatim
!>
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NNB)
!>          The values of the blocksize NB.
!> \endverbatim
!>
!> \param[in] NNS
!> \verbatim
!>          NNS is INTEGER
!>          The number of values of NRHS contained in the vector NSVAL.
!> \endverbatim
!>
!> \param[in] NSVAL
!> \verbatim
!>          NSVAL is INTEGER array, dimension (NNS)
!>          The values of the number of right hand sides NRHS.
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
!> \param[in] TSTERR
!> \verbatim
!>          TSTERR is LOGICAL
!>          Flag that indicates whether error exits are to be tested.
!> \endverbatim
!>
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
!>                      (max(NMAX,2*NSMAX))
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (NMAX)
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
   SUBROUTINE SCHKPO( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, &
                      THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, &
                      XACT, WORK, RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NMAX, NN, NNB, NNS, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
   REAL               A( * ), AFAC( * ), AINV( * ), B( * ), &
                      RWORK( * ), WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO
   PARAMETER          ( ZERO = 0.0E+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 9 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 8 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, &
                      IZERO, K, KL, KU, LDA, MODE, N, NB, NERRS, &
                      NFAIL, NIMAT, NRHS, NRUN
   REAL               ANORM, CNDNUM, RCOND, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   REAL               SGET06, SLANSY
   EXTERNAL           SGET06, SLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, SERRPO, SGET04, SLACPY, &
                      SLARHS, SLATB4, SLATMS, SPOCON, SPORFS, SPOT01, &
                      SPOT02, SPOT03, SPOT05, SPOTRF, SPOTRI, SPOTRS, &
                      XLAENV
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, NUNIT
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Single precision'
   PATH( 2: 3 ) = 'PO'
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
!     Test the error exits
!
   IF( TSTERR ) &
      CALL SERRPO( PATH, NOUT )
   INFOT = 0
   CALL XLAENV( 2, 2 )
!
!     Do for each value of N in NVAL
!
   DO IN = 1, NN
      N = NVAL( IN )
      LDA = MAX( N, 1 )
      XTYPE = 'N'
      NIMAT = NTYPES
      IF( N <= 0 ) &
         NIMAT = 1
!
      IZERO = 0
      DO IMAT = 1, NIMAT
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 110
!
!           Skip types 3, 4, or 5 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 5
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 110
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
!
!              Set up parameters with SLATB4 and generate a test matrix
!              with SLATMS.
!
            CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
            SRNAMT = 'SLATMS'
            CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, &
                         INFO )
!
!              Check error code from SLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 100
            END IF
!
!              For types 3-5, zero one row and column of the matrix to
!              test that INFO is returned correctly.
!
            IF( ZEROT ) THEN
               IF( IMAT == 3 ) THEN
                  IZERO = 1
               ELSE IF( IMAT == 4 ) THEN
                  IZERO = N
               ELSE
                  IZERO = N / 2 + 1
               END IF
               IOFF = ( IZERO-1 )*LDA
!
!                 Set row and column IZERO of A to 0.
!
               IF( IUPLO == 1 ) THEN
                  DO I = 1, IZERO - 1
                     A( IOFF+I ) = ZERO
                  ENDDO
                  IOFF = IOFF + IZERO
                  DO I = IZERO, N
                     A( IOFF ) = ZERO
                     IOFF = IOFF + LDA
                  ENDDO
               ELSE
                  IOFF = IZERO
                  DO I = 1, IZERO - 1
                     A( IOFF ) = ZERO
                     IOFF = IOFF + LDA
                  ENDDO
                  IOFF = IOFF - IZERO
                  DO I = IZERO, N
                     A( IOFF+I ) = ZERO
                  ENDDO
               END IF
            ELSE
               IZERO = 0
            END IF
!
!              Do for each value of NB in NBVAL
!
            DO INB = 1, NNB
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
!
!                 Compute the L*L' or U'*U factorization of the matrix.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'SPOTRF'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SPOTRF( UPLO, N, AFAC, LDA, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SPOTRF : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from SPOTRF.
!
               IF( INFO /= IZERO ) THEN
                  CALL ALAERH( PATH, 'SPOTRF', INFO, IZERO, UPLO, N, &
                               N, -1, -1, NB, IMAT, NFAIL, NERRS, &
                               NOUT )
                  GO TO 90
               END IF
!
!                 Skip the tests if INFO is not 0.
!
               IF( INFO /= 0 ) &
                  GO TO 90
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               CALL SPOT01( UPLO, N, A, LDA, AINV, LDA, RWORK, &
                            RESULT( 1 ) )
!
!+    TEST 2
!                 Form the inverse and compute the residual.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'SPOTRI'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SPOTRI( UPLO, N, AINV, LDA, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SPOTRI : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from SPOTRI.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'SPOTRI', INFO, 0, UPLO, N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               CALL SPOT03( UPLO, N, A, LDA, AINV, LDA, WORK, LDA, &
                            RWORK, RCONDC, RESULT( 2 ) )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 1, 2
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, &
                        RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
               ENDDO
               NRUN = NRUN + 2
!
!                 Skip the rest of the tests unless this is the first
!                 blocksize.
!
               IF( INB /= 1 ) &
                  GO TO 90
!
               DO IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
!
!+    TEST 3
!                 Solve and compute residual for A * X = B .
!
                  SRNAMT = 'SLARHS'
                  CALL SLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                               NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'SPOTRS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SPOTRS( UPLO, N, NRHS, AFAC, LDA, X, LDA, &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SPOTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                 Check error code from SPOTRS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'SPOTRS', INFO, 0, UPLO, N, &
                                  N, -1, -1, NRHS, IMAT, NFAIL, &
                                  NERRS, NOUT )
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  CALL SPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, &
                               LDA, RWORK, RESULT( 3 ) )
!
!+    TEST 4
!                 Check solution from generated exact solution.
!
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 4 ) )
!
!+    TESTS 5, 6, and 7
!                 Use iterative refinement to improve the solution.
!
                  SRNAMT = 'SPORFS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SPORFS( UPLO, N, NRHS, A, LDA, AFAC, LDA, B, &
                               LDA, X, LDA, RWORK, RWORK( NRHS+1 ), &
                               WORK, IWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SPORFS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                 Check error code from SPORFS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'SPORFS', INFO, 0, UPLO, N, &
                                  N, -1, -1, NRHS, IMAT, NFAIL, &
                                  NERRS, NOUT )
!
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 5 ) )
                  CALL SPOT05( UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, &
                               XACT, LDA, RWORK, RWORK( NRHS+1 ), &
                               RESULT( 6 ) )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = 3, 7
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, &
                           IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                  ENDDO
                  NRUN = NRUN + 5
               ENDDO
!
!+    TEST 8
!                 Get an estimate of RCOND = 1/CNDNUM.
!
               ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
               SRNAMT = 'SPOCON'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SPOCON( UPLO, N, AFAC, LDA, ANORM, RCOND, WORK, &
                            IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SPOCON : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from SPOCON.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'SPOCON', INFO, 0, UPLO, N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               RESULT( 8 ) = SGET06( RCOND, RCONDC )
!
!                 Print the test ratio if it is  >=  THRESH.
!
               IF( RESULT( 8 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 8, &
                     RESULT( 8 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
90          CONTINUE
            ENDDO
  100       CONTINUE
            ENDDO
  110    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', &
         I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ',', 10X, ' type ', I2, &
         ', test(', I2, ') =', G12.5 )
   RETURN
!
!     End of SCHKPO
!
END
                                                                                                                                                                                                                                                                                                            




