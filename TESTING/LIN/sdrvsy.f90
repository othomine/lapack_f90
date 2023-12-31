!> \brief \b SDRVSY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
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
!> SDRVSY tests the driver routines SSYSV and -SVX.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand side vectors to be generated for
!>          each linear system.
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
!>          B is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (NMAX*max(2,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX+2*NRHS)
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
   SUBROUTINE SDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, &
                      A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, &
                      NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NMAX, NN, NOUT, NRHS
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NVAL( * )
   REAL               A( * ), AFAC( * ), AINV( * ), B( * ), &
                      RWORK( * ), WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTYPES, NTESTS
   PARAMETER          ( NTYPES = 10, NTESTS = 6 )
   INTEGER            NFACT
   PARAMETER          ( NFACT = 2 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, FACT, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, &
                      IZERO, J, K, K1, KL, KU, LDA, LWORK, MODE, N, &
                      NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT
   REAL               AINVNM, ANORM, CNDNUM, RCOND, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          FACTS( NFACT ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   REAL               SGET06, SLANSY
   EXTERNAL           SGET06, SLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, SERRVX, SGET04, SLACPY, &
                      SLARHS, SLASET, SLATB4, SLATMS, SPOT02, SPOT05, &
                      SSYSV, SSYSVX, SSYT01, SSYTRF, SSYTRI2, XLAENV
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
   INTRINSIC          MAX, MIN
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Single precision'
   PATH( 2: 3 ) = 'SY'
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
   LWORK = MAX( 2*NMAX, NMAX*NRHS )
!
!     Test the error exits
!
   IF( TSTERR ) &
      CALL SERRVX( PATH, NOUT )
   INFOT = 0
!
!     Set the block size and minimum block size for testing.
!
   NB = 1
   NBMIN = 2
   CALL XLAENV( 1, NB )
   CALL XLAENV( 2, NBMIN )
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
      DO IMAT = 1, NIMAT
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 170
!
!           Skip types 3, 4, 5, or 6 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 6
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 170
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
               GO TO 160
            END IF
!
!              For types 3-6, zero one or more rows and columns of the
!              matrix to test that INFO is returned correctly.
!
            IF( ZEROT ) THEN
               IF( IMAT == 3 ) THEN
                  IZERO = 1
               ELSE IF( IMAT == 4 ) THEN
                  IZERO = N
               ELSE
                  IZERO = N / 2 + 1
               END IF
!
               IF( IMAT < 6 ) THEN
!
!                    Set row and column IZERO to zero.
!
                  IF( IUPLO == 1 ) THEN
                     IOFF = ( IZERO-1 )*LDA
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
                  IOFF = 0
                  IF( IUPLO == 1 ) THEN
!
!                       Set the first IZERO rows and columns to zero.
!
                     DO J = 1, N
                        I2 = MIN( J, IZERO )
                        DO I = 1, I2
                           A( IOFF+I ) = ZERO
                        ENDDO
                        IOFF = IOFF + LDA
                     ENDDO
                  ELSE
!
!                       Set the last IZERO rows and columns to zero.
!
                     DO J = 1, N
                        I1 = MAX( J, IZERO )
                        DO I = I1, N
                           A( IOFF+I ) = ZERO
                        ENDDO
                        IOFF = IOFF + LDA
                     ENDDO
                  END IF
               END IF
            ELSE
               IZERO = 0
            END IF
!
            DO IFACT = 1, NFACT
!
!                 Do first for FACT = 'F', then for other values.
!
               FACT = FACTS( IFACT )
!
!                 Compute the condition number for comparison with
!                 the value returned by SSYSVX.
!
               IF( ZEROT ) THEN
                  IF( IFACT == 1 ) &
                     GO TO 150
                  RCONDC = ZERO
!
               ELSE IF( IFACT == 1 ) THEN
!
!                    Compute the 1-norm of A.
!
                  ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
!
!                    Factor the matrix A.
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
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SSYTRF( UPLO, N, AFAC, LDA, IWORK, WORK, &
                               LWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SSYTRF : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compute inv(A) and take its norm.
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
                  LWORK = (N+NB+1)*(NB+3)
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SSYTRI2( UPLO, N, AINV, LDA, IWORK, WORK, &
                               LWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SSYTRI2 : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  AINVNM = SLANSY( '1', UPLO, N, AINV, LDA, RWORK )
!
!                    Compute the 1-norm condition number of A.
!
                  IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                     RCONDC = ONE
                  ELSE
                     RCONDC = ( ONE / ANORM ) / AINVNM
                  END IF
               END IF
!
!                 Form an exact solution and set the right hand side.
!
               SRNAMT = 'SLARHS'
               CALL SLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                            NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, &
                            INFO )
               XTYPE = 'C'
!
!                 --- Test SSYSV  ---
!
               IF( IFACT == 2 ) THEN
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
!                    Factor the matrix and solve the system using SSYSV.
!
                  SRNAMT = 'SSYSV '
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SSYSV( UPLO, N, NRHS, AFAC, LDA, IWORK, X, &
                              LDA, WORK, LWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SSYSV : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Adjust the expected value of INFO to account for
!                    pivoting.
!
                  K = IZERO
                  IF( K > 0 ) THEN
  100                   CONTINUE
                     IF( IWORK( K ) < 0 ) THEN
                        IF( IWORK( K ) /= -K ) THEN
                           K = -IWORK( K )
                           GO TO 100
                        END IF
                     ELSE IF( IWORK( K ) /= K ) THEN
                        K = IWORK( K )
                        GO TO 100
                     END IF
                  END IF
!
!                    Check error code from SSYSV .
!
                  IF( INFO /= K ) THEN
                     CALL ALAERH( PATH, 'SSYSV ', INFO, K, UPLO, N, &
                                  N, -1, -1, NRHS, IMAT, NFAIL, &
                                  NERRS, NOUT )
                     GO TO 120
                  ELSE IF( INFO /= 0 ) THEN
                     GO TO 120
                  END IF
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                  CALL SSYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK, &
                               AINV, LDA, RWORK, RESULT( 1 ) )
!
!                    Compute residual of the computed solution.
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
                               LDA, RWORK, RESULT( 2 ) )
!
!                    Check solution from generated exact solution.
!
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 3 ) )
                  NT = 3
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = 1, NT
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )'SSYSV ', UPLO, N, &
                           IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                     ENDDO
                  NRUN = NRUN + NT
  120                CONTINUE
               END IF
!
!                 --- Test SSYSVX ---
!
               IF( IFACT == 2 )  THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLASET( UPLO, N, N, ZERO, ZERO, AFAC, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               ENDIF
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLASET( 'Full', N, NRHS, ZERO, ZERO, X, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Solve the system and compute the condition number and
!                 error bounds using SSYSVX.
!
               SRNAMT = 'SSYSVX'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SSYSVX( FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, &
                            IWORK, B, LDA, X, LDA, RCOND, RWORK, &
                            RWORK( NRHS+1 ), WORK, LWORK, &
                            IWORK( N+1 ), INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SSYSVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Adjust the expected value of INFO to account for
!                 pivoting.
!
               K = IZERO
               IF( K > 0 ) THEN
  130                CONTINUE
                  IF( IWORK( K ) < 0 ) THEN
                     IF( IWORK( K ) /= -K ) THEN
                        K = -IWORK( K )
                        GO TO 130
                     END IF
                  ELSE IF( IWORK( K ) /= K ) THEN
                     K = IWORK( K )
                     GO TO 130
                  END IF
               END IF
!
!                 Check the error code from SSYSVX.
!
               IF( INFO /= K ) THEN
                  CALL ALAERH( PATH, 'SSYSVX', INFO, K, FACT // UPLO, &
                               N, N, -1, -1, NRHS, IMAT, NFAIL, &
                               NERRS, NOUT )
                  GO TO 150
               END IF
!
               IF( INFO == 0 ) THEN
                  IF( IFACT >= 2 ) THEN
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                     CALL SSYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK, &
                                  AINV, LDA, RWORK( 2*NRHS+1 ), &
                                  RESULT( 1 ) )
                     K1 = 1
                  ELSE
                     K1 = 2
                  END IF
!
!                    Compute residual of the computed solution.
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
                               LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )
!
!                    Check solution from generated exact solution.
!
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 3 ) )
!
!                    Check the error bounds from iterative refinement.
!
                  CALL SPOT05( UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, &
                               XACT, LDA, RWORK, RWORK( NRHS+1 ), &
                               RESULT( 4 ) )
               ELSE
                  K1 = 6
               END IF
!
!                 Compare RCOND from SSYSVX with the computed value
!                 in RCONDC.
!
               RESULT( 6 ) = SGET06( RCOND, RCONDC )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = K1, 6
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALADHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 )'SSYSVX', FACT, UPLO, &
                        N, IMAT, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
                  ENDDO
               NRUN = NRUN + 7 - K1
!
  150          CONTINUE
               ENDDO
!
  160       CONTINUE
            ENDDO
  170    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, &
         ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5, &
         ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
   RETURN
!
!     End of SDRVSY
!
END
                                                                                                                                                                                                                                                                                                            




