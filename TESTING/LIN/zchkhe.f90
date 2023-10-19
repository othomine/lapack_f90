!> \brief \b ZCHKHE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKHE( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                          THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
!                          XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKHE tests ZHETRF, -TRI2, -TRS, -TRS2, -RFS, and -CON.
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
!>          THRESH is DOUBLE PRECISION
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
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(NMAX,2*NSMAX))
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
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZCHKHE( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, &
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
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D+0 )
   COMPLEX*16         CZERO
   PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 10 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 9 )
!     ..
!     .. Local Scalars ..
   LOGICAL            TRFCON, ZEROT
   CHARACTER          DIST, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, &
                      IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, &
                      N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT
   DOUBLE PRECISION   ANORM, CNDNUM, RCOND, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DGET06, ZLANHE
   EXTERNAL           DGET06, ZLANHE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, XLAENV, ZERRHE, ZGET04, &
                      ZHECON, ZHERFS, ZHET01, ZHETRF, ZHETRI2, &
                      ZHETRS, ZLACPY, ZLAIPD, ZLARHS, ZLATB4, ZLATMS, &
                      ZPOT02, ZPOT03, ZPOT05
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
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
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Zomplex precision'
   PATH( 2: 3 ) = 'HE'
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
      CALL ZERRHE( PATH, NOUT )
   INFOT = 0
!
!     Set the minimum block size for which the block routine should
!     be used, which will be later returned by ILAENV
!
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
!              Set up parameters with ZLATB4 for the matrix generator
!              based on the type of matrix to be generated.
!
            CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
!              Generate a matrix with ZLATMS.
!
            SRNAMT = 'ZLATMS'
            CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, &
                         INFO )
!
!              Check error code from ZLATMS and handle error.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
!                 Skip all tests for this generated matrix
!
               GO TO 160
            END IF
!
!              For types 3-6, zero one or more rows and columns of
!              the matrix to test that INFO is returned correctly.
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
                        A( IOFF+I ) = CZERO
                     ENDDO
                     IOFF = IOFF + IZERO
                     DO I = IZERO, N
                        A( IOFF ) = CZERO
                        IOFF = IOFF + LDA
                     ENDDO
                  ELSE
                     IOFF = IZERO
                     DO I = 1, IZERO - 1
                        A( IOFF ) = CZERO
                        IOFF = IOFF + LDA
                     ENDDO
                     IOFF = IOFF - IZERO
                     DO I = IZERO, N
                        A( IOFF+I ) = CZERO
                     ENDDO
                  END IF
               ELSE
                  IF( IUPLO == 1 ) THEN
!
!                       Set the first IZERO rows and columns to zero.
!
                     IOFF = 0
                     DO J = 1, N
                        I2 = MIN( J, IZERO )
                        DO I = 1, I2
                           A( IOFF+I ) = CZERO
                        ENDDO
                        IOFF = IOFF + LDA
                     ENDDO
                  ELSE
!
!                       Set the last IZERO rows and columns to zero.
!
                     IOFF = 0
                     DO J = 1, N
                        I1 = MAX( J, IZERO )
                        DO I = I1, N
                           A( IOFF+I ) = CZERO
                        ENDDO
                        IOFF = IOFF + LDA
                     ENDDO
                  END IF
               END IF
            ELSE
               IZERO = 0
            END IF
!
!              End generate test matrix A.
!
!
!              Set the imaginary part of the diagonals.
!
            CALL ZLAIPD( N, A, LDA+1, 0 )
!
!              Do for each value of NB in NBVAL
!
            DO INB = 1, NNB
!
!                 Set the optimal blocksize, which will be later
!                 returned by ILAENV.
!
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
!
!                 Copy the test matrix A into matrix AFAC which
!                 will be factorized in place. This is needed to
!                 preserve the test matrix A for subsequent tests.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Compute the L*D*L**T or U*D*U**T factorization of the
!                 matrix. IWORK stores details of the interchanges and
!                 the block structure of D. AINV is a work array for
!                 block factorization, LWORK is the length of AINV.
!
               LWORK = MAX( 2, NB )*LDA
               SRNAMT = 'ZHETRF'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHETRF( UPLO, N, AFAC, LDA, IWORK, AINV, LWORK, &
                            INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHETRF : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Adjust the expected value of INFO to account for
!                 pivoting.
!
               K = IZERO
               IF( K > 0 ) THEN
  100                CONTINUE
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
!                 Check error code from ZHETRF and handle error.
!
               IF( INFO /= K ) &
                  CALL ALAERH( PATH, 'ZHETRF', INFO, K, UPLO, N, N, &
                               -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
!
!                 Set the condition estimate flag if the INFO is not 0.
!
               IF( INFO /= 0 ) THEN
                  TRFCON = .TRUE.
               ELSE
                  TRFCON = .FALSE.
               END IF
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
               CALL ZHET01( UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, &
                            LDA, RWORK, RESULT( 1 ) )
               NT = 1
!
!+    TEST 2
!                 Form the inverse and compute the residual.
!
               IF( INB == 1 .AND. .NOT.TRFCON ) THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  SRNAMT = 'ZHETRI2'
                  LWORK = (N+NB+1)*(NB+3)
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHETRI2( UPLO, N, AINV, LDA, IWORK, WORK, &
                               LWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHETRI2 : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from ZHETRI and handle error.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'ZHETRI', INFO, -1, UPLO, N, &
                                  N, -1, -1, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
!                    Compute the residual for a symmetric matrix times
!                    its inverse.
!
                  CALL ZPOT03( UPLO, N, A, LDA, AINV, LDA, WORK, LDA, &
                               RWORK, RCONDC, RESULT( 2 ) )
                  NT = 2
               END IF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 1, NT
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, &
                        RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
                  ENDDO
               NRUN = NRUN + NT
!
!                 Skip the other tests if this is not the first block
!                 size.
!
               IF( INB > 1 ) &
                  GO TO 150
!
!                 Do only the condition estimate if INFO is not 0.
!
               IF( TRFCON ) THEN
                  RCONDC = ZERO
                  GO TO 140
               END IF
!
!                 Do for each value of NRHS in NSVAL.
!
               DO IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
!
!+    TEST 3 (Using TRS)
!                 Solve and compute residual for  A * X = B.
!
!                    Choose a set of NRHS random solution vectors
!                    stored in XACT and set up the right hand side B
!
                  SRNAMT = 'ZLARHS'
                  CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                               NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'ZHETRS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHETRS( UPLO, N, NRHS, AFAC, LDA, IWORK, X, &
                               LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHETRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from ZHETRS and handle error.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'ZHETRS', INFO, 0, UPLO, N, &
                                  N, -1, -1, NRHS, IMAT, NFAIL, &
                                  NERRS, NOUT )
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compute the residual for the solution
!
                  CALL ZPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, &
                               LDA, RWORK, RESULT( 3 ) )
!
!+    TEST 4 (Using TRS2)
!                 Solve and compute residual for  A * X = B.
!
!                    Choose a set of NRHS random solution vectors
!                    stored in XACT and set up the right hand side B
!
                  SRNAMT = 'ZLARHS'
                  CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                               NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'ZHETRS2'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHETRS2( UPLO, N, NRHS, AFAC, LDA, IWORK, X, &
                               LDA, WORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHETRS2 : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from ZHETRS2 and handle error.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'ZHETRS2', INFO, 0, UPLO, N, &
                                  N, -1, -1, NRHS, IMAT, NFAIL, &
                                  NERRS, NOUT )
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compute the residual for the solution
!
                  CALL ZPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, &
                               LDA, RWORK, RESULT( 4 ) )
!
!+    TEST 5
!                 Check solution from generated exact solution.
!
                  CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 5 ) )
!
!+    TESTS 6, 7, and 8
!                 Use iterative refinement to improve the solution.
!
                  SRNAMT = 'ZHERFS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZHERFS( UPLO, N, NRHS, A, LDA, AFAC, LDA, &
                               IWORK, B, LDA, X, LDA, RWORK, &
                               RWORK( NRHS+1 ), WORK, &
                               RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZHERFS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from ZHERFS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'ZHERFS', INFO, 0, UPLO, N, &
                                  N, -1, -1, NRHS, IMAT, NFAIL, &
                                  NERRS, NOUT )
!
                  CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 6 ) )
                  CALL ZPOT05( UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, &
                               XACT, LDA, RWORK, RWORK( NRHS+1 ), &
                               RESULT( 7 ) )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = 3, 8
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, &
                           IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                     ENDDO
                  NRUN = NRUN + 6
!
!                 End do for each value of NRHS in NSVAL.
!
                  ENDDO
!
!+    TEST 9
!                 Get an estimate of RCOND = 1/CNDNUM.
!
  140             CONTINUE
               ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
               SRNAMT = 'ZHECON'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZHECON( UPLO, N, AFAC, LDA, IWORK, ANORM, RCOND, &
                            WORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZHECON : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from ZHECON and handle error.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'ZHECON', INFO, 0, UPLO, N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               RESULT( 9 ) = DGET06( RCOND, RCONDC )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               IF( RESULT( 9 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 9, &
                     RESULT( 9 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
  150          CONTINUE
               ENDDO
  160       CONTINUE
            ENDDO
  170    CONTINUE
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
!     End of ZCHKHE
!
END
                                                                                                                                                                                                                                                                                                            




