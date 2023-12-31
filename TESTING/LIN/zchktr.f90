!> \brief \b ZCHKTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKTR( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                          THRESH, TSTERR, NMAX, A, AINV, B, X, XACT,
!                          WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NBVAL( * ), NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AINV( * ), B( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKTR tests ZTRTRI, -TRS, -RFS, and -CON, and ZLATRS(3)
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
!>          The values of the matrix column dimension N.
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
!>          The leading dimension of the work arrays.
!>          NMAX >= the maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
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
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      (max(NMAX,2*NSMAX))
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
   SUBROUTINE ZCHKTR( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, &
                      THRESH, TSTERR, NMAX, A, AINV, B, X, XACT, &
                      WORK, RWORK, NOUT )
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
   INTEGER            NBVAL( * ), NSVAL( * ), NVAL( * )
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( * ), AINV( * ), B( * ), WORK( * ), X( * ), &
                      XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NTYPE1, NTYPES
   PARAMETER          ( NTYPE1 = 10, NTYPES = 18 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 10 )
   INTEGER            NTRAN
   PARAMETER          ( NTRAN = 3 )
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
   CHARACTER          DIAG, NORM, TRANS, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IDIAG, IMAT, IN, INB, INFO, IRHS, ITRAN, &
                      IUPLO, K, LDA, N, NB, NERRS, NFAIL, NRHS, NRUN
   DOUBLE PRECISION   AINVNM, ANORM, BIGNUM, DUMMY, RCOND, RCONDC, &
                      RCONDI, RCONDO, RES, SCALE, DLAMCH
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          TRANSS( NTRAN ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS ), RWORK2( 2*NMAX ), &
                      SCALE3( 2 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   ZLANTR
   EXTERNAL           LSAME, ZLANTR
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, DLAMCH, XLAENV, ZCOPY, &
                      ZDSCAL, ZERRTR, ZGET04, ZLACPY, ZLARHS, ZLATRS, &
                      ZLATRS3, ZLATTR, ZTRCON, ZTRRFS, ZTRT01, &
                      ZTRT02, ZTRT03, ZTRT05, ZTRT06, ZTRTRI, ZTRTRS
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
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' / , TRANSS / 'N', 'T', 'C' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Zomplex precision'
   PATH( 2: 3 ) = 'TR'
   BIGNUM = DLAMCH('Overflow') / DLAMCH('Precision')
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
      CALL ZERRTR( PATH, NOUT )
   INFOT = 0
!
   DO IN = 1, NN
!
!        Do for each value of N in NVAL
!
      N = NVAL( IN )
      LDA = MAX( 1, N )
      XTYPE = 'N'
!
      DO IMAT = 1, NTYPE1
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 80
!
         DO IUPLO = 1, 2
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
            UPLO = UPLOS( IUPLO )
!
!              Call ZLATTR to generate a triangular test matrix.
!
            SRNAMT = 'ZLATTR'
            CALL ZLATTR( IMAT, UPLO, 'No transpose', DIAG, ISEED, N, &
                         A, LDA, X, WORK, RWORK, INFO )
!
!              Set IDIAG = 1 for non-unit matrices, 2 for unit.
!
            IF( LSAME( DIAG, 'N' ) ) THEN
               IDIAG = 1
            ELSE
               IDIAG = 2
            END IF
!
            DO INB = 1, NNB
!
!                 Do for each blocksize in NBVAL
!
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
!
!+    TEST 1
!                 Form the inverse of A.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( UPLO, N, N, A, LDA, AINV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'ZTRTRI'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZTRTRI( UPLO, DIAG, N, AINV, LDA, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZTRTRI : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from ZTRTRI.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'ZTRTRI', INFO, 0, UPLO // DIAG, &
                               N, N, -1, -1, NB, IMAT, NFAIL, NERRS, &
                               NOUT )
!
!                 Compute the infinity-norm condition number of A.
!
               ANORM = ZLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK )
               AINVNM = ZLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, &
                        RWORK )
               IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDI = ONE
               ELSE
                  RCONDI = ( ONE / ANORM ) / AINVNM
               END IF
!
!                 Compute the residual for the triangular matrix times
!                 its inverse.  Also compute the 1-norm condition number
!                 of A.
!
               CALL ZTRT01( UPLO, DIAG, N, A, LDA, AINV, LDA, RCONDO, &
                            RWORK, RESULT( 1 ) )
!                 Print the test ratio if it is  >=  THRESH.
!
               IF( RESULT( 1 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9999 )UPLO, DIAG, N, NB, IMAT, &
                     1, RESULT( 1 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
!
!                 Skip remaining tests if not the first block size.
!
               IF( INB /= 1 ) &
                  GO TO 60
!
               DO IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
                  XTYPE = 'N'
!
                  DO ITRAN = 1, NTRAN
!
!                    Do for op(A) = A, A**T, or A**H.
!
                     TRANS = TRANSS( ITRAN )
                     IF( ITRAN == 1 ) THEN
                        NORM = 'O'
                        RCONDC = RCONDO
                     ELSE
                        NORM = 'I'
                        RCONDC = RCONDI
                     END IF
!
!+    TEST 2
!                       Solve and compute residual for op(A)*x = b.
!
                     SRNAMT = 'ZLARHS'
                     CALL ZLARHS( PATH, XTYPE, UPLO, TRANS, N, N, 0, &
                                  IDIAG, NRHS, A, LDA, XACT, LDA, B, &
                                  LDA, ISEED, INFO )
                     XTYPE = 'C'
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
                     SRNAMT = 'ZTRTRS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, &
                                  X, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : ZTRTRS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from ZTRTRS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'ZTRTRS', INFO, 0, &
                                     UPLO // TRANS // DIAG, N, N, -1, &
                                     -1, NRHS, IMAT, NFAIL, NERRS, &
                                     NOUT )
!
!                       This line is needed on a Sun SPARCstation.
!
                     IF( N > 0 ) &
                        DUMMY = DBLE( A( 1 ) )
!
                     CALL ZTRT02( UPLO, TRANS, DIAG, N, NRHS, A, LDA, &
                                  X, LDA, B, LDA, WORK, RWORK, &
                                  RESULT( 2 ) )
!
!+    TEST 3
!                       Check solution from generated exact solution.
!
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 3 ) )
!
!+    TESTS 4, 5, and 6
!                       Use iterative refinement to improve the solution
!                       and compute error bounds.
!
                     SRNAMT = 'ZTRRFS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, &
                                  B, LDA, X, LDA, RWORK, &
                                  RWORK( NRHS+1 ), WORK, &
                                  RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : ZTRRFS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from ZTRRFS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'ZTRRFS', INFO, 0, &
                                     UPLO // TRANS // DIAG, N, N, -1, &
                                     -1, NRHS, IMAT, NFAIL, NERRS, &
                                     NOUT )
!
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 4 ) )
                     CALL ZTRT05( UPLO, TRANS, DIAG, N, NRHS, A, LDA, &
                                  B, LDA, X, LDA, XACT, LDA, RWORK, &
                                  RWORK( NRHS+1 ), RESULT( 5 ) )
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                     DO K = 2, 6
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9998 )UPLO, TRANS, &
                              DIAG, N, NRHS, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
                     ENDDO
                     NRUN = NRUN + 5
                  ENDDO
               ENDDO
!
!+    TEST 7
!                       Get an estimate of RCOND = 1/CNDNUM.
!
               DO ITRAN = 1, 2
                  IF( ITRAN == 1 ) THEN
                     NORM = 'O'
                     RCONDC = RCONDO
                  ELSE
                     NORM = 'I'
                     RCONDC = RCONDI
                  END IF
                  SRNAMT = 'ZTRCON'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, &
                               WORK, RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZTRCON : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                       Check error code from ZTRCON.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'ZTRCON', INFO, 0, &
                                  NORM // UPLO // DIAG, N, N, -1, -1, &
                                  -1, IMAT, NFAIL, NERRS, NOUT )
!
                  CALL ZTRT06( RCOND, RCONDC, UPLO, DIAG, N, A, LDA, &
                               RWORK, RESULT( 7 ) )
!
!                    Print the test ratio if it is  >=  THRESH.
!
                  IF( RESULT( 7 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )NORM, UPLO, N, IMAT, &
                        7, RESULT( 7 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
               ENDDO
60          CONTINUE
            ENDDO
         ENDDO
80    CONTINUE
      ENDDO
!
!        Use pathological test matrices to test ZLATRS.
!
      DO IMAT = NTYPE1 + 1, NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 110
!
         DO IUPLO = 1, 2
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
            UPLO = UPLOS( IUPLO )
            DO ITRAN = 1, NTRAN
!
!                 Do for op(A) = A, A**T, and A**H.
!
               TRANS = TRANSS( ITRAN )
!
!                 Call ZLATTR to generate a triangular test matrix.
!
               SRNAMT = 'ZLATTR'
               CALL ZLATTR( IMAT, UPLO, TRANS, DIAG, ISEED, N, A, &
                            LDA, X, WORK, RWORK, INFO )
!
!+    TEST 8
!                 Solve the system op(A)*x = b.
!
               SRNAMT = 'ZLATRS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZCOPY( N, X, 1, B, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLATRS( UPLO, TRANS, DIAG, 'N', N, A, LDA, B, &
                            SCALE, RWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLATRS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from ZLATRS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'ZLATRS', INFO, 0, &
                               UPLO // TRANS // DIAG // 'N', N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, &
                            RWORK, ONE, B, LDA, X, LDA, WORK, &
                            RESULT( 8 ) )
!
!+    TEST 9
!                 Solve op(A)*X = b again with NORMIN = 'Y'.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZCOPY( N, X, 1, B( N+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLATRS( UPLO, TRANS, DIAG, 'Y', N, A, LDA, &
                            B( N+1 ), SCALE, RWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLATRS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from ZLATRS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'ZLATRS', INFO, 0, &
                               UPLO // TRANS // DIAG // 'Y', N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, &
                            RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, &
                            RESULT( 9 ) )
!
!+    TEST 10
!                 Solve op(A)*X = B
!
               SRNAMT = 'ZLATRS3'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZCOPY( N, X, 1, B, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZCOPY( N, X, 1, B( N+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZDSCAL( N, BIGNUM, B( N+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZDSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLATRS3( UPLO, TRANS, DIAG, 'N', N, 2, A, LDA, &
                             B, MAX(1, N), SCALE3, RWORK, RWORK2, &
                             2*NMAX, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLATRS3 : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from ZLATRS3.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'ZLATRS3', INFO, 0, &
                               UPLO // TRANS // DIAG // 'N', N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
               CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, &
                            SCALE3( 1 ), RWORK, ONE, B( 1 ), LDA, &
                            X, LDA, WORK, RESULT( 10 ) )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZDSCAL( N, BIGNUM, X, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZDSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, &
                            SCALE3( 2 ), RWORK, ONE, B( N+1 ), LDA, &
                            X, LDA, WORK, RES )
               RESULT( 10 ) = MAX( RESULT( 10 ), RES )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               IF( RESULT( 8 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9996 )'ZLATRS', UPLO, TRANS, &
                     DIAG, 'N', N, IMAT, 8, RESULT( 8 )
                  NFAIL = NFAIL + 1
               END IF
               IF( RESULT( 9 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9996 )'ZLATRS', UPLO, TRANS, &
                     DIAG, 'Y', N, IMAT, 9, RESULT( 9 )
                  NFAIL = NFAIL + 1
               END IF
               IF( RESULT( 10 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9996 )'ZLATRS3', UPLO, TRANS, &
                     DIAG, 'N', N, IMAT, 10, RESULT( 10 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 3
            ENDDO
            ENDDO
  110    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' UPLO=''', A1, ''', DIAG=''', A1, ''', N=', I5, ', NB=', &
         I4, ', type ', I2, ', test(', I2, ')= ', G12.5 )
 9998 FORMAT( ' UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, &
         ''', N=', I5, ', NB=', I4, ', type ', I2, ', test(', &
         I2, ')= ', G12.5 )
 9997 FORMAT( ' NORM=''', A1, ''', UPLO =''', A1, ''', N=', I5, ',', &
         11X, ' type ', I2, ', test(', I2, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''', ''', &
         A1, ''',', I5, ', ... ), type ', I2, ', test(', I2, ')=', &
         G12.5 )
   RETURN
!
!     End of ZCHKTR
!
END
                                                                                                                                                                                                                                                                                                            




