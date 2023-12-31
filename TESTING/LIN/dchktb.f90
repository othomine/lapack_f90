!> \brief \b DCHKTB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, AB, AINV, B, X, XACT, WORK, RWORK, IWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   AB( * ), AINV( * ), B( * ), RWORK( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKTB tests DTBTRS, -RFS, and -CON, and DLATBS.
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
!> \param[out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      NMAX, AB, AINV, B, X, XACT, WORK, RWORK, IWORK, &
                      NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NMAX, NN, NNS, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
   DOUBLE PRECISION   AB( * ), AINV( * ), B( * ), RWORK( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NTYPE1, NTYPES
   PARAMETER          ( NTYPE1 = 9, NTYPES = 17 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 8 )
   INTEGER            NTRAN
   PARAMETER          ( NTRAN = 3 )
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   CHARACTER          DIAG, NORM, TRANS, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IDIAG, IK, IMAT, IN, INFO, IRHS, ITRAN, &
                      IUPLO, J, K, KD, LDA, LDAB, N, NERRS, NFAIL, &
                      NIMAT, NIMAT2, NK, NRHS, NRUN
   DOUBLE PRECISION   AINVNM, ANORM, RCOND, RCONDC, RCONDI, RCONDO, &
                      SCALE
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          TRANSS( NTRAN ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLANTB, DLANTR
   EXTERNAL           LSAME, DLANTB, DLANTR
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, DCOPY, DERRTR, DGET04, &
                      DLACPY, DLARHS, DLASET, DLATBS, DLATTB, DTBCON, &
                      DTBRFS, DTBSV, DTBT02, DTBT03, DTBT05, DTBT06, &
                      DTBTRS
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
   INTRINSIC          MAX, MIN
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' / , TRANSS / 'N', 'T', 'C' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Double precision'
   PATH( 2: 3 ) = 'TB'
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
      CALL DERRTR( PATH, NOUT )
   INFOT = 0
!
   DO IN = 1, NN
!
!        Do for each value of N in NVAL
!
      N = NVAL( IN )
      LDA = MAX( 1, N )
      XTYPE = 'N'
      NIMAT = NTYPE1
      NIMAT2 = NTYPES
      IF( N <= 0 ) THEN
         NIMAT = 1
         NIMAT2 = NTYPE1 + 1
      END IF
!
      NK = MIN( N+1, 4 )
      DO IK = 1, NK
!
!           Do for KD = 0, N, (3N-1)/4, and (N+1)/4. This order makes
!           it easier to skip redundant values for small values of N.
!
         IF( IK == 1 ) THEN
            KD = 0
         ELSE IF( IK == 2 ) THEN
            KD = MAX( N, 0 )
         ELSE IF( IK == 3 ) THEN
            KD = ( 3*N-1 ) / 4
         ELSE IF( IK == 4 ) THEN
            KD = ( N+1 ) / 4
         END IF
         LDAB = KD + 1
!
         DO IMAT = 1, NIMAT
!
!              Do the tests only if DOTYPE( IMAT ) is true.
!
            IF( .NOT.DOTYPE( IMAT ) ) &
               GO TO 90
!
            DO IUPLO = 1, 2
!
!                 Do first for UPLO = 'U', then for UPLO = 'L'
!
               UPLO = UPLOS( IUPLO )
!
!                 Call DLATTB to generate a triangular test matrix.
!
               SRNAMT = 'DLATTB'
               CALL DLATTB( IMAT, UPLO, 'No transpose', DIAG, ISEED, &
                            N, KD, AB, LDAB, X, WORK, INFO )
!
!                 Set IDIAG = 1 for non-unit matrices, 2 for unit.
!
               IF( LSAME( DIAG, 'N' ) ) THEN
                  IDIAG = 1
               ELSE
                  IDIAG = 2
               END IF
!
!                 Form the inverse of A so we can get a good estimate
!                 of RCONDC = 1/(norm(A) * norm(inv(A))).
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLASET( 'Full', N, N, ZERO, ONE, AINV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( LSAME( UPLO, 'U' ) ) THEN
                  DO J = 1, N
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DTBSV( UPLO, 'No transpose', DIAG, J, KD, &
                                 AB, LDAB, AINV( ( J-1 )*LDA+1 ), 1 )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DTBSV : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDDO
               ELSE
                  DO J = 1, N
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DTBSV( UPLO, 'No transpose', DIAG, N-J+1, &
                                 KD, AB( ( J-1 )*LDAB+1 ), LDAB, &
                                 AINV( ( J-1 )*LDA+J ), 1 )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DTBSV : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDDO
               END IF
!
!                 Compute the 1-norm condition number of A.
!
               ANORM = DLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, &
                       RWORK )
               AINVNM = DLANTR( '1', UPLO, DIAG, N, N, AINV, LDA, &
                        RWORK )
               IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDO = ONE
               ELSE
                  RCONDO = ( ONE / ANORM ) / AINVNM
               END IF
!
!                 Compute the infinity-norm condition number of A.
!
               ANORM = DLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, &
                       RWORK )
               AINVNM = DLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, &
                        RWORK )
               IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDI = ONE
               ELSE
                  RCONDI = ( ONE / ANORM ) / AINVNM
               END IF
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
!+    TEST 1
!                    Solve and compute residual for op(A)*x = b.
!
                     SRNAMT = 'DLARHS'
                     CALL DLARHS( PATH, XTYPE, UPLO, TRANS, N, N, KD, &
                                  IDIAG, NRHS, AB, LDAB, XACT, LDA, &
                                  B, LDA, ISEED, INFO )
                     XTYPE = 'C'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
                     SRNAMT = 'DTBTRS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
                                  LDAB, X, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DTBTRS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                    Check error code from DTBTRS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'DTBTRS', INFO, 0, &
                                     UPLO // TRANS // DIAG, N, N, KD, &
                                     KD, NRHS, IMAT, NFAIL, NERRS, &
                                     NOUT )
!
                     CALL DTBT02( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
                                  LDAB, X, LDA, B, LDA, WORK, &
                                  RESULT( 1 ) )
!
!+    TEST 2
!                    Check solution from generated exact solution.
!
                     CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 2 ) )
!
!+    TESTS 3, 4, and 5
!                    Use iterative refinement to improve the solution
!                    and compute error bounds.
!
                     SRNAMT = 'DTBRFS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
                                  LDAB, B, LDA, X, LDA, RWORK, &
                                  RWORK( NRHS+1 ), WORK, IWORK, &
                                  INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DTBRFS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                    Check error code from DTBRFS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'DTBRFS', INFO, 0, &
                                     UPLO // TRANS // DIAG, N, N, KD, &
                                     KD, NRHS, IMAT, NFAIL, NERRS, &
                                     NOUT )
!
                     CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 3 ) )
                     CALL DTBT05( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
                                  LDAB, B, LDA, X, LDA, XACT, LDA, &
                                  RWORK, RWORK( NRHS+1 ), &
                                  RESULT( 4 ) )
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                     DO K = 1, 5
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )UPLO, TRANS, &
                              DIAG, N, KD, NRHS, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
                     ENDDO
                     NRUN = NRUN + 5
                  ENDDO
               ENDDO
!
!+    TEST 6
!                    Get an estimate of RCOND = 1/CNDNUM.
!
               DO ITRAN = 1, 2
                  IF( ITRAN == 1 ) THEN
                     NORM = 'O'
                     RCONDC = RCONDO
                  ELSE
                     NORM = 'I'
                     RCONDC = RCONDI
                  END IF
                  SRNAMT = 'DTBCON'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, &
                               RCOND, WORK, IWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DTBCON : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from DTBCON.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'DTBCON', INFO, 0, &
                                  NORM // UPLO // DIAG, N, N, KD, KD, &
                                  -1, IMAT, NFAIL, NERRS, NOUT )
!
                  CALL DTBT06( RCOND, RCONDC, UPLO, DIAG, N, KD, AB, &
                               LDAB, RWORK, RESULT( 6 ) )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  IF( RESULT( 6 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 ) 'DTBCON', NORM, UPLO, &
                        DIAG, N, KD, IMAT, 6, RESULT( 6 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
               ENDDO
            ENDDO
90       CONTINUE
         ENDDO
!
!           Use pathological test matrices to test DLATBS.
!
         DO IMAT = NTYPE1 + 1, NIMAT2
!
!              Do the tests only if DOTYPE( IMAT ) is true.
!
            IF( .NOT.DOTYPE( IMAT ) ) &
               GO TO 120
!
            DO IUPLO = 1, 2
!
!                 Do first for UPLO = 'U', then for UPLO = 'L'
!
               UPLO = UPLOS( IUPLO )
               DO ITRAN = 1, NTRAN
!
!                    Do for op(A) = A, A**T, and A**H.
!
                  TRANS = TRANSS( ITRAN )
!
!                    Call DLATTB to generate a triangular test matrix.
!
                  SRNAMT = 'DLATTB'
                  CALL DLATTB( IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, &
                               AB, LDAB, X, WORK, INFO )
!
!+    TEST 7
!                    Solve the system op(A)*x = b
!
                  SRNAMT = 'DLATBS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DCOPY( N, X, 1, B, 1 )
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
                  CALL DLATBS( UPLO, TRANS, DIAG, 'N', N, KD, AB, &
                               LDAB, B, SCALE, RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLATBS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from DLATBS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'DLATBS', INFO, 0, &
                                  UPLO // TRANS // DIAG // 'N', N, N, &
                                  KD, KD, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
                  CALL DTBT03( UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, &
                               SCALE, RWORK, ONE, B, LDA, X, LDA, &
                               WORK, RESULT( 7 ) )
!
!+    TEST 8
!                    Solve op(A)*x = b again with NORMIN = 'Y'.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DCOPY( N, X, 1, B, 1 )
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
                  CALL DLATBS( UPLO, TRANS, DIAG, 'Y', N, KD, AB, &
                               LDAB, B, SCALE, RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLATBS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from DLATBS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'DLATBS', INFO, 0, &
                                  UPLO // TRANS // DIAG // 'Y', N, N, &
                                  KD, KD, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
                  CALL DTBT03( UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, &
                               SCALE, RWORK, ONE, B, LDA, X, LDA, &
                               WORK, RESULT( 8 ) )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  IF( RESULT( 7 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )'DLATBS', UPLO, TRANS, &
                        DIAG, 'N', N, KD, IMAT, 7, RESULT( 7 )
                     NFAIL = NFAIL + 1
                  END IF
                  IF( RESULT( 8 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )'DLATBS', UPLO, TRANS, &
                        DIAG, 'Y', N, KD, IMAT, 8, RESULT( 8 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 2
                  ENDDO
               ENDDO
  120       CONTINUE
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' UPLO=''', A1, ''', TRANS=''', A1, ''', &
         DIAG=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I5, &
         ', type ', I2, ', test(', I2, ')=', G12.5 )
 9998 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''',', &
         I5, ',', I5, ',  ... ), type ', I2, ', test(', I2, ')=', &
         G12.5 )
 9997 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''', ''', &
         A1, ''',', I5, ',', I5, ', ...  ),  type ', I2, ', test(', &
         I1, ')=', G12.5 )
   RETURN
!
!     End of DCHKTB
!
END
                                                                                                                                                                                                                                                                                                            




