!> \brief \b CCHKTB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, AB, AINV, B, X, XACT, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NSVAL( * ), NVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            AB( * ), AINV( * ), B( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKTB tests CTBTRS, -RFS, and -CON, and CLATBS.
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
!>          The leading dimension of the work arrays.
!>          NMAX >= the maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
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
!> \ingroup complex_lin
!
!  =====================================================================
   SUBROUTINE CCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      NMAX, AB, AINV, B, X, XACT, WORK, RWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NMAX, NN, NNS, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            NSVAL( * ), NVAL( * )
   REAL               RWORK( * )
   COMPLEX            AB( * ), AINV( * ), B( * ), WORK( * ), X( * ), &
                      XACT( * )
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
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   CHARACTER          DIAG, NORM, TRANS, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IDIAG, IK, IMAT, IN, INFO, IRHS, ITRAN, &
                      IUPLO, J, K, KD, LDA, LDAB, N, NERRS, NFAIL, &
                      NIMAT, NIMAT2, NK, NRHS, NRUN
   REAL               AINVNM, ANORM, RCOND, RCONDC, RCONDI, RCONDO, &
                      SCALE
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          TRANSS( NTRAN ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANTB, CLANTR
   EXTERNAL           LSAME, CLANTB, CLANTR
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, CCOPY, CERRTR, CGET04, &
                      CLACPY, CLARHS, CLASET, CLATBS, CLATTB, CTBCON, &
                      CTBRFS, CTBSV, CTBT02, CTBT03, CTBT05, CTBT06, &
                      CTBTRS
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
   INTRINSIC          CMPLX, MAX, MIN
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' / , TRANSS / 'N', 'T', 'C' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Complex precision'
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
      CALL CERRTR( PATH, NOUT )
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
!                 Call CLATTB to generate a triangular test matrix.
!
               SRNAMT = 'CLATTB'
               CALL CLATTB( IMAT, UPLO, 'No transpose', DIAG, ISEED, &
                            N, KD, AB, LDAB, X, WORK, RWORK, INFO )
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
               CALL CLASET( 'Full', N, N, CMPLX( ZERO ), &
                            CMPLX( ONE ), AINV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( LSAME( UPLO, 'U' ) ) THEN
                  DO J = 1, N
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CTBSV( UPLO, 'No transpose', DIAG, J, KD, &
                                 AB, LDAB, AINV( ( J-1 )*LDA+1 ), 1 )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CTBSV : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDDO
               ELSE
                  DO J = 1, N
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CTBSV( UPLO, 'No transpose', DIAG, N-J+1, &
                                 KD, AB( ( J-1 )*LDAB+1 ), LDAB, &
                                 AINV( ( J-1 )*LDA+J ), 1 )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CTBSV : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDDO
               END IF
!
!                 Compute the 1-norm condition number of A.
!
               ANORM = CLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, &
                       RWORK )
               AINVNM = CLANTR( '1', UPLO, DIAG, N, N, AINV, LDA, &
                        RWORK )
               IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDO = ONE
               ELSE
                  RCONDO = ( ONE / ANORM ) / AINVNM
               END IF
!
!                 Compute the infinity-norm condition number of A.
!
               ANORM = CLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, &
                       RWORK )
               AINVNM = CLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, &
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
                     SRNAMT = 'CLARHS'
                     CALL CLARHS( PATH, XTYPE, UPLO, TRANS, N, N, KD, &
                                  IDIAG, NRHS, AB, LDAB, XACT, LDA, &
                                  B, LDA, ISEED, INFO )
                     XTYPE = 'C'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
                     SRNAMT = 'CTBTRS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
                                  LDAB, X, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CTBTRS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                    Check error code from CTBTRS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'CTBTRS', INFO, 0, &
                                     UPLO // TRANS // DIAG, N, N, KD, &
                                     KD, NRHS, IMAT, NFAIL, NERRS, &
                                     NOUT )
!
                     CALL CTBT02( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
                                  LDAB, X, LDA, B, LDA, WORK, RWORK, &
                                  RESULT( 1 ) )
!
!+    TEST 2
!                    Check solution from generated exact solution.
!
                     CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 2 ) )
!
!+    TESTS 3, 4, and 5
!                    Use iterative refinement to improve the solution
!                    and compute error bounds.
!
                     SRNAMT = 'CTBRFS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL CTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
                                  LDAB, B, LDA, X, LDA, RWORK, &
                                  RWORK( NRHS+1 ), WORK, &
                                  RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : CTBRFS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                    Check error code from CTBRFS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'CTBRFS', INFO, 0, &
                                     UPLO // TRANS // DIAG, N, N, KD, &
                                     KD, NRHS, IMAT, NFAIL, NERRS, &
                                     NOUT )
!
                     CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 3 ) )
                     CALL CTBT05( UPLO, TRANS, DIAG, N, KD, NRHS, AB, &
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
                  SRNAMT = 'CTBCON'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, &
                               RCOND, WORK, RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CTBCON : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from CTBCON.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'CTBCON', INFO, 0, &
                                  NORM // UPLO // DIAG, N, N, KD, KD, &
                                  -1, IMAT, NFAIL, NERRS, NOUT )
!
                  CALL CTBT06( RCOND, RCONDC, UPLO, DIAG, N, KD, AB, &
                               LDAB, RWORK, RESULT( 6 ) )
!
!                    Print the test ratio if it is  >=  THRESH.
!
                  IF( RESULT( 6 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 ) 'CTBCON', NORM, UPLO, &
                        DIAG, N, KD, IMAT, 6, RESULT( 6 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
               ENDDO
            ENDDO
90       CONTINUE
         ENDDO
!
!           Use pathological test matrices to test CLATBS.
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
!                    Call CLATTB to generate a triangular test matrix.
!
                  SRNAMT = 'CLATTB'
                  CALL CLATTB( IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, &
                               AB, LDAB, X, WORK, RWORK, INFO )
!
!+    TEST 7
!                    Solve the system op(A)*x = b
!
                  SRNAMT = 'CLATBS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CCOPY( N, X, 1, B, 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLATBS( UPLO, TRANS, DIAG, 'N', N, KD, AB, &
                               LDAB, B, SCALE, RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLATBS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from CLATBS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'CLATBS', INFO, 0, &
                                  UPLO // TRANS // DIAG // 'N', N, N, &
                                  KD, KD, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
                  CALL CTBT03( UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, &
                               SCALE, RWORK, ONE, B, LDA, X, LDA, &
                               WORK, RESULT( 7 ) )
!
!+    TEST 8
!                    Solve op(A)*x = b again with NORMIN = 'Y'.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CCOPY( N, X, 1, B, 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CCOPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL CLATBS( UPLO, TRANS, DIAG, 'Y', N, KD, AB, &
                               LDAB, B, SCALE, RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : CLATBS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from CLATBS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'CLATBS', INFO, 0, &
                                  UPLO // TRANS // DIAG // 'Y', N, N, &
                                  KD, KD, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
                  CALL CTBT03( UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, &
                               SCALE, RWORK, ONE, B, LDA, X, LDA, &
                               WORK, RESULT( 8 ) )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  IF( RESULT( 7 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )'CLATBS', UPLO, TRANS, &
                        DIAG, 'N', N, KD, IMAT, 7, RESULT( 7 )
                     NFAIL = NFAIL + 1
                  END IF
                  IF( RESULT( 8 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )'CLATBS', UPLO, TRANS, &
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
!     End of CCHKTB
!
END
                                                                                                                                                                                                                                                                                                            



