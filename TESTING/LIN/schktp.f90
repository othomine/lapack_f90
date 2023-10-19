!> \brief \b SCHKTP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKTP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, AP, AINVP, B, X, XACT, WORK, RWORK,
!                          IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
!       REAL               AINVP( * ), AP( * ), B( * ), RWORK( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKTP tests STPTRI, -TRS, -RFS, and -CON, and SLATPS
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
!>          The leading dimension of the work arrays.  NMAX >= the
!>          maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] AP
!> \verbatim
!>          AP is REAL array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AINVP
!> \verbatim
!>          AINVP is REAL array, dimension
!>                      (NMAX*(NMAX+1)/2)
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
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (NMAX)
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
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SCHKTP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      NMAX, AP, AINVP, B, X, XACT, WORK, RWORK, &
                      IWORK, NOUT )
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
   INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
   REAL               AINVP( * ), AP( * ), B( * ), RWORK( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NTYPE1, NTYPES
   PARAMETER          ( NTYPE1 = 10, NTYPES = 18 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 9 )
   INTEGER            NTRAN
   PARAMETER          ( NTRAN = 3 )
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   CHARACTER          DIAG, NORM, TRANS, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IDIAG, IMAT, IN, INFO, IRHS, ITRAN, IUPLO, &
                      K, LAP, LDA, N, NERRS, NFAIL, NRHS, NRUN
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
   REAL               SLANTP
   EXTERNAL           LSAME, SLANTP
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, SCOPY, SERRTR, SGET04, &
                      SLACPY, SLARHS, SLATPS, SLATTP, STPCON, STPRFS, &
                      STPT01, STPT02, STPT03, STPT05, STPT06, STPTRI, &
                      STPTRS
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
   PATH( 1: 1 ) = 'Single precision'
   PATH( 2: 3 ) = 'TP'
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
      CALL SERRTR( PATH, NOUT )
   INFOT = 0
!
   DO IN = 1, NN
!
!        Do for each value of N in NVAL
!
      N = NVAL( IN )
      LDA = MAX( 1, N )
      LAP = LDA*( LDA+1 ) / 2
      XTYPE = 'N'
!
      DO IMAT = 1, NTYPE1
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 70
!
         DO IUPLO = 1, 2
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
            UPLO = UPLOS( IUPLO )
!
!              Call SLATTP to generate a triangular test matrix.
!
            SRNAMT = 'SLATTP'
            CALL SLATTP( IMAT, UPLO, 'No transpose', DIAG, ISEED, N, &
                         AP, X, WORK, INFO )
!
!              Set IDIAG = 1 for non-unit matrices, 2 for unit.
!
            IF( LSAME( DIAG, 'N' ) ) THEN
               IDIAG = 1
            ELSE
               IDIAG = 2
            END IF
!
!+    TEST 1
!              Form the inverse of A.
!
            IF( N > 0 )  THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SCOPY( LAP, AP, 1, AINVP, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            ENDIF
            SRNAMT = 'STPTRI'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL STPTRI( UPLO, DIAG, N, AINVP, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : STPTRI : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from STPTRI.
!
            IF( INFO /= 0 ) &
               CALL ALAERH( PATH, 'STPTRI', INFO, 0, UPLO // DIAG, N, &
                            N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
!              Compute the infinity-norm condition number of A.
!
            ANORM = SLANTP( 'I', UPLO, DIAG, N, AP, RWORK )
            AINVNM = SLANTP( 'I', UPLO, DIAG, N, AINVP, RWORK )
            IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
               RCONDI = ONE
            ELSE
               RCONDI = ( ONE / ANORM ) / AINVNM
            END IF
!
!              Compute the residual for the triangular matrix times its
!              inverse.  Also compute the 1-norm condition number of A.
!
            CALL STPT01( UPLO, DIAG, N, AP, AINVP, RCONDO, RWORK, &
                         RESULT( 1 ) )
!
!              Print the test ratio if it is  >=  THRESH.
!
            IF( RESULT( 1 ) >= THRESH ) THEN
               IF( NFAIL == 0 .AND. NERRS == 0 ) &
                  CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )UPLO, DIAG, N, IMAT, 1, &
                  RESULT( 1 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
!
            DO IRHS = 1, NNS
               NRHS = NSVAL( IRHS )
               XTYPE = 'N'
!
               DO ITRAN = 1, NTRAN
!
!                 Do for op(A) = A, A**T, or A**H.
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
!                 Solve and compute residual for op(A)*x = b.
!
                  SRNAMT = 'SLARHS'
                  CALL SLARHS( PATH, XTYPE, UPLO, TRANS, N, N, 0, &
                               IDIAG, NRHS, AP, LAP, XACT, LDA, B, &
                               LDA, ISEED, INFO )
                  XTYPE = 'C'
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
                  SRNAMT = 'STPTRS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL STPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, X, &
                               LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : STPTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                 Check error code from STPTRS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'STPTRS', INFO, 0, &
                                  UPLO // TRANS // DIAG, N, N, -1, &
                                  -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
                  CALL STPT02( UPLO, TRANS, DIAG, N, NRHS, AP, X, &
                               LDA, B, LDA, WORK, RESULT( 2 ) )
!
!+    TEST 3
!                 Check solution from generated exact solution.
!
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 3 ) )
!
!+    TESTS 4, 5, and 6
!                 Use iterative refinement to improve the solution and
!                 compute error bounds.
!
                  SRNAMT = 'STPRFS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL STPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, &
                               LDA, X, LDA, RWORK, RWORK( NRHS+1 ), &
                               WORK, IWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : STPRFS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                 Check error code from STPRFS.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'STPRFS', INFO, 0, &
                                  UPLO // TRANS // DIAG, N, N, -1, &
                                  -1, NRHS, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 4 ) )
                  CALL STPT05( UPLO, TRANS, DIAG, N, NRHS, AP, B, &
                               LDA, X, LDA, XACT, LDA, RWORK, &
                               RWORK( NRHS+1 ), RESULT( 5 ) )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = 2, 6
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9998 )UPLO, TRANS, DIAG, &
                           N, NRHS, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                  ENDDO
                  NRUN = NRUN + 5
               ENDDO
            ENDDO
!
!+    TEST 7
!                 Get an estimate of RCOND = 1/CNDNUM.
!
            DO ITRAN = 1, 2
               IF( ITRAN == 1 ) THEN
                  NORM = 'O'
                  RCONDC = RCONDO
               ELSE
                  NORM = 'I'
                  RCONDC = RCONDI
               END IF
!
               SRNAMT = 'STPCON'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL STPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, &
                            IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : STPCON : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from STPCON.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'STPCON', INFO, 0, &
                               NORM // UPLO // DIAG, N, N, -1, -1, &
                               -1, IMAT, NFAIL, NERRS, NOUT )
!
               CALL STPT06( RCOND, RCONDC, UPLO, DIAG, N, AP, RWORK, &
                            RESULT( 7 ) )
!
!                 Print the test ratio if it is  >=  THRESH.
!
               IF( RESULT( 7 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9997 ) 'STPCON', NORM, UPLO, &
                     DIAG, N, IMAT, 7, RESULT( 7 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
            ENDDO
         ENDDO
70    CONTINUE
      ENDDO
!
!        Use pathological test matrices to test SLATPS.
!
      DO IMAT = NTYPE1 + 1, NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 100
!
         DO IUPLO = 1, 2
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
            UPLO = UPLOS( IUPLO )
            DO ITRAN = 1, NTRAN
!
!                 Do for op(A) = A, A**T, or A**H.
!
               TRANS = TRANSS( ITRAN )
!
!                 Call SLATTP to generate a triangular test matrix.
!
               SRNAMT = 'SLATTP'
               CALL SLATTP( IMAT, UPLO, TRANS, DIAG, ISEED, N, AP, X, &
                            WORK, INFO )
!
!+    TEST 8
!                 Solve the system op(A)*x = b.
!
               SRNAMT = 'SLATPS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SCOPY( N, X, 1, B, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLATPS( UPLO, TRANS, DIAG, 'N', N, AP, B, SCALE, &
                            RWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLATPS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from SLATPS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'SLATPS', INFO, 0, &
                               UPLO // TRANS // DIAG // 'N', N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               CALL STPT03( UPLO, TRANS, DIAG, N, 1, AP, SCALE, &
                            RWORK, ONE, B, LDA, X, LDA, WORK, &
                            RESULT( 8 ) )
!
!+    TEST 9
!                 Solve op(A)*x = b again with NORMIN = 'Y'.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SCOPY( N, X, 1, B( N+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLATPS( UPLO, TRANS, DIAG, 'Y', N, AP, B( N+1 ), &
                            SCALE, RWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLATPS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from SLATPS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'SLATPS', INFO, 0, &
                               UPLO // TRANS // DIAG // 'Y', N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               CALL STPT03( UPLO, TRANS, DIAG, N, 1, AP, SCALE, &
                            RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, &
                            RESULT( 9 ) )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               IF( RESULT( 8 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9996 )'SLATPS', UPLO, TRANS, &
                     DIAG, 'N', N, IMAT, 8, RESULT( 8 )
                  NFAIL = NFAIL + 1
               END IF
               IF( RESULT( 9 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9996 )'SLATPS', UPLO, TRANS, &
                     DIAG, 'Y', N, IMAT, 9, RESULT( 9 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 2
            ENDDO
         ENDDO
  100    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' UPLO=''', A1, ''', DIAG=''', A1, ''', N=', I5, &
         ', type ', I2, ', test(', I2, ')= ', G12.5 )
 9998 FORMAT( ' UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, &
         ''', N=', I5, ''', NRHS=', I5, ', type ', I2, ', test(', &
         I2, ')= ', G12.5 )
 9997 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''',', &
         I5, ', ... ), type ', I2, ', test(', I2, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''', ''', &
         A1, ''',', I5, ', ... ), type ', I2, ', test(', I2, ')=', &
         G12.5 )
   RETURN
!
!     End of SCHKTP
!
END
                                                                                                                                                                                                                                                                                                            




