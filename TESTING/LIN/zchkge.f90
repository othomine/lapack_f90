!> \brief \b ZCHKGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS,
!                          NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B,
!                          X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NMAX, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ),
!      $                   NVAL( * )
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
!> ZCHKGE tests ZGETRF, -TRI, -TRS, -RFS, and -CON.
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
!>          The maximum value permitted for M or N, used in dimensioning
!>          the work arrays.
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
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      (max(2*NMAX,2*NSMAX+NWORK))
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
   SUBROUTINE ZCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS, &
                      NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, &
                      X, XACT, WORK, RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NM, NMAX, NN, NNB, NNS, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), &
                      NVAL( * )
   DOUBLE PRECISION   RWORK( * )
   COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 11 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 8 )
   INTEGER            NTRAN
   PARAMETER          ( NTRAN = 3 )
!     ..
!     .. Local Scalars ..
   LOGICAL            TRFCON, ZEROT
   CHARACTER          DIST, NORM, TRANS, TYPE, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN, &
                      IZERO, K, KL, KU, LDA, LWORK, M, MODE, N, NB, &
                      NERRS, NFAIL, NIMAT, NRHS, NRUN, NT
   DOUBLE PRECISION   AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, DUMMY, &
                      RCOND, RCONDC, RCONDI, RCONDO
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          TRANSS( NTRAN )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DGET06, ZLANGE
   EXTERNAL           DGET06, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, XLAENV, ZERRGE, ZGECON, &
                      ZGERFS, ZGET01, ZGET02, ZGET03, ZGET04, ZGET07, &
                      ZGETRF, ZGETRI, ZGETRS, ZLACPY, ZLARHS, ZLASET, &
                      ZLATB4, ZLATMS
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DCMPLX, MAX, MIN
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
   DATA               ISEEDY / 1988, 1989, 1990, 1991 / , &
                      TRANSS / 'N', 'T', 'C' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Zomplex precision'
   PATH( 2: 3 ) = 'GE'
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
!     Test the error exits
!
   CALL XLAENV( 1, 1 )
   IF( TSTERR ) &
      CALL ZERRGE( PATH, NOUT )
   INFOT = 0
   CALL XLAENV( 2, 2 )
!
!     Do for each value of M in MVAL
!
   DO IM = 1, NM
      M = MVAL( IM )
      LDA = MAX( 1, M )
!
!        Do for each value of N in NVAL
!
      DO IN = 1, NN
         N = NVAL( IN )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( M <= 0 .OR. N <= 0 ) &
            NIMAT = 1
!
         DO IMAT = 1, NIMAT
!
!              Do the tests only if DOTYPE( IMAT ) is true.
!
            IF( .NOT.DOTYPE( IMAT ) ) &
               GO TO 100
!
!              Skip types 5, 6, or 7 if the matrix size is too small.
!
            ZEROT = IMAT >= 5 .AND. IMAT <= 7
            IF( ZEROT .AND. N < IMAT-4 ) &
               GO TO 100
!
!              Set up parameters with ZLATB4 and generate a test matrix
!              with ZLATMS.
!
            CALL ZLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
            SRNAMT = 'ZLATMS'
            CALL ZLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, &
                         WORK, INFO )
!
!              Check error code from ZLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 100
            END IF
!
!              For types 5-7, zero one or more columns of the matrix to
!              test that INFO is returned correctly.
!
            IF( ZEROT ) THEN
               IF( IMAT == 5 ) THEN
                  IZERO = 1
               ELSE IF( IMAT == 6 ) THEN
                  IZERO = MIN( M, N )
               ELSE
                  IZERO = MIN( M, N ) / 2 + 1
               END IF
               IOFF = ( IZERO-1 )*LDA
               IF( IMAT < 7 ) THEN
                  DO I = 1, M
                     A( IOFF+I ) = ZERO
                  ENDDO
               ELSE
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLASET( 'Full', M, N-IZERO+1, DCMPLX( ZERO ), &
                               DCMPLX( ZERO ), A( IOFF+1 ), LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               END IF
            ELSE
               IZERO = 0
            END IF
!
!              These lines, if used in place of the calls in the DO 60
!              loop, cause the code to bomb on a Sun SPARCstation.
!
!               ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
!               ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
!
!              Do for each blocksize in NBVAL
!
            DO INB = 1, NNB
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
!
!                 Compute the LU factorization of the matrix.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( 'Full', M, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'ZGETRF'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZGETRF( M, N, AFAC, LDA, IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZGETRF : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from ZGETRF.
!
               IF( INFO /= IZERO ) &
                  CALL ALAERH( PATH, 'ZGETRF', INFO, IZERO, ' ', M, &
                               N, -1, -1, NB, IMAT, NFAIL, NERRS, &
                               NOUT )
               TRFCON = .FALSE.
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLACPY( 'Full', M, N, AFAC, LDA, AINV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               CALL ZGET01( M, N, A, LDA, AINV, LDA, IWORK, RWORK, &
                            RESULT( 1 ) )
               NT = 1
!
!+    TEST 2
!                 Form the inverse if the factorization was successful
!                 and compute the residual.
!
               IF( M == N .AND. INFO == 0 ) THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLACPY( 'Full', N, N, AFAC, LDA, AINV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  SRNAMT = 'ZGETRI'
                  NRHS = NSVAL( 1 )
                  LWORK = NMAX*MAX( 3, NRHS )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZGETRI( N, AINV, LDA, IWORK, WORK, LWORK, &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZGETRI : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from ZGETRI.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'ZGETRI', INFO, 0, ' ', N, N, &
                                  -1, -1, NB, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
!                    Compute the residual for the matrix times its
!                    inverse.  Also compute the 1-norm condition number
!                    of A.
!
                  CALL ZGET03( N, A, LDA, AINV, LDA, WORK, LDA, &
                               RWORK, RCONDO, RESULT( 2 ) )
                  ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
!
!                    Compute the infinity-norm condition number of A.
!
                  ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
                  AINVNM = ZLANGE( 'I', N, N, AINV, LDA, RWORK )
                  IF( ANORMI <= ZERO .OR. AINVNM <= ZERO ) THEN
                     RCONDI = ONE
                  ELSE
                     RCONDI = ( ONE / ANORMI ) / AINVNM
                  END IF
                  NT = 2
               ELSE
!
!                    Do only the condition estimate if INFO > 0.
!
                  TRFCON = .TRUE.
                  ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
                  ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
                  RCONDO = ZERO
                  RCONDI = ZERO
               END IF
!
!                 Print information about the tests so far that did not
!                 pass the threshold.
!
               DO K = 1, NT
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9999 )M, N, NB, IMAT, K, &
                        RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
               ENDDO
               NRUN = NRUN + NT
!
!                 Skip the remaining tests if this is not the first
!                 block size or if M .ne. N.  Skip the solve tests if
!                 the matrix is singular.
!
               IF( INB > 1 .OR. M /= N ) &
                  GO TO 90
               IF( TRFCON ) &
                  GO TO 70
!
               DO IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
                  XTYPE = 'N'
!
                  DO ITRAN = 1, NTRAN
                     TRANS = TRANSS( ITRAN )
                     IF( ITRAN == 1 ) THEN
                        RCONDC = RCONDO
                     ELSE
                        RCONDC = RCONDI
                     END IF
!
!+    TEST 3
!                       Solve and compute residual for A * X = B.
!
                     SRNAMT = 'ZLARHS'
                     CALL ZLARHS( PATH, XTYPE, ' ', TRANS, N, N, KL, &
                                  KU, NRHS, A, LDA, XACT, LDA, B, &
                                  LDA, ISEED, INFO )
                     XTYPE = 'C'
!
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
                     SRNAMT = 'ZGETRS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL ZGETRS( TRANS, N, NRHS, AFAC, LDA, IWORK, &
                                  X, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : ZGETRS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from ZGETRS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'ZGETRS', INFO, 0, TRANS, &
                                     N, N, -1, -1, NRHS, IMAT, NFAIL, &
                                     NERRS, NOUT )
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL ZGET02( TRANS, N, N, NRHS, A, LDA, X, LDA, &
                                  WORK, LDA, RWORK, RESULT( 3 ) )
!
!+    TEST 4
!                       Check solution from generated exact solution.
!
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 4 ) )
!
!+    TESTS 5, 6, and 7
!                       Use iterative refinement to improve the
!                       solution.
!
                     SRNAMT = 'ZGERFS'
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL ZGERFS( TRANS, N, NRHS, A, LDA, AFAC, LDA, &
                                  IWORK, B, LDA, X, LDA, RWORK, &
                                  RWORK( NRHS+1 ), WORK, &
                                  RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : ZGERFS : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from ZGERFS.
!
                     IF( INFO /= 0 ) &
                        CALL ALAERH( PATH, 'ZGERFS', INFO, 0, TRANS, &
                                     N, N, -1, -1, NRHS, IMAT, NFAIL, &
                                     NERRS, NOUT )
!
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 5 ) )
                     CALL ZGET07( TRANS, N, NRHS, A, LDA, B, LDA, X, &
                                  LDA, XACT, LDA, RWORK, .TRUE., &
                                  RWORK( NRHS+1 ), RESULT( 6 ) )
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                     DO K = 3, 7
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, &
                              IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
                     ENDDO
                     NRUN = NRUN + 5
                  ENDDO
               ENDDO
!
!+    TEST 8
!                    Get an estimate of RCOND = 1/CNDNUM.
!
70             CONTINUE
               DO ITRAN = 1, 2
                  IF( ITRAN == 1 ) THEN
                     ANORM = ANORMO
                     RCONDC = RCONDO
                     NORM = 'O'
                  ELSE
                     ANORM = ANORMI
                     RCONDC = RCONDI
                     NORM = 'I'
                  END IF
                  SRNAMT = 'ZGECON'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZGECON( NORM, N, AFAC, LDA, ANORM, RCOND, &
                               WORK, RWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZGECON : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                       Check error code from ZGECON.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( PATH, 'ZGECON', INFO, 0, NORM, N, &
                                  N, -1, -1, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
!                       This line is needed on a Sun SPARCstation.
!
                  DUMMY = RCOND
!
                  RESULT( 8 ) = DGET06( RCOND, RCONDC )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  IF( RESULT( 8 ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 8, &
                        RESULT( 8 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
               ENDDO
90          CONTINUE
            ENDDO
  100       CONTINUE
            ENDDO
!
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' M = ', I5, ', N =', I5, ', NB =', I4, ', type ', I2, &
         ', test(', I2, ') =', G12.5 )
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2, &
         ', test(', I2, ') =', G12.5 )
   RETURN
!
!     End of ZCHKGE
!
END
                                                                                                                                                                                                                                                                                                            




