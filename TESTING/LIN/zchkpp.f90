!> \brief \b ZCHKPP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKPP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NSVAL( * ), NVAL( * )
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
!> ZCHKPP tests ZPPTRF, -TRI, -TRS, -RFS, and -CON
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
!>          A is COMPLEX*16 array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension
!>                      (NMAX*(NMAX+1)/2)
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
   SUBROUTINE ZCHKPP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, &
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
   INTEGER            NSVAL( * ), NVAL( * )
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
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 9 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 8 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, PACKIT, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IMAT, IN, INFO, IOFF, IRHS, IUPLO, IZERO, K, &
                      KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NPP, &
                      NRHS, NRUN
   DOUBLE PRECISION   ANORM, CNDNUM, RCOND, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          PACKS( 2 ), UPLOS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DGET06, ZLANHP
   EXTERNAL           DGET06, ZLANHP
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, ZCOPY, ZERRPO, ZGET04, &
                      ZLACPY, ZLAIPD, ZLARHS, ZLATB4, ZLATMS, ZPPCON, &
                      ZPPRFS, ZPPT01, ZPPT02, ZPPT03, ZPPT05, ZPPTRF, &
                      ZPPTRI, ZPPTRS
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
   DATA               UPLOS / 'U', 'L' / , PACKS / 'C', 'R' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Zomplex precision'
   PATH( 2: 3 ) = 'PP'
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
      CALL ZERRPO( PATH, NOUT )
   INFOT = 0
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
            GO TO 100
!
!           Skip types 3, 4, or 5 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 5
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 100
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
            PACKIT = PACKS( IUPLO )
!
!              Set up parameters with ZLATB4 and generate a test matrix
!              with ZLATMS.
!
            CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
            SRNAMT = 'ZLATMS'
            CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, &
                         INFO )
!
!              Check error code from ZLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 90
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
!
!                 Set row and column IZERO of A to 0.
!
               IF( IUPLO == 1 ) THEN
                  IOFF = ( IZERO-1 )*IZERO / 2
                  DO I = 1, IZERO - 1
                     A( IOFF+I ) = ZERO
                  ENDDO
                  IOFF = IOFF + IZERO
                  DO I = IZERO, N
                     A( IOFF ) = ZERO
                     IOFF = IOFF + I
                  ENDDO
               ELSE
                  IOFF = IZERO
                  DO I = 1, IZERO - 1
                     A( IOFF ) = ZERO
                     IOFF = IOFF + N - I
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
!              Set the imaginary part of the diagonals.
!
            IF( IUPLO == 1 ) THEN
               CALL ZLAIPD( N, A, 2, 1 )
            ELSE
               CALL ZLAIPD( N, A, N, -1 )
            END IF
!
!              Compute the L*L' or U'*U factorization of the matrix.
!
            NPP = N*( N+1 ) / 2
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZCOPY( NPP, A, 1, AFAC, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            SRNAMT = 'ZPPTRF'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZPPTRF( UPLO, N, AFAC, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZPPTRF : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from ZPPTRF.
!
            IF( INFO /= IZERO ) THEN
               CALL ALAERH( PATH, 'ZPPTRF', INFO, IZERO, UPLO, N, N, &
                            -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 90
            END IF
!
!              Skip the tests if INFO is not 0.
!
            IF( INFO /= 0 ) &
               GO TO 90
!
!+    TEST 1
!              Reconstruct matrix from factors and compute residual.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZCOPY( NPP, AFAC, 1, AINV, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            CALL ZPPT01( UPLO, N, A, AINV, RWORK, RESULT( 1 ) )
!
!+    TEST 2
!              Form the inverse and compute the residual.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZCOPY( NPP, AFAC, 1, AINV, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            SRNAMT = 'ZPPTRI'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZPPTRI( UPLO, N, AINV, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZPPTRI : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from ZPPTRI.
!
            IF( INFO /= 0 ) &
               CALL ALAERH( PATH, 'ZPPTRI', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
            CALL ZPPT03( UPLO, N, A, AINV, WORK, LDA, RWORK, RCONDC, &
                         RESULT( 2 ) )
!
!              Print information about the tests that did not pass
!              the threshold.
!
            DO K = 1, 2
               IF( RESULT( K ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, K, &
                     RESULT( K )
                  NFAIL = NFAIL + 1
               END IF
            ENDDO
            NRUN = NRUN + 2
!
            DO IRHS = 1, NNS
               NRHS = NSVAL( IRHS )
!
!+    TEST 3
!              Solve and compute residual for  A * X = B.
!
               SRNAMT = 'ZLARHS'
               CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                            NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, &
                            INFO )
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
               SRNAMT = 'ZPPTRS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZPPTRS( UPLO, N, NRHS, AFAC, X, LDA, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZPPTRS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!              Check error code from ZPPTRS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'ZPPTRS', INFO, 0, UPLO, N, N, &
                               -1, -1, NRHS, IMAT, NFAIL, NERRS, &
                               NOUT )
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
               CALL ZPPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, &
                            RWORK, RESULT( 3 ) )
!
!+    TEST 4
!              Check solution from generated exact solution.
!
               CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 4 ) )
!
!+    TESTS 5, 6, and 7
!              Use iterative refinement to improve the solution.
!
               SRNAMT = 'ZPPRFS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZPPRFS( UPLO, N, NRHS, A, AFAC, B, LDA, X, LDA, &
                            RWORK, RWORK( NRHS+1 ), WORK, &
                            RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZPPRFS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!              Check error code from ZPPRFS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'ZPPRFS', INFO, 0, UPLO, N, N, &
                               -1, -1, NRHS, IMAT, NFAIL, NERRS, &
                               NOUT )
!
               CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 5 ) )
               CALL ZPPT05( UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, &
                            LDA, RWORK, RWORK( NRHS+1 ), &
                            RESULT( 6 ) )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 3, 7
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, &
                        K, RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
               ENDDO
               NRUN = NRUN + 5
            ENDDO
!
!+    TEST 8
!              Get an estimate of RCOND = 1/CNDNUM.
!
            ANORM = ZLANHP( '1', UPLO, N, A, RWORK )
            SRNAMT = 'ZPPCON'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZPPCON( UPLO, N, AFAC, ANORM, RCOND, WORK, RWORK, &
                         INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZPPCON : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from ZPPCON.
!
            IF( INFO /= 0 ) &
               CALL ALAERH( PATH, 'ZPPCON', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
            RESULT( 8 ) = DGET06( RCOND, RCONDC )
!
!              Print the test ratio if greater than or equal to THRESH.
!
            IF( RESULT( 8 ) >= THRESH ) THEN
               IF( NFAIL == 0 .AND. NERRS == 0 ) &
                  CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, 8, &
                  RESULT( 8 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
!
90       CONTINUE
         ENDDO
  100    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', type ', I2, ', test ', &
         I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test(', I2, ') =', G12.5 )
   RETURN
!
!     End of ZCHKPP
!
END
                                                                                                                                                                                                                                                                                                            




