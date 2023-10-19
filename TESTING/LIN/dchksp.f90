!> \brief \b DCHKSP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKSP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK,
!                          IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKSP tests DSPTRF, -TRI, -TRS, -RFS, and -CON
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
!>          A is DOUBLE PRECISION array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension
!>                      (NMAX*(NMAX+1)/2)
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
!>                      (NMAX*max(2,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array,
!>                                 dimension (NMAX+2*NSMAX)
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DCHKSP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, &
                      IWORK, NOUT )
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
   DOUBLE PRECISION   A( * ), AFAC( * ), AINV( * ), B( * ), &
                      RWORK( * ), WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 10 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 8 )
!     ..
!     .. Local Scalars ..
   LOGICAL            TRFCON, ZEROT
   CHARACTER          DIST, PACKIT, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, I1, I2, IMAT, IN, INFO, IOFF, IRHS, IUPLO, &
                      IZERO, J, K, KL, KU, LDA, MODE, N, NERRS, &
                      NFAIL, NIMAT, NPP, NRHS, NRUN, NT
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
   LOGICAL            LSAME
   DOUBLE PRECISION   DGET06, DLANSP
   EXTERNAL           LSAME, DGET06, DLANSP
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, DCOPY, DERRSY, DGET04, &
                      DLACPY, DLARHS, DLATB4, DLATMS, DPPT02, DPPT03, &
                      DPPT05, DSPCON, DSPRFS, DSPT01, DSPTRF, DSPTRI, &
                      DSPTRS
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
   PATH( 1: 1 ) = 'Double precision'
   PATH( 2: 3 ) = 'SP'
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
      CALL DERRSY( PATH, NOUT )
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
      IZERO = 0
      DO IMAT = 1, NIMAT
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( .NOT.DOTYPE( IMAT ) ) &
            GO TO 160
!
!           Skip types 3, 4, 5, or 6 if the matrix size is too small.
!
         ZEROT = IMAT >= 3 .AND. IMAT <= 6
         IF( ZEROT .AND. N < IMAT-2 ) &
            GO TO 160
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
            IF( LSAME( UPLO, 'U' ) ) THEN
               PACKIT = 'C'
            ELSE
               PACKIT = 'R'
            END IF
!
!              Set up parameters with DLATB4 and generate a test matrix
!              with DLATMS.
!
            CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                         CNDNUM, DIST )
!
            SRNAMT = 'DLATMS'
            CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                         CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, &
                         INFO )
!
!              Check error code from DLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 150
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
                        IOFF = IOFF + J
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
                        IOFF = IOFF + N - J
                     ENDDO
                  END IF
               END IF
            ELSE
               IZERO = 0
            END IF
!
!              Compute the L*D*L' or U*D*U' factorization of the matrix.
!
            NPP = N*( N+1 ) / 2
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DCOPY( NPP, A, 1, AFAC, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            SRNAMT = 'DSPTRF'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSPTRF( UPLO, N, AFAC, IWORK, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSPTRF : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Adjust the expected value of INFO to account for
!              pivoting.
!
            K = IZERO
            IF( K > 0 ) THEN
  100             CONTINUE
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
!              Check error code from DSPTRF.
!
            IF( INFO /= K ) &
               CALL ALAERH( PATH, 'DSPTRF', INFO, K, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
            IF( INFO /= 0 ) THEN
               TRFCON = .TRUE.
            ELSE
               TRFCON = .FALSE.
            END IF
!
!+    TEST 1
!              Reconstruct matrix from factors and compute residual.
!
            CALL DSPT01( UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK, &
                         RESULT( 1 ) )
            NT = 1
!
!+    TEST 2
!              Form the inverse and compute the residual.
!
            IF( .NOT.TRFCON ) THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( NPP, AFAC, 1, AINV, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'DSPTRI'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DSPTRI( UPLO, N, AINV, IWORK, WORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DSPTRI : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!              Check error code from DSPTRI.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'DSPTRI', INFO, 0, UPLO, N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
               CALL DPPT03( UPLO, N, A, AINV, WORK, LDA, RWORK, &
                            RCONDC, RESULT( 2 ) )
               NT = 2
            END IF
!
!              Print information about the tests that did not pass
!              the threshold.
!
            DO K = 1, NT
               IF( RESULT( K ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, K, &
                     RESULT( K )
                  NFAIL = NFAIL + 1
               END IF
               ENDDO
            NRUN = NRUN + NT
!
!              Do only the condition estimate if INFO is not 0.
!
            IF( TRFCON ) THEN
               RCONDC = ZERO
               GO TO 140
            END IF
!
            DO IRHS = 1, NNS
               NRHS = NSVAL( IRHS )
!
!+    TEST 3
!              Solve and compute residual for  A * X = B.
!
               SRNAMT = 'DLARHS'
               CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                            NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, &
                            INFO )
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
               SRNAMT = 'DSPTRS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DSPTRS( UPLO, N, NRHS, AFAC, IWORK, X, LDA, &
                            INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DSPTRS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!              Check error code from DSPTRS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'DSPTRS', INFO, 0, UPLO, N, N, &
                               -1, -1, NRHS, IMAT, NFAIL, NERRS, &
                               NOUT )
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               CALL DPPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, &
                            RWORK, RESULT( 3 ) )
!
!+    TEST 4
!              Check solution from generated exact solution.
!
               CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 4 ) )
!
!+    TESTS 5, 6, and 7
!              Use iterative refinement to improve the solution.
!
               SRNAMT = 'DSPRFS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DSPRFS( UPLO, N, NRHS, A, AFAC, IWORK, B, LDA, X, &
                            LDA, RWORK, RWORK( NRHS+1 ), WORK, &
                            IWORK( N+1 ), INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DSPRFS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!              Check error code from DSPRFS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'DSPRFS', INFO, 0, UPLO, N, N, &
                               -1, -1, NRHS, IMAT, NFAIL, NERRS, &
                               NOUT )
!
               CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 5 ) )
               CALL DPPT05( UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, &
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
  140          CONTINUE
            ANORM = DLANSP( '1', UPLO, N, A, RWORK )
            SRNAMT = 'DSPCON'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSPCON( UPLO, N, AFAC, IWORK, ANORM, RCOND, WORK, &
                         IWORK( N+1 ), INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSPCON : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from DSPCON.
!
            IF( INFO /= 0 ) &
               CALL ALAERH( PATH, 'DSPCON', INFO, 0, UPLO, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
            RESULT( 8 ) = DGET06( RCOND, RCONDC )
!
!              Print the test ratio if it is  >=  THRESH.
!
            IF( RESULT( 8 ) >= THRESH ) THEN
               IF( NFAIL == 0 .AND. NERRS == 0 ) &
                  CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, 8, &
                  RESULT( 8 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
  150       CONTINUE
            ENDDO
  160    CONTINUE
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
!     End of DCHKSP
!
END
                                                                                                                                                                                                                                                                                                            




