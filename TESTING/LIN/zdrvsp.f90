!> \brief \b ZDRVSP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVSP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
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
!> ZDRVSP tests the driver routines ZSPSV and -SVX.
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
!>          B is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(2,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (NMAX+2*NRHS)
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
   SUBROUTINE ZDRVSP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, &
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
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NVAL( * )
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
   INTEGER            NTYPES, NTESTS
   PARAMETER          ( NTYPES = 11, NTESTS = 6 )
   INTEGER            NFACT
   PARAMETER          ( NFACT = 2 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, FACT, PACKIT, TYPE, UPLO, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, &
                      IZERO, J, K, K1, KL, KU, LDA, MODE, N, NB, &
                      NBMIN, NERRS, NFAIL, NIMAT, NPP, NRUN, NT
   DOUBLE PRECISION   AINVNM, ANORM, CNDNUM, RCOND, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          FACTS( NFACT )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DGET06, ZLANSP
   EXTERNAL           DGET06, ZLANSP
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, XLAENV, ZCOPY, ZERRVX, &
                      ZGET04, ZLACPY, ZLARHS, ZLASET, ZLATB4, ZLATMS, &
                      ZLATSP, ZPPT05, ZSPSV, ZSPSVX, ZSPT01, ZSPT02, &
                      ZSPTRF, ZSPTRI
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
   INTRINSIC          DCMPLX, MAX, MIN
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               FACTS / 'F', 'N' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Zomplex precision'
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
      CALL ZERRVX( PATH, NOUT )
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
      NPP = N*( N+1 ) / 2
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
            IF( IUPLO == 1 ) THEN
               UPLO = 'U'
               PACKIT = 'C'
            ELSE
               UPLO = 'L'
               PACKIT = 'R'
            END IF
!
            IF( IMAT /= NTYPES ) THEN
!
!                 Set up parameters with ZLATB4 and generate a test
!                 matrix with ZLATMS.
!
               CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, &
                            MODE, CNDNUM, DIST )
!
               SRNAMT = 'ZLATMS'
               CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, &
                            CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, &
                            WORK, INFO )
!
!                 Check error code from ZLATMS.
!
               IF( INFO /= 0 ) THEN
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, &
                               -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 160
               END IF
!
!                 For types 3-6, zero one or more rows and columns of
!                 the matrix to test that INFO is returned correctly.
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
!                       Set row and column IZERO to zero.
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
                     IF( IUPLO == 1 ) THEN
!
!                          Set the first IZERO rows and columns to zero.
!
                        IOFF = 0
                        DO J = 1, N
                           I2 = MIN( J, IZERO )
                           DO I = 1, I2
                              A( IOFF+I ) = ZERO
                           ENDDO
                           IOFF = IOFF + J
                        ENDDO
                     ELSE
!
!                          Set the last IZERO rows and columns to zero.
!
                        IOFF = 0
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
            ELSE
!
!                 Use a special block diagonal matrix to test alternate
!                 code for the 2-by-2 blocks.
!
               CALL ZLATSP( UPLO, N, A, ISEED )
            END IF
!
            DO IFACT = 1, NFACT
!
!                 Do first for FACT = 'F', then for other values.
!
               FACT = FACTS( IFACT )
!
!                 Compute the condition number for comparison with
!                 the value returned by ZSPSVX.
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
                  ANORM = ZLANSP( '1', UPLO, N, A, RWORK )
!
!                    Factor the matrix A.
!
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
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZSPTRF( UPLO, N, AFAC, IWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZSPTRF : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compute inv(A) and take its norm.
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
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZSPTRI( UPLO, N, AINV, IWORK, WORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZSPTRI : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  AINVNM = ZLANSP( '1', UPLO, N, AINV, RWORK )
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
               SRNAMT = 'ZLARHS'
               CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, &
                            NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, &
                            INFO )
               XTYPE = 'C'
!
!                 --- Test ZSPSV  ---
!
               IF( IFACT == 2 ) THEN
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
!                    Factor the matrix and solve the system using ZSPSV.
!
                  SRNAMT = 'ZSPSV '
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZSPSV( UPLO, N, NRHS, AFAC, IWORK, X, LDA, &
                              INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZSPSV : ',&
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
!                    Check error code from ZSPSV .
!
                  IF( INFO /= K ) THEN
                     CALL ALAERH( PATH, 'ZSPSV ', INFO, K, UPLO, N, &
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
                  CALL ZSPT01( UPLO, N, A, AFAC, IWORK, AINV, LDA, &
                               RWORK, RESULT( 1 ) )
!
!                    Compute residual of the computed solution.
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
                  CALL ZSPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, &
                               RWORK, RESULT( 2 ) )
!
!                    Check solution from generated exact solution.
!
                  CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
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
                        WRITE( NOUT, FMT = 9999 )'ZSPSV ', UPLO, N, &
                           IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                     ENDDO
                  NRUN = NRUN + NT
  120                CONTINUE
               END IF
!
!                 --- Test ZSPSVX ---
!
               IF( IFACT == 2 .AND. NPP > 0 )  THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZLASET( 'Full', NPP, 1, DCMPLX( ZERO ), &
                               DCMPLX( ZERO ), AFAC, NPP )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               ENDIF
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZLASET( 'Full', N, NRHS, DCMPLX( ZERO ), &
                            DCMPLX( ZERO ), X, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Solve the system and compute the condition number and
!                 error bounds using ZSPSVX.
!
               SRNAMT = 'ZSPSVX'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZSPSVX( FACT, UPLO, N, NRHS, A, AFAC, IWORK, B, &
                            LDA, X, LDA, RCOND, RWORK, &
                            RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), &
                            INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZSPSVX : ',&
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
!                 Check the error code from ZSPSVX.
!
               IF( INFO /= K ) THEN
                  CALL ALAERH( PATH, 'ZSPSVX', INFO, K, FACT // UPLO, &
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
                     CALL ZSPT01( UPLO, N, A, AFAC, IWORK, AINV, LDA, &
                                  RWORK( 2*NRHS+1 ), RESULT( 1 ) )
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
                  CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  CALL ZSPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, &
                               RWORK( 2*NRHS+1 ), RESULT( 2 ) )
!
!                    Check solution from generated exact solution.
!
                  CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 3 ) )
!
!                    Check the error bounds from iterative refinement.
!
                  CALL ZPPT05( UPLO, N, NRHS, A, B, LDA, X, LDA, &
                               XACT, LDA, RWORK, RWORK( NRHS+1 ), &
                               RESULT( 4 ) )
               ELSE
                  K1 = 6
               END IF
!
!                 Compare RCOND from ZSPSVX with the computed value
!                 in RCONDC.
!
               RESULT( 6 ) = DGET06( RCOND, RCONDC )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = K1, 6
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALADHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 )'ZSPSVX', FACT, UPLO, &
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
!     End of ZDRVSP
!
END
                                                                                                                                                                                                                                                                                                            




