!> \brief \b DCHKPT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          A, D, E, B, X, XACT, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   A( * ), B( * ), D( * ), E( * ), RWORK( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKPT tests DPTTRF, -TRS, -RFS, and -CON
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
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (NMAX*2)
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
   SUBROUTINE DCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      A, D, E, B, X, XACT, WORK, RWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NN, NNS, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            NSVAL( * ), NVAL( * )
   DOUBLE PRECISION   A( * ), B( * ), D( * ), E( * ), RWORK( * ), &
                      WORK( * ), X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 12 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 7 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, TYPE
   CHARACTER*3        PATH
   INTEGER            I, IA, IMAT, IN, INFO, IRHS, IX, IZERO, J, K, &
                      KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, &
                      NRHS, NRUN
   DOUBLE PRECISION   AINVNM, ANORM, COND, DMAX, RCOND, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS ), Z( 3 )
!     ..
!     .. External Functions ..
   INTEGER            IDAMAX
   DOUBLE PRECISION   DASUM, DGET06, DLANST
   EXTERNAL           IDAMAX, DASUM, DGET06, DLANST
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, DCOPY, DERRGT, DGET04, &
                      DLACPY, DLAPTM, DLARNV, DLATB4, DLATMS, DPTCON, &
                      DPTRFS, DPTT01, DPTT02, DPTT05, DPTTRF, DPTTRS, &
                      DSCAL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX
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
   DATA               ISEEDY / 0, 0, 0, 1 /
!     ..
!     .. Executable Statements ..
!
   PATH( 1: 1 ) = 'Double precision'
   PATH( 2: 3 ) = 'PT'
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
      CALL DERRGT( PATH, NOUT )
   INFOT = 0
!
   DO IN = 1, NN
!
!        Do for each value of N in NVAL.
!
      N = NVAL( IN )
      LDA = MAX( 1, N )
      NIMAT = NTYPES
      IF( N <= 0 ) &
         NIMAT = 1
!
      DO IMAT = 1, NIMAT
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
         IF( N > 0 .AND. .NOT.DOTYPE( IMAT ) ) &
            GO TO 100
!
!           Set up parameters with DLATB4.
!
         CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                      COND, DIST )
!
         ZEROT = IMAT >= 8 .AND. IMAT <= 10
         IF( IMAT <= 6 ) THEN
!
!              Type 1-6:  generate a symmetric tridiagonal matrix of
!              known condition number in lower triangular band storage.
!
            SRNAMT = 'DLATMS'
            CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, &
                         ANORM, KL, KU, 'B', A, 2, WORK, INFO )
!
!              Check the error code from DLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'DLATMS', INFO, 0, ' ', N, N, KL, &
                            KU, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 100
            END IF
            IZERO = 0
!
!              Copy the matrix to D and E.
!
            IA = 1
            DO I = 1, N - 1
               D( I ) = A( IA )
               E( I ) = A( IA+1 )
               IA = IA + 2
            ENDDO
            IF( N > 0 ) &
               D( N ) = A( IA )
         ELSE
!
!              Type 7-12:  generate a diagonally dominant matrix with
!              unknown condition number in the vectors D and E.
!
            IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) THEN
!
!                 Let D and E have values from [-1,1].
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLARNV( 2, ISEED, N, D )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLARNV( 2, ISEED, N-1, E )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Make the tridiagonal matrix diagonally dominant.
!
               IF( N == 1 ) THEN
                  D( 1 ) = ABS( D( 1 ) )
               ELSE
                  D( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                  D( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
                  DO I = 2, N - 1
                     D( I ) = ABS( D( I ) ) + ABS( E( I ) ) + &
                              ABS( E( I-1 ) )
                  ENDDO
               END IF
!
!                 Scale D and E so the maximum element is ANORM.
!
               IX = IDAMAX( N, D, 1 )
               DMAX = D( IX )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DSCAL( N, ANORM / DMAX, D, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DSCAL( N-1, ANORM / DMAX, E, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
            ELSE IF( IZERO > 0 ) THEN
!
!                 Reuse the last matrix by copying back the zeroed out
!                 elements.
!
               IF( IZERO == 1 ) THEN
                  D( 1 ) = Z( 2 )
                  IF( N > 1 ) &
                     E( 1 ) = Z( 3 )
               ELSE IF( IZERO == N ) THEN
                  E( N-1 ) = Z( 1 )
                  D( N ) = Z( 2 )
               ELSE
                  E( IZERO-1 ) = Z( 1 )
                  D( IZERO ) = Z( 2 )
                  E( IZERO ) = Z( 3 )
               END IF
            END IF
!
!              For types 8-10, set one row and column of the matrix to
!              zero.
!
            IZERO = 0
            IF( IMAT == 8 ) THEN
               IZERO = 1
               Z( 2 ) = D( 1 )
               D( 1 ) = ZERO
               IF( N > 1 ) THEN
                  Z( 3 ) = E( 1 )
                  E( 1 ) = ZERO
               END IF
            ELSE IF( IMAT == 9 ) THEN
               IZERO = N
               IF( N > 1 ) THEN
                  Z( 1 ) = E( N-1 )
                  E( N-1 ) = ZERO
               END IF
               Z( 2 ) = D( N )
               D( N ) = ZERO
            ELSE IF( IMAT == 10 ) THEN
               IZERO = ( N+1 ) / 2
               IF( IZERO > 1 ) THEN
                  Z( 1 ) = E( IZERO-1 )
                  E( IZERO-1 ) = ZERO
                  Z( 3 ) = E( IZERO )
                  E( IZERO ) = ZERO
               END IF
               Z( 2 ) = D( IZERO )
               D( IZERO ) = ZERO
            END IF
         END IF
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DCOPY( N, D, 1, D( N+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( N > 1 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DCOPY( N-1, E, 1, E( N+1 ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
!+    TEST 1
!           Factor A as L*D*L' and compute the ratio
!              norm(L*D*L' - A) / (n * norm(A) * EPS )
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DPTTRF( N, D( N+1 ), E( N+1 ), INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DPTTRF : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from DPTTRF.
!
         IF( INFO /= IZERO ) THEN
            CALL ALAERH( PATH, 'DPTTRF', INFO, IZERO, ' ', N, N, -1, &
                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
            GO TO 100
         END IF
!
         IF( INFO > 0 ) THEN
            RCONDC = ZERO
            GO TO 90
         END IF
!
         CALL DPTT01( N, D, E, D( N+1 ), E( N+1 ), WORK, &
                      RESULT( 1 ) )
!
!           Print the test ratio if greater than or equal to THRESH.
!
         IF( RESULT( 1 ) >= THRESH ) THEN
            IF( NFAIL == 0 .AND. NERRS == 0 ) &
               CALL ALAHD( NOUT, PATH )
            WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 )
            NFAIL = NFAIL + 1
         END IF
         NRUN = NRUN + 1
!
!           Compute RCONDC = 1 / (norm(A) * norm(inv(A))
!
!           Compute norm(A).
!
         ANORM = DLANST( '1', N, D, E )
!
!           Use DPTTRS to solve for one column at a time of inv(A),
!           computing the maximum column sum as we go.
!
         AINVNM = ZERO
         DO I = 1, N
            DO J = 1, N
               X( J ) = ZERO
            ENDDO
            X( I ) = ONE
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DPTTRS( N, 1, D( N+1 ), E( N+1 ), X, LDA, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DPTTRS : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            AINVNM = MAX( AINVNM, DASUM( N, X, 1 ) )
         ENDDO
         RCONDC = ONE / MAX( ONE, ANORM*AINVNM )
!
         DO IRHS = 1, NNS
            NRHS = NSVAL( IRHS )
!
!           Generate NRHS random solution vectors.
!
            IX = 1
            DO J = 1, NRHS
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLARNV( 2, ISEED, N, XACT( IX ) )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IX = IX + LDA
            ENDDO
!
!           Set the right hand side.
!
            CALL DLAPTM( N, NRHS, ONE, D, E, XACT, LDA, ZERO, B, &
                         LDA )
!
!+    TEST 2
!           Solve A*x = b and compute the residual.
!
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
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DPTTRS( N, NRHS, D( N+1 ), E( N+1 ), X, LDA, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DPTTRS : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!           Check error code from DPTTRS.
!
            IF( INFO /= 0 ) &
               CALL ALAERH( PATH, 'DPTTRS', INFO, 0, ' ', N, N, -1, &
                            -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
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
            CALL DPTT02( N, NRHS, D, E, X, LDA, WORK, LDA, &
                         RESULT( 2 ) )
!
!+    TEST 3
!           Check solution from generated exact solution.
!
            CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                         RESULT( 3 ) )
!
!+    TESTS 4, 5, and 6
!           Use iterative refinement to improve the solution.
!
            SRNAMT = 'DPTRFS'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DPTRFS( N, NRHS, D, E, D( N+1 ), E( N+1 ), B, LDA, &
                         X, LDA, RWORK, RWORK( NRHS+1 ), WORK, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DPTRFS : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!           Check error code from DPTRFS.
!
            IF( INFO /= 0 ) &
               CALL ALAERH( PATH, 'DPTRFS', INFO, 0, ' ', N, N, -1, &
                            -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
!
            CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                         RESULT( 4 ) )
            CALL DPTT05( N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA, &
                         RWORK, RWORK( NRHS+1 ), RESULT( 5 ) )
!
!           Print information about the tests that did not pass the
!           threshold.
!
            DO K = 2, 6
               IF( RESULT( K ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALAHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9998 )N, NRHS, IMAT, K, &
                     RESULT( K )
                  NFAIL = NFAIL + 1
               END IF
            ENDDO
            NRUN = NRUN + 5
         ENDDO
!
!+    TEST 7
!           Estimate the reciprocal of the condition number of the
!           matrix.
!
90       CONTINUE
         SRNAMT = 'DPTCON'
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DPTCON( N, D( N+1 ), E( N+1 ), ANORM, RCOND, RWORK, &
                      INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DPTCON : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from DPTCON.
!
         IF( INFO /= 0 ) &
            CALL ALAERH( PATH, 'DPTCON', INFO, 0, ' ', N, N, -1, -1, &
                         -1, IMAT, NFAIL, NERRS, NOUT )
!
         RESULT( 7 ) = DGET06( RCOND, RCONDC )
!
!           Print the test ratio if greater than or equal to THRESH.
!
         IF( RESULT( 7 ) >= THRESH ) THEN
            IF( NFAIL == 0 .AND. NERRS == 0 ) &
               CALL ALAHD( NOUT, PATH )
            WRITE( NOUT, FMT = 9999 )N, IMAT, 7, RESULT( 7 )
            NFAIL = NFAIL + 1
         END IF
         NRUN = NRUN + 1
  100    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( ' N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', &
         G12.5 )
 9998 FORMAT( ' N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, &
         ') = ', G12.5 )
   RETURN
!
!     End of DCHKPT
!
END
                                                                                                                                                                                                                                                                                                            




