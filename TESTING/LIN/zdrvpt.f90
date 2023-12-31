!> \brief \b ZDRVPT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVPT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, D,
!                          E, B, X, XACT, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NVAL( * )
!       DOUBLE PRECISION   D( * ), RWORK( * )
!       COMPLEX*16         A( * ), B( * ), E( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDRVPT tests ZPTSV and -SVX.
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
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (NMAX*2)
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
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (NMAX+2*NRHS)
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
   SUBROUTINE ZDRVPT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, D, &
                      E, B, X, XACT, WORK, RWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NN, NOUT, NRHS
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            NVAL( * )
   DOUBLE PRECISION   D( * ), RWORK( * )
   COMPLEX*16         A( * ), B( * ), E( * ), WORK( * ), X( * ), &
                      XACT( * )
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
   PARAMETER          ( NTESTS = 6 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   CHARACTER          DIST, FACT, TYPE
   CHARACTER*3        PATH
   INTEGER            I, IA, IFACT, IMAT, IN, INFO, IX, IZERO, J, K, &
                      K1, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, &
                      NRUN, NT
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
   DOUBLE PRECISION   DGET06, DZASUM, ZLANHT
   EXTERNAL           IDAMAX, DGET06, DZASUM, ZLANHT
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, DCOPY, DLARNV, DSCAL, &
                      ZCOPY, ZDSCAL, ZERRVX, ZGET04, ZLACPY, ZLAPTM, &
                      ZLARNV, ZLASET, ZLATB4, ZLATMS, ZPTSV, ZPTSVX, &
                      ZPTT01, ZPTT02, ZPTT05, ZPTTRF, ZPTTRS
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DCMPLX, MAX
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
   PATH( 1: 1 ) = 'Zomplex precision'
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
      CALL ZERRVX( PATH, NOUT )
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
            GO TO 110
!
!           Set up parameters with ZLATB4.
!
         CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                      COND, DIST )
!
         ZEROT = IMAT >= 8 .AND. IMAT <= 10
         IF( IMAT <= 6 ) THEN
!
!              Type 1-6:  generate a symmetric tridiagonal matrix of
!              known condition number in lower triangular band storage.
!
            SRNAMT = 'ZLATMS'
            CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, &
                         ANORM, KL, KU, 'B', A, 2, WORK, INFO )
!
!              Check the error code from ZLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', N, N, KL, &
                            KU, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 110
            END IF
            IZERO = 0
!
!              Copy the matrix to D and E.
!
            IA = 1
            DO I = 1, N - 1
               D( I ) = DBLE( A( IA ) )
               E( I ) = A( IA+1 )
               IA = IA + 2
            ENDDO
            IF( N > 0 ) &
               D( N ) = DBLE( A( IA ) )
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
               CALL ZLARNV( 2, ISEED, N-1, E )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
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
               IF( N > 1 )  THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL ZDSCAL( N-1, ANORM / DMAX, E, 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZDSCAL : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               ENDIF
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
                  Z( 3 ) = DBLE( E( 1 ) )
                  E( 1 ) = ZERO
               END IF
            ELSE IF( IMAT == 9 ) THEN
               IZERO = N
               IF( N > 1 ) THEN
                  Z( 1 ) = DBLE( E( N-1 ) )
                  E( N-1 ) = ZERO
               END IF
               Z( 2 ) = D( N )
               D( N ) = ZERO
            ELSE IF( IMAT == 10 ) THEN
               IZERO = ( N+1 ) / 2
               IF( IZERO > 1 ) THEN
                  Z( 1 ) = DBLE( E( IZERO-1 ) )
                  E( IZERO-1 ) = ZERO
                  Z( 3 ) = DBLE( E( IZERO ) )
                  E( IZERO ) = ZERO
               END IF
               Z( 2 ) = D( IZERO )
               D( IZERO ) = ZERO
            END IF
         END IF
!
!           Generate NRHS random solution vectors.
!
         IX = 1
         DO J = 1, NRHS
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLARNV( 2, ISEED, N, XACT( IX ) )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IX = IX + LDA
         ENDDO
!
!           Set the right hand side.
!
         CALL ZLAPTM( 'Lower', N, NRHS, ONE, D, E, XACT, LDA, ZERO, &
                      B, LDA )
!
         DO IFACT = 1, 2
            IF( IFACT == 1 ) THEN
               FACT = 'F'
            ELSE
               FACT = 'N'
            END IF
!
!              Compute the condition number for comparison with
!              the value returned by ZPTSVX.
!
            IF( ZEROT ) THEN
               IF( IFACT == 1 ) &
                  GO TO 100
               RCONDC = ZERO
!
            ELSE IF( IFACT == 1 ) THEN
!
!                 Compute the 1-norm of A.
!
               ANORM = ZLANHT( '1', N, D, E )
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
                  CALL ZCOPY( N-1, E, 1, E( N+1 ), 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               ENDIF
!
!                 Factor the matrix A.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZPTTRF( N, D( N+1 ), E( N+1 ), INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZPTTRF : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Use ZPTTRS to solve for one column at a time of
!                 inv(A), computing the maximum column sum as we go.
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
                  CALL ZPTTRS( 'Lower', N, 1, D( N+1 ), E( N+1 ), X, &
                               LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZPTTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  AINVNM = MAX( AINVNM, DZASUM( N, X, 1 ) )
               ENDDO
!
!                 Compute the 1-norm condition number of A.
!
               IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDC = ONE
               ELSE
                  RCONDC = ( ONE / ANORM ) / AINVNM
               END IF
            END IF
!
            IF( IFACT == 2 ) THEN
!
!                 --- Test ZPTSV --
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
                  CALL ZCOPY( N-1, E, 1, E( N+1 ), 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               ENDIF
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
!                 Factor A as L*D*L' and solve the system A*X = B.
!
               SRNAMT = 'ZPTSV '
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL ZPTSV( N, NRHS, D( N+1 ), E( N+1 ), X, LDA, &
                           INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : ZPTSV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from ZPTSV .
!
               IF( INFO /= IZERO ) &
                  CALL ALAERH( PATH, 'ZPTSV ', INFO, IZERO, ' ', N, &
                               N, 1, 1, NRHS, IMAT, NFAIL, NERRS, &
                               NOUT )
               NT = 0
               IF( IZERO == 0 ) THEN
!
!                    Check the factorization by computing the ratio
!                       norm(L*D*L' - A) / (n * norm(A) * EPS )
!
                  CALL ZPTT01( N, D, E, D( N+1 ), E( N+1 ), WORK, &
                               RESULT( 1 ) )
!
!                    Compute the residual in the solution.
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
                  CALL ZPTT02( 'Lower', N, NRHS, D, E, X, LDA, WORK, &
                               LDA, RESULT( 2 ) )
!
!                    Check solution from generated exact solution.
!
                  CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 3 ) )
                  NT = 3
               END IF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 1, NT
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALADHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9999 )'ZPTSV ', N, IMAT, K, &
                        RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
               ENDDO
               NRUN = NRUN + NT
            END IF
!
!              --- Test ZPTSVX ---
!
            IF( IFACT > 1 ) THEN
!
!                 Initialize D( N+1:2*N ) and E( N+1:2*N ) to zero.
!
               DO I = 1, N - 1
                  D( N+I ) = ZERO
                  E( N+I ) = ZERO
               ENDDO
               IF( N > 0 ) &
                  D( N+N ) = ZERO
            END IF
!
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
!              Solve the system and compute the condition number and
!              error bounds using ZPTSVX.
!
            SRNAMT = 'ZPTSVX'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZPTSVX( FACT, N, NRHS, D, E, D( N+1 ), E( N+1 ), B, &
                         LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), &
                         WORK, RWORK( 2*NRHS+1 ), INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZPTSVX : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check the error code from ZPTSVX.
!
            IF( INFO /= IZERO ) &
               CALL ALAERH( PATH, 'ZPTSVX', INFO, IZERO, FACT, N, N, &
                            1, 1, NRHS, IMAT, NFAIL, NERRS, NOUT )
            IF( IZERO == 0 ) THEN
               IF( IFACT == 2 ) THEN
!
!                    Check the factorization by computing the ratio
!                       norm(L*D*L' - A) / (n * norm(A) * EPS )
!
                  K1 = 1
                  CALL ZPTT01( N, D, E, D( N+1 ), E( N+1 ), WORK, &
                               RESULT( 1 ) )
               ELSE
                  K1 = 2
               END IF
!
!                 Compute the residual in the solution.
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
               CALL ZPTT02( 'Lower', N, NRHS, D, E, X, LDA, WORK, &
                            LDA, RESULT( 2 ) )
!
!                 Check solution from generated exact solution.
!
               CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 3 ) )
!
!                 Check error bounds from iterative refinement.
!
               CALL ZPTT05( N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA, &
                            RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
            ELSE
               K1 = 6
            END IF
!
!              Check the reciprocal of the condition number.
!
            RESULT( 6 ) = DGET06( RCOND, RCONDC )
!
!              Print information about the tests that did not pass
!              the threshold.
!
            DO K = K1, 6
               IF( RESULT( K ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALADHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9998 )'ZPTSVX', FACT, N, IMAT, &
                     K, RESULT( K )
                  NFAIL = NFAIL + 1
               END IF
            ENDDO
            NRUN = NRUN + 7 - K1
  100       CONTINUE
            ENDDO
  110    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test ', I2, &
         ', ratio = ', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', N =', I5, ', type ', I2, &
         ', test ', I2, ', ratio = ', G12.5 )
   RETURN
!
!     End of ZDRVPT
!
END
                                                                                                                                                                                                                                                                                                            




