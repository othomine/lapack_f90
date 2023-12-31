!> \brief \b DCHKGT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKGT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          A, AF, B, X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ),
!      $                   X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKGT tests DGTTRF, -TRS, -RFS, and -CON
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
!>          A is DOUBLE PRECISION array, dimension (NMAX*4)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (NMAX*4)
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
   SUBROUTINE DCHKGT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, &
                      A, AF, B, X, XACT, WORK, RWORK, IWORK, NOUT )
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
   INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
   DOUBLE PRECISION   A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ), &
                      X( * ), XACT( * )
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
   LOGICAL            TRFCON, ZEROT
   CHARACTER          DIST, NORM, TRANS, TYPE
   CHARACTER*3        PATH
   INTEGER            I, IMAT, IN, INFO, IRHS, ITRAN, IX, IZERO, J, &
                      K, KL, KOFF, KU, LDA, M, MODE, N, NERRS, NFAIL, &
                      NIMAT, NRHS, NRUN
   DOUBLE PRECISION   AINVNM, ANORM, COND, RCOND, RCONDC, RCONDI, &
                      RCONDO
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          TRANSS( 3 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS ), Z( 3 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DASUM, DGET06, DLANGT
   EXTERNAL           DASUM, DGET06, DLANGT
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAERH, ALAHD, ALASUM, DCOPY, DERRGE, DGET04, &
                      DGTCON, DGTRFS, DGTT01, DGTT02, DGTT05, DGTTRF, &
                      DGTTRS, DLACPY, DLAGTM, DLARNV, DLATB4, DLATMS, &
                      DSCAL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
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
   DATA               ISEEDY / 0, 0, 0, 1 / , TRANSS / 'N', 'T', &
                      'C' /
!     ..
!     .. Executable Statements ..
!
   PATH( 1: 1 ) = 'Double precision'
   PATH( 2: 3 ) = 'GT'
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
      CALL DERRGE( PATH, NOUT )
   INFOT = 0
!
   DO IN = 1, NN
!
!        Do for each value of N in NVAL.
!
      N = NVAL( IN )
      M = MAX( N-1, 0 )
      LDA = MAX( 1, N )
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
!           Set up parameters with DLATB4.
!
         CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                      COND, DIST )
!
         ZEROT = IMAT >= 8 .AND. IMAT <= 10
         IF( IMAT <= 6 ) THEN
!
!              Types 1-6:  generate matrices of known condition number.
!
            KOFF = MAX( 2-KU, 3-MAX( 1, N ) )
            SRNAMT = 'DLATMS'
            CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, &
                         ANORM, KL, KU, 'Z', AF( KOFF ), 3, WORK, &
                         INFO )
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
            IF( N > 1 ) THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( N-1, AF( 4 ), 3, A, 1 )
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
               CALL DCOPY( N-1, AF( 3 ), 3, A( N+M+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            END IF
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DCOPY( N, AF( 2 ), 3, A( M+1 ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ELSE
!
!              Types 7-12:  generate tridiagonal matrices with
!              unknown condition numbers.
!
            IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) THEN
!
!                 Generate a matrix with elements from [-1,1].
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLARNV( 2, ISEED, N+2*M, A )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( ANORM /= ONE )  THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DSCAL( N+2*M, ANORM, A, 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
               ENDIF
            ELSE IF( IZERO > 0 ) THEN
!
!                 Reuse the last matrix by copying back the zeroed out
!                 elements.
!
               IF( IZERO == 1 ) THEN
                  A( N ) = Z( 2 )
                  IF( N > 1 ) &
                     A( 1 ) = Z( 3 )
               ELSE IF( IZERO == N ) THEN
                  A( 3*N-2 ) = Z( 1 )
                  A( 2*N-1 ) = Z( 2 )
               ELSE
                  A( 2*N-2+IZERO ) = Z( 1 )
                  A( N-1+IZERO ) = Z( 2 )
                  A( IZERO ) = Z( 3 )
               END IF
            END IF
!
!              If IMAT > 7, set one column of the matrix to 0.
!
            IF( .NOT.ZEROT ) THEN
               IZERO = 0
            ELSE IF( IMAT == 8 ) THEN
               IZERO = 1
               Z( 2 ) = A( N )
               A( N ) = ZERO
               IF( N > 1 ) THEN
                  Z( 3 ) = A( 1 )
                  A( 1 ) = ZERO
               END IF
            ELSE IF( IMAT == 9 ) THEN
               IZERO = N
               Z( 1 ) = A( 3*N-2 )
               Z( 2 ) = A( 2*N-1 )
               A( 3*N-2 ) = ZERO
               A( 2*N-1 ) = ZERO
            ELSE
               IZERO = ( N+1 ) / 2
               DO I = IZERO, N - 1
                  A( 2*N-2+I ) = ZERO
                  A( N-1+I ) = ZERO
                  A( I ) = ZERO
               ENDDO
               A( 3*N-2 ) = ZERO
               A( 2*N-1 ) = ZERO
            END IF
         END IF
!
!+    TEST 1
!           Factor A as L*U and compute the ratio
!              norm(L*U - A) / (n * norm(A) * EPS )
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DCOPY( N+2*M, A, 1, AF, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         SRNAMT = 'DGTTRF'
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGTTRF( N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), &
                      IWORK, INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGTTRF : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Check error code from DGTTRF.
!
         IF( INFO /= IZERO ) &
            CALL ALAERH( PATH, 'DGTTRF', INFO, IZERO, ' ', N, N, 1, &
                         1, -1, IMAT, NFAIL, NERRS, NOUT )
         TRFCON = INFO /= 0
!
         CALL DGTT01( N, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), &
                      AF( N+M+1 ), AF( N+2*M+1 ), IWORK, WORK, LDA, &
                      RWORK, RESULT( 1 ) )
!
!           Print the test ratio if it is  >=  THRESH.
!
         IF( RESULT( 1 ) >= THRESH ) THEN
            IF( NFAIL == 0 .AND. NERRS == 0 ) &
               CALL ALAHD( NOUT, PATH )
            WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 )
            NFAIL = NFAIL + 1
         END IF
         NRUN = NRUN + 1
!
         DO ITRAN = 1, 2
            TRANS = TRANSS( ITRAN )
            IF( ITRAN == 1 ) THEN
               NORM = 'O'
            ELSE
               NORM = 'I'
            END IF
            ANORM = DLANGT( NORM, N, A, A( M+1 ), A( N+M+1 ) )
!
            IF( .NOT.TRFCON ) THEN
!
!                 Use DGTTRS to solve for one column at a time of inv(A)
!                 or inv(A^T), computing the maximum column sum as we
!                 go.
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
                  CALL DGTTRS( TRANS, N, 1, AF, AF( M+1 ), &
                               AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, &
                               LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DGTTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  AINVNM = MAX( AINVNM, DASUM( N, X, 1 ) )
               ENDDO
!
!                 Compute RCONDC = 1 / (norm(A) * norm(inv(A))
!
               IF( ANORM <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDC = ONE
               ELSE
                  RCONDC = ( ONE / ANORM ) / AINVNM
               END IF
               IF( ITRAN == 1 ) THEN
                  RCONDO = RCONDC
               ELSE
                  RCONDI = RCONDC
               END IF
            ELSE
               RCONDC = ZERO
            END IF
!
!+    TEST 7
!              Estimate the reciprocal of the condition number of the
!              matrix.
!
            SRNAMT = 'DGTCON'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGTCON( NORM, N, AF, AF( M+1 ), AF( N+M+1 ), &
                         AF( N+2*M+1 ), IWORK, ANORM, RCOND, WORK, &
                         IWORK( N+1 ), INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGTCON : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
!              Check error code from DGTCON.
!
            IF( INFO /= 0 ) &
               CALL ALAERH( PATH, 'DGTCON', INFO, 0, NORM, N, N, -1, &
                            -1, -1, IMAT, NFAIL, NERRS, NOUT )
!
            RESULT( 7 ) = DGET06( RCOND, RCONDC )
!
!              Print the test ratio if it is  >=  THRESH.
!
            IF( RESULT( 7 ) >= THRESH ) THEN
               IF( NFAIL == 0 .AND. NERRS == 0 ) &
                  CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 7, &
                  RESULT( 7 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
         ENDDO
!
!           Skip the remaining tests if the matrix is singular.
!
         IF( TRFCON ) &
            GO TO 100
!
         DO IRHS = 1, NNS
            NRHS = NSVAL( IRHS )
!
!              Generate NRHS random solution vectors.
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
            DO ITRAN = 1, 3
               TRANS = TRANSS( ITRAN )
               IF( ITRAN == 1 ) THEN
                  RCONDC = RCONDO
               ELSE
                  RCONDC = RCONDI
               END IF
!
!                 Set the right hand side.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLAGTM( TRANS, N, NRHS, ONE, A, A( M+1 ), &
                            A( N+M+1 ), XACT, LDA, ZERO, B, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLAGTM : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!+    TEST 2
!                 Solve op(A) * X = B and compute the residual.
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
               SRNAMT = 'DGTTRS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DGTTRS( TRANS, N, NRHS, AF, AF( M+1 ), &
                            AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, &
                            LDA, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DGTTRS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from DGTTRS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'DGTTRS', INFO, 0, TRANS, N, N, &
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
               CALL DGTT02( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), &
                            X, LDA, WORK, LDA, RESULT( 2 ) )
!
!+    TEST 3
!                 Check solution from generated exact solution.
!
               CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 3 ) )
!
!+    TESTS 4, 5, and 6
!                 Use iterative refinement to improve the solution.
!
               SRNAMT = 'DGTRFS'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DGTRFS( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), &
                            AF, AF( M+1 ), AF( N+M+1 ), &
                            AF( N+2*M+1 ), IWORK, B, LDA, X, LDA, &
                            RWORK, RWORK( NRHS+1 ), WORK, &
                            IWORK( N+1 ), INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DGTRFS : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from DGTRFS.
!
               IF( INFO /= 0 ) &
                  CALL ALAERH( PATH, 'DGTRFS', INFO, 0, TRANS, N, N, &
                               -1, -1, NRHS, IMAT, NFAIL, NERRS, &
                               NOUT )
!
               CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 4 ) )
               CALL DGTT05( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), &
                            B, LDA, X, LDA, XACT, LDA, RWORK, &
                            RWORK( NRHS+1 ), RESULT( 5 ) )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = 2, 6
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, &
                        K, RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
               ENDDO
               NRUN = NRUN + 5
            ENDDO
         ENDDO
!
  100    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 12X, 'N =', I5, ',', 10X, ' type ', I2, ', test(', I2, &
         ') = ', G12.5 )
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', &
         I2, ', test(', I2, ') = ', G12.5 )
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2, &
         ', test(', I2, ') = ', G12.5 )
   RETURN
!
!     End of DCHKGT
!
END
                                                                                                                                                                                                                                                                                                            




