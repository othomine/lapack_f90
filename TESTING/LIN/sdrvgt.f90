!> \brief \b SDRVGT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVGT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, AF,
!                          B, X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       REAL               A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ),
!      $                   X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SDRVGT tests SGTSV and -SVX.
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
!>          The number of right hand sides, NRHS >= 0.
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
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (NMAX*4)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is REAL array, dimension (NMAX*4)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
!>                      (max(NMAX,2*NRHS))
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
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SDRVGT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, AF, &
                      B, X, XACT, WORK, RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NN, NOUT, NRHS
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NVAL( * )
   REAL               A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ), &
                      X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 12 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 6 )
!     ..
!     .. Local Scalars ..
   LOGICAL            TRFCON, ZEROT
   CHARACTER          DIST, FACT, TRANS, TYPE
   CHARACTER*3        PATH
   INTEGER            I, IFACT, IMAT, IN, INFO, ITRAN, IX, IZERO, J, &
                      K, K1, KL, KOFF, KU, LDA, M, MODE, N, NERRS, &
                      NFAIL, NIMAT, NRUN, NT
   REAL               AINVNM, ANORM, ANORMI, ANORMO, COND, RCOND, &
                      RCONDC, RCONDI, RCONDO
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          TRANSS( 3 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS ), Z( 3 )
!     ..
!     .. External Functions ..
   REAL               SASUM, SGET06, SLANGT
   EXTERNAL           SASUM, SGET06, SLANGT
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, SCOPY, SERRVX, SGET04, &
                      SGTSV, SGTSVX, SGTT01, SGTT02, SGTT05, SGTTRF, &
                      SGTTRS, SLACPY, SLAGTM, SLARNV, SLASET, SLATB4, &
                      SLATMS, SSCAL
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
   PATH( 1: 1 ) = 'Single precision'
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
      CALL SERRVX( PATH, NOUT )
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
            GO TO 130
!
!           Set up parameters with SLATB4.
!
         CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                      COND, DIST )
!
         ZEROT = IMAT >= 8 .AND. IMAT <= 10
         IF( IMAT <= 6 ) THEN
!
!              Types 1-6:  generate matrices of known condition number.
!
            KOFF = MAX( 2-KU, 3-MAX( 1, N ) )
            SRNAMT = 'SLATMS'
            CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, &
                         ANORM, KL, KU, 'Z', AF( KOFF ), 3, WORK, &
                         INFO )
!
!              Check the error code from SLATMS.
!
            IF( INFO /= 0 ) THEN
               CALL ALAERH( PATH, 'SLATMS', INFO, 0, ' ', N, N, KL, &
                            KU, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 130
            END IF
            IZERO = 0
!
            IF( N > 1 ) THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SCOPY( N-1, AF( 4 ), 3, A, 1 )
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
               CALL SCOPY( N-1, AF( 3 ), 3, A( N+M+1 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            END IF
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( N, AF( 2 ), 3, A( M+1 ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
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
               CALL SLARNV( 2, ISEED, N+2*M, A )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLARNV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IF( ANORM /= ONE )  THEN
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SSCAL( N+2*M, ANORM, A, 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SSCAL : ',&
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
         DO IFACT = 1, 2
            IF( IFACT == 1 ) THEN
               FACT = 'F'
            ELSE
               FACT = 'N'
            END IF
!
!              Compute the condition number for comparison with
!              the value returned by SGTSVX.
!
            IF( ZEROT ) THEN
               IF( IFACT == 1 ) &
                  GO TO 120
               RCONDO = ZERO
               RCONDI = ZERO
!
            ELSE IF( IFACT == 1 ) THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SCOPY( N+2*M, A, 1, AF, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Compute the 1-norm and infinity-norm of A.
!
               ANORMO = SLANGT( '1', N, A, A( M+1 ), A( N+M+1 ) )
               ANORMI = SLANGT( 'I', N, A, A( M+1 ), A( N+M+1 ) )
!
!                 Factor the matrix A.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGTTRF( N, AF, AF( M+1 ), AF( N+M+1 ), &
                            AF( N+2*M+1 ), IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGTTRF : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Use SGTTRS to solve for one column at a time of
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
                  CALL SGTTRS( 'No transpose', N, 1, AF, AF( M+1 ), &
                               AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, &
                               LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGTTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  AINVNM = MAX( AINVNM, SASUM( N, X, 1 ) )
               ENDDO
!
!                 Compute the 1-norm condition number of A.
!
               IF( ANORMO <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDO = ONE
               ELSE
                  RCONDO = ( ONE / ANORMO ) / AINVNM
               END IF
!
!                 Use SGTTRS to solve for one column at a time of
!                 inv(A'), computing the maximum column sum as we go.
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
                  CALL SGTTRS( 'Transpose', N, 1, AF, AF( M+1 ), &
                               AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, &
                               LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGTTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  AINVNM = MAX( AINVNM, SASUM( N, X, 1 ) )
               ENDDO
!
!                 Compute the infinity-norm condition number of A.
!
               IF( ANORMI <= ZERO .OR. AINVNM <= ZERO ) THEN
                  RCONDI = ONE
               ELSE
                  RCONDI = ( ONE / ANORMI ) / AINVNM
               END IF
            END IF
!
            DO ITRAN = 1, 3
               TRANS = TRANSS( ITRAN )
               IF( ITRAN == 1 ) THEN
                  RCONDC = RCONDO
               ELSE
                  RCONDC = RCONDI
               END IF
!
!                 Generate NRHS random solution vectors.
!
               IX = 1
               DO J = 1, NRHS
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLARNV( 2, ISEED, N, XACT( IX ) )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLARNV : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IX = IX + LDA
               ENDDO
!
!                 Set the right hand side.
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLAGTM( TRANS, N, NRHS, ONE, A, A( M+1 ), &
                            A( N+M+1 ), XACT, LDA, ZERO, B, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLAGTM : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
               IF( IFACT == 2 .AND. ITRAN == 1 ) THEN
!
!                    --- Test SGTSV  ---
!
!                    Solve the system using Gaussian elimination with
!                    partial pivoting.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SCOPY( N+2*M, A, 1, AF, 1 )
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
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'SGTSV '
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SGTSV( N, NRHS, AF, AF( M+1 ), AF( N+M+1 ), X, &
                              LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGTSV : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from SGTSV .
!
                  IF( INFO /= IZERO ) &
                     CALL ALAERH( PATH, 'SGTSV ', INFO, IZERO, ' ', &
                                  N, N, 1, 1, NRHS, IMAT, NFAIL, &
                                  NERRS, NOUT )
                  NT = 1
                  IF( IZERO == 0 ) THEN
!
!                       Check residual of computed solution.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL SGTT02( TRANS, N, NRHS, A, A( M+1 ), &
                                  A( N+M+1 ), X, LDA, WORK, LDA, &
                                  RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                                  RESULT( 3 ) )
                     NT = 3
                  END IF
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  DO K = 2, NT
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )'SGTSV ', N, IMAT, &
                           K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                  ENDDO
                  NRUN = NRUN + NT - 1
               END IF
!
!                 --- Test SGTSVX ---
!
               IF( IFACT > 1 ) THEN
!
!                    Initialize AF to zero.
!
                  DO I = 1, 3*N - 2
                     AF( I ) = ZERO
                  ENDDO
               END IF
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLASET( 'Full', N, NRHS, ZERO, ZERO, X, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Solve the system and compute the condition number and
!                 error bounds using SGTSVX.
!
               SRNAMT = 'SGTSVX'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGTSVX( FACT, TRANS, N, NRHS, A, A( M+1 ), &
                            A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), &
                            AF( N+2*M+1 ), IWORK, B, LDA, X, LDA, &
                            RCOND, RWORK, RWORK( NRHS+1 ), WORK, &
                            IWORK( N+1 ), INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGTSVX : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check the error code from SGTSVX.
!
               IF( INFO /= IZERO ) &
                  CALL ALAERH( PATH, 'SGTSVX', INFO, IZERO, &
                               FACT // TRANS, N, N, 1, 1, NRHS, IMAT, &
                               NFAIL, NERRS, NOUT )
!
               IF( IFACT >= 2 ) THEN
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                  CALL SGTT01( N, A, A( M+1 ), A( N+M+1 ), AF, &
                               AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), &
                               IWORK, WORK, LDA, RWORK, RESULT( 1 ) )
                  K1 = 1
               ELSE
                  K1 = 2
               END IF
!
               IF( INFO == 0 ) THEN
                  TRFCON = .FALSE.
!
!                    Check residual of computed solution.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  CALL SGTT02( TRANS, N, NRHS, A, A( M+1 ), &
                               A( N+M+1 ), X, LDA, WORK, LDA, &
                               RESULT( 2 ) )
!
!                    Check solution from generated exact solution.
!
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                               RESULT( 3 ) )
!
!                    Check the error bounds from iterative refinement.
!
                  CALL SGTT05( TRANS, N, NRHS, A, A( M+1 ), &
                               A( N+M+1 ), B, LDA, X, LDA, XACT, LDA, &
                               RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                  NT = 5
               END IF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
               DO K = K1, NT
                  IF( RESULT( K ) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) &
                        CALL ALADHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 )'SGTSVX', FACT, TRANS, &
                        N, IMAT, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
                  ENDDO
!
!                 Check the reciprocal of the condition number.
!
               RESULT( 6 ) = SGET06( RCOND, RCONDC )
               IF( RESULT( 6 ) >= THRESH ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) &
                     CALL ALADHD( NOUT, PATH )
                  WRITE( NOUT, FMT = 9998 )'SGTSVX', FACT, TRANS, N, &
                     IMAT, K, RESULT( K )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + NT - K1 + 2
!
               ENDDO
  120       CONTINUE
            ENDDO
  130    CONTINUE
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test ', I2, &
         ', ratio = ', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N =', &
         I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
   RETURN
!
!     End of SDRVGT
!
END
                                                                                                                                                                                                                                                                                                            




