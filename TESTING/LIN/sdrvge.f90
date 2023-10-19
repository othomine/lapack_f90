!> \brief \b SDRVGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       REAL               A( * ), AFAC( * ), ASAV( * ), B( * ),
!      $                   BSAV( * ), RWORK( * ), S( * ), WORK( * ),
!      $                   X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SDRVGE tests the driver routines SGESV and -SVX.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand side vectors to be generated for
!>          each linear system.
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
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is REAL array, dimension (NMAX*NRHS)
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
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (2*NMAX)
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
!>          RWORK is REAL array, dimension (2*NRHS+NMAX)
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
   SUBROUTINE SDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, &
                      A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, &
                      RWORK, IWORK, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            TSTERR
   INTEGER            NMAX, NN, NOUT, NRHS
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            IWORK( * ), NVAL( * )
   REAL               A( * ), AFAC( * ), ASAV( * ), B( * ), &
                      BSAV( * ), RWORK( * ), S( * ), WORK( * ), &
                      X( * ), XACT( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTYPES
   PARAMETER          ( NTYPES = 11 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 7 )
   INTEGER            NTRAN
   PARAMETER          ( NTRAN = 3 )
!     ..
!     .. Local Scalars ..
   LOGICAL            EQUIL, NOFACT, PREFAC, TRFCON, ZEROT
   CHARACTER          DIST, EQUED, FACT, TRANS, TYPE, XTYPE
   CHARACTER*3        PATH
   INTEGER            I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, ITRAN, &
                      IZERO, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, &
                      NBMIN, NERRS, NFACT, NFAIL, NIMAT, NRUN, NT
   REAL               AINVNM, AMAX, ANORM, ANORMI, ANORMO, CNDNUM, &
                      COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, &
                      ROLDI, ROLDO, ROWCND, RPVGRW
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SGET06, SLAMCH, SLANGE, SLANTR
   EXTERNAL           LSAME, SGET06, SLAMCH, SLANGE, SLANTR
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, SERRVX, SGEEQU, SGESV, &
                      SGESVX, SGET01, SGET02, SGET04, SGET07, SGETRF, &
                      SGETRI, SLACPY, SLAQGE, SLARHS, SLASET, SLATB4, &
                      SLATMS, XLAENV
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
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               TRANSS / 'N', 'T', 'C' /
   DATA               FACTS / 'F', 'N', 'E' /
   DATA               EQUEDS / 'N', 'R', 'C', 'B' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   PATH( 1: 1 ) = 'Single precision'
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
   IF( TSTERR ) &
      CALL SERRVX( PATH, NOUT )
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
            GO TO 80
!
!           Skip types 5, 6, or 7 if the matrix size is too small.
!
         ZEROT = IMAT >= 5 .AND. IMAT <= 7
         IF( ZEROT .AND. N < IMAT-4 ) &
            GO TO 80
!
!           Set up parameters with SLATB4 and generate a test matrix
!           with SLATMS.
!
         CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, &
                      CNDNUM, DIST )
         RCONDC = ONE / CNDNUM
!
         SRNAMT = 'SLATMS'
         CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, &
                      ANORM, KL, KU, 'No packing', A, LDA, WORK, &
                      INFO )
!
!           Check error code from SLATMS.
!
         IF( INFO /= 0 ) THEN
            CALL ALAERH( PATH, 'SLATMS', INFO, 0, ' ', N, N, -1, -1, &
                         -1, IMAT, NFAIL, NERRS, NOUT )
            GO TO 80
         END IF
!
!           For types 5-7, zero one or more columns of the matrix to
!           test that INFO is returned correctly.
!
         IF( ZEROT ) THEN
            IF( IMAT == 5 ) THEN
               IZERO = 1
            ELSE IF( IMAT == 6 ) THEN
               IZERO = N
            ELSE
               IZERO = N / 2 + 1
            END IF
            IOFF = ( IZERO-1 )*LDA
            IF( IMAT < 7 ) THEN
               DO I = 1, N
                  A( IOFF+I ) = ZERO
               ENDDO
            ELSE
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLASET( 'Full', N, N-IZERO+1, ZERO, ZERO, &
                            A( IOFF+1 ), LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            END IF
         ELSE
            IZERO = 0
         END IF
!
!           Save a copy of the matrix A in ASAV.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLACPY( 'Full', N, N, A, LDA, ASAV, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         DO IEQUED = 1, 4
            EQUED = EQUEDS( IEQUED )
            IF( IEQUED == 1 ) THEN
               NFACT = 3
            ELSE
               NFACT = 1
            END IF
!
            DO IFACT = 1, NFACT
               FACT = FACTS( IFACT )
               PREFAC = LSAME( FACT, 'F' )
               NOFACT = LSAME( FACT, 'N' )
               EQUIL = LSAME( FACT, 'E' )
!
               IF( ZEROT ) THEN
                  IF( PREFAC ) &
                     GO TO 60
                  RCONDO = ZERO
                  RCONDI = ZERO
!
               ELSE IF( .NOT.NOFACT ) THEN
!
!                    Compute the condition number for comparison with
!                    the value returned by SGESVX (FACT = 'N' reuses
!                    the condition number from the previous iteration
!                    with FACT = 'F').
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, N, ASAV, LDA, AFAC, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  IF( EQUIL .OR. IEQUED > 1 ) THEN
!
!                       Compute row and column scale factors to
!                       equilibrate the matrix A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SGEEQU( N, N, AFAC, LDA, S, S( N+1 ), &
                                  ROWCND, COLCND, AMAX, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SGEEQU : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     IF( INFO == 0 .AND. N > 0 ) THEN
                        IF( LSAME( EQUED, 'R' ) ) THEN
                           ROWCND = ZERO
                           COLCND = ONE
                        ELSE IF( LSAME( EQUED, 'C' ) ) THEN
                           ROWCND = ONE
                           COLCND = ZERO
                        ELSE IF( LSAME( EQUED, 'B' ) ) THEN
                           ROWCND = ZERO
                           COLCND = ZERO
                        END IF
!
!                          Equilibrate the matrix.
!
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL SLAQGE( N, N, AFAC, LDA, S, S( N+1 ), &
                                     ROWCND, COLCND, AMAX, EQUED )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : SLAQGE : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                     END IF
                  END IF
!
!                    Save the condition number of the non-equilibrated
!                    system for use in SGET04.
!
                  IF( EQUIL ) THEN
                     ROLDO = RCONDO
                     ROLDI = RCONDI
                  END IF
!
!                    Compute the 1-norm and infinity-norm of A.
!
                  ANORMO = SLANGE( '1', N, N, AFAC, LDA, RWORK )
                  ANORMI = SLANGE( 'I', N, N, AFAC, LDA, RWORK )
!
!                    Factor the matrix A.
!
                  SRNAMT = 'SGETRF'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SGETRF( N, N, AFAC, LDA, IWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGETRF : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Form the inverse of A.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, N, AFAC, LDA, A, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  LWORK = NMAX*MAX( 3, NRHS )
                  SRNAMT = 'SGETRI'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SGETRI( N, A, LDA, IWORK, WORK, LWORK, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGETRI : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compute the 1-norm condition number of A.
!
                  AINVNM = SLANGE( '1', N, N, A, LDA, RWORK )
                  IF( ANORMO <= ZERO .OR. AINVNM <= ZERO ) THEN
                     RCONDO = ONE
                  ELSE
                     RCONDO = ( ONE / ANORMO ) / AINVNM
                  END IF
!
!                    Compute the infinity-norm condition number of A.
!
                  AINVNM = SLANGE( 'I', N, N, A, LDA, RWORK )
                  IF( ANORMI <= ZERO .OR. AINVNM <= ZERO ) THEN
                     RCONDI = ONE
                  ELSE
                     RCONDI = ( ONE / ANORMI ) / AINVNM
                  END IF
               END IF
!
               DO ITRAN = 1, NTRAN
!
!                    Do for each value of TRANS.
!
                  TRANS = TRANSS( ITRAN )
                  IF( ITRAN == 1 ) THEN
                     RCONDC = RCONDO
                  ELSE
                     RCONDC = RCONDI
                  END IF
!
!                    Restore the matrix A.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Form an exact solution and set the right hand side.
!
                  SRNAMT = 'SLARHS'
                  CALL SLARHS( PATH, XTYPE, 'Full', TRANS, N, N, KL, &
                               KU, NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
                  XTYPE = 'C'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  IF( NOFACT .AND. ITRAN == 1 ) THEN
!
!                       --- Test SGESV  ---
!
!                       Compute the LU factorization of the matrix and
!                       solve the system.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'Full', N, N, A, LDA, AFAC, LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
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
                     SRNAMT = 'SGESV '
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SGESV( N, NRHS, AFAC, LDA, IWORK, X, LDA, &
                                 INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SGESV : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Check error code from SGESV .
!
                     IF( INFO /= IZERO ) &
                        CALL ALAERH( PATH, 'SGESV ', INFO, IZERO, &
                                     ' ', N, N, -1, -1, NRHS, IMAT, &
                                     NFAIL, NERRS, NOUT )
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                     CALL SGET01( N, N, A, LDA, AFAC, LDA, IWORK, &
                                  RWORK, RESULT( 1 ) )
                     NT = 1
                     IF( IZERO == 0 ) THEN
!
!                          Compute residual of the computed solution.
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
                        CALL SGET02( 'No transpose', N, N, NRHS, A, &
                                     LDA, X, LDA, WORK, LDA, RWORK, &
                                     RESULT( 2 ) )
!
!                          Check solution from generated exact solution.
!
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     RCONDC, RESULT( 3 ) )
                        NT = 3
                     END IF
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                     DO K = 1, NT
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALADHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )'SGESV ', N, &
                              IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
                     ENDDO
                     NRUN = NRUN + NT
                  END IF
!
!                    --- Test SGESVX ---
!
                  IF( .NOT.PREFAC )  THEN
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLASET( 'Full', N, N, ZERO, ZERO, AFAC, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ENDIF
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
                  IF( IEQUED > 1 .AND. N > 0 ) THEN
!
!                       Equilibrate the matrix if FACT = 'F' and
!                       EQUED = 'R', 'C', or 'B'.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLAQGE( N, N, A, LDA, S, S( N+1 ), ROWCND, &
                                  COLCND, AMAX, EQUED )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLAQGE : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
!                    Solve the system and compute the condition number
!                    and error bounds using SGESVX.
!
                  SRNAMT = 'SGESVX'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SGESVX( FACT, TRANS, N, NRHS, A, LDA, AFAC, &
                               LDA, IWORK, EQUED, S, S( N+1 ), B, &
                               LDA, X, LDA, RCOND, RWORK, &
                               RWORK( NRHS+1 ), WORK, IWORK( N+1 ), &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGESVX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check the error code from SGESVX.
!
                  IF( INFO /= IZERO ) &
                     CALL ALAERH( PATH, 'SGESVX', INFO, IZERO, &
                                  FACT // TRANS, N, N, -1, -1, NRHS, &
                                  IMAT, NFAIL, NERRS, NOUT )
!
!                    Compare WORK(1) from SGESVX with the computed
!                    reciprocal pivot growth factor RPVGRW
!
                  IF( INFO /= 0 .AND. INFO <= N) THEN
                     RPVGRW = SLANTR( 'M', 'U', 'N', INFO, INFO, &
                              AFAC, LDA, WORK )
                     IF( RPVGRW == ZERO ) THEN
                        RPVGRW = ONE
                     ELSE
                        RPVGRW = SLANGE( 'M', N, INFO, A, LDA, &
                                 WORK ) / RPVGRW
                     END IF
                  ELSE
                     RPVGRW = SLANTR( 'M', 'U', 'N', N, N, AFAC, LDA, &
                              WORK )
                     IF( RPVGRW == ZERO ) THEN
                        RPVGRW = ONE
                     ELSE
                        RPVGRW = SLANGE( 'M', N, N, A, LDA, WORK ) / &
                                 RPVGRW
                     END IF
                  END IF
                  RESULT( 7 ) = ABS( RPVGRW-WORK( 1 ) ) / &
                                MAX( WORK( 1 ), RPVGRW ) / &
                                SLAMCH( 'E' )
!
                  IF( .NOT.PREFAC ) THEN
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                     CALL SGET01( N, N, A, LDA, AFAC, LDA, IWORK, &
                                  RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                     K1 = 1
                  ELSE
                     K1 = 2
                  END IF
!
                  IF( INFO == 0 ) THEN
                     TRFCON = .FALSE.
!
!                       Compute residual of the computed solution.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, &
                                  LDA )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     CALL SGET02( TRANS, N, N, NRHS, ASAV, LDA, X, &
                                  LDA, WORK, LDA, RWORK( 2*NRHS+1 ), &
                                  RESULT( 2 ) )
!
!                       Check solution from generated exact solution.
!
                     IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, &
                         'N' ) ) ) THEN
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     RCONDC, RESULT( 3 ) )
                     ELSE
                        IF( ITRAN == 1 ) THEN
                           ROLDC = ROLDO
                        ELSE
                           ROLDC = ROLDI
                        END IF
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, &
                                     ROLDC, RESULT( 3 ) )
                     END IF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                     CALL SGET07( TRANS, N, NRHS, ASAV, LDA, B, LDA, &
                                  X, LDA, XACT, LDA, RWORK, .TRUE., &
                                  RWORK( NRHS+1 ), RESULT( 4 ) )
                  ELSE
                     TRFCON = .TRUE.
                  END IF
!
!                    Compare RCOND from SGESVX with the computed value
!                    in RCONDC.
!
                  RESULT( 6 ) = SGET06( RCOND, RCONDC )
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                  IF( .NOT.TRFCON ) THEN
                     DO K = K1, NTESTS
                        IF( RESULT( K ) >= THRESH ) THEN
                           IF( NFAIL == 0 .AND. NERRS == 0 ) &
                              CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'SGESVX', &
                                 FACT, TRANS, N, EQUED, IMAT, K, &
                                 RESULT( K )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'SGESVX', &
                                 FACT, TRANS, N, IMAT, K, RESULT( K )
                           END IF
                           NFAIL = NFAIL + 1
                        END IF
                     ENDDO
                     NRUN = NRUN + NTESTS - K1 + 1
                  ELSE
                     IF( RESULT( 1 ) >= THRESH .AND. .NOT.PREFAC ) &
                          THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        IF( PREFAC ) THEN
                           WRITE( NOUT, FMT = 9997 )'SGESVX', FACT, &
                              TRANS, N, EQUED, IMAT, 1, RESULT( 1 )
                        ELSE
                           WRITE( NOUT, FMT = 9998 )'SGESVX', FACT, &
                              TRANS, N, IMAT, 1, RESULT( 1 )
                        END IF
                        NFAIL = NFAIL + 1
                        NRUN = NRUN + 1
                     END IF
                     IF( RESULT( 6 ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        IF( PREFAC ) THEN
                           WRITE( NOUT, FMT = 9997 )'SGESVX', FACT, &
                              TRANS, N, EQUED, IMAT, 6, RESULT( 6 )
                        ELSE
                           WRITE( NOUT, FMT = 9998 )'SGESVX', FACT, &
                              TRANS, N, IMAT, 6, RESULT( 6 )
                        END IF
                        NFAIL = NFAIL + 1
                        NRUN = NRUN + 1
                     END IF
                     IF( RESULT( 7 ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, PATH )
                        IF( PREFAC ) THEN
                           WRITE( NOUT, FMT = 9997 )'SGESVX', FACT, &
                              TRANS, N, EQUED, IMAT, 7, RESULT( 7 )
                        ELSE
                           WRITE( NOUT, FMT = 9998 )'SGESVX', FACT, &
                              TRANS, N, IMAT, 7, RESULT( 7 )
                        END IF
                        NFAIL = NFAIL + 1
                        NRUN = NRUN + 1
                     END IF
!
                  END IF
!
               ENDDO
60          CONTINUE
            ENDDO
         ENDDO
80    CONTINUE
      ENDDO
   ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test(', I2, ') =', &
         G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, &
         ', type ', I2, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, &
         ', EQUED=''', A1, ''', type ', I2, ', test(', I1, ')=', &
         G12.5 )
   RETURN
!
!     End of SDRVGE
!
END
                                                                                                                                                                                                                                                                                                            




