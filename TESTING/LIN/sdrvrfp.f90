!> \brief \b SDRVRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL,
!      +              THRESH, A, ASAV, AFAC, AINV, B,
!      +              BSAV, XACT, X, ARF, ARFINV,
!      +              S_WORK_SLATMS, S_WORK_SPOT01, S_TEMP_SPOT02,
!      +              S_TEMP_SPOT03, S_WORK_SLANSY,
!      +              S_WORK_SPOT02, S_WORK_SPOT03 )
!
!       .. Scalar Arguments ..
!       INTEGER            NN, NNS, NNT, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
!       REAL               A( * )
!       REAL               AINV( * )
!       REAL               ASAV( * )
!       REAL               B( * )
!       REAL               BSAV( * )
!       REAL               AFAC( * )
!       REAL               ARF( * )
!       REAL               ARFINV( * )
!       REAL               XACT( * )
!       REAL               X( * )
!       REAL               S_WORK_SLATMS( * )
!       REAL               S_WORK_SPOT01( * )
!       REAL               S_TEMP_SPOT02( * )
!       REAL               S_TEMP_SPOT03( * )
!       REAL               S_WORK_SLANSY( * )
!       REAL               S_WORK_SPOT02( * )
!       REAL               S_WORK_SPOT03( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SDRVRFP tests the LAPACK RFP routines:
!>     SPFTRF, SPFTRS, and SPFTRI.
!>
!> This testing routine follow the same tests as DDRVPO (test for the full
!> format Symmetric Positive Definite solver).
!>
!> The tests are performed in Full Format, conversion back and forth from
!> full format to RFP format are performed using the routines STRTTF and
!> STFTTR.
!>
!> First, a specific matrix A of size N is created. There is nine types of
!> different matrixes possible.
!>  1. Diagonal                        6. Random, CNDNUM = sqrt(0.1/EPS)
!>  2. Random, CNDNUM = 2              7. Random, CNDNUM = 0.1/EPS
!> *3. First row and column zero       8. Scaled near underflow
!> *4. Last row and column zero        9. Scaled near overflow
!> *5. Middle row and column zero
!> (* - tests error exits from SPFTRF, no test ratios are computed)
!> A solution XACT of size N-by-NRHS is created and the associated right
!> hand side B as well. Then SPFTRF is called to compute L (or U), the
!> Cholesky factor of A. Then L (or U) is used to solve the linear system
!> of equations AX = B. This gives X. Then L (or U) is used to compute the
!> inverse of A, AINV. The following four tests are then performed:
!> (1) norm( L*L' - A ) / ( N * norm(A) * EPS ) or
!>     norm( U'*U - A ) / ( N * norm(A) * EPS ),
!> (2) norm(B - A*X) / ( norm(A) * norm(X) * EPS ),
!> (3) norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
!> (4) ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS ),
!> where EPS is the machine precision, RCOND the condition number of A, and
!> norm( . ) the 1-norm for (1,2,3) and the inf-norm for (4).
!> Errors occur when INFO parameter is not as expected. Failures occur when
!> a test ratios is greater than THRES.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>                The unit number for output.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>                The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>                The values of the matrix dimension N.
!> \endverbatim
!>
!> \param[in] NNS
!> \verbatim
!>          NNS is INTEGER
!>                The number of values of NRHS contained in the vector NSVAL.
!> \endverbatim
!>
!> \param[in] NSVAL
!> \verbatim
!>          NSVAL is INTEGER array, dimension (NNS)
!>                The values of the number of right-hand sides NRHS.
!> \endverbatim
!>
!> \param[in] NNT
!> \verbatim
!>          NNT is INTEGER
!>                The number of values of MATRIX TYPE contained in the vector NTVAL.
!> \endverbatim
!>
!> \param[in] NTVAL
!> \verbatim
!>          NTVAL is INTEGER array, dimension (NNT)
!>                The values of matrix type (between 0 and 9 for PO/PP/PF matrices).
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>                The threshold value for the test ratios.  A result is
!>                included in the output file if RESULT >= THRESH.  To have
!>                every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is REAL array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] ARF
!> \verbatim
!>          ARF is REAL array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] ARFINV
!> \verbatim
!>          ARFINV is REAL array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] S_WORK_SLATMS
!> \verbatim
!>          S_WORK_SLATMS is REAL array, dimension ( 3*NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_SPOT01
!> \verbatim
!>          S_WORK_SPOT01 is REAL array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] S_TEMP_SPOT02
!> \verbatim
!>          S_TEMP_SPOT02 is REAL array, dimension ( NMAX*MAXRHS )
!> \endverbatim
!>
!> \param[out] S_TEMP_SPOT03
!> \verbatim
!>          S_TEMP_SPOT03 is REAL array, dimension ( NMAX*NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_SLANSY
!> \verbatim
!>          S_WORK_SLANSY is REAL array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_SPOT02
!> \verbatim
!>          S_WORK_SPOT02 is REAL array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_SPOT03
!> \verbatim
!>          S_WORK_SPOT03 is REAL array, dimension ( NMAX )
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
   SUBROUTINE SDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, &
                 THRESH, A, ASAV, AFAC, AINV, B, &
                 BSAV, XACT, X, ARF, ARFINV, &
                 S_WORK_SLATMS, S_WORK_SPOT01, S_TEMP_SPOT02, &
                 S_TEMP_SPOT03, S_WORK_SLANSY, &
                 S_WORK_SPOT02, S_WORK_SPOT03 )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NN, NNS, NNT, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
   REAL               A( * )
   REAL               AINV( * )
   REAL               ASAV( * )
   REAL               B( * )
   REAL               BSAV( * )
   REAL               AFAC( * )
   REAL               ARF( * )
   REAL               ARFINV( * )
   REAL               XACT( * )
   REAL               X( * )
   REAL               S_WORK_SLATMS( * )
   REAL               S_WORK_SPOT01( * )
   REAL               S_TEMP_SPOT02( * )
   REAL               S_TEMP_SPOT03( * )
   REAL               S_WORK_SLANSY( * )
   REAL               S_WORK_SPOT02( * )
   REAL               S_WORK_SPOT03( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 4 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZEROT
   INTEGER            I, INFO, IUPLO, LDA, LDB, IMAT, NERRS, NFAIL, &
                      NRHS, NRUN, IZERO, IOFF, K, NT, N, IFORM, IIN, &
                      IIT, IIS
   CHARACTER          DIST, CTYPE, UPLO, CFORM
   INTEGER            KL, KU, MODE
   REAL               ANORM, AINVNM, CNDNUM, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 ), FORMS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   REAL               SLANSY
   EXTERNAL           SLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, SGET04, STFTTR, SLACPY, &
                      SLARHS, SLATB4, SLATMS, SPFTRI, SPFTRF, SPFTRS, &
                      SPOT01, SPOT02, SPOT03, SPOTRI, SPOTRF, STRTTF
!     ..
!     .. Scalars in Common ..
   CHARACTER*32       SRNAMT
!     ..
!     .. Common blocks ..
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS / 'U', 'L' /
   DATA               FORMS / 'N', 'T' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
   DO IIN = 1, NN
!
      N = NVAL( IIN )
      LDA = MAX( N, 1 )
      LDB = MAX( N, 1 )
!
      DO IIS = 1, NNS
!
         NRHS = NSVAL( IIS )
!
         DO IIT = 1, NNT
!
            IMAT = NTVAL( IIT )
!
!              If N == 0, only consider the first type
!
            IF( N == 0 .AND. IIT >= 1 ) GO TO 120
!
!              Skip types 3, 4, or 5 if the matrix size is too small.
!
            IF( IMAT == 4 .AND. N <= 1 ) GO TO 120
            IF( IMAT == 5 .AND. N <= 2 ) GO TO 120
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
            DO IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
!
!                 Do first for CFORM = 'N', then for CFORM = 'C'
!
               DO IFORM = 1, 2
                  CFORM = FORMS( IFORM )
!
!                    Set up parameters with SLATB4 and generate a test
!                    matrix with SLATMS.
!
                  CALL SLATB4( 'SPO', IMAT, N, N, CTYPE, KL, KU, &
                               ANORM, MODE, CNDNUM, DIST )
!
                  SRNAMT = 'SLATMS'
                  CALL SLATMS( N, N, DIST, ISEED, CTYPE, &
                               S_WORK_SLATMS, &
                               MODE, CNDNUM, ANORM, KL, KU, UPLO, A, &
                               LDA, S_WORK_SLATMS, INFO )
!
!                    Check error code from SLATMS.
!
                  IF( INFO /= 0 ) THEN
                     CALL ALAERH( 'SPF', 'SLATMS', INFO, 0, UPLO, N, &
                                  N, -1, -1, -1, IIT, NFAIL, NERRS, &
                                  NOUT )
                     GO TO 100
                  END IF
!
!                    For types 3-5, zero one row and column of the matrix to
!                    test that INFO is returned correctly.
!
                  ZEROT = IMAT >= 3 .AND. IMAT <= 5
                  IF( ZEROT ) THEN
                     IF( IIT == 3 ) THEN
                        IZERO = 1
                     ELSE IF( IIT == 4 ) THEN
                        IZERO = N
                     ELSE
                        IZERO = N / 2 + 1
                     END IF
                     IOFF = ( IZERO-1 )*LDA
!
!                       Set row and column IZERO of A to 0.
!
                     IF( IUPLO == 1 ) THEN
                        DO I = 1, IZERO - 1
                           A( IOFF+I ) = ZERO
                        ENDDO
                        IOFF = IOFF + IZERO
                        DO I = IZERO, N
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
                        ENDDO
                     ELSE
                        IOFF = IZERO
                        DO I = 1, IZERO - 1
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
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
!                    Save a copy of the matrix A in ASAV.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( UPLO, N, N, A, LDA, ASAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compute the condition number of A (RCONDC).
!
                  IF( ZEROT ) THEN
                     RCONDC = ZERO
                  ELSE
!
!                       Compute the 1-norm of A.
!
                     ANORM = SLANSY( '1', UPLO, N, A, LDA, &
                            S_WORK_SLANSY )
!
!                       Factor the matrix A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SPOTRF( UPLO, N, A, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SPOTRF : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Form the inverse of A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SPOTRI( UPLO, N, A, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SPOTRI : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif

                     IF ( N  /=  0 ) THEN
!
!                          Compute the 1-norm condition number of A.
!
                        AINVNM = SLANSY( '1', UPLO, N, A, LDA, &
                              S_WORK_SLANSY )
                        RCONDC = ( ONE / ANORM ) / AINVNM
!
!                          Restore the matrix A.
!
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL SLACPY( UPLO, N, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                     END IF
!
                  END IF
!
!                    Form an exact solution and set the right hand side.
!
                  SRNAMT = 'SLARHS'
                  CALL SLARHS( 'SPO', 'N', UPLO, ' ', N, N, KL, KU, &
                               NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
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
!                    Compute the L*L' or U'*U factorization of the
!                    matrix and solve the system.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
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
                  CALL SLACPY( 'Full', N, NRHS, B, LDB, X, LDB )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'STRTTF'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL STRTTF( CFORM, UPLO, N, AFAC, LDA, ARF, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : STRTTF : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  SRNAMT = 'SPFTRF'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SPFTRF( CFORM, UPLO, N, ARF, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SPFTRF : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from SPFTRF.
!
                  IF( INFO /= IZERO ) THEN
!
!                       LANGOU: there is a small hick here: IZERO should
!                       always be INFO however if INFO is ZERO, ALAERH does not
!                       complain.
!
                      CALL ALAERH( 'SPF', 'SPFSV ', INFO, IZERO, &
                                   UPLO, N, N, -1, -1, NRHS, IIT, &
                                   NFAIL, NERRS, NOUT )
                      GO TO 100
                   END IF
!
!                    Skip the tests if INFO is not 0.
!
                  IF( INFO /= 0 ) THEN
                     GO TO 100
                  END IF
!
                  SRNAMT = 'SPFTRS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SPFTRS( CFORM, UPLO, N, NRHS, ARF, X, LDB, &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SPFTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'STFTTR'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL STFTTR( CFORM, UPLO, N, ARF, AFAC, LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : STFTTR : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( UPLO, N, N, AFAC, LDA, ASAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  CALL SPOT01( UPLO, N, A, LDA, AFAC, LDA, &
                                S_WORK_SPOT01, RESULT( 1 ) )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( UPLO, N, N, ASAV, LDA, AFAC, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Form the inverse and compute the residual.
!
                  IF(MOD(N,2) == 0)THEN
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'A', N+1, N/2, ARF, N+1, ARFINV, &
                                  N+1 )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ELSE
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL SLACPY( 'A', N, (N+1)/2, ARF, N, ARFINV, &
                                  N )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
                  SRNAMT = 'SPFTRI'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SPFTRI( CFORM, UPLO, N, ARFINV , INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SPFTRI : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'STFTTR'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL STFTTR( CFORM, UPLO, N, ARFINV, AINV, LDA, &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : STFTTR : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from SPFTRI.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( 'SPO', 'SPFTRI', INFO, 0, UPLO, N, &
                                  N, -1, -1, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
                  CALL SPOT03( UPLO, N, A, LDA, AINV, LDA, &
                               S_TEMP_SPOT03, LDA, S_WORK_SPOT03, &
                               RCONDC, RESULT( 2 ) )
!
!                    Compute residual of the computed solution.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, &
                               S_TEMP_SPOT02, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  CALL SPOT02( UPLO, N, NRHS, A, LDA, X, LDA, &
                               S_TEMP_SPOT02, LDA, S_WORK_SPOT02, &
                               RESULT( 3 ) )
!
!                    Check solution from generated exact solution.

                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 4 ) )
                  NT = 4
!
!                    Print information about the tests that did not
!                    pass the threshold.
!
                  DO K = 1, NT
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, 'SPF' )
                        WRITE( NOUT, FMT = 9999 )'SPFSV ', UPLO, &
                               N, IIT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
                  ENDDO
                  NRUN = NRUN + NT
  100             CONTINUE
                  ENDDO
               ENDDO
  120       CONTINUE
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   CALL ALASVM( 'SPF', NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, &
         ', test(', I1, ')=', G12.5 )
!
   RETURN
!
!     End of SDRVRFP
!
END
                                                                                                                                                                                                                                                                                                            




