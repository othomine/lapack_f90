!> \brief \b DDRVRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL,
!      +              THRESH, A, ASAV, AFAC, AINV, B,
!      +              BSAV, XACT, X, ARF, ARFINV,
!      +              D_WORK_DLATMS, D_WORK_DPOT01, D_TEMP_DPOT02,
!      +              D_TEMP_DPOT03, D_WORK_DLANSY,
!      +              D_WORK_DPOT02, D_WORK_DPOT03 )
!
!       .. Scalar Arguments ..
!       INTEGER            NN, NNS, NNT, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
!       DOUBLE PRECISION   A( * )
!       DOUBLE PRECISION   AINV( * )
!       DOUBLE PRECISION   ASAV( * )
!       DOUBLE PRECISION   B( * )
!       DOUBLE PRECISION   BSAV( * )
!       DOUBLE PRECISION   AFAC( * )
!       DOUBLE PRECISION   ARF( * )
!       DOUBLE PRECISION   ARFINV( * )
!       DOUBLE PRECISION   XACT( * )
!       DOUBLE PRECISION   X( * )
!       DOUBLE PRECISION   D_WORK_DLATMS( * )
!       DOUBLE PRECISION   D_WORK_DPOT01( * )
!       DOUBLE PRECISION   D_TEMP_DPOT02( * )
!       DOUBLE PRECISION   D_TEMP_DPOT03( * )
!       DOUBLE PRECISION   D_WORK_DLANSY( * )
!       DOUBLE PRECISION   D_WORK_DPOT02( * )
!       DOUBLE PRECISION   D_WORK_DPOT03( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRVRFP tests the LAPACK RFP routines:
!>     DPFTRF, DPFTRS, and DPFTRI.
!>
!> This testing routine follow the same tests as DDRVPO (test for the full
!> format Symmetric Positive Definite solver).
!>
!> The tests are performed in Full Format, conversion back and forth from
!> full format to RFP format are performed using the routines DTRTTF and
!> DTFTTR.
!>
!> First, a specific matrix A of size N is created. There is nine types of
!> different matrixes possible.
!>  1. Diagonal                        6. Random, CNDNUM = sqrt(0.1/EPS)
!>  2. Random, CNDNUM = 2              7. Random, CNDNUM = 0.1/EPS
!> *3. First row and column zero       8. Scaled near underflow
!> *4. Last row and column zero        9. Scaled near overflow
!> *5. Middle row and column zero
!> (* - tests error exits from DPFTRF, no test ratios are computed)
!> A solution XACT of size N-by-NRHS is created and the associated right
!> hand side B as well. Then DPFTRF is called to compute L (or U), the
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
!>          THRESH is DOUBLE PRECISION
!>                The threshold value for the test ratios.  A result is
!>                included in the output file if RESULT >= THRESH.  To have
!>                every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] ARF
!> \verbatim
!>          ARF is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] ARFINV
!> \verbatim
!>          ARFINV is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] D_WORK_DLATMS
!> \verbatim
!>          D_WORK_DLATMS is DOUBLE PRECISION array, dimension ( 3*NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DPOT01
!> \verbatim
!>          D_WORK_DPOT01 is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_TEMP_DPOT02
!> \verbatim
!>          D_TEMP_DPOT02 is DOUBLE PRECISION array, dimension ( NMAX*MAXRHS )
!> \endverbatim
!>
!> \param[out] D_TEMP_DPOT03
!> \verbatim
!>          D_TEMP_DPOT03 is DOUBLE PRECISION array, dimension ( NMAX*NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DLANSY
!> \verbatim
!>          D_WORK_DLANSY is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DPOT02
!> \verbatim
!>          D_WORK_DPOT02 is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DPOT03
!> \verbatim
!>          D_WORK_DPOT03 is DOUBLE PRECISION array, dimension ( NMAX )
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
   SUBROUTINE DDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, &
                 THRESH, A, ASAV, AFAC, AINV, B, &
                 BSAV, XACT, X, ARF, ARFINV, &
                 D_WORK_DLATMS, D_WORK_DPOT01, D_TEMP_DPOT02, &
                 D_TEMP_DPOT03, D_WORK_DLANSY, &
                 D_WORK_DPOT02, D_WORK_DPOT03 )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NN, NNS, NNT, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
   DOUBLE PRECISION   A( * )
   DOUBLE PRECISION   AINV( * )
   DOUBLE PRECISION   ASAV( * )
   DOUBLE PRECISION   B( * )
   DOUBLE PRECISION   BSAV( * )
   DOUBLE PRECISION   AFAC( * )
   DOUBLE PRECISION   ARF( * )
   DOUBLE PRECISION   ARFINV( * )
   DOUBLE PRECISION   XACT( * )
   DOUBLE PRECISION   X( * )
   DOUBLE PRECISION   D_WORK_DLATMS( * )
   DOUBLE PRECISION   D_WORK_DPOT01( * )
   DOUBLE PRECISION   D_TEMP_DPOT02( * )
   DOUBLE PRECISION   D_TEMP_DPOT03( * )
   DOUBLE PRECISION   D_WORK_DLANSY( * )
   DOUBLE PRECISION   D_WORK_DPOT02( * )
   DOUBLE PRECISION   D_WORK_DPOT03( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
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
   DOUBLE PRECISION   ANORM, AINVNM, CNDNUM, RCONDC
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 ), FORMS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLANSY
   EXTERNAL           DLANSY
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALADHD, ALAERH, ALASVM, DGET04, DTFTTR, DLACPY, &
                      DLARHS, DLATB4, DLATMS, DPFTRI, DPFTRF, DPFTRS, &
                      DPOT01, DPOT02, DPOT03, DPOTRI, DPOTRF, DTRTTF
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
!                    Set up parameters with DLATB4 and generate a test
!                    matrix with DLATMS.
!
                  CALL DLATB4( 'DPO', IMAT, N, N, CTYPE, KL, KU, &
                               ANORM, MODE, CNDNUM, DIST )
!
                  SRNAMT = 'DLATMS'
                  CALL DLATMS( N, N, DIST, ISEED, CTYPE, &
                               D_WORK_DLATMS, &
                               MODE, CNDNUM, ANORM, KL, KU, UPLO, A, &
                               LDA, D_WORK_DLATMS, INFO )
!
!                    Check error code from DLATMS.
!
                  IF( INFO /= 0 ) THEN
                     CALL ALAERH( 'DPF', 'DLATMS', INFO, 0, UPLO, N, &
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
                  CALL DLACPY( UPLO, N, N, A, LDA, ASAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
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
                     ANORM = DLANSY( '1', UPLO, N, A, LDA, &
                            D_WORK_DLANSY )
!
!                       Factor the matrix A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DPOTRF( UPLO, N, A, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DPOTRF : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
!
!                       Form the inverse of A.
!
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DPOTRI( UPLO, N, A, LDA, INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DPOTRI : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif

                     IF ( N  /=  0 ) THEN
!
!                          Compute the 1-norm condition number of A.
!
                        AINVNM = DLANSY( '1', UPLO, N, A, LDA, &
                              D_WORK_DLANSY )
                        RCONDC = ( ONE / ANORM ) / AINVNM
!
!                          Restore the matrix A.
!
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL DLACPY( UPLO, N, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                     END IF
!
                  END IF
!
!                    Form an exact solution and set the right hand side.
!
                  SRNAMT = 'DLARHS'
                  CALL DLARHS( 'DPO', 'N', UPLO, ' ', N, N, KL, KU, &
                               NRHS, A, LDA, XACT, LDA, B, LDA, &
                               ISEED, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
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
                  CALL DLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
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
                  CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDB )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'DTRTTF'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DTRTTF( CFORM, UPLO, N, AFAC, LDA, ARF, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DTRTTF : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  SRNAMT = 'DPFTRF'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DPFTRF( CFORM, UPLO, N, ARF, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DPFTRF : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from DPFTRF.
!
                  IF( INFO /= IZERO ) THEN
!
!                       LANGOU: there is a small hick here: IZERO should
!                       always be INFO however if INFO is ZERO, ALAERH does not
!                       complain.
!
                      CALL ALAERH( 'DPF', 'DPFSV ', INFO, IZERO, &
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
                  SRNAMT = 'DPFTRS'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DPFTRS( CFORM, UPLO, N, NRHS, ARF, X, LDB, &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DPFTRS : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'DTFTTR'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DTFTTR( CFORM, UPLO, N, ARF, AFAC, LDA, INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DTFTTR : ',&
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
                  CALL DLACPY( UPLO, N, N, AFAC, LDA, ASAV, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  CALL DPOT01( UPLO, N, A, LDA, AFAC, LDA, &
                                D_WORK_DPOT01, RESULT( 1 ) )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DLACPY( UPLO, N, N, ASAV, LDA, AFAC, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
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
                     CALL DLACPY( 'A', N+1, N/2, ARF, N+1, ARFINV, &
                                  N+1 )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  ELSE
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DLACPY( 'A', N, (N+1)/2, ARF, N, ARFINV, &
                                  N )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                  END IF
!
                  SRNAMT = 'DPFTRI'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DPFTRI( CFORM, UPLO, N, ARFINV , INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DPFTRI : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
                  SRNAMT = 'DTFTTR'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DTFTTR( CFORM, UPLO, N, ARFINV, AINV, LDA, &
                               INFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DTFTTR : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Check error code from DPFTRI.
!
                  IF( INFO /= 0 ) &
                     CALL ALAERH( 'DPO', 'DPFTRI', INFO, 0, UPLO, N, &
                                  N, -1, -1, -1, IMAT, NFAIL, NERRS, &
                                  NOUT )
!
                  CALL DPOT03( UPLO, N, A, LDA, AINV, LDA, &
                               D_TEMP_DPOT03, LDA, D_WORK_DPOT03, &
                               RCONDC, RESULT( 2 ) )
!
!                    Compute residual of the computed solution.
!
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL DLACPY( 'Full', N, NRHS, B, LDA, &
                               D_TEMP_DPOT02, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  CALL DPOT02( UPLO, N, NRHS, A, LDA, X, LDA, &
                               D_TEMP_DPOT02, LDA, D_WORK_DPOT02, &
                               RESULT( 3 ) )
!
!                    Check solution from generated exact solution.

                  CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, &
                            RESULT( 4 ) )
                  NT = 4
!
!                    Print information about the tests that did not
!                    pass the threshold.
!
                  DO K = 1, NT
                     IF( RESULT( K ) >= THRESH ) THEN
                        IF( NFAIL == 0 .AND. NERRS == 0 ) &
                           CALL ALADHD( NOUT, 'DPF' )
                        WRITE( NOUT, FMT = 9999 )'DPFSV ', UPLO, &
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
   CALL ALASVM( 'DPF', NOUT, NFAIL, NRUN, NERRS )
!
 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, &
         ', test(', I1, ')=', G12.5 )
!
   RETURN
!
!     End of DDRVRFP
!
END
                                                                                                                                                                                                                                                                                                            




