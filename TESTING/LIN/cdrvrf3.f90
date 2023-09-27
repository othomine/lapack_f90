!> \brief \b CDRVRF3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CDRVRF3( NOUT, NN, NVAL, THRESH, A, LDA, ARF, B1, B2,
!      +                    S_WORK_CLANGE, C_WORK_CGEQRF, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, NN, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       REAL               S_WORK_CLANGE( * )
!       COMPLEX            A( LDA, * ), ARF( * ), B1( LDA, * ),
!      +                   B2( LDA, * )
!       COMPLEX            C_WORK_CGEQRF( * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CDRVRF3 tests the LAPACK RFP routines:
!>     CTFSM
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
!>          A is COMPLEX*16 array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>                The leading dimension of the array A.  LDA >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] ARF
!> \verbatim
!>          ARF is COMPLEX array, dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] B1
!> \verbatim
!>          B1 is COMPLEX array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[out] B2
!> \verbatim
!>          B2 is COMPLEX array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[out] S_WORK_CLANGE
!> \verbatim
!>          S_WORK_CLANGE is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] C_WORK_CGEQRF
!> \verbatim
!>          C_WORK_CGEQRF is COMPLEX array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (NMAX)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex_lin
!
!  =====================================================================
   SUBROUTINE CDRVRF3( NOUT, NN, NVAL, THRESH, A, LDA, ARF, B1, B2, &
                       S_WORK_CLANGE, C_WORK_CGEQRF, TAU )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, NN, NOUT
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            NVAL( NN )
   REAL               S_WORK_CLANGE( * )
   COMPLEX            A( LDA, * ), ARF( * ), B1( LDA, * ), &
                      B2( LDA, * )
   COMPLEX            C_WORK_CGEQRF( * ), TAU( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Parameters ..
   COMPLEX            ZERO, ONE
   PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) , &
                        ONE  = ( 1.0E+0, 0.0E+0 ) )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 1 )
!     ..
!     .. Local Scalars ..
   CHARACTER          UPLO, CFORM, DIAG, TRANS, SIDE
   INTEGER            I, IFORM, IIM, IIN, INFO, IUPLO, J, M, N, NA, &
                      NFAIL, NRUN, ISIDE, IDIAG, IALPHA, ITRANS
   COMPLEX            ALPHA
   REAL               EPS
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 ), FORMS( 2 ), TRANSS( 2 ), &
                      DIAGS( 2 ), SIDES( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   REAL               RESULT( NTESTS )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLAMCH, CLANGE
   COMPLEX            CLARND
   EXTERNAL           SLAMCH, CLARND, CLANGE, LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CTRTTF, CGEQRF, CGEQLF, CTFSM, CTRSM
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, SQRT
!     ..
!     .. Scalars in Common ..
   CHARACTER*32       SRNAMT
!     ..
!     .. Common blocks ..
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               ISEEDY / 1988, 1989, 1990, 1991 /
   DATA               UPLOS  / 'U', 'L' /
   DATA               FORMS  / 'N', 'C' /
   DATA               SIDES  / 'L', 'R' /
   DATA               TRANSS / 'N', 'C' /
   DATA               DIAGS  / 'N', 'U' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   NRUN = 0
   NFAIL = 0
   INFO = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
   EPS = SLAMCH( 'Precision' )
!
   DO IIM = 1, NN
!
      M = NVAL( IIM )
!
      DO IIN = 1, NN
!
         N = NVAL( IIN )
!
         DO IFORM = 1, 2
!
            CFORM = FORMS( IFORM )
!
            DO IUPLO = 1, 2
!
               UPLO = UPLOS( IUPLO )
!
               DO ISIDE = 1, 2
!
                  SIDE = SIDES( ISIDE )
!
                  DO ITRANS = 1, 2
!
                     TRANS = TRANSS( ITRANS )
!
                     DO IDIAG = 1, 2
!
                        DIAG = DIAGS( IDIAG )
!
                        DO IALPHA = 1, 3
!
                           IF ( IALPHA == 1 ) THEN
                              ALPHA = ZERO
                           ELSE IF ( IALPHA == 2 ) THEN
                              ALPHA = ONE
                           ELSE
                              ALPHA = CLARND( 4, ISEED )
                           END IF
!
!                             All the parameters are set:
!                                CFORM, SIDE, UPLO, TRANS, DIAG, M, N,
!                                and ALPHA
!                             READY TO TEST!
!
                           NRUN = NRUN + 1
!
                           IF ( ISIDE == 1 ) THEN
!
!                                The case ISIDE == 1 is when SIDE == 'L'
!                                -> A is M-by-M ( B is M-by-N )
!
                              NA = M
!
                           ELSE
!
!                                The case ISIDE == 2 is when SIDE == 'R'
!                                -> A is N-by-N ( B is M-by-N )
!
                              NA = N
!
                           END IF
!
!                             Generate A our NA--by--NA triangular
!                             matrix.
!                             Our test is based on forward error so we
!                             do want A to be well conditioned! To get
!                             a well-conditioned triangular matrix, we
!                             take the R factor of the QR/LQ factorization
!                             of a random matrix.
!
                           DO J = 1, NA
                              DO I = 1, NA
                                 A( I, J ) = CLARND( 4, ISEED )
                              END DO
                           END DO
!
                           IF ( IUPLO == 1 ) THEN
!
!                                The case IUPLO == 1 is when SIDE == 'U'
!                                -> QR factorization.
!
                              SRNAMT = 'CGEQRF'
                              CALL CGEQRF( NA, NA, A, LDA, TAU, &
                                           C_WORK_CGEQRF, LDA, &
                                           INFO )
!
!                                Forcing main diagonal of test matrix to
!                                be unit makes it ill-conditioned for
!                                some test cases
!
                              IF ( LSAME( DIAG, 'U' ) ) THEN
                                 DO J = 1, NA
                                    DO I = 1, J
                                       A( I, J ) = A( I, J ) / &
                                               ( 2.0 * A( J, J ) )
                                    END DO
                                 END DO
                              END IF
!
                           ELSE
!
!                                The case IUPLO == 2 is when SIDE == 'L'
!                                -> QL factorization.
!
                              SRNAMT = 'CGELQF'
                              CALL CGELQF( NA, NA, A, LDA, TAU, &
                                           C_WORK_CGEQRF, LDA, &
                                           INFO )
!
!                                Forcing main diagonal of test matrix to
!                                be unit makes it ill-conditioned for
!                                some test cases
!
                              IF ( LSAME( DIAG, 'U' ) ) THEN
                                 DO I = 1, NA
                                    DO J = 1, I
                                       A( I, J ) = A( I, J ) / &
                                               ( 2.0 * A( I, I ) )
                                    END DO
                                 END DO
                              END IF
!
                           END IF
!
!                             After the QR factorization, the diagonal
!                             of A is made of real numbers, we multiply
!                             by a random complex number of absolute
!                             value 1.0E+00.
!
                           DO J = 1, NA
                              A( J, J ) = A( J, J ) * &
                                      CLARND( 5, ISEED )
                           END DO
!
!                             Store a copy of A in RFP format (in ARF).
!
                           SRNAMT = 'CTRTTF'
                           CALL CTRTTF( CFORM, UPLO, NA, A, LDA, ARF, &
                                        INFO )
!
!                             Generate B1 our M--by--N right-hand side
!                             and store a copy in B2.
!
                           DO J = 1, N
                              DO I = 1, M
                                 B1( I, J ) = CLARND( 4, ISEED )
                                 B2( I, J ) = B1( I, J )
                              END DO
                           END DO
!
!                             Solve op( A ) X = B or X op( A ) = B
!                             with CTRSM
!
                           SRNAMT = 'CTRSM'
                           CALL CTRSM( SIDE, UPLO, TRANS, DIAG, M, N, &
                                  ALPHA, A, LDA, B1, LDA )
!
!                             Solve op( A ) X = B or X op( A ) = B
!                             with CTFSM
!
                           SRNAMT = 'CTFSM'
                           CALL CTFSM( CFORM, SIDE, UPLO, TRANS, &
                                       DIAG, M, N, ALPHA, ARF, B2, &
                                       LDA )
!
!                             Check that the result agrees.
!
                           DO J = 1, N
                              DO I = 1, M
                                 B1( I, J ) = B2( I, J ) - B1( I, J )
                              END DO
                           END DO
!
                           RESULT( 1 ) = CLANGE( 'I', M, N, B1, LDA, &
                                               S_WORK_CLANGE )
!
                           RESULT( 1 ) = RESULT( 1 ) / SQRT( EPS ) &
                                       / MAX ( MAX( M, N ), 1 )
!
                           IF( RESULT( 1 ) >= THRESH ) THEN
                              IF( NFAIL == 0 ) THEN
                                 WRITE( NOUT, * )
                                 WRITE( NOUT, FMT = 9999 )
                              END IF
                              WRITE( NOUT, FMT = 9997 ) 'CTFSM', &
                                 CFORM, SIDE, UPLO, TRANS, DIAG, M, &
                                 N, RESULT( 1 )
                              NFAIL = NFAIL + 1
                           END IF
!
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   IF ( NFAIL == 0 ) THEN
      WRITE( NOUT, FMT = 9996 ) 'CTFSM', NRUN
   ELSE
      WRITE( NOUT, FMT = 9995 ) 'CTFSM', NFAIL, NRUN
   END IF
!
 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing CTFSM &
            ***')
 9997 FORMAT( 1X, '     Failure in ',A5,', CFORM=''',A1,''',', &
    ' SIDE=''',A1,''',',' UPLO=''',A1,''',',' TRANS=''',A1,''',', &
    ' DIAG=''',A1,''',',' M=',I3,', N =', I3,', test=',G12.5)
 9996 FORMAT( 1X, 'All tests for ',A5,' auxiliary routine passed the ', &
           'threshold ( ',I5,' tests run)')
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5, &
           ' tests failed to pass the threshold')
!
   RETURN
!
!     End of CDRVRF3
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
