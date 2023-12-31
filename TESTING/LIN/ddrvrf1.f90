!> \brief \b DDRVRF1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVRF1( NOUT, NN, NVAL, THRESH, A, LDA, ARF, WORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, NN, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       DOUBLE PRECISION   A( LDA, * ), ARF( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRVRF1 tests the LAPACK RFP routines:
!>     DLANSF
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
!>          A is DOUBLE PRECISION array, dimension (LDA,NMAX)
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
!>          ARF is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension ( NMAX )
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
   SUBROUTINE DDRVRF1( NOUT, NN, NVAL, THRESH, A, LDA, ARF, WORK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, NN, NOUT
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   INTEGER            NVAL( NN )
   DOUBLE PRECISION   A( LDA, * ), ARF( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Parameters ..
   DOUBLE PRECISION   ONE
   PARAMETER          ( ONE = 1.0D+0 )
   INTEGER            NTESTS
   PARAMETER          ( NTESTS = 1 )
!     ..
!     .. Local Scalars ..
   CHARACTER          UPLO, CFORM, NORM
   INTEGER            I, IFORM, IIN, IIT, INFO, INORM, IUPLO, J, N, &
                      NERRS, NFAIL, NRUN
   DOUBLE PRECISION   EPS, LARGE, NORMA, NORMARF, SMALL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 ), FORMS( 2 ), NORMS( 4 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
   DOUBLE PRECISION   RESULT( NTESTS )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLANSY, DLANSF, DLARND
   EXTERNAL           DLAMCH, DLANSY, DLANSF, DLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           DTRTTF
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
   DATA               NORMS / 'M', '1', 'I', 'F' /
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
   NRUN = 0
   NFAIL = 0
   NERRS = 0
   INFO = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
   EPS = DLAMCH( 'Precision' )
   SMALL = DLAMCH( 'Safe minimum' )
   LARGE = ONE / SMALL
   SMALL = SMALL * LDA * LDA
   LARGE = LARGE / LDA / LDA
!
   DO IIN = 1, NN
!
      N = NVAL( IIN )
!
         DO IIT = 1, 3
!           Nothing to do for N=0
         IF ( N  ==  0 ) EXIT
!
!           IIT = 1 : random matrix
!           IIT = 2 : random matrix scaled near underflow
!           IIT = 3 : random matrix scaled near overflow
!
         DO J = 1, N
            DO I = 1, N
               A( I, J) = DLARND( 2, ISEED )
            END DO
         END DO
!
         IF ( IIT == 2 ) THEN
            DO J = 1, N
               DO I = 1, N
                  A( I, J) = A( I, J ) * LARGE
               END DO
            END DO
         END IF
!
         IF ( IIT == 3 ) THEN
            DO J = 1, N
               DO I = 1, N
                  A( I, J) = A( I, J) * SMALL
               END DO
            END DO
         END IF
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO IUPLO = 1, 2
!
            UPLO = UPLOS( IUPLO )
!
!              Do first for CFORM = 'N', then for CFORM = 'C'
!
            DO IFORM = 1, 2
!
               CFORM = FORMS( IFORM )
!
               SRNAMT = 'DTRTTF'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DTRTTF( CFORM, UPLO, N, A, LDA, ARF, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DTRTTF : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Check error code from DTRTTF
!
               IF( INFO /= 0 ) THEN
                  IF( NFAIL == 0 .AND. NERRS == 0 ) THEN
                     WRITE( NOUT, * )
                     WRITE( NOUT, FMT = 9999 )
                  END IF
                  WRITE( NOUT, FMT = 9998 ) SRNAMT, UPLO, CFORM, N
                  NERRS = NERRS + 1
                  GO TO 100
               END IF
!
               DO INORM = 1, 4
!
!                    Check all four norms: 'M', '1', 'I', 'F'
!
                  NORM = NORMS( INORM )
                  NORMARF = DLANSF( NORM, CFORM, UPLO, N, ARF, WORK )
                  NORMA = DLANSY( NORM, UPLO, N, A, LDA, WORK )
!
                  RESULT(1) = ( NORMA - NORMARF ) / NORMA / EPS
                  NRUN = NRUN + 1
!
                  IF( RESULT(1) >= THRESH ) THEN
                     IF( NFAIL == 0 .AND. NERRS == 0 ) THEN
                        WRITE( NOUT, * )
                        WRITE( NOUT, FMT = 9999 )
                     END IF
                     WRITE( NOUT, FMT = 9997 ) 'DLANSF', &
                         N, IIT, UPLO, CFORM, NORM, RESULT(1)
                     NFAIL = NFAIL + 1
                  END IF
               ENDDO
  100          CONTINUE
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   IF ( NFAIL == 0 ) THEN
      WRITE( NOUT, FMT = 9996 ) 'DLANSF', NRUN
   ELSE
      WRITE( NOUT, FMT = 9995 ) 'DLANSF', NFAIL, NRUN
   END IF
   IF ( NERRS /= 0 ) THEN
      WRITE( NOUT, FMT = 9994 ) NERRS, 'DLANSF'
   END IF
!
 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing DLANSF &
            ***')
 9998 FORMAT( 1X, '     Error in ',A6,' with UPLO=''',A1,''', FORM=''', &
           A1,''', N=',I5)
 9997 FORMAT( 1X, '     Failure in ',A6,' N=',I5,' TYPE=',I5,' UPLO=''', &
           A1, ''', FORM =''',A1,''', NORM=''',A1,''', test=',G12.5)
 9996 FORMAT( 1X, 'All tests for ',A6,' auxiliary routine passed the ', &
           'threshold ( ',I5,' tests run)')
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5, &
           ' tests failed to pass the threshold')
 9994 FORMAT( 26X, I5,' error message recorded (',A6,')')
!
   RETURN
!
!     End of DDRVRF1
!
END
                                                                                                                                                                                                                                                                                                            




