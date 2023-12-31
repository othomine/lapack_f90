!> \brief \b DDRVRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVRF2( NOUT, NN, NVAL, A, LDA, ARF, AP, ASAV  )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, NN, NOUT
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       DOUBLE PRECISION   A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRVRF2 tests the LAPACK RFP conversion routines.
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
!> \param[out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is DOUBLE PRECISION array, dimension (LDA,NMAX)
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
   SUBROUTINE DDRVRF2( NOUT, NN, NVAL, A, LDA, ARF, AP, ASAV  )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, NN, NOUT
!     ..
!     .. Array Arguments ..
   INTEGER            NVAL( NN )
   DOUBLE PRECISION   A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            LOWER, OK1, OK2
   CHARACTER          UPLO, CFORM
   INTEGER            I, IFORM, IIN, INFO, IUPLO, J, N, &
                      NERRS, NRUN
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          UPLOS( 2 ), FORMS( 2 )
   INTEGER            ISEED( 4 ), ISEEDY( 4 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLARND
   EXTERNAL           DLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           DTFTTR, DTFTTP, DTRTTF, DTRTTP, DTPTTR, DTPTTF
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
   NERRS = 0
   INFO = 0
   DO I = 1, 4
      ISEED( I ) = ISEEDY( I )
   ENDDO
!
   DO IIN = 1, NN
!
      N = NVAL( IIN )
!
!        Do first for UPLO = 'U', then for UPLO = 'L'
!
      DO IUPLO = 1, 2
!
         UPLO = UPLOS( IUPLO )
         LOWER = .TRUE.
         IF ( IUPLO == 1 ) LOWER = .FALSE.
!
!           Do first for CFORM = 'N', then for CFORM = 'T'
!
         DO IFORM = 1, 2
!
            CFORM = FORMS( IFORM )
!
            NRUN = NRUN + 1
!
            DO J = 1, N
               DO I = 1, N
                  A( I, J) = DLARND( 2, ISEED )
               END DO
            END DO
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
            SRNAMT = 'DTFTTP'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DTFTTP( CFORM, UPLO, N, ARF, AP, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DTFTTP : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            SRNAMT = 'DTPTTR'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DTPTTR( UPLO, N, AP, ASAV, LDA, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DTPTTR : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            OK1 = .TRUE.
            IF ( LOWER ) THEN
               DO J = 1, N
                  DO I = J, N
                     IF ( A(I,J) /= ASAV(I,J) ) THEN
                        OK1 = .FALSE.
                     END IF
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = 1, J
                     IF ( A(I,J) /= ASAV(I,J) ) THEN
                        OK1 = .FALSE.
                     END IF
                  END DO
               END DO
            END IF
!
            NRUN = NRUN + 1
!
            SRNAMT = 'DTRTTP'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DTRTTP( UPLO, N, A, LDA, AP, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DTRTTP : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            SRNAMT = 'DTPTTF'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DTPTTF( CFORM, UPLO, N, AP, ARF, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DTPTTF : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            SRNAMT = 'DTFTTR'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DTFTTR( CFORM, UPLO, N, ARF, ASAV, LDA, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DTFTTR : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            OK2 = .TRUE.
            IF ( LOWER ) THEN
               DO J = 1, N
                  DO I = J, N
                     IF ( A(I,J) /= ASAV(I,J) ) THEN
                        OK2 = .FALSE.
                     END IF
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = 1, J
                     IF ( A(I,J) /= ASAV(I,J) ) THEN
                        OK2 = .FALSE.
                     END IF
                  END DO
               END DO
            END IF
!
            IF (( .NOT.OK1 ).OR.( .NOT.OK2 )) THEN
               IF( NERRS == 0 ) THEN
                  WRITE( NOUT, * )
                  WRITE( NOUT, FMT = 9999 )
               END IF
               WRITE( NOUT, FMT = 9998 ) N, UPLO, CFORM
               NERRS = NERRS + 1
            END IF
!
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
   IF ( NERRS == 0 ) THEN
      WRITE( NOUT, FMT = 9997 ) NRUN
   ELSE
      WRITE( NOUT, FMT = 9996 ) NERRS, NRUN
   END IF
!
 9999 FORMAT( 1X, ' *** Error(s) while testing the RFP conversion', &
            ' routines ***')
 9998 FORMAT( 1X, '     Error in RFP,conversion routines N=',I5, &
           ' UPLO=''', A1, ''', FORM =''',A1,'''')
 9997 FORMAT( 1X, 'All tests for the RFP conversion routines passed ( ', &
           I5,' tests run)')
 9996 FORMAT( 1X, 'RFP conversion routines: ',I5,' out of ',I5, &
           ' error message recorded')
!
   RETURN
!
!     End of DDRVRF2
!
END
                                                                                                                                                                                                                                                                                                            




