!> \brief \b ZQRT13
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, M, N, SCALE
!       DOUBLE PRECISION   NORMA
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZQRT13 generates a full-rank matrix that may be scaled to have large
!> or small norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SCALE
!> \verbatim
!>          SCALE is INTEGER
!>          SCALE = 1: normally scaled matrix
!>          SCALE = 2: matrix scaled up
!>          SCALE = 3: matrix scaled down
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of A.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[out] NORMA
!> \verbatim
!>          NORMA is DOUBLE PRECISION
!>          The one-norm of A.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is integer array, dimension (4)
!>          Seed for random number generator
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
   SUBROUTINE ZQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, M, N, SCALE
   DOUBLE PRECISION   NORMA
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE
   PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
   INTEGER            INFO, J
   DOUBLE PRECISION   BIGNUM, SMLNUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
   EXTERNAL           DLAMCH, DZASUM, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZLARNV, ZLASCL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, DCMPLX, SIGN
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   DUMMY( 1 )
!     ..
!     .. Executable Statements ..
!
   IF( M <= 0 .OR. N <= 0 ) &
      RETURN
!
!     benign matrix
!
   DO J = 1, N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLARNV( 2, ISEED, M, A( 1, J ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( J <= M ) THEN
         A( J, J ) = A( J, J ) + DCMPLX( SIGN( DZASUM( M, A( 1, J ), &
                     1 ), DBLE( A( J, J ) ) ) )
      END IF
   ENDDO
!
!     scaled versions
!
   IF( SCALE /= 1 ) THEN
      NORMA = ZLANGE( 'Max', M, N, A, LDA, DUMMY )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SMLNUM / DLAMCH( 'Epsilon' )
      BIGNUM = ONE / SMLNUM
!
      IF( SCALE == 2 ) THEN
!
!           matrix scaled up
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLASCL( 'General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, &
                      INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLASCL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ELSE IF( SCALE == 3 ) THEN
!
!           matrix scaled down
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLASCL( 'General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, &
                      INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLASCL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      END IF
   END IF
!
   NORMA = ZLANGE( 'One-norm', M, N, A, LDA, DUMMY )
   RETURN
!
!     End of ZQRT13
!
END
                                                                                                                                                                                                                                                                                                            




