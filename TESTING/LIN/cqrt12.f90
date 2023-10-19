!> \brief \b CQRT12
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CQRT12( M, N, A, LDA, S, WORK, LWORK,
!                        RWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * ), S( * )
!       COMPLEX            A( LDA, * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CQRT12 computes the singular values `svlues' of the upper trapezoid
!> of A(1:M,1:N) and returns the ratio
!>
!>      || s - svlues||/(||svlues||*eps*max(M,N))
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The M-by-N matrix A. Only the upper trapezoid is referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (min(M,N))
!>          The singular values of the matrix A.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK. LWORK >= M*N + 2*min(M,N) +
!>          max(M,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (4*min(M,N))
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
!> \ingroup complex_lin
!
!  =====================================================================
   REAL             FUNCTION CQRT12( M, N, A, LDA, S, WORK, LWORK, &
                    RWORK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   REAL               RWORK( * ), S( * )
   COMPLEX            A( LDA, * ), WORK( LWORK )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, ISCL, J, MN
   REAL               ANRM, BIGNUM, NRMSVL, SMLNUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   REAL               DUMMY( 1 )
!     ..
!     .. External Functions ..
   REAL               CLANGE, SASUM, SLAMCH, SNRM2
   EXTERNAL           CLANGE, SASUM, SLAMCH, SNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEBD2, CLASCL, CLASET, SAXPY, SBDSQR, SLASCL, &
                      XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          CMPLX, MAX, MIN, REAL
!     ..
!     .. Executable Statements ..
!
   CQRT12 = ZERO
!
!     Test that enough workspace is supplied
!
   IF( LWORK < M*N+2*MIN( M, N )+MAX( M, N ) ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'CQRT12', 7 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      RETURN
   END IF
!
!     Quick return if possible
!
   MN = MIN( M, N )
   IF( MN <= ZERO ) &
      RETURN
!
   NRMSVL = SNRM2( MN, S, 1 )
!
!     Copy upper triangle of A into work
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), WORK, M )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DO J = 1, N
      DO I = 1, MIN( J, M )
         WORK( ( J-1 )*M+I ) = A( I, J )
      ENDDO
   ENDDO
!
!     Get machine parameters
!
   SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
   BIGNUM = ONE / SMLNUM
!
!     Scale work if max entry outside range [SMLNUM,BIGNUM]
!
   ANRM = CLANGE( 'M', M, N, WORK, M, DUMMY )
   ISCL = 0
   IF( ANRM > ZERO .AND. ANRM < SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CLASCL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      ISCL = 1
   ELSE IF( ANRM > BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CLASCL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      ISCL = 1
   END IF
!
   IF( ANRM /= ZERO ) THEN
!
!        Compute SVD of work
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL CGEBD2( M, N, WORK, M, RWORK( 1 ), RWORK( MN+1 ), &
                   WORK( M*N+1 ), WORK( M*N+MN+1 ), &
                   WORK( M*N+2*MN+1 ), INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : CGEBD2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SBDSQR( 'Upper', MN, 0, 0, 0, RWORK( 1 ), RWORK( MN+1 ), &
                   DUMMY, MN, DUMMY, 1, DUMMY, MN, RWORK( 2*MN+1 ), &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SBDSQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      IF( ISCL == 1 ) THEN
         IF( ANRM > BIGNUM ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MN, 1, RWORK( 1 ), &
                         MN, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         END IF
         IF( ANRM < SMLNUM ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MN, 1, RWORK( 1 ), &
                         MN, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         END IF
      END IF
!
   ELSE
!
      DO I = 1, MN
         RWORK( I ) = ZERO
      ENDDO
   END IF
!
!     Compare s and singular values of work
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SAXPY( MN, -ONE, S, 1, RWORK( 1 ), 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SAXPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   CQRT12 = SASUM( MN, RWORK( 1 ), 1 ) / &
            ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )
   IF( NRMSVL /= ZERO ) &
      CQRT12 = CQRT12 / NRMSVL
!
   RETURN
!
!     End of CQRT12
!
END
                                                                                                                                                                                                                                                                                                            



