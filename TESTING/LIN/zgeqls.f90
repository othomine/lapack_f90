!> \brief \b ZGEQLS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEQLS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Solve the least squares problem
!>     min || A*X - B ||
!> using the QL factorization
!>     A = Q*L
!> computed by ZGEQLF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  M >= N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          Details of the QL factorization of the original matrix A as
!>          returned by ZGEQLF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= M.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N)
!>          Details of the orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the m-by-nrhs right hand side matrix B.
!>          On exit, the n-by-nrhs solution matrix X, stored in rows
!>          m-n+1:m.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= M.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK must be at least NRHS,
!>          and should be at least NRHS*NB, where NB is the block size
!>          for this environment.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
   SUBROUTINE ZGEQLS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, &
                      INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Array Arguments ..
   COMPLEX*16         A( LDA, * ), B( LDB, * ), TAU( * ), &
                      WORK( LWORK )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX*16         ONE
   PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZTRSM, ZUNMQL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 .OR. N > M ) THEN
      INFO = -2
   ELSE IF( NRHS < 0 ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -5
   ELSE IF( LDB < MAX( 1, M ) ) THEN
      INFO = -8
   ELSE IF( LWORK < 1 .OR. LWORK < NRHS .AND. M > 0 .AND. N > 0 ) &
             THEN
      INFO = -10
   END IF
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'ZGEQLS', -INFO )
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
   IF( N == 0 .OR. NRHS == 0 .OR. M == 0 ) &
      RETURN
!
!     B := Q' * B
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZUNMQL( 'Left', 'Conjugate transpose', M, NRHS, N, A, LDA, &
                TAU, B, LDB, WORK, LWORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZUNMQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Solve L*X = B(m-n+1:m,:)
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, &
               ONE, A( M-N+1, 1 ), LDA, B( M-N+1, 1 ), LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZTRSM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   RETURN
!
!     End of ZGEQLS
!
END
                                                                                                                                                                                                                                                                                                            




