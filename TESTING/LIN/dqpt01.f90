!> \brief \b DQPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DQPT01( M, N, K, A, AF, LDA, TAU, JPVT,
!                        WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQPT01 tests the QR-factorization with pivoting of a matrix A.  The
!> array AF contains the (possibly partial) QR-factorization of A, where
!> the upper triangle of AF(1:k,1:k) is a partial triangular factor,
!> the entries below the diagonal in the first k columns are the
!> Householder vectors, and the rest of AF contains a partially updated
!> matrix.
!>
!> This function returns ||A*P - Q*R||/(||norm(A)||*eps*M)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and AF.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and AF.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of columns of AF that have been reduced
!>          to upper triangular form.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          The original matrix A.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          The (possibly partial) output of DGEQPF.  The upper triangle
!>          of AF(1:k,1:k) is a partial triangular factor, the entries
!>          below the diagonal in the first k columns are the Householder
!>          vectors, and the rest of AF contains a partially updated
!>          matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A and AF.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          Details of the Householder transformations as returned by
!>          DGEQPF.
!> \endverbatim
!>
!> \param[in] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          Pivot information as returned by DGEQPF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= M*N+N.
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
   DOUBLE PRECISION FUNCTION DQPT01( M, N, K, A, AF, LDA, TAU, JPVT, &
                    WORK, LWORK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            JPVT( * )
   DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), TAU( * ), &
                      WORK( LWORK )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, J
   DOUBLE PRECISION   NORMA
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   RWORK( 1 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLANGE
   EXTERNAL           DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DAXPY, DCOPY, DORMQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Executable Statements ..
!
   DQPT01 = ZERO
!
!     Test if there is enough workspace
!
   IF( LWORK < M*N+N ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DQPT01', 10 )
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
   IF( M <= 0 .OR. N <= 0 ) &
      RETURN
!
   NORMA = DLANGE( 'One-norm', M, N, A, LDA, RWORK )
!
   DO J = 1, K
      DO I = 1, MIN( J, M )
         WORK( ( J-1 )*M+I ) = AF( I, J )
      ENDDO
      DO I = J + 1, M
         WORK( ( J-1 )*M+I ) = ZERO
      ENDDO
   ENDDO
   DO J = K + 1, N
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDDO
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DORMQR( 'Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, &
                M, WORK( M*N+1 ), LWORK-M*N, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DORMQR : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO J = 1, N
!
!        Compare i-th column of QR and jpvt(i)-th column of A
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DAXPY( M, -ONE, A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), &
                  1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DAXPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDDO
!
   DQPT01 = DLANGE( 'One-norm', M, N, WORK, M, RWORK ) / &
            ( DBLE( MAX( M, N ) )*DLAMCH( 'Epsilon' ) )
   IF( NORMA /= ZERO ) &
      DQPT01 = DQPT01 / NORMA
!
   RETURN
!
!     End of DQPT01
!
END
                                                                                                                                                                                                                                                                                                            




