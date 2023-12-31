!> \brief \b CUNGRQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNGRQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungrq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungrq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungrq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNGRQ generates an M-by-N complex matrix Q with orthonormal rows,
!> which is defined as the last M rows of a product of K elementary
!> reflectors of order N
!>
!>       Q  =  H(1)**H H(2)**H . . . H(k)**H
!>
!> as returned by CGERQF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. N >= M.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. M >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the (m-k+i)-th row must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by CGERQF in the last k rows of its array argument
!>          A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGERQF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,M).
!>          For optimum performance LWORK >= M*NB, where NB is the
!>          optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
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
!> \ingroup ungrq
!
!  =====================================================================
   SUBROUTINE CUNGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY
   INTEGER            I, IB, II, IINFO, IWS, J, KK, L, LDWORK, &
                      LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLARFB, CLARFT, CUNGR2, XERBLA
!     ..
!     .. External Functions ..
   INTEGER            ILAENV
   EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   LQUERY = ( LWORK == -1 )
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < M ) THEN
      INFO = -2
   ELSE IF( K < 0 .OR. K > M ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -5
   END IF
!
   IF( INFO == 0 ) THEN
      IF( M <= 0 ) THEN
         LWKOPT = 1
      ELSE
         NB = ILAENV( 1, 'CUNGRQ', ' ', M, N, K, -1 )
         LWKOPT = M*NB
      END IF
      WORK( 1 ) = LWKOPT
!
      IF( LWORK < MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CUNGRQ', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M <= 0 ) RETURN
!
   NBMIN = 2
   NX = 0
   IWS = M
   IF( NB > 1 .AND. NB < K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
      NX = MAX( 0, ILAENV( 3, 'CUNGRQ', ' ', M, N, K, -1 ) )
      IF( NX < K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
         LDWORK = M
         IWS = LDWORK*NB
         IF( LWORK < IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'CUNGRQ', ' ', M, N, K, -1 ) )
         END IF
      END IF
   END IF
!
   IF( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN
!
!        Use blocked code after the first block.
!        The last kk rows are handled by the block method.
!
      KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
!
!        Set A(1:m-kk,n-kk+1:n) to zero.
!
      A(1:M-KK,N-KK+1:N) = (0.0E+0,0.0E+0)
   ELSE
      KK = 0
   END IF
!
!     Use unblocked code for the first or only block.
!
   CALL CUNGR2( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
!
   IF( KK > 0 ) THEN
!
!        Use blocked code
!
      DO I = K - KK + 1, K, NB
         IB = MIN( NB, K-I+1 )
         II = M - K + I
         IF( II > 1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL CLARFT( 'Backward', 'Rowwise', N-K+I+IB-1, IB, &
                         A( II, 1 ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H**H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
!
            CALL CLARFB( 'Right', 'Conjugate transpose', 'Backward', &
                         'Rowwise', II-1, N-K+I+IB-1, IB, A( II, 1 ), &
                         LDA, WORK, LDWORK, A, LDA, WORK( IB+1 ), &
                         LDWORK )
         END IF
!
!           Apply H**H to columns 1:n-k+i+ib-1 of current block
!
         CALL CUNGR2( IB, N-K+I+IB-1, IB, A( II, 1 ), LDA, TAU( I ), &
                      WORK, IINFO )
!
!           Set columns n-k+i+ib:n of current block to zero
!
         A(II:II+IB-1,N-K+I+IB:N) = (0.0E+0,0.0E+0)
      ENDDO
   END IF
!
   WORK( 1 ) = IWS
   RETURN
!
!     End of CUNGRQ
!
END
