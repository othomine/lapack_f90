!> \brief \b CUNGR2 generates all or part of the unitary matrix Q from an RQ factorization determined by cgerqf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNGR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
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
!> CUNGR2 generates an m by n complex matrix Q with orthonormal rows,
!> which is defined as the last m rows of a product of k elementary
!> reflectors of order n
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
!>          On exit, the m-by-n matrix Q.
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
!>          WORK is COMPLEX array, dimension (M)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
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
!> \ingroup ungr2
!
!  =====================================================================
   SUBROUTINE CUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, II, J, L
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLARF, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < M ) THEN
      INFO = -2
   ELSE IF( K < 0 .OR. K > M ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -5
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CUNGR2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M <= 0 ) RETURN
!
   IF( K < M ) THEN
!
!        Initialise rows 1:m-k to rows of the unit matrix
!
      DO J = 1, N
         A(1:M-K, J ) = (0.0E+0,0.0E+0)
         IF( J > N-M .AND. J <= N-K ) A( M-N+J, J ) = (1.0E+0,0.0E+0)
      ENDDO
   END IF
!
   DO I = 1, K
      II = M - K + I
!
!        Apply H(i)**H to A(1:m-k+i,1:n-k+i) from the right
!
      A(II,1:N-M+II-1) = CONJG(A(II,1:N-M+II-1))
      A(II,N-M+II) = (1.0E+0,0.0E+0)
      CALL CLARF( 'Right', II-1, N-M+II, A( II, 1 ), LDA, &
                  CONJG( TAU( I ) ), A, LDA, WORK )
      A(II,1:N-M+II-1) = CONJG(-TAU( I )*A(II,1:N-M+II-1))
      A( II, N-M+II ) = (1.0E+0,0.0E+0) - CONJG( TAU( I ) )
!
!        Set A(m-k+i,n-k+i+1:n) to zero
!
      A( II,N-M+II+1:N) = (0.0E+0,0.0E+0)
   ENDDO
   RETURN
!
!     End of CUNGR2
!
END
