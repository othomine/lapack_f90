!> \brief \b CUNG2L generates all or part of the unitary matrix Q from a QL factorization determined by cgeqlf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNG2L + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cung2l.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cung2l.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cung2l.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
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
!> CUNG2L generates an m by n complex matrix Q with orthonormal columns,
!> which is defined as the last n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by CGEQLF.
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
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the (n-k+i)-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by CGEQLF in the last k columns of its array
!>          argument A.
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
!>          reflector H(i), as returned by CGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
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
!> \ingroup ung2l
!
!  =====================================================================
   SUBROUTINE CUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
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
   EXTERNAL           CLARF, CSCAL, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 .OR. N > M ) THEN
      INFO = -2
   ELSE IF( K < 0 .OR. K > N ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -5
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CUNG2L', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N <= 0 ) RETURN
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
   A(1:M,1:N-K) = (0.0E+0,0.0E+0)
   DO J = 1, N - K
      A(M-N+J,J) = (1.0E+0,0.0E+0)
   ENDDO
!
   DO I = 1, K
      II = N - K + I
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
      A( M-N+II, II ) = (1.0E+0,0.0E+0)
      CALL CLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, &
                  LDA, WORK )
      A(1:M-N+II-1,II) = -TAU( I )*A(1:M-N+II-1,II)
      A( M-N+II, II ) = (1.0E+0,0.0E+0) - TAU( I )
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
      A(M-N+II+1:M,II) = (0.0E+0,0.0E+0)
   ENDDO
   RETURN
!
!     End of CUNG2L
!
END
