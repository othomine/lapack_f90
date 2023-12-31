!> \brief \b DGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGETF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETF2 computes an LU factorization of a general m-by-n matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the right-looking Level 2 BLAS version of the algorithm.
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
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the m by n matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!>               has been completed, but the factor U is exactly
!>               singular, and division by zero will occur if it is used
!>               to solve a system of equations.
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
!> \ingroup getf2
!
!  =====================================================================
   SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION   SFMIN
   INTEGER            I, J, JP
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   INTEGER            IDAMAX
   EXTERNAL           DLAMCH, IDAMAX
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -4
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DGETF2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M == 0 .OR. N == 0 ) RETURN
!
!     Compute machine safe minimum
!
   SFMIN = DLAMCH('S')
!
   DO J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
      JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
      IPIV( J ) = JP
      IF( A( JP, J ) /= 0.0D0 ) THEN
!
!           Apply the interchange to columns 1:N.
!
         IF( JP /= J ) CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
         IF( J < M ) A(J+1:M,J) = A(J+1:M,J) / A( J, J )
!
      ELSE IF( INFO == 0 ) THEN
!
         INFO = J
      END IF
!
      IF( J < MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
         CALL DGER( M-J, N-J, -1.0D0, A( J+1, J ), 1, A( J, J+1 ), LDA, &
                    A( J+1, J+1 ), LDA )
      END IF
   ENDDO
   RETURN
!
!     End of DGETF2
!
END
