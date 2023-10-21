!> \brief \b CGETRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGETRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGETRF computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the right-looking Level 3 BLAS version of the algorithm.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix to be factored.
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and division by zero will occur if it is used
!>                to solve a system of equations.
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
!> \ingroup getrf
!
!  =====================================================================
   SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
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
   COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CGETRF2, CLASWP, CTRSM, XERBLA
!     ..
!     .. External Functions ..
   INTEGER            ILAENV
   EXTERNAL           ILAENV
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
      CALL XERBLA( 'CGETRF', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M == 0 .OR. N == 0 ) RETURN
!
!     Determine the block size for this environment.
!
   NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 )
   IF( NB <= 1 .OR. NB >= MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
      CALL CGETRF2( M, N, A, LDA, IPIV, INFO )
   ELSE
!
!        Use blocked code.
!
      DO J = 1, MIN( M, N ), NB
         JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
         CALL CGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
         IF( INFO == 0 .AND. IINFO > 0 ) &
            INFO = IINFO + J - 1
         DO I = J, MIN( M, J+JB-1 )
            IPIV( I ) = J - 1 + IPIV( I )
         ENDDO
!
!           Apply interchanges to columns 1:J-1.
!
         CALL CLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
         IF( J+JB <= N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
            CALL CLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                         IPIV, 1 )
!
!              Compute block row of U.
!
            CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                        N-J-JB+1, (1.0E+0,0.0E+0), A( J, J ), LDA, A( J, J+JB ), &
                        LDA )
            IF( J+JB <= M ) THEN
!
!                 Update trailing submatrix.
!
               CALL CGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                           N-J-JB+1, JB, -(1.0E+0,0.0E+0), A( J+JB, J ), LDA, &
                           A( J, J+JB ), LDA, (1.0E+0,0.0E+0), A( J+JB, J+JB ), &
                           LDA )
            END IF
         END IF
      ENDDO
   END IF
   RETURN
!
!     End of CGETRF
!
END
