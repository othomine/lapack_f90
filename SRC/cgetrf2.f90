!> \brief \b CGETRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE CGETRF2( M, N, A, LDA, IPIV, INFO )
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
!> CGETRF2 computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the recursive version of the algorithm. It divides
!> the matrix into four submatrices:
!>
!>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
!>    A = [ -----|----- ]  with n1 = min(m,n)/2
!>        [  A21 | A22  ]       n2 = n-n1
!>
!>                                       [ A11 ]
!> The subroutine calls itself to factor [ --- ],
!>                                       [ A12 ]
!>                 [ A12 ]
!> do the swaps on [ --- ], solve A12, update A22,
!>                 [ A22 ]
!>
!> then calls itself to factor A22 and do the swaps on A21.
!>
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
!> \ingroup getrf2
!
!  =====================================================================
   RECURSIVE SUBROUTINE CGETRF2( M, N, A, LDA, IPIV, INFO )
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
   REAL               SFMIN
   COMPLEX            TEMP
   INTEGER            I, IINFO, N1, N2
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   INTEGER            ICAMAX
   EXTERNAL           SLAMCH, ICAMAX
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CLASWP, CTRSM, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
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
      CALL XERBLA( 'CGETRF2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M == 0 .OR. N == 0 ) RETURN

   IF ( M == 1 ) THEN
!
!        Use unblocked code for one row case
!        Just need to handle IPIV and INFO
!
      IPIV( 1 ) = 1
      IF ( A(1,1) == (0.0E+0,0.0E+0) ) INFO = 1
!
   ELSE IF( N == 1 ) THEN
!
!        Use unblocked code for one column case
!
!
!        Compute machine safe minimum
!
      SFMIN = SLAMCH('S')
!
!        Find pivot and test for singularity
!
      I = maxloc(ABS(REAL(A(1:M,1))) + ABS(AIMAG(A(1:M,1))),1)
      IPIV( 1 ) = I
      IF( A( I, 1 ) /= (0.0E+0,0.0E+0) ) THEN
!
!           Apply the interchange
!
         IF( I /= 1 ) THEN
            TEMP = A( 1, 1 )
            A( 1, 1 ) = A( I, 1 )
            A( I, 1 ) = TEMP
         END IF
!
!           Compute elements 2:M of the column
!
         IF( ABS(A( 1, 1 ))  >=  SFMIN ) THEN
            A(2:M,1) = A(2:M,1)/A(1,1)
         ELSE
            A(2:M,1) = A(2:M,1)/A(1,1)
         END IF
!
      ELSE
         INFO = 1
      END IF
!
   ELSE
!
!        Use recursive code
!
      N1 = MIN( M, N ) / 2
      N2 = N-N1
!
!               [ A11 ]
!        Factor [ --- ]
!               [ A21 ]
!
      CALL CGETRF2( M, N1, A, LDA, IPIV, IINFO )

      IF ( INFO == 0 .AND. IINFO > 0 ) INFO = IINFO
!
!                              [ A12 ]
!        Apply interchanges to [ --- ]
!                              [ A22 ]
!
      CALL CLASWP( N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 )
!
!        Solve A12
!
      CALL CTRSM( 'L', 'L', 'N', 'U', N1, N2, (1.0E+0,0.0E+0), A, LDA, &
                  A( 1, N1+1 ), LDA )
!
!        Update A22
!
      CALL CGEMM( 'N', 'N', M-N1, N2, N1, -(1.0E+0,0.0E+0), A( N1+1, 1 ), LDA, &
                  A( 1, N1+1 ), LDA, (1.0E+0,0.0E+0), A( N1+1, N1+1 ), LDA )
!
!        Factor A22
!
      CALL CGETRF2( M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ), &
                    IINFO )
!
!        Adjust INFO and the pivot indices
!
      IF ( INFO == 0 .AND. IINFO > 0 ) INFO = IINFO + N1
      DO I = N1+1, MIN( M, N )
         IPIV( I ) = IPIV( I ) + N1
      ENDDO
!
!        Apply interchanges to A21
!
      CALL CLASWP( N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 )
!
   END IF
   RETURN
!
!     End of CGETRF2
!
END
