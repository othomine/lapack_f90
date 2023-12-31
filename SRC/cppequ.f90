!> \brief \b CPPEQU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPPEQU + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppequ.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppequ.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppequ.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       REAL               AMAX, SCOND
!       ..
!       .. Array Arguments ..
!       REAL               S( * )
!       COMPLEX            AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPPEQU computes row and column scalings intended to equilibrate a
!> Hermitian positive definite matrix A in packed storage and reduce
!> its condition number (with respect to the two-norm).  S contains the
!> scale factors, S(i)=1/sqrt(A(i,i)), chosen so that the scaled matrix
!> B with elements B(i,j)=S(i)*A(i,j)*S(j) has ones on the diagonal.
!> This choice of S puts the condition number of B within a factor N of
!> the smallest possible condition number over all possible diagonal
!> scalings.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          The upper or lower triangle of the Hermitian matrix A, packed
!>          columnwise in a linear array.  The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          If INFO = 0, S contains the scale factors for A.
!> \endverbatim
!>
!> \param[out] SCOND
!> \verbatim
!>          SCOND is REAL
!>          If INFO = 0, S contains the ratio of the smallest S(i) to
!>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
!>          large nor too small, it is not worth scaling by S.
!> \endverbatim
!>
!> \param[out] AMAX
!> \verbatim
!>          AMAX is REAL
!>          Absolute value of largest matrix element.  If AMAX is very
!>          close to overflow or very close to underflow, the matrix
!>          should be scaled.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.
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
!> \ingroup ppequ
!
!  =====================================================================
   SUBROUTINE CPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, N
   REAL               AMAX, SCOND
!     ..
!     .. Array Arguments ..
   REAL               S( * )
   COMPLEX            AP( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I, JJ
   REAL               SMIN
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CPPEQU', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) THEN
      SCOND = 1.0E+0
      AMAX = 0.0E+0
      RETURN
   END IF
!
!     Initialize SMIN and AMAX.
!
   S( 1 ) = REAL( AP( 1 ) )
   SMIN = S( 1 )
   AMAX = S( 1 )
!
   IF( UPPER ) THEN
!
!        UPLO = 'U':  Upper triangle of A is stored.
!        Find the minimum and maximum diagonal elements.
!
      JJ = 1
      DO I = 2, N
         JJ = JJ + I
         S( I ) = REAL( AP( JJ ) )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
      ENDDO
!
   ELSE
!
!        UPLO = 'L':  Lower triangle of A is stored.
!        Find the minimum and maximum diagonal elements.
!
      JJ = 1
      DO I = 2, N
         JJ = JJ + N - I + 2
         S( I ) = REAL( AP( JJ ) )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
      ENDDO
   END IF
!
   IF( SMIN <= 0.0E+0 ) THEN
!
!        Find the first non-positive diagonal element and return.
!
      DO I = 1, N
         IF( S( I ) <= 0.0E+0 ) THEN
            INFO = I
            RETURN
         END IF
      ENDDO
   ELSE
!
!        Set the scale factors to the reciprocals
!        of the diagonal elements.
!
      S(1:N) = 1.0E+0 / SQRT( S(1:N) )
!
!        Compute SCOND = min(S(I)) / max(S(I))
!
      SCOND = SQRT( SMIN ) / SQRT( AMAX )
   END IF
   RETURN
!
!     End of CPPEQU
!
END
