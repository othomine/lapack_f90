!> \brief \b CPOTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPOTF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPOTF2( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPOTF2 computes the Cholesky factorization of a complex Hermitian
!> positive definite matrix A.
!>
!> The factorization has the form
!>    A = U**H * U ,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          n by n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**H *U  or A = L*L**H.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, the leading principal minor of order k
!>               is not positive, and the factorization could not be
!>               completed.
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
!> \ingroup potf2
!
!  =====================================================================
   SUBROUTINE CPOTF2( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            J
   REAL               AJJ
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, SISNAN
   COMPLEX            CDOTC
   EXTERNAL           LSAME, CDOTC, SISNAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMV, XERBLA
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
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CPOTF2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
   IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U**H *U.
!
      DO J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
         AJJ = REAL( REAL( A( J, J ) ) - CDOTC( J-1, A( 1, J ), 1, &
               A( 1, J ), 1 ) )
         IF( AJJ <= 0.0E+0.OR.SISNAN( AJJ ) ) THEN
            A( J, J ) = AJJ
            GO TO 30
         END IF
         AJJ = SQRT( AJJ )
         A( J, J ) = AJJ
!
!           Compute elements J+1:N of row J.
!
         IF( J < N ) THEN
            A(1:J-1,J) = CONJG(A(1:J-1,J))
            CALL CGEMV( 'Transpose', J-1, N-J, -(1.0E+0,0.0E+0), A( 1, J+1 ), &
                        LDA, A( 1, J ), 1, (1.0E+0,0.0E+0), A( J, J+1 ), LDA )
            A(1:J-1,J) = CONJG(A(1:J-1,J))
            A(J,J+1:N) = A(J,J+1:N) / AJJ
         END IF
      ENDDO
   ELSE
!
!        Compute the Cholesky factorization A = L*L**H.
!
      DO J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
         AJJ = REAL( REAL( A( J, J ) ) - CDOTC( J-1, A( J, 1 ), LDA, &
               A( J, 1 ), LDA ) )
         IF( AJJ <= 0.0E+0.OR.SISNAN( AJJ ) ) THEN
            A( J, J ) = AJJ
            GO TO 30
         END IF
         AJJ = SQRT( AJJ )
         A( J, J ) = AJJ
!
!           Compute elements J+1:N of column J.
!
         IF( J < N ) THEN
            A(J,1:J-1) = CONJG(A(J,1:J-1))
            CALL CGEMV( 'No transpose', N-J, J-1, -(1.0E+0,0.0E+0), A( J+1, 1 ), &
                        LDA, A( J, 1 ), LDA, (1.0E+0,0.0E+0), A( J+1, J ), 1 )
            A(J,1:J-1) = CONJG(A(J,1:J-1))
            A(J+1:N,J) = A(J+1:N,J) / AJJ
         END IF
      ENDDO
   END IF
   GO TO 40
!
30 CONTINUE
   INFO = J
!
40 CONTINUE
   RETURN
!
!     End of CPOTF2
!
END
