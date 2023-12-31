!> \brief \b CPOTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPOTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPOTRF( UPLO, N, A, LDA, INFO )
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
!> CPOTRF computes the Cholesky factorization of a complex Hermitian
!> positive definite matrix A.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the block version of the algorithm, calling Level 3 BLAS.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**H*U or A = L*L**H.
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading principal minor of order i
!>                is not positive, and the factorization could not be
!>                completed.
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
!> \ingroup potrf
!
!  =====================================================================
   SUBROUTINE CPOTRF( UPLO, N, A, LDA, INFO )
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
   INTEGER            J, JB, NB
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CHERK, CPOTRF2, CTRSM, XERBLA
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
      CALL XERBLA( 'CPOTRF', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Determine the block size for this environment.
!
   NB = ILAENV( 1, 'CPOTRF', UPLO, N, -1, -1, -1 )
   IF( NB <= 1 .OR. NB >= N ) THEN
!
!        Use unblocked code.
!
      CALL CPOTRF2( UPLO, N, A, LDA, INFO )
   ELSE
!
!        Use blocked code.
!
      IF( UPPER ) THEN
!
!           Compute the Cholesky factorization A = U**H *U.
!
         DO J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            JB = MIN( NB, N-J+1 )
            CALL CHERK( 'Upper', 'Conjugate transpose', JB, J-1, &
                        -1.0E+0, A( 1, J ), LDA, 1.0E+0, A( J, J ), LDA )
            CALL CPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
            IF( INFO /= 0 ) &
               GO TO 30
            IF( J+JB <= N ) THEN
!
!                 Compute the current block row.
!
               CALL CGEMM( 'Conjugate transpose', 'No transpose', JB, &
                           N-J-JB+1, J-1, -(1.0E+0,0.0E+0), A( 1, J ), LDA, &
                           A( 1, J+JB ), LDA, (1.0E+0,0.0E+0), A( J, J+JB ), &
                           LDA )
               CALL CTRSM( 'Left', 'Upper', 'Conjugate transpose', &
                           'Non-unit', JB, N-J-JB+1, (1.0E+0,0.0E+0), A( J, J ), &
                           LDA, A( J, J+JB ), LDA )
            END IF
         ENDDO
!
      ELSE
!
!           Compute the Cholesky factorization A = L*L**H.
!
         DO J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            JB = MIN( NB, N-J+1 )
            CALL CHERK( 'Lower', 'No transpose', JB, J-1, -1.0E+0, &
                        A( J, 1 ), LDA, 1.0E+0, A( J, J ), LDA )
            CALL CPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
            IF( INFO /= 0 ) &
               GO TO 30
            IF( J+JB <= N ) THEN
!
!                 Compute the current block column.
!
               CALL CGEMM( 'No transpose', 'Conjugate transpose', &
                           N-J-JB+1, JB, J-1, -(1.0E+0,0.0E+0), A( J+JB, 1 ), &
                           LDA, A( J, 1 ), LDA, (1.0E+0,0.0E+0), A( J+JB, J ), &
                           LDA )
               CALL CTRSM( 'Right', 'Lower', 'Conjugate transpose', &
                           'Non-unit', N-J-JB+1, JB, (1.0E+0,0.0E+0), A( J, J ), &
                           LDA, A( J+JB, J ), LDA )
            END IF
         ENDDO
      END IF
   END IF
   GO TO 40
!
30 CONTINUE
   INFO = INFO + J - 1
!
40 CONTINUE
   RETURN
!
!     End of CPOTRF
!
END
