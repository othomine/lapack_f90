!> \brief \b CTRTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTRTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrtri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrtri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrtri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
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
!> CTRTRI computes the inverse of a complex upper or lower triangular
!> matrix A.
!>
!> This is the Level 3 BLAS version of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          = 'N':  A is non-unit triangular;
!>          = 'U':  A is unit triangular.
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
!>          On entry, the triangular matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of the array A contains
!>          the upper triangular matrix, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of the array A contains
!>          the lower triangular matrix, and the strictly upper
!>          triangular part of A is not referenced.  If DIAG = 'U', the
!>          diagonal elements of A are also not referenced and are
!>          assumed to be 1.
!>          On exit, the (triangular) inverse of the original matrix, in
!>          the same storage format.
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
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!>               matrix is singular and its inverse can not be computed.
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
!> \ingroup trtri
!
!  =====================================================================
   SUBROUTINE CTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, UPLO
   INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            NOUNIT, UPPER
   INTEGER            J, JB, NB, NN
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
   EXTERNAL           CTRMM, CTRSM, CTRTI2, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   NOUNIT = LSAME( DIAG, 'N' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -5
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTRTRI', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Check for singularity if non-unit.
!
   IF( NOUNIT ) THEN
      DO INFO = 1, N
         IF( A( INFO, INFO ) == (0.0E+0,0.0E+0) ) RETURN
      ENDDO
      INFO = 0
   END IF
!
!     Determine the block size for this environment.
!
   NB = ILAENV( 1, 'CTRTRI', UPLO // DIAG, N, -1, -1, -1 )
   IF( NB <= 1 .OR. NB >= N ) THEN
!
!        Use unblocked code
!
      CALL CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
   ELSE
!
!        Use blocked code
!
      IF( UPPER ) THEN
!
!           Compute inverse of upper triangular matrix
!
         DO J = 1, N, NB
            JB = MIN( NB, N-J+1 )
!
!              Compute rows 1:j-1 of current block column
!
            CALL CTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1, &
                        JB, (1.0E+0,0.0E+0), A, LDA, A( 1, J ), LDA )
            CALL CTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1, &
                        JB, -(1.0E+0,0.0E+0), A( J, J ), LDA, A( 1, J ), LDA )
!
!              Compute inverse of current diagonal block
!
            CALL CTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
         ENDDO
      ELSE
!
!           Compute inverse of lower triangular matrix
!
         NN = ( ( N-1 ) / NB )*NB + 1
         DO J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
            IF( J+JB <= N ) THEN
!
!                 Compute rows j+jb:n of current block column
!
               CALL CTRMM( 'Left', 'Lower', 'No transpose', DIAG, &
                           N-J-JB+1, JB, (1.0E+0,0.0E+0), A( J+JB, J+JB ), LDA, &
                           A( J+JB, J ), LDA )
               CALL CTRSM( 'Right', 'Lower', 'No transpose', DIAG, &
                           N-J-JB+1, JB, -(1.0E+0,0.0E+0), A( J, J ), LDA, &
                           A( J+JB, J ), LDA )
            END IF
!
!              Compute inverse of current diagonal block
!
            CALL CTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
         ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of CTRTRI
!
END
