!> \brief \b CHEGST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHEGST + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegst.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegst.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegst.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHEGST reduces a complex Hermitian-definite generalized
!> eigenproblem to standard form.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
!>
!> B must have been previously factorized as U**H*U or L*L**H by CPOTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
!>          = 2 or 3: compute U*A*U**H or L**H*A*L.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored and B is factored as
!>                  U**H*U;
!>          = 'L':  Lower triangle of A is stored and B is factored as
!>                  L*L**H.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
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
!>          On exit, if INFO = 0, the transformed matrix, stored in the
!>          same format as A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          The triangular factor from the Cholesky factorization of B,
!>          as returned by CPOTRF.
!>          B is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup hegst
!
!  =====================================================================
   SUBROUTINE CHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            K, KB, NB
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHEGS2, CHEMM, CHER2K, CTRMM, CTRSM, XERBLA
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( ITYPE < 1 .OR. ITYPE > 3 ) THEN
      INFO = -1
   ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -5
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -7
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CHEGST', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Determine the block size for this environment.
!
   NB = ILAENV( 1, 'CHEGST', UPLO, N, -1, -1, -1 )
!
   IF( NB <= 1 .OR. NB >= N ) THEN
!
!        Use unblocked code
!
      CALL CHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
   ELSE
!
!        Use blocked code
!
      IF( ITYPE == 1 ) THEN
         IF( UPPER ) THEN
!
!              Compute inv(U**H)*A*inv(U)
!
            DO K = 1, N, NB
               KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(k:n,k:n)
!
               CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                            B( K, K ), LDB, INFO )
               IF( K+KB <= N ) THEN
                  CALL CTRSM( 'Left', UPLO, 'Conjugate transpose', &
                              'Non-unit', KB, N-K-KB+1, (1.0E+0,0.0E+0), &
                              B( K, K ), LDB, A( K, K+KB ), LDA )
                  CALL CHEMM( 'Left', UPLO, KB, N-K-KB+1, -(0.5E+0,0.0E+0), &
                              A( K, K ), LDA, B( K, K+KB ), LDB, &
                              (1.0E+0,0.0E+0), A( K, K+KB ), LDA )
                  CALL CHER2K( UPLO, 'Conjugate transpose', N-K-KB+1, &
                               KB, -(1.0E+0,0.0E+0), A( K, K+KB ), LDA, &
                               B( K, K+KB ), LDB, 1.0E+0, &
                               A( K+KB, K+KB ), LDA )
                  CALL CHEMM( 'Left', UPLO, KB, N-K-KB+1, -(0.5E+0,0.0E+0), &
                              A( K, K ), LDA, B( K, K+KB ), LDB, &
                              (1.0E+0,0.0E+0), A( K, K+KB ), LDA )
                  CALL CTRSM( 'Right', UPLO, 'No transpose', &
                              'Non-unit', KB, N-K-KB+1, (1.0E+0,0.0E+0), &
                              B( K+KB, K+KB ), LDB, A( K, K+KB ), &
                              LDA )
               END IF
            ENDDO
         ELSE
!
!              Compute inv(L)*A*inv(L**H)
!
            DO K = 1, N, NB
               KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(k:n,k:n)
!
               CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                            B( K, K ), LDB, INFO )
               IF( K+KB <= N ) THEN
                  CALL CTRSM( 'Right', UPLO, 'Conjugate transpose', &
                              'Non-unit', N-K-KB+1, KB, (1.0E+0,0.0E+0), &
                              B( K, K ), LDB, A( K+KB, K ), LDA )
                  CALL CHEMM( 'Right', UPLO, N-K-KB+1, KB, -(0.5E+0,0.0E+0), &
                              A( K, K ), LDA, B( K+KB, K ), LDB, &
                              (1.0E+0,0.0E+0), A( K+KB, K ), LDA )
                  CALL CHER2K( UPLO, 'No transpose', N-K-KB+1, KB, &
                               -(1.0E+0,0.0E+0), A( K+KB, K ), LDA, &
                               B( K+KB, K ), LDB, 1.0E+0, &
                               A( K+KB, K+KB ), LDA )
                  CALL CHEMM( 'Right', UPLO, N-K-KB+1, KB, -(0.5E+0,0.0E+0), &
                              A( K, K ), LDA, B( K+KB, K ), LDB, &
                              (1.0E+0,0.0E+0), A( K+KB, K ), LDA )
                  CALL CTRSM( 'Left', UPLO, 'No transpose', &
                              'Non-unit', N-K-KB+1, KB, (1.0E+0,0.0E+0), &
                              B( K+KB, K+KB ), LDB, A( K+KB, K ), &
                              LDA )
               END IF
            ENDDO
         END IF
      ELSE
         IF( UPPER ) THEN
!
!              Compute U*A*U**H
!
            DO K = 1, N, NB
               KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
!
               CALL CTRMM( 'Left', UPLO, 'No transpose', 'Non-unit', &
                           K-1, KB, (1.0E+0,0.0E+0), B, LDB, A( 1, K ), LDA )
               CALL CHEMM( 'Right', UPLO, K-1, KB, (0.5E+0,0.0E+0), A( K, K ), &
                           LDA, B( 1, K ), LDB, (1.0E+0,0.0E+0), A( 1, K ), &
                           LDA )
               CALL CHER2K( UPLO, 'No transpose', K-1, KB, (1.0E+0,0.0E+0), &
                            A( 1, K ), LDA, B( 1, K ), LDB, 1.0E+0, A, &
                            LDA )
               CALL CHEMM( 'Right', UPLO, K-1, KB, (0.5E+0,0.0E+0), A( K, K ), &
                           LDA, B( 1, K ), LDB, (1.0E+0,0.0E+0), A( 1, K ), &
                           LDA )
               CALL CTRMM( 'Right', UPLO, 'Conjugate transpose', &
                           'Non-unit', K-1, KB, (1.0E+0,0.0E+0), B( K, K ), LDB, &
                           A( 1, K ), LDA )
               CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                            B( K, K ), LDB, INFO )
            ENDDO
         ELSE
!
!              Compute L**H*A*L
!
            DO K = 1, N, NB
               KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
!
               CALL CTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', &
                           KB, K-1, (1.0E+0,0.0E+0), B, LDB, A( K, 1 ), LDA )
               CALL CHEMM( 'Left', UPLO, KB, K-1, (0.5E+0,0.0E+0), A( K, K ), &
                           LDA, B( K, 1 ), LDB, (1.0E+0,0.0E+0), A( K, 1 ), &
                           LDA )
               CALL CHER2K( UPLO, 'Conjugate transpose', K-1, KB, &
                            (1.0E+0,0.0E+0), A( K, 1 ), LDA, B( K, 1 ), LDB, &
                            1.0E+0, A, LDA )
               CALL CHEMM( 'Left', UPLO, KB, K-1, (0.5E+0,0.0E+0), A( K, K ), &
                           LDA, B( K, 1 ), LDB, (1.0E+0,0.0E+0), A( K, 1 ), &
                           LDA )
               CALL CTRMM( 'Left', UPLO, 'Conjugate transpose', &
                           'Non-unit', KB, K-1, (1.0E+0,0.0E+0), B( K, K ), LDB, &
                           A( K, 1 ), LDA )
               CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                            B( K, K ), LDB, INFO )
            ENDDO
         END IF
      END IF
   END IF
   RETURN
!
!     End of CHEGST
!
END
