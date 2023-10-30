!> \brief \b CHEGS2 reduces a Hermitian definite generalized eigenproblem to standard form, using the factorization results obtained from cpotrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHEGS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
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
!> CHEGS2 reduces a complex Hermitian-definite generalized
!> eigenproblem to standard form.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H *A*L.
!>
!> B must have been previously factorized as U**H *U or L*L**H by ZPOTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
!>          = 2 or 3: compute U*A*U**H or L**H *A*L.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored, and how B has been factorized.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
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
!>          n by n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of A contains the lower
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
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup hegs2
!
!  =====================================================================
   SUBROUTINE CHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
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
   INTEGER            K
   REAL               AKK, BKK
   COMPLEX            CT
!     ..
!     .. External Subroutines ..
   EXTERNAL           CAXPY, CHER2, CTRMV, CTRSV, &
                      XERBLA
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
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
      CALL XERBLA( 'CHEGS2', -INFO )
      RETURN
   END IF
!
   IF( ITYPE == 1 ) THEN
      IF( UPPER ) THEN
!
!           Compute inv(U**H)*A*inv(U)
!
         DO K = 1, N
!
!              Update the upper triangle of A(k:n,k:n)
!
            AKK = REAL( A( K, K ) )
            BKK = REAL( B( K, K ) )
            AKK = AKK / BKK**2
            A( K, K ) = AKK
            IF( K < N ) THEN
               A(K,K+1:N) = A(K,K+1:N)/BKK
               CT = -0.5E+0*AKK
               A(K,K+1:N) = CONJG(A(K,K+1:N))
               B(K,K+1:N) = CONJG(B(K,K+1:N))
               A(K,K+1:N) = A(K,K+1:N) + CT*B(K,K+1:N)
               CALL CHER2( UPLO, N-K, -(1.0E+0,0.0E+0), A( K, K+1 ), LDA, &
                           B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )
               A(K,K+1:N) = A(K,K+1:N) + CT*B(K,K+1:N)
               B(K,K+1:N) = CONJG(B(K,K+1:N))
               CALL CTRSV( UPLO, 'Conjugate transpose', 'Non-unit', &
                           N-K, B( K+1, K+1 ), LDB, A( K, K+1 ), &
                           LDA )
               A(K,K+1:N) = CONJG(A(K,K+1:N))
            END IF
         ENDDO
      ELSE
!
!           Compute inv(L)*A*inv(L**H)
!
         DO K = 1, N
!
!              Update the lower triangle of A(k:n,k:n)
!
            AKK = REAL( A( K, K ) )
            BKK = REAL( B( K, K ) )
            AKK = AKK / BKK**2
            A( K, K ) = AKK
            IF( K < N ) THEN
               A(K+1:N,K) = A(K+1:N,K)/BKK
               CT = -0.5E+0*AKK
               A(K+1:N,K) = A(K+1:N,K) + CT*B(K+1:N,K)
               CALL CHER2( UPLO, N-K, -(1.0E+0,0.0E+0), A( K+1, K ), 1, &
                           B( K+1, K ), 1, A( K+1, K+1 ), LDA )
               A(K+1:N,K) = A(K+1:N,K) + CT*B(K+1:N,K)
               CALL CTRSV( UPLO, 'No transpose', 'Non-unit', N-K, &
                           B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
            END IF
         ENDDO
      END IF
   ELSE
      IF( UPPER ) THEN
!
!           Compute U*A*U**H
!
         DO K = 1, N
!
!              Update the upper triangle of A(1:k,1:k)
!
            AKK = REAL( A( K, K ) )
            BKK = REAL( B( K, K ) )
            CALL CTRMV( UPLO, 'No transpose', 'Non-unit', K-1, B, &
                        LDB, A( 1, K ), 1 )
            CT = 0.5E+0*AKK
            A(1:K-1,K) = A(1:K-1,K) + CT*B(1:K-1,K)
            CALL CHER2( UPLO, K-1, (1.0E+0,0.0E+0), A( 1, K ), 1, B( 1, K ), 1, &
                        A, LDA )
            A(1:K-1,K) = BKK*(A(1:K-1,K) + CT*B(1:K-1,K))
            A( K, K ) = AKK*BKK**2
         ENDDO
      ELSE
!
!           Compute L**H *A*L
!
         DO K = 1, N
!
!              Update the lower triangle of A(1:k,1:k)
!
            AKK = REAL( A( K, K ) )
            BKK = REAL( B( K, K ) )
            A(K,1:K-1) = CONJG(A(K,1:K-1))
            CALL CTRMV( UPLO, 'Conjugate transpose', 'Non-unit', K-1, &
                        B, LDB, A( K, 1 ), LDA )
            CT = 0.5E+0*AKK
            B(K,1:K-1) = CONJG(B(K,1:K-1))
            A(K,1:K-1) = A(K,1:K-1) + CT*B(K,1:K-1)
            CALL CHER2( UPLO, K-1, (1.0E+0,0.0E+0), A( K, 1 ), LDA, B( K, 1 ), &
                        LDB, A, LDA )
            A(K,1:K-1) = A(K,1:K-1) + CT*B(K,1:K-1)
            B(K,1:K-1) = CONJG(B(K,1:K-1))
            A(K,1:K-1) = BKK*A(K,1:K-1)
            A(K,1:K-1) = CONJG(A(K,1:K-1))
            A( K, K ) = AKK*BKK**2
         ENDDO
      END IF
   END IF
   RETURN
!
!     End of CHEGS2
!
END
