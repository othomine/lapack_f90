!> \brief \b ZLAVSY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAVSY( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B,
!                          LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAVSY performs one of the matrix-vector operations
!>    x := A*x  or  x := A'*x,
!> where x is an N element vector and  A is one of the factors
!> from the block U*D*U' or L*D*L' factorization computed by ZSYTRF.
!>
!> If TRANS = 'N', multiplies by U  or U * D  (or L  or L * D)
!> If TRANS = 'T', multiplies by U' or D * U' (or L' or D * L')
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the factor stored in A is upper or lower
!>          triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation to be performed:
!>          = 'N':  x := A*x
!>          = 'T':  x := A'*x
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the diagonal blocks are unit
!>          matrices.  If the diagonal blocks are assumed to be unit,
!>          then A = U or A = L, otherwise A = U*D or A = L*D.
!>          = 'U':  Diagonal blocks are assumed to be unit matrices.
!>          = 'N':  Diagonal blocks are assumed to be non-unit matrices.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of vectors
!>          x to be multiplied by A.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by ZSYTRF.
!>          Stored as a 2-D triangular matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D,
!>          as determined by ZSYTRF.
!>
!>          If UPLO = 'U':
!>               If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>               were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>               (If IPIV( k ) = k, no interchange was done).
!>
!>               If IPIV(k) = IPIV(k-1) < 0, then rows and
!>               columns k-1 and -IPIV(k) were interchanged,
!>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>               If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>               were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>               (If IPIV( k ) = k, no interchange was done).
!>
!>               If IPIV(k) = IPIV(k+1) < 0, then rows and
!>               columns k+1 and -IPIV(k) were interchanged,
!>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, B contains NRHS vectors of length N.
!>          On exit, B is overwritten with the product A * B.
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
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_lin
!
!  =====================================================================
   SUBROUTINE ZLAVSY( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B, &
                      LDB, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, TRANS, UPLO
   INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX*16         CONE
   PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            NOUNIT
   INTEGER            J, K, KP
   COMPLEX*16         D11, D12, D21, D22, T1, T2
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZGEMV, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) ) &
             THEN
      INFO = -2
   ELSE IF( .NOT.LSAME( DIAG, 'U' ) .AND. .NOT.LSAME( DIAG, 'N' ) ) &
             THEN
      INFO = -3
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -6
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -9
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'ZLAVSY ', -INFO )
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF( N == 0 ) &
      RETURN
!
   NOUNIT = LSAME( DIAG, 'N' )
!------------------------------------------
!
!     Compute  B := A * B  (No transpose)
!
!------------------------------------------
   IF( LSAME( TRANS, 'N' ) ) THEN
!
!        Compute  B := U*B
!        where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Loop forward applying the transformations.
!
         K = 1
10       CONTINUE
         IF( K > N ) &
            GO TO 30
         IF( IPIV( K ) > 0 ) THEN
!
!              1 x 1 pivot block
!
!              Multiply by the diagonal element if forming U * D.
!
            IF( NOUNIT ) &
               CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
!
!              Multiply by  P(K) * inv(U(K))  if K > 1.
!
            IF( K > 1 ) THEN
!
!                 Apply the transformation.
!
               CALL ZGERU( K-1, NRHS, CONE, A( 1, K ), 1, B( K, 1 ), &
                           LDB, B( 1, 1 ), LDB )
!
!                 Interchange if P(K) != I.
!
               KP = IPIV( K )
               IF( KP /= K ) &
                  CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
            K = K + 1
         ELSE
!
!              2 x 2 pivot block
!
!              Multiply by the diagonal block if forming U * D.
!
            IF( NOUNIT ) THEN
               D11 = A( K, K )
               D22 = A( K+1, K+1 )
               D12 = A( K, K+1 )
               D21 = D12
               DO J = 1, NRHS
                  T1 = B( K, J )
                  T2 = B( K+1, J )
                  B( K, J ) = D11*T1 + D12*T2
                  B( K+1, J ) = D21*T1 + D22*T2
               ENDDO
            END IF
!
!              Multiply by  P(K) * inv(U(K))  if K > 1.
!
            IF( K > 1 ) THEN
!
!                 Apply the transformations.
!
               CALL ZGERU( K-1, NRHS, CONE, A( 1, K ), 1, B( K, 1 ), &
                           LDB, B( 1, 1 ), LDB )
               CALL ZGERU( K-1, NRHS, CONE, A( 1, K+1 ), 1, &
                           B( K+1, 1 ), LDB, B( 1, 1 ), LDB )
!
!                 Interchange if P(K) != I.
!
               KP = ABS( IPIV( K ) )
               IF( KP /= K ) &
                  CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
            K = K + 2
         END IF
         GO TO 10
30       CONTINUE
!
!        Compute  B := L*B
!        where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
!
      ELSE
!
!           Loop backward applying the transformations to B.
!
         K = N
40       CONTINUE
         IF( K < 1 ) &
            GO TO 60
!
!           Test the pivot index.  If greater than zero, a 1 x 1
!           pivot was used, otherwise a 2 x 2 pivot was used.
!
         IF( IPIV( K ) > 0 ) THEN
!
!              1 x 1 pivot block:
!
!              Multiply by the diagonal element if forming L * D.
!
            IF( NOUNIT ) &
               CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
!
!              Multiply by  P(K) * inv(L(K))  if K < N.
!
            IF( K /= N ) THEN
               KP = IPIV( K )
!
!                 Apply the transformation.
!
               CALL ZGERU( N-K, NRHS, CONE, A( K+1, K ), 1, &
                            B( K, 1 ), LDB, B( K+1, 1 ), LDB )
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
               IF( KP /= K ) &
                  CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
            K = K - 1
!
         ELSE
!
!              2 x 2 pivot block:
!
!              Multiply by the diagonal block if forming L * D.
!
            IF( NOUNIT ) THEN
               D11 = A( K-1, K-1 )
               D22 = A( K, K )
               D21 = A( K, K-1 )
               D12 = D21
               DO J = 1, NRHS
                  T1 = B( K-1, J )
                  T2 = B( K, J )
                  B( K-1, J ) = D11*T1 + D12*T2
                  B( K, J ) = D21*T1 + D22*T2
               ENDDO
            END IF
!
!              Multiply by  P(K) * inv(L(K))  if K < N.
!
            IF( K /= N ) THEN
!
!                 Apply the transformation.
!
               CALL ZGERU( N-K, NRHS, CONE, A( K+1, K ), 1, &
                           B( K, 1 ), LDB, B( K+1, 1 ), LDB )
               CALL ZGERU( N-K, NRHS, CONE, A( K+1, K-1 ), 1, &
                           B( K-1, 1 ), LDB, B( K+1, 1 ), LDB )
!
!                 Interchange if a permutation was applied at the
!                 K-th step of the factorization.
!
               KP = ABS( IPIV( K ) )
               IF( KP /= K ) &
                  CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
            K = K - 2
         END IF
         GO TO 40
60       CONTINUE
      END IF
!----------------------------------------
!
!     Compute  B := A' * B  (transpose)
!
!----------------------------------------
   ELSE IF( LSAME( TRANS, 'T' ) ) THEN
!
!        Form  B := U'*B
!        where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
!        and   U' = inv(U'(1))*P(1)* ... *inv(U'(m))*P(m)
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Loop backward applying the transformations.
!
         K = N
70       CONTINUE
         IF( K < 1 ) &
            GO TO 90
!
!           1 x 1 pivot block.
!
         IF( IPIV( K ) > 0 ) THEN
            IF( K > 1 ) THEN
!
!                 Interchange if P(K) != I.
!
               KP = IPIV( K )
               IF( KP /= K ) &
                  CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!                 Apply the transformation
!
               CALL ZGEMV( 'Transpose', K-1, NRHS, CONE, B, LDB, &
                           A( 1, K ), 1, CONE, B( K, 1 ), LDB )
            END IF
            IF( NOUNIT ) &
               CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
            K = K - 1
!
!           2 x 2 pivot block.
!
         ELSE
            IF( K > 2 ) THEN
!
!                 Interchange if P(K) != I.
!
               KP = ABS( IPIV( K ) )
               IF( KP /= K-1 ) &
                  CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), &
                              LDB )
!
!                 Apply the transformations
!
               CALL ZGEMV( 'Transpose', K-2, NRHS, CONE, B, LDB, &
                           A( 1, K ), 1, CONE, B( K, 1 ), LDB )
               CALL ZGEMV( 'Transpose', K-2, NRHS, CONE, B, LDB, &
                           A( 1, K-1 ), 1, CONE, B( K-1, 1 ), LDB )
            END IF
!
!              Multiply by the diagonal block if non-unit.
!
            IF( NOUNIT ) THEN
               D11 = A( K-1, K-1 )
               D22 = A( K, K )
               D12 = A( K-1, K )
               D21 = D12
               DO J = 1, NRHS
                  T1 = B( K-1, J )
                  T2 = B( K, J )
                  B( K-1, J ) = D11*T1 + D12*T2
                  B( K, J ) = D21*T1 + D22*T2
               ENDDO
            END IF
            K = K - 2
         END IF
         GO TO 70
90       CONTINUE
!
!        Form  B := L'*B
!        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
!        and   L' = inv(L'(m))*P(m)* ... *inv(L'(1))*P(1)
!
      ELSE
!
!           Loop forward applying the L-transformations.
!
         K = 1
  100       CONTINUE
         IF( K > N ) &
            GO TO 120
!
!           1 x 1 pivot block
!
         IF( IPIV( K ) > 0 ) THEN
            IF( K < N ) THEN
!
!                 Interchange if P(K) != I.
!
               KP = IPIV( K )
               IF( KP /= K ) &
                  CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!                 Apply the transformation
!
               CALL ZGEMV( 'Transpose', N-K, NRHS, CONE, B( K+1, 1 ), &
                          LDB, A( K+1, K ), 1, CONE, B( K, 1 ), LDB )
            END IF
            IF( NOUNIT ) &
               CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
            K = K + 1
!
!           2 x 2 pivot block.
!
         ELSE
            IF( K < N-1 ) THEN
!
!              Interchange if P(K) != I.
!
               KP = ABS( IPIV( K ) )
               IF( KP /= K+1 ) &
                  CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), &
                              LDB )
!
!                 Apply the transformation
!
               CALL ZGEMV( 'Transpose', N-K-1, NRHS, CONE, &
                           B( K+2, 1 ), LDB, A( K+2, K+1 ), 1, CONE, &
                           B( K+1, 1 ), LDB )
               CALL ZGEMV( 'Transpose', N-K-1, NRHS, CONE, &
                           B( K+2, 1 ), LDB, A( K+2, K ), 1, CONE, &
                           B( K, 1 ), LDB )
            END IF
!
!              Multiply by the diagonal block if non-unit.
!
            IF( NOUNIT ) THEN
               D11 = A( K, K )
               D22 = A( K+1, K+1 )
               D21 = A( K+1, K )
               D12 = D21
               DO J = 1, NRHS
                  T1 = B( K, J )
                  T2 = B( K+1, J )
                  B( K, J ) = D11*T1 + D12*T2
                  B( K+1, J ) = D21*T1 + D22*T2
                  ENDDO
            END IF
            K = K + 2
         END IF
         GO TO 100
  120       CONTINUE
      END IF
   END IF
   RETURN
!
!     End of ZLAVSY
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
