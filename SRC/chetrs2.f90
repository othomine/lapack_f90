!> \brief \b CHETRS2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
!                           WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRS2 solves a system of linear equations A*X = B with a complex
!> Hermitian matrix A using the factorization A = U*D*U**H or
!> A = L*D*L**H computed by CHETRF and converted by CSYCONV.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**H;
!>          = 'L':  Lower triangular, form is A = L*D*L**H.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by CHETRF.
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
!>          Details of the interchanges and the block structure of D
!>          as determined by CHETRF.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
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
!> \ingroup hetrs2
!
!  =====================================================================
   SUBROUTINE CHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I, IINFO, J, K, KP
   REAL               S
   COMPLEX            AK, AKM1, AKM1K, BK, BKM1, DENOM, B_TMP( NRHS )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CSYCONV, CTRSM, XERBLA
!     ..
!     .. Executable Statements ..
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( NRHS < 0 ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -5
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -8
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CHETRS2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 .OR. NRHS == 0 ) RETURN
!
!     Convert A
!
   CALL CSYCONV( UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO )
!
   IF( UPPER ) THEN
!
!        Solve A*X = B, where A = U*D*U**H.
!
!       P**T * B
     K=N
     DO WHILE ( K  >=  1 )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_TMP(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K-1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
         KP = -IPIV( K )
         IF( KP == -IPIV( K-1 ) ) THEN
            B_TMP(1:NRHS) = B(K-1,1:NRHS)
            B(K-1,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K-2
      END IF
     END DO
!
!  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
!
     CALL CTRSM('L','U','N','U',N,NRHS,1.0E+0,A,LDA,B,LDB)
!
!  Compute D \ B -> B   [ D \ (U \P**T * B) ]
!
      I=N
      DO WHILE ( I  >=  1 )
         IF( IPIV(I)  >  0 ) THEN
           S = REAL( 1.0E+0 ) / REAL( A( I, I ) )
           B(I,1:NRHS) = S*B(I,1:NRHS)
         ELSEIF ( I  >  1) THEN
            IF ( IPIV(I-1)  ==  IPIV(I) ) THEN
               AKM1K = WORK(I)
               AKM1 = A( I-1, I-1 ) / AKM1K
               AK = A( I, I ) / CONJG( AKM1K )
               DENOM = AKM1*AK - 1.0E+0
               DO J = 1, NRHS
                  BKM1 = B( I-1, J ) / AKM1K
                  BK = B( I, J ) / CONJG( AKM1K )
                  B( I-1, J ) = ( AK*BKM1-BK ) / DENOM
                  B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM
                 ENDDO
            I = I - 1
            ENDIF
         ENDIF
         I = I - 1
      END DO
!
!      Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ]
!
      CALL CTRSM('L','U','C','U',N,NRHS,1.0E+0,A,LDA,B,LDB)
!
!       P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ]
!
     K=1
     DO WHILE ( K  <=  N )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_TMP(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K+1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
         KP = -IPIV( K )
         IF( K  <  N .AND. KP == -IPIV( K+1 ) ) THEN
            B_TMP(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K+2
      ENDIF
     END DO
!
   ELSE
!
!        Solve A*X = B, where A = L*D*L**H.
!
!       P**T * B
     K=1
     DO WHILE ( K  <=  N )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_TMP(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K+1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K and -IPIV(K+1).
         KP = -IPIV( K+1 )
         IF( KP == -IPIV( K ) ) THEN
            B_TMP(1:NRHS) = B(K+1,1:NRHS)
            B(K+1,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K+2
      ENDIF
     END DO
!
!  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
!
     CALL CTRSM('L','L','N','U',N,NRHS,1.0E+0,A,LDA,B,LDB)
!
!  Compute D \ B -> B   [ D \ (L \P**T * B) ]
!
      I=1
      DO WHILE ( I  <=  N )
         IF( IPIV(I)  >  0 ) THEN
           S = REAL( 1.0E+0 ) / REAL( A( I, I ) )
           B(I,1:NRHS) = S*B(I,1:NRHS)
         ELSE
               AKM1K = WORK(I)
               AKM1 = A( I, I ) / CONJG( AKM1K )
               AK = A( I+1, I+1 ) / AKM1K
               DENOM = AKM1*AK - 1.0E+0
               DO J = 1, NRHS
                  BKM1 = B( I, J ) / CONJG( AKM1K )
                  BK = B( I+1, J ) / AKM1K
                  B( I, J ) = ( AK*BKM1-BK ) / DENOM
                  B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
                 ENDDO
               I = I + 1
         ENDIF
         I = I + 1
      END DO
!
!  Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ]
!
     CALL CTRSM('L','L','C','U',N,NRHS,1.0E+0,A,LDA,B,LDB)
!
!       P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ]
!
     K=N
     DO WHILE ( K  >=  1 )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_TMP(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K-1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
         KP = -IPIV( K )
         IF( K > 1 .AND. KP == -IPIV( K-1 ) ) THEN
            B_TMP(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_TMP(1:NRHS)
         ENDIF
         K=K-2
      ENDIF
     END DO
!
   END IF
!
!     Revert A
!
   CALL CSYCONV( UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO )
!
   RETURN
!
!     End of CHETRS2
!
END

