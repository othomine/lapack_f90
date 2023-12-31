!> \brief \b CSYTRS2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYTRS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
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
!> CSYTRS2 solves a system of linear equations A*X = B with a complex
!> symmetric matrix A using the factorization A = U*D*U**T or
!> A = L*D*L**T computed by CSYTRF and converted by CSYCONV.
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
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by CSYTRF.
!>          Note that A is input / output. This might be counter-intuitive,
!>          and one may think that A is input only. A is input / output. This
!>          is because, at the start of the subroutine, we permute A in a
!>          "better" form and then we permute A back to its original form at
!>          the end.
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
!>          as determined by CSYTRF.
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
   SUBROUTINE CSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
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
!     .. Local Array ..
   COMPLEX            B_tmp(NRHS)
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I, IINFO, J, K, KP
   COMPLEX            AK, AKM1, AKM1K, BK, BKM1, DENOM
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
      CALL XERBLA( 'CSYTRS2', -INFO )
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
!        Solve A*X = B, where A = U*D*U**T.
!
!       P**T * B
     K=N
     DO WHILE ( K  >=  1 )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_tmp(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
         ENDIF
         K=K-1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
         KP = -IPIV( K )
         IF( KP == -IPIV( K-1 ) ) THEN
            B_tmp(1:NRHS) = B(K-1,1:NRHS)
            B(K-1,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
         ENDIF
         K=K-2
      END IF
     END DO
!
!  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
!
     CALL CTRSM('L','U','N','U',N,NRHS,(1.0E+0,0.0E+0),A,LDA,B,LDB)
!
!  Compute D \ B -> B   [ D \ (U \P**T * B) ]
!
      I=N
      DO WHILE ( I  >=  1 )
         IF( IPIV(I)  >  0 ) THEN
           B(I,1:NRHS) = B(I,1:NRHS) / A( I, I )
         ELSEIF ( I  >  1) THEN
            IF ( IPIV(I-1)  ==  IPIV(I) ) THEN
               AKM1K = (1.0E+0,0.0E+0) / WORK(I)
               AKM1 = A( I-1, I-1 ) * AKM1K
               AK = A( I, I ) * AKM1K
               DENOM = (1.0E+0,0.0E+0) / (AKM1*AK - (1.0E+0,0.0E+0 ))
               DO J = 1, NRHS
                  BKM1 = B( I-1, J ) * AKM1K
                  BK = B( I, J ) * AKM1K
                  B( I-1, J ) = ( AK*BKM1-BK ) * DENOM
                  B( I, J ) = ( AKM1*BK-BKM1 ) * DENOM
               ENDDO
            I = I - 1
            ENDIF
         ENDIF
         I = I - 1
      END DO
!
!      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
!
      CALL CTRSM('L','U','T','U',N,NRHS,(1.0E+0,0.0E+0 ),A,LDA,B,LDB)
!
!       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
!
     K=1
     DO WHILE ( K  <=  N )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_tmp(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
         ENDIF
         K=K+1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
         KP = -IPIV( K )
         IF( K  <  N .AND. KP == -IPIV( K+1 ) ) THEN
            B_tmp(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
         ENDIF
         K=K+2
      ENDIF
     END DO
!
   ELSE
!
!        Solve A*X = B, where A = L*D*L**T.
!
!       P**T * B
     K=1
     DO WHILE ( K  <=  N )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_tmp(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
         ENDIF
         K=K+1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K and -IPIV(K+1).
         KP = -IPIV( K+1 )
         IF( KP == -IPIV( K ) ) THEN
            B_tmp(1:NRHS) = B(K+1,1:NRHS)
            B(K+1,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
         ENDIF
         K=K+2
      ENDIF
     END DO
!
!  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
!
     CALL CTRSM('L','L','N','U',N,NRHS,(1.0E+0,0.0E+0 ),A,LDA,B,LDB)
!
!  Compute D \ B -> B   [ D \ (L \P**T * B) ]
!
      I=1
      DO WHILE ( I  <=  N )
         IF( IPIV(I)  >  0 ) THEN
           B(I,1:NRHS) = B(I,1:NRHS) / A( I, I )
         ELSE
               AKM1K = (1.0E+0,0.0E+0 ) / WORK(I)
               AKM1 = A( I, I ) * AKM1K
               AK = A( I+1, I+1 ) * AKM1K
               DENOM = (1.0E+0,0.0E+0 ) / (AKM1*AK - (1.0E+0,0.0E+0 ))
               DO J = 1, NRHS
                  BKM1 = B( I, J ) * AKM1K
                  BK = B( I+1, J ) * AKM1K
                  B( I, J ) = ( AK*BKM1-BK ) * DENOM
                  B( I+1, J ) = ( AKM1*BK-BKM1 ) * DENOM
               ENDDO
               I = I + 1
         ENDIF
         I = I + 1
      END DO
!
!  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
!
     CALL CTRSM('L','L','T','U',N,NRHS,(1.0E+0,0.0E+0 ),A,LDA,B,LDB)
!
!       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
!
     K=N
     DO WHILE ( K  >=  1 )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal block
!           Interchange rows K and IPIV(K).
         KP = IPIV( K )
         IF( KP /= K ) THEN
            B_tmp(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
         ENDIF
         K=K-1
      ELSE
!           2 x 2 diagonal block
!           Interchange rows K-1 and -IPIV(K).
         KP = -IPIV( K )
         IF( K > 1 .AND. KP == -IPIV( K-1 ) ) THEN
            B_tmp(1:NRHS) = B(K,1:NRHS)
            B(K,1:NRHS) = B(KP,1:NRHS)
            B(KP,1:NRHS) = B_tmp(1:NRHS)
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
!     End of CSYTRS2
!
END
