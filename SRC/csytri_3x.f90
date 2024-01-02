!> \brief \b CSYTRI_3X
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYTRI_3X + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytri_3x.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytri_3x.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytri_3x.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYTRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ),  E( * ), WORK( N+NB+1, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> CSYTRI_3X computes the inverse of a complex symmetric indefinite
!> matrix A using the factorization computed by CSYTRF_RK or CSYTRF_BK:
!>
!>     A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**T (or L**T) is the transpose of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is symmetric and block
!> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are
!>          stored as an upper or lower triangular matrix.
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
!>          On entry, diagonal of the block diagonal matrix D and
!>          factors U or L as computed by CSYTRF_RK and CSYTRF_BK:
!>            a) ONLY diagonal elements of the symmetric block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                should be provided on entry in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!>
!>          On exit, if INFO = 0, the symmetric inverse of the original
!>          matrix.
!>             If UPLO = 'U': the upper triangular part of the inverse
!>             is formed and the part of A below the diagonal is not
!>             referenced;
!>             If UPLO = 'L': the lower triangular part of the inverse
!>             is formed and the part of A above the diagonal is not
!>             referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX array, dimension (N)
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the symmetric block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) not referenced;
!>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) not referenced.
!>
!>          NOTE: For 1-by-1 diagonal block D(k), where
!>          1 <= k <= N, the element E(k) is not referenced in both
!>          UPLO = 'U' or UPLO = 'L' cases.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by CSYTRF_RK or CSYTRF_BK.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N+NB+1,NB+3).
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Block size.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
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
!> \ingroup hetri_3x
!
!> \par Contributors:
!  ==================
!> \verbatim
!>
!>  June 2017,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE CSYTRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, N, NB
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   COMPLEX            A( LDA, * ), E( * ), WORK( N+NB+1, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            CUT, I, ICOUNT, INVD, IP, K, NNB, J, U11
   COMPLEX            AK, AKKP1, AKP1, UoD, UoT, U01_I_J, U01_IP1_J, &
                      U11_I_J, U11_IP1_J
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CSYSWAPR, CTRTRI, CTRMM, XERBLA
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
!
!     Quick return if possible
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSYTRI_3X', -INFO )
      RETURN
   END IF
   IF( N == 0 ) RETURN
!
!     Workspace got Non-diag elements of D
!
   WORK(1:N, 1 ) = E(1:N)
!
!     Check that the diagonal matrix D is nonsingular.
!
   IF( UPPER ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
      DO INFO = N, 1, -1
         IF( IPIV( INFO ) > 0 .AND. A( INFO, INFO ) == (0.0E+0,0.0E+0) ) RETURN
      END DO
   ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
      DO INFO = 1, N
         IF( IPIV( INFO ) > 0 .AND. A( INFO, INFO ) == (0.0E+0,0.0E+0) ) RETURN
      END DO
   END IF
!
   INFO = 0
!
!     Splitting Workspace
!     U01 is a block ( N, NB+1 )
!     The first element of U01 is in WORK( 1, 1 )
!     U11 is a block ( NB+1, NB+1 )
!     The first element of U11 is in WORK( N+1, 1 )
!
   U11 = N
!
!     INVD is a block ( N, 2 )
!     The first element of INVD is in WORK( 1, INVD )
!
   INVD = NB + 2

   IF( UPPER ) THEN
!
!        Begin Upper
!
!        invA = P * inv(U**T) * inv(D) * inv(U) * P**T.
!
      CALL CTRTRI( UPLO, 'U', N, A, LDA, INFO )
!
!        inv(D) and inv(D) * inv(U)
!
      K = 1
      DO WHILE( K <= N )
         IF( IPIV( K ) > 0 ) THEN
!              1 x 1 diagonal NNB
            WORK( K, INVD ) = (1.0E+0,0.0E+0) /  A( K, K )
            WORK( K, INVD+1 ) = (0.0E+0,0.0E+0)
         ELSE
!              2 x 2 diagonal NNB
            uoT = (1.0E+0,0.0E+0) /  WORK( K+1, 1 )
            AK = A( K, K ) * uoT
            AKP1 = A( K+1, K+1 ) * uoT
            AKKP1 = WORK( K+1, 1 )  * uoT
            uoD = uoT / ( AK*AKP1-(1.0E+0,0.0E+0) )
            WORK( K, INVD ) = AKP1 * uoD
            WORK( K+1, INVD+1 ) = AK * uoD
            WORK( K, INVD+1 ) = -AKKP1 * uoD
            WORK( K+1, INVD ) = WORK( K, INVD+1 )
            K = K + 1
         END IF
         K = K + 1
      END DO
!
!        inv(U**T) = (inv(U))**T
!
!        inv(U**T) * inv(D) * inv(U)
!
      CUT = N
      DO WHILE( CUT > 0 )
         NNB = NB
         IF( CUT <= NNB ) THEN
            NNB = CUT
         ELSE
!              count negative elements,
            ICOUNT = count(IPIV(CUT+1-NNB:CUT) < 0)
!              need a even number for a clear cut
            IF( MOD( ICOUNT, 2 ) == 1 ) NNB = NNB + 1
         END IF

         CUT = CUT - NNB
!
!           U01 Block
!
         DO I = 1, CUT
            WORK(I,1:NNB) = A(I,CUT+1:CUT+NNB)
         END DO
!
!           U11 Block
!
         DO I = 1, NNB
            WORK( U11+I,1:I-1) = (0.0E+0,0.0E+0)
            WORK( U11+I, I ) = (1.0E+0,0.0E+0)
            WORK( U11+I, I+1:NNB ) = A( CUT+I, CUT+I+1:CUT+NNB )
         END DO
!
!           invD * U01
!
         I = 1
         DO WHILE( I <= CUT )
            IF( IPIV( I ) > 0 ) THEN
               WORK(I,1:NNB) = WORK( I, INVD ) * WORK( I,1:NNB)
            ELSE
               DO J = 1, NNB
                  U01_I_J = WORK( I, J )
                  U01_IP1_J = WORK( I+1, J )
                  WORK( I, J ) = WORK( I, INVD ) * U01_I_J + WORK( I, INVD+1 ) * U01_IP1_J
                  WORK( I+1, J ) = WORK( I+1, INVD ) * U01_I_J + WORK( I+1, INVD+1 ) * U01_IP1_J
               END DO
               I = I + 1
            END IF
            I = I + 1
         END DO
!
!           invD1 * U11
!
         I = 1
         DO WHILE ( I <= NNB )
            IF( IPIV( CUT+I ) > 0 ) THEN
               WORK( U11+I,I:NNB) = WORK(CUT+I,INVD) * WORK(U11+I,I:NNB)
            ELSE
               DO J = I, NNB
                  U11_I_J = WORK(U11+I,J)
                  U11_IP1_J = WORK(U11+I+1,J)
                  WORK( U11+I, J ) = WORK(CUT+I,INVD) * WORK(U11+I,J) &
                               + WORK(CUT+I,INVD+1) * WORK(U11+I+1,J)
                  WORK( U11+I+1, J ) = WORK(CUT+I+1,INVD) * U11_I_J &
                                  + WORK(CUT+I+1,INVD+1) * U11_IP1_J
               END DO
               I = I + 1
            END IF
            I = I + 1
         END DO
!
!           U11**T * invD1 * U11 -> U11
!
         CALL CTRMM( 'L', 'U', 'T', 'U', NNB, NNB, &
                    (1.0E+0,0.0E+0), A( CUT+1, CUT+1 ), LDA, WORK( U11+1, 1 ), &
                    N+NB+1 )
!
         DO I = 1, NNB
            A( CUT+I, CUT+I:CUT+NNB ) = WORK( U11+I, I:NNB )
         END DO
!
!           U01**T * invD * U01 -> A( CUT+I, CUT+J )
!
         CALL CGEMM( 'T', 'N', NNB, NNB, CUT, (1.0E+0,0.0E+0), A( 1, CUT+1 ), &
                     LDA, WORK, N+NB+1, (0.0E+0,0.0E+0), WORK(U11+1,1), &
                     N+NB+1 )

!
!           U11 =  U11**T * invD1 * U11 + U01**T * invD * U01
!
         DO I = 1, NNB
            A( CUT+I, CUT+I:CUT+NNB ) = A( CUT+I, CUT+I:CUT+NNB ) + WORK(U11+I,I:NNB)
         END DO
!
!           U01 =  U00**T * invD0 * U01
!
         CALL CTRMM( 'L', UPLO, 'T', 'U', CUT, NNB, &
                     (1.0E+0,0.0E+0), A, LDA, WORK, N+NB+1 )

!
!           Update U01
!
         A(1:CUT, CUT+1:CUT+NNB ) = WORK(1:CUT,1:NNB)
!
!           Next Block
!
      END DO
!
!        Apply PERMUTATIONS P and P**T:
!        P * inv(U**T) * inv(D) * inv(U) * P**T.
!        Interchange rows and columns I and IPIV(I) in reverse order
!        from the formation order of IPIV vector for Upper case.
!
!        ( We can use a loop over IPIV with increment 1,
!        since the ABS value of IPIV(I) represents the row (column)
!        index of the interchange with row (column) i in both 1x1
!        and 2x2 pivot cases, i.e. we don't need separate code branches
!        for 1x1 and 2x2 pivot cases )
!
      DO I = 1, N
          IP = ABS( IPIV( I ) )
          IF (I  <  IP) THEN
             CALL CSYSWAPR( UPLO, N, A, LDA, I ,IP )
          ELSEIF (I  >  IP) THEN
             CALL CSYSWAPR( UPLO, N, A, LDA, IP ,I )
          END IF
      END DO
!
   ELSE
!
!        Begin Lower
!
!        inv A = P * inv(L**T) * inv(D) * inv(L) * P**T.
!
      CALL CTRTRI( UPLO, 'U', N, A, LDA, INFO )
!
!        inv(D) and inv(D) * inv(L)
!
      K = N
      DO WHILE ( K  >=  1 )
         IF( IPIV( K ) > 0 ) THEN
!              1 x 1 diagonal NNB
            WORK( K, INVD ) = (1.0E+0,0.0E+0) / A( K, K )
            WORK( K, INVD+1 ) = (0.0E+0,0.0E+0)
         ELSE
!              2 x 2 diagonal NNB
            uoT = (1.0E+0,0.0E+0) / WORK( K-1, 1 )
            AK = A( K-1, K-1 ) * uoT
            AKP1 = A( K, K ) * uoT
            AKKP1 = WORK( K-1, 1 ) * uoT
            uoD = uoT/( AK*AKP1-(1.0E+0,0.0E+0) )
            WORK( K-1, INVD ) = AKP1 * uoD
            WORK( K, INVD ) = AK * uoD
            WORK( K, INVD+1 ) = -AKKP1 * uoD
            WORK( K-1, INVD+1 ) = WORK( K, INVD+1 )
            K = K - 1
         END IF
         K = K - 1
      END DO
!
!        inv(L**T) = (inv(L))**T
!
!        inv(L**T) * inv(D) * inv(L)
!
      CUT = 0
      DO WHILE( CUT < N )
         NNB = NB
         IF( (CUT + NNB) > N ) THEN
            NNB = N - CUT
         ELSE
!              count negative elements,
            ICOUNT = COUNT(IPIV(CUT + 1:CUT+NNB) < 0)
!              need a even number for a clear cut
            IF( MOD( ICOUNT, 2 ) == 1 ) NNB = NNB + 1
         END IF
!
!           L21 Block
!
         WORK(1:N-CUT-NNB,1:NNB) = A( CUT+NNB+1:N, CUT+1:CUT+NNB )
!
!           L11 Block
!
         DO I = 1, NNB
            WORK( U11+I,1:I-1) = A( CUT+I, CUT+1:CUT+I-1 )
            WORK( U11+I, I) = (1.0E+0,0.0E+0)
            WORK( U11+I,I+1:NNB) = (0.0E+0,0.0E+0)
         END DO
!
!           invD*L21
!
         I = N-CUT-NNB
         DO WHILE( I >= 1 )
            IF( IPIV( CUT+NNB+I ) > 0 ) THEN
               WORK( I,1:NNB) = WORK( CUT+NNB+I, INVD) * WORK( I,1:NNB)
            ELSE
               DO J = 1, NNB
                  U01_I_J = WORK(I,J)
                  U01_IP1_J = WORK(I-1,J)
                  WORK(I,J)=WORK(CUT+NNB+I,INVD)*U01_I_J+ &
                           WORK(CUT+NNB+I,INVD+1)*U01_IP1_J
                  WORK(I-1,J)=WORK(CUT+NNB+I-1,INVD+1)*U01_I_J+ &
                           WORK(CUT+NNB+I-1,INVD)*U01_IP1_J
               END DO
               I = I - 1
            END IF
            I = I - 1
         END DO
!
!           invD1*L11
!
         I = NNB
         DO WHILE( I >= 1 )
            IF( IPIV( CUT+I ) > 0 ) THEN
               WORK( U11+I,1:NNB) = WORK( CUT+I, INVD)*WORK(U11+I,1:NNB)
            ELSE
               DO J = 1, NNB
                  U11_I_J = WORK( U11+I, J )
                  U11_IP1_J = WORK( U11+I-1, J )
                  WORK( U11+I, J ) = WORK(CUT+I,INVD) * WORK(U11+I,J) &
                                   + WORK(CUT+I,INVD+1) * U11_IP1_J
                  WORK( U11+I-1, J ) = WORK(CUT+I-1,INVD+1) * U11_I_J &
                                     + WORK(CUT+I-1,INVD) * U11_IP1_J
               END DO
               I = I - 1
            END IF
            I = I - 1
         END DO
!
!           L11**T * invD1 * L11 -> L11
!
         CALL CTRMM( 'L', UPLO, 'T', 'U', NNB, NNB, (1.0E+0,0.0E+0), &
                      A( CUT+1, CUT+1 ), LDA, WORK( U11+1, 1 ), &
                      N+NB+1 )

!
         DO I = 1, NNB
            A( CUT+I, CUT+1:CUT+I ) = WORK( U11+I,1:I)
         END DO
!
         IF( (CUT+NNB) < N ) THEN
!
!              L21**T * invD2*L21 -> A( CUT+I, CUT+J )
!
            CALL CGEMM( 'T', 'N', NNB, NNB, N-NNB-CUT, (1.0E+0,0.0E+0), &
                        A( CUT+NNB+1, CUT+1 ), LDA, WORK, N+NB+1, &
                        (0.0E+0,0.0E+0), WORK( U11+1, 1 ), N+NB+1 )

!
!              L11 =  L11**T * invD1 * L11 + U01**T * invD * U01
!
            DO I = 1, NNB
               A( CUT+I, CUT+1:CUT+I ) = A( CUT+I, CUT+1:CUT+I )+WORK(U11+I,1:I)
            END DO
!
!              L01 =  L22**T * invD2 * L21
!
            CALL CTRMM( 'L', UPLO, 'T', 'U', N-NNB-CUT, NNB, (1.0E+0,0.0E+0), &
                        A( CUT+NNB+1, CUT+NNB+1 ), LDA, WORK, &
                        N+NB+1 )
!
!              Update L21
!
            A(CUT+NNB+1:N,CUT+1:CUT+NNB) = WORK(1:N-CUT-NNB,1:NNB)
!
         ELSE
!
!              L11 =  L11**T * invD1 * L11
!
            DO I = 1, NNB
               A( CUT+I, CUT+1:CUT+I ) = WORK( U11+I,1:I)
            END DO
         END IF
!
!           Next Block
!
         CUT = CUT + NNB
!
      END DO
!
!        Apply PERMUTATIONS P and P**T:
!        P * inv(L**T) * inv(D) * inv(L) * P**T.
!        Interchange rows and columns I and IPIV(I) in reverse order
!        from the formation order of IPIV vector for Lower case.
!
!        ( We can use a loop over IPIV with increment -1,
!        since the ABS value of IPIV(I) represents the row (column)
!        index of the interchange with row (column) i in both 1x1
!        and 2x2 pivot cases, i.e. we don't need separate code branches
!        for 1x1 and 2x2 pivot cases )
!
      DO I = N, 1, -1
          IP = ABS( IPIV( I ) )
          IF (I  <  IP) THEN
             CALL CSYSWAPR( UPLO, N, A, LDA, I ,IP )
          ELSEIF (I  >  IP) THEN
             CALL CSYSWAPR( UPLO, N, A, LDA, IP ,I )
          END IF
      END DO
!
   END IF
!
   RETURN
!
!     End of CSYTRI_3X
!
   END
