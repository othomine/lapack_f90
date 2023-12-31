!> \brief \b CSYTRF_AA_2STAGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYTRF_AA_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrf_aa_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrf_aa_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrf_aa_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!      SUBROUTINE CSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV,
!                                   IPIV2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, LDA, LTB, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IPIV2( * )
!       COMPLEX            A( LDA, * ), TB( * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYTRF_AA_2STAGE computes the factorization of a complex symmetric matrix A
!> using the Aasen's algorithm.  The form of the factorization is
!>
!>    A = U**T*T*U  or  A = L*T*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and T is a complex symmetric band matrix with the
!> bandwidth of NB (NB is internally selected and stored in TB( 1 ), and T is
!> LU factorized with partial pivoting).
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
!>          On entry, the hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, L is stored below (or above) the subdiagonal blocks,
!>          when UPLO  is 'L' (or 'U').
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TB
!> \verbatim
!>          TB is COMPLEX array, dimension (LTB)
!>          On exit, details of the LU factorization of the band matrix.
!> \endverbatim
!>
!> \param[in] LTB
!> \verbatim
!>          LTB is INTEGER
!>          The size of the array TB. LTB >= 4*N, internally
!>          used to select NB such that LTB >= (3*NB+1)*N.
!>
!>          If LTB = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of LTB,
!>          returns this value as the first entry of TB, and
!>          no error message related to LTB is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          On exit, it contains the details of the interchanges, i.e.,
!>          the row and column k of A were interchanged with the
!>          row and column IPIV(k).
!> \endverbatim
!>
!> \param[out] IPIV2
!> \verbatim
!>          IPIV2 is INTEGER array, dimension (N)
!>          On exit, it contains the details of the interchanges, i.e.,
!>          the row and column k of T were interchanged with the
!>          row and column IPIV(k).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX workspace of size LWORK
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The size of WORK. LWORK >= N, internally used to select NB
!>          such that LWORK >= N*NB.
!>
!>          If LWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the WORK array,
!>          returns this value as the first entry of the WORK array, and
!>          no error message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, band LU factorization failed on i-th column
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
!> \ingroup hetrf_aa_2stage
!
!  =====================================================================
   SUBROUTINE CSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, &
                                IPIV2, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
   IMPLICIT NONE
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            N, LDA, LTB, LWORK, INFO
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * ), IPIV2( * )
   COMPLEX            A( LDA, * ), TB( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   LOGICAL            UPPER, TQUERY, WQUERY
   INTEGER            I, J, K, I1, I2, TD
   INTEGER            LDTB, NB, KB, JB, NT, IINFO
   COMPLEX            PIV
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
   EXTERNAL           CCOPY, CGBTRF, CGEMM, CGETRF, CLACPY, &
                      CLASET, CTRSM, CSWAP, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   WQUERY = ( LWORK == -1 )
   TQUERY = ( LTB == -1 )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   ELSE IF ( LTB  <  4*N .AND. .NOT.TQUERY ) THEN
      INFO = -6
   ELSE IF ( LWORK  <  N .AND. .NOT.WQUERY ) THEN
      INFO = -10
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSYTRF_AA_2STAGE', -INFO )
      RETURN
   END IF
!
!     Answer the query
!
   NB = ILAENV( 1, 'CSYTRF_AA_2STAGE', UPLO, N, -1, -1, -1 )
   IF( INFO == 0 ) THEN
      IF( TQUERY ) TB( 1 ) = (3*NB+1)*N
      IF( WQUERY ) WORK( 1 ) = N*NB
   END IF
   IF( TQUERY .OR. WQUERY ) RETURN
!
!     Quick return
!
   IF ( N == 0 ) RETURN
!
!     Determine the number of the block size
!
   LDTB = LTB/N
   IF( LDTB  <  3*NB+1 ) NB = (LDTB-1)/3
   IF( LWORK  <  NB*N ) NB = LWORK/N
!
!     Determine the number of the block columns
!
   NT = (N+NB-1)/NB
   TD = 2*NB
   KB = MIN(NB, N)
!
!     Initialize vectors/matrices
!
   DO J = 1, KB
      IPIV( J ) = J
   END DO
!
!     Save NB
!
   TB( 1 ) = NB
!
   IF( UPPER ) THEN
!
!        .....................................................
!        Factorize A as U**T*D*U using the upper triangle of A
!        .....................................................
!
      DO J = 0, NT-1
!
!           Generate Jth column of W and H
!
         KB = MIN(NB, N-J*NB)
         DO I = 1, J-1
            IF( I == 1 ) THEN
!                  H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J)
               IF( I  ==  (J-1) ) THEN
                  JB = NB+KB
               ELSE
                  JB = 2*NB
               END IF
               CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                       NB, KB, JB, &
                       (1.0E+0,0.0E+0),  TB( TD+1 + (I*NB)*LDTB ), LDTB-1, &
                              A( (I-1)*NB+1, J*NB+1 ), LDA, &
                       (0.0E+0,0.0E+0), WORK( I*NB+1 ), N )
            ELSE
!                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
               IF( I  ==  J-1) THEN
                  JB = 2*NB+KB
               ELSE
                  JB = 3*NB
               END IF
               CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                       NB, KB, JB, &
                       (1.0E+0,0.0E+0),  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), &
                          LDTB-1, &
                              A( (I-2)*NB+1, J*NB+1 ), LDA, &
                       (0.0E+0,0.0E+0), WORK( I*NB+1 ), N )
            END IF
         END DO
!
!           Compute T(J,J)
!
         CALL CLACPY( 'Upper', KB, KB, A( J*NB+1, J*NB+1 ), LDA, &
                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
         IF( J > 1 ) THEN
!              T(J,J) = U(1:J,J)'*H(1:J)
            CALL CGEMM( 'Transpose', 'NoTranspose', &
                    KB, KB, (J-1)*NB, &
                   -(1.0E+0,0.0E+0), A( 1, J*NB+1 ), LDA, &
                          WORK( NB+1 ), N, &
                    (1.0E+0,0.0E+0), TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
!              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
            CALL CGEMM( 'Transpose', 'NoTranspose', &
                    KB, NB, KB, &
                    (1.0E+0,0.0E+0),  A( (J-1)*NB+1, J*NB+1 ), LDA, &
                           TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, &
                    (0.0E+0,0.0E+0), WORK( 1 ), N )
            CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                    KB, KB, NB, &
                   -(1.0E+0,0.0E+0), WORK( 1 ), N, &
                          A( (J-2)*NB+1, J*NB+1 ), LDA, &
                    (1.0E+0,0.0E+0), TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
         END IF
!
!           Expand T(J,J) into full format
!
         DO I = 1, KB
            DO K = I+1, KB
               TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB ) &
                   = TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB )
            END DO
         END DO
         IF( J > 0 ) THEN
!               CALL CHEGST( 1, 'Upper', KB,
!     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1,
!     $                      A( (J-1)*NB+1, J*NB+1 ), LDA, IINFO )
            CALL CTRSM( 'L', 'U', 'T', 'N', KB, KB, (1.0E+0,0.0E+0), &
                        A( (J-1)*NB+1, J*NB+1 ), LDA, &
                        TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            CALL CTRSM( 'R', 'U', 'N', 'N', KB, KB, (1.0E+0,0.0E+0), &
                        A( (J-1)*NB+1, J*NB+1 ), LDA, &
                        TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
         END IF
!
         IF( J < NT-1 ) THEN
            IF( J > 0 ) THEN
!
!                 Compute H(J,J)
!
               IF( J == 1 ) THEN
                  CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                          KB, KB, KB, &
                          (1.0E+0,0.0E+0),  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, &
                                 A( (J-1)*NB+1, J*NB+1 ), LDA, &
                          (0.0E+0,0.0E+0), WORK( J*NB+1 ), N )
               ELSE
                  CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                         KB, KB, NB+KB, &
                         (1.0E+0,0.0E+0), TB( TD+NB+1 + ((J-1)*NB)*LDTB ), &
                            LDTB-1, &
                                A( (J-2)*NB+1, J*NB+1 ), LDA, &
                         (0.0E+0,0.0E+0), WORK( J*NB+1 ), N )
               END IF
!
!                 Update with the previous column
!
               CALL CGEMM( 'Transpose', 'NoTranspose', &
                       NB, N-(J+1)*NB, J*NB, &
                       -(1.0E+0,0.0E+0), WORK( NB+1 ), N, &
                              A( 1, (J+1)*NB+1 ), LDA, &
                        (1.0E+0,0.0E+0), A( J*NB+1, (J+1)*NB+1 ), LDA )
            END IF
!
!              Copy panel to workspace to call CGETRF
!
            DO K = 1, NB
                CALL CCOPY( N-(J+1)*NB, &
                            A( J*NB+K, (J+1)*NB+1 ), LDA, &
                            WORK( 1+(K-1)*N ), 1 )
            END DO
!
!              Factorize panel
!
            CALL CGETRF( N-(J+1)*NB, NB, &
                         WORK, N, &
                         IPIV( (J+1)*NB+1 ), IINFO )
!               IF (IINFO /= 0 .AND. INFO == 0) THEN
!                  INFO = IINFO+(J+1)*NB
!               END IF
!
!              Copy panel back
!
            DO K = 1, NB
                CALL CCOPY( N-(J+1)*NB, &
                            WORK( 1+(K-1)*N ), 1, &
                            A( J*NB+K, (J+1)*NB+1 ), LDA )
            END DO
!
!              Compute T(J+1, J), zero out for GEMM update
!
            KB = MIN(NB, N-(J+1)*NB)
            CALL CLASET( 'Full', KB, NB, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), &
                         TB( TD+NB+1 + (J*NB)*LDTB), LDTB-1 )
            CALL CLACPY( 'Upper', KB, NB, &
                         WORK, N, &
                         TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
            IF( J > 0 ) THEN
               CALL CTRSM( 'R', 'U', 'N', 'U', KB, NB, (1.0E+0,0.0E+0), &
                           A( (J-1)*NB+1, J*NB+1 ), LDA, &
                           TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
            END IF
!
!              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
!              updates
!
            DO K = 1, NB
               DO I = 1, KB
                  TB( TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB ) &
                     = TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB )
               END DO
            END DO
            CALL CLASET( 'Lower', KB, NB, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), &
                         A( J*NB+1, (J+1)*NB+1), LDA )
!
!              Apply pivots to trailing submatrix of A
!
            DO K = 1, KB
!                 > Adjust ipiv
               IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB
!
               I1 = (J+1)*NB+K
               I2 = IPIV( (J+1)*NB+K )
               IF( I1 /= I2 ) THEN
!                    > Apply pivots to previous columns of L
                  CALL CSWAP( K-1, A( (J+1)*NB+1, I1 ), 1, A( (J+1)*NB+1, I2 ), 1 )
!                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                  IF( I2 > (I1+1) ) &
                     CALL CSWAP( I2-I1-1, A( I1, I1+1 ), LDA, A( I1+1, I2 ), 1 )
!                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                  IF( I2 < N ) &
                     CALL CSWAP( N-I2, A( I1, I2+1 ), LDA, A( I2, I2+1 ), LDA )
!                    > Swap A(I1, I1) with A(I2, I2)
                  PIV = A( I1, I1 )
                  A( I1, I1 ) = A( I2, I2 )
                  A( I2, I2 ) = PIV
!                    > Apply pivots to previous columns of L
                  IF( J > 0 ) THEN
                     CALL CSWAP( J*NB, A( 1, I1 ), 1, A( 1, I2 ), 1 )
                  END IF
               ENDIF
            END DO
         END IF
      END DO
   ELSE
!
!        .....................................................
!        Factorize A as L*D*L**T using the lower triangle of A
!        .....................................................
!
      DO J = 0, NT-1
!
!           Generate Jth column of W and H
!
         KB = MIN(NB, N-J*NB)
         DO I = 1, J-1
            IF( I == 1 ) THEN
!                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)'
               IF( I  ==  (J-1) ) THEN
                  JB = NB+KB
               ELSE
                  JB = 2*NB
               END IF
               CALL CGEMM( 'NoTranspose', 'Transpose', &
                       NB, KB, JB, &
                       (1.0E+0,0.0E+0), TB( TD+1 + (I*NB)*LDTB ), LDTB-1, &
                             A( J*NB+1, (I-1)*NB+1 ), LDA, &
                       (0.0E+0,0.0E+0), WORK( I*NB+1 ), N )
            ELSE
!                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
               IF( I  ==  (J-1) ) THEN
                  JB = 2*NB+KB
               ELSE
                  JB = 3*NB
               END IF
               CALL CGEMM( 'NoTranspose', 'Transpose', &
                       NB, KB, JB, &
                       (1.0E+0,0.0E+0),  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), &
                          LDTB-1, &
                              A( J*NB+1, (I-2)*NB+1 ), LDA, &
                       (0.0E+0,0.0E+0), WORK( I*NB+1 ), N )
            END IF
         END DO
!
!           Compute T(J,J)
!
         CALL CLACPY( 'Lower', KB, KB, A( J*NB+1, J*NB+1 ), LDA, &
                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
         IF( J > 1 ) THEN
!              T(J,J) = L(J,1:J)*H(1:J)
            CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                    KB, KB, (J-1)*NB, &
                   -(1.0E+0,0.0E+0), A( J*NB+1, 1 ), LDA, &
                          WORK( NB+1 ), N, &
                    (1.0E+0,0.0E+0), TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
!              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
            CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                    KB, NB, KB, &
                    (1.0E+0,0.0E+0),  A( J*NB+1, (J-1)*NB+1 ), LDA, &
                           TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, &
                    (0.0E+0,0.0E+0), WORK( 1 ), N )
            CALL CGEMM( 'NoTranspose', 'Transpose', &
                    KB, KB, NB, &
                   -(1.0E+0,0.0E+0), WORK( 1 ), N, &
                          A( J*NB+1, (J-2)*NB+1 ), LDA, &
                    (1.0E+0,0.0E+0), TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
         END IF
!
!           Expand T(J,J) into full format
!
         DO I = 1, KB
            DO K = I+1, KB
               TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB ) &
                   = TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB )
            END DO
         END DO
         IF( J > 0 ) THEN
!               CALL CHEGST( 1, 'Lower', KB,
!     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1,
!     $                      A( J*NB+1, (J-1)*NB+1 ), LDA, IINFO )
            CALL CTRSM( 'L', 'L', 'N', 'N', KB, KB, (1.0E+0,0.0E+0), &
                        A( J*NB+1, (J-1)*NB+1 ), LDA, &
                        TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            CALL CTRSM( 'R', 'L', 'T', 'N', KB, KB, (1.0E+0,0.0E+0), &
                        A( J*NB+1, (J-1)*NB+1 ), LDA, &
                        TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
         END IF
!
!           Symmetrize T(J,J)
!
         DO I = 1, KB
            DO K = I+1, KB
               TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB ) &
                   = TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB )
            END DO
         END DO
!
         IF( J < NT-1 ) THEN
            IF( J > 0 ) THEN
!
!                 Compute H(J,J)
!
               IF( J == 1 ) THEN
                  CALL CGEMM( 'NoTranspose', 'Transpose', &
                          KB, KB, KB, &
                          (1.0E+0,0.0E+0),  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, &
                                 A( J*NB+1, (J-1)*NB+1 ), LDA, &
                          (0.0E+0,0.0E+0), WORK( J*NB+1 ), N )
               ELSE
                  CALL CGEMM( 'NoTranspose', 'Transpose', &
                         KB, KB, NB+KB, &
                         (1.0E+0,0.0E+0), TB( TD+NB+1 + ((J-1)*NB)*LDTB ), &
                            LDTB-1, &
                                A( J*NB+1, (J-2)*NB+1 ), LDA, &
                         (0.0E+0,0.0E+0), WORK( J*NB+1 ), N )
               END IF
!
!                 Update with the previous column
!
               CALL CGEMM( 'NoTranspose', 'NoTranspose', &
                       N-(J+1)*NB, NB, J*NB, &
                       -(1.0E+0,0.0E+0), A( (J+1)*NB+1, 1 ), LDA, &
                              WORK( NB+1 ), N, &
                        (1.0E+0,0.0E+0), A( (J+1)*NB+1, J*NB+1 ), LDA )
            END IF
!
!              Factorize panel
!
            CALL CGETRF( N-(J+1)*NB, NB, &
                         A( (J+1)*NB+1, J*NB+1 ), LDA, &
                         IPIV( (J+1)*NB+1 ), IINFO )
!               IF (IINFO /= 0 .AND. INFO == 0) THEN
!                  INFO = IINFO+(J+1)*NB
!               END IF
!
!              Compute T(J+1, J), zero out for GEMM update
!
            KB = MIN(NB, N-(J+1)*NB)
            CALL CLASET( 'Full', KB, NB, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), &
                         TB( TD+NB+1 + (J*NB)*LDTB), LDTB-1 )
            CALL CLACPY( 'Upper', KB, NB, &
                         A( (J+1)*NB+1, J*NB+1 ), LDA, &
                         TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
            IF( J > 0 ) THEN
               CALL CTRSM( 'R', 'L', 'T', 'U', KB, NB, (1.0E+0,0.0E+0), &
                           A( J*NB+1, (J-1)*NB+1 ), LDA, &
                           TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
            END IF
!
!              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
!              updates
!
            DO K = 1, NB
               DO I = 1, KB
                  TB( TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB ) = &
                     TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB )
               END DO
            END DO
            CALL CLASET( 'Upper', KB, NB, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), &
                         A( (J+1)*NB+1, J*NB+1 ), LDA )
!
!              Apply pivots to trailing submatrix of A
!
            DO K = 1, KB
!                 > Adjust ipiv
               IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB
!
               I1 = (J+1)*NB+K
               I2 = IPIV( (J+1)*NB+K )
               IF( I1 /= I2 ) THEN
!                    > Apply pivots to previous columns of L
                  CALL CSWAP( K-1, A( I1, (J+1)*NB+1 ), LDA, A( I2, (J+1)*NB+1 ), LDA )
!                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                  IF( I2 > (I1+1) ) &
                     CALL CSWAP( I2-I1-1, A( I1+1, I1 ), 1, A( I2, I1+1 ), LDA )
!                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                  IF( I2 < N ) &
                     CALL CSWAP( N-I2, A( I2+1, I1 ), 1, A( I2+1, I2 ), 1 )
!                    > Swap A(I1, I1) with A(I2, I2)
                  PIV = A( I1, I1 )
                  A( I1, I1 ) = A( I2, I2 )
                  A( I2, I2 ) = PIV
!                    > Apply pivots to previous columns of L
                  IF( J > 0 ) THEN
                     CALL CSWAP( J*NB, A( I1, 1 ), LDA, A( I2, 1 ), LDA )
                  END IF
               ENDIF
            END DO
!
!              Apply pivots to previous columns of L
!
!               CALL CLASWP( J*NB, A( 1, 1 ), LDA,
!     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
         END IF
      END DO
   END IF
!
!     Factor the band matrix
   CALL CGBTRF( N, N, NB, NB, TB, LDTB, IPIV2, INFO )
!
   RETURN
!
!     End of CSYTRF_AA_2STAGE
!
END
