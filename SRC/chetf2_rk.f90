!> \brief \b CHETF2_RK computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETF2_RK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetf2_rk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetf2_rk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetf2_rk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), E ( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> CHETF2_RK computes the factorization of a complex Hermitian matrix A
!> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
!>
!>    A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**H (or L**H) is the conjugate of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is Hermitian and block
!> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
!> For more information see Further Details section.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored:
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
!>          On entry, the Hermitian matrix A.
!>            If UPLO = 'U': the leading N-by-N upper triangular part
!>            of A contains the upper triangular part of the matrix A,
!>            and the strictly lower triangular part of A is not
!>            referenced.
!>
!>            If UPLO = 'L': the leading N-by-N lower triangular part
!>            of A contains the lower triangular part of the matrix A,
!>            and the strictly upper triangular part of A is not
!>            referenced.
!>
!>          On exit, contains:
!>            a) ONLY diagonal elements of the Hermitian block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                are stored on exit in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is COMPLEX array, dimension (N)
!>          On exit, contains the superdiagonal (or subdiagonal)
!>          elements of the Hermitian block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
!>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
!>
!>          NOTE: For 1-by-1 diagonal block D(k), where
!>          1 <= k <= N, the element E(k) is set to 0 in both
!>          UPLO = 'U' or UPLO = 'L' cases.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          IPIV describes the permutation matrix P in the factorization
!>          of matrix A as follows. The absolute value of IPIV(k)
!>          represents the index of row and column that were
!>          interchanged with the k-th row and column. The value of UPLO
!>          describes the order in which the interchanges were applied.
!>          Also, the sign of IPIV represents the block structure of
!>          the Hermitian block diagonal matrix D with 1-by-1 or 2-by-2
!>          diagonal blocks which correspond to 1 or 2 interchanges
!>          at each factorization step. For more info see Further
!>          Details section.
!>
!>          If UPLO = 'U',
!>          ( in factorization order, k decreases from N to 1 ):
!>            a) A single positive entry IPIV(k) > 0 means:
!>               D(k,k) is a 1-by-1 diagonal block.
!>               If IPIV(k) != k, rows and columns k and IPIV(k) were
!>               interchanged in the matrix A(1:N,1:N);
!>               If IPIV(k) = k, no interchange occurred.
!>
!>            b) A pair of consecutive negative entries
!>               IPIV(k) < 0 and IPIV(k-1) < 0 means:
!>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>               (NOTE: negative entries in IPIV appear ONLY in pairs).
!>               1) If -IPIV(k) != k, rows and columns
!>                  k and -IPIV(k) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k) = k, no interchange occurred.
!>               2) If -IPIV(k-1) != k-1, rows and columns
!>                  k-1 and -IPIV(k-1) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k-1) = k-1, no interchange occurred.
!>
!>            c) In both cases a) and b), always ABS( IPIV(k) ) <= k.
!>
!>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
!>
!>          If UPLO = 'L',
!>          ( in factorization order, k increases from 1 to N ):
!>            a) A single positive entry IPIV(k) > 0 means:
!>               D(k,k) is a 1-by-1 diagonal block.
!>               If IPIV(k) != k, rows and columns k and IPIV(k) were
!>               interchanged in the matrix A(1:N,1:N).
!>               If IPIV(k) = k, no interchange occurred.
!>
!>            b) A pair of consecutive negative entries
!>               IPIV(k) < 0 and IPIV(k+1) < 0 means:
!>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!>               (NOTE: negative entries in IPIV appear ONLY in pairs).
!>               1) If -IPIV(k) != k, rows and columns
!>                  k and -IPIV(k) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k) = k, no interchange occurred.
!>               2) If -IPIV(k+1) != k+1, rows and columns
!>                  k-1 and -IPIV(k-1) were interchanged
!>                  in the matrix A(1:N,1:N).
!>                  If -IPIV(k+1) = k+1, no interchange occurred.
!>
!>            c) In both cases a) and b), always ABS( IPIV(k) ) >= k.
!>
!>            d) NOTE: Any entry IPIV(k) is always NONZERO on output.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>
!>          < 0: If INFO = -k, the k-th argument had an illegal value
!>
!>          > 0: If INFO = k, the matrix A is singular, because:
!>                 If UPLO = 'U': column k in the upper
!>                 triangular part of A contains all zeros.
!>                 If UPLO = 'L': column k in the lower
!>                 triangular part of A contains all zeros.
!>
!>               Therefore D(k,k) is exactly zero, and superdiagonal
!>               elements of column k of U (or subdiagonal elements of
!>               column k of L ) are all zeros. The factorization has
!>               been completed, but the block diagonal matrix D is
!>               exactly singular, and division by zero will occur if
!>               it is used to solve a system of equations.
!>
!>               NOTE: INFO only stores the first occurrence of
!>               a singularity, any subsequent occurrence of singularity
!>               is not stored in INFO even though the factorization
!>               always completes.
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
!> \ingroup hetf2_rk
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> TODO: put further details
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  December 2016,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!>  01-01-96 - Based on modifications by
!>    J. Lewis, Boeing Computer Services Company
!>    A. Petitet, Computer Science Dept.,
!>                Univ. of Tenn., Knoxville abd , USA
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE CHETF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )
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
   INTEGER            IPIV( * )
   COMPLEX            A( LDA, * ), E( * )
!     ..
!
!  ======================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            DONE, UPPER
   INTEGER            I, II, IMAX, ITEMP, J, JMAX, K, KK, KP, KSTEP, &
                      P
   REAL               ABSAKK, COLMAX, D, D11, D22, R1, STEMP, &
                      ROWMAX, TT, SFMIN
!
!     Initialize ALPHA for use in choosing pivot block size.
!
   REAL, PARAMETER :: ALPHA = ( 1.0E+0+SQRT( 17.0E+0 ) ) / 8.0E+0
   COMPLEX            D12, D21, T, WK, WKM1, WKP1, Z, A_TMP( LDA ), AT_TMP( N )
!     ..
!     .. External Functions ..
!
   LOGICAL            LSAME
   INTEGER            ICAMAX
   REAL               SLAMCH, SLAPY2, CABS1
   EXTERNAL           LSAME, ICAMAX, SLAMCH, SLAPY2, CABS1
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, CHER
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
      CALL XERBLA( 'CHETF2_RK', -INFO )
      RETURN
   END IF
!
!     Compute machine safe minimum
!
   SFMIN = SLAMCH( 'S' )
!
   IF( UPPER ) THEN
!
!        Factorize A as U*D*U**H using the upper triangle of A
!
!        Initialize the first entry of array E, where superdiagonal
!        elements of D are stored
!
      E( 1 ) = (0.0E+0,0.0E+0)
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
      K = N
10    CONTINUE
!
!        If K < 1, exit from loop
!
      IF( K < 1 ) GO TO 34
      KSTEP = 1
      P = K
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
      ABSAKK = ABS( REAL( A( K, K ) ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
      IF( K > 1 ) THEN
         IMAX = ICAMAX( K-1, A( 1, K ), 1 )
         COLMAX = CABS1( A( IMAX, K ) )
      ELSE
         COLMAX = 0.0E+0
      END IF
!
      IF( ( MAX( ABSAKK, COLMAX ) == 0.0E+0 ) ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
         IF( INFO == 0 ) INFO = K
         KP = K
         A( K, K ) = REAL( A( K, K ) )
!
!           Set E( K ) to zero
!
         IF( K > 1 ) E( K ) = (0.0E+0,0.0E+0)
!
      ELSE
!
!           ============================================================
!
!           BEGIN pivot search
!
!           Case(1)
!           Equivalent to testing for ABSAKK >= ALPHA*COLMAX
!           (used to handle NaN and Inf)
!
         IF( ABSAKK >= ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
            KP = K
!
         ELSE
!
            DONE = .FALSE.
!
!              Loop until pivot found
!
12          CONTINUE
!
!                 BEGIN pivot search loop body
!
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
               IF( IMAX /= K ) THEN
                  JMAX = IMAX + ICAMAX( K-IMAX, A( IMAX, IMAX+1 ), &
                                        LDA )
                  ROWMAX = CABS1( A( IMAX, JMAX ) )
               ELSE
                  ROWMAX = 0.0E+0
               END IF
!
               IF( IMAX > 1 ) THEN
                  ITEMP = ICAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  STEMP = CABS1( A( ITEMP, IMAX ) )
                  IF( STEMP > ROWMAX ) THEN
                     ROWMAX = STEMP
                     JMAX = ITEMP
                  END IF
               END IF
!
!                 Case(2)
!                 Equivalent to testing for
!                 ABS( REAL( W( IMAX,KW-1 ) ) ) >= ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
               IF( .NOT.( ABS( REAL( A( IMAX, IMAX ) ) ) &
                           < ALPHA*ROWMAX ) ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                  KP = IMAX
                  DONE = .TRUE.
!
!                 Case(3)
!                 Equivalent to testing for ROWMAX == COLMAX,
!                 (used to handle NaN and Inf)
!
               ELSE IF( ( P == JMAX ) .OR. ( ROWMAX <= COLMAX ) ) THEN
!
!                    interchange rows and columns K-1 and IMAX,
!                    use 2-by-2 pivot block
!
                  KP = IMAX
                  KSTEP = 2
                  DONE = .TRUE.
!
!                 Case(4)
               ELSE
!
!                    Pivot not found: set params and repeat
!
                  P = IMAX
                  COLMAX = ROWMAX
                  IMAX = JMAX
               END IF
!
!                 END pivot search loop body
!
            IF( .NOT.DONE ) GOTO 12
!
         END IF
!
!           END pivot search
!
!           ============================================================
!
!           KK is the column of A where pivoting step stopped
!
         KK = K - KSTEP + 1
!
!           For only a 2x2 pivot, interchange rows and columns K and P
!           in the leading submatrix A(1:k,1:k)
!
         IF( ( KSTEP == 2 ) .AND. ( P /= K ) ) THEN
!              (1) Swap columnar parts
            IF( P > 1 ) THEN
               A_TMP(1:P-1) = A(1:P-1,K)
               A(1:P-1,K) = A(1:P-1,P)
               A(1:P-1,P) = A_TMP(1:P-1)
            ENDIF
!              (2) Swap and conjugate middle parts
            A_TMP(1:K-1-P) = CONJG(A(P+1:K-1,K) )
            A(P+1:K-1, K ) = CONJG(A(P,P+1:K-1) )
            A(P,P+1:K-1) = A_TMP(1:K-1-P)
!              (3) Swap and conjugate corner elements at row-col intersection
            A( P, K ) = CONJG( A( P, K ) )
!              (4) Swap diagonal elements at row-col intersection
            R1 = REAL( A( K, K ) )
            A( K, K ) = REAL( A( P, P ) )
            A( P, P ) = R1
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
            IF( K < N ) THEN
               AT_TMP(1:N-K) = A(K,K+1:N)
               A(K,K+1:N) = A(P,K+1:N)
               A(P,K+1:N) = AT_TMP(1:N-K)
            ENDIF
!
         END IF
!
!           For both 1x1 and 2x2 pivots, interchange rows and
!           columns KK and KP in the leading submatrix A(1:k,1:k)
!
         IF( KP /= KK ) THEN
!              (1) Swap columnar parts
            IF( KP > 1 ) THEN
               A_TMP(1:KP-1) = A(1:KP-1,KK)
               A(1:KP-1,KK) = A(1:KP-1,KP)
               A(1:KP-1,KP) = A_TMP(1:KP-1)
            ENDIF
!              (2) Swap and conjugate middle parts
            A_TMP(1:KK-KP-1) = CONJG(A(KP+1:KK-1,KK))
            A(KP+1:KK-1,KK) = CONJG(A(KP,KP+1:KK-1))
            A(KP,KP+1:KK-1) = A_TMP(1:KK-KP-1)
!              (3) Swap and conjugate corner elements at row-col intersection
            A( KP, KK ) = CONJG( A( KP, KK ) )
!              (4) Swap diagonal elements at row-col intersection
            R1 = REAL( A( KK, KK ) )
            A( KK, KK ) = REAL( A( KP, KP ) )
            A( KP, KP ) = R1
!
            IF( KSTEP == 2 ) THEN
!                 (*) Make sure that diagonal element of pivot is real
               A( K, K ) = REAL( A( K, K ) )
!                 (5) Swap row elements
               T = A( K-1, K )
               A( K-1, K ) = A( KP, K )
               A( KP, K ) = T
            END IF
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
            IF( K < N ) THEN
              AT_TMP(1:N-K) = A(KK,K+1:N)
              A(KK,K+1:N) = A(KP,K+1:N)
              A(KP,K+1:N) = AT_TMP(1:N-K)
            ENDIF
!
         ELSE
!              (*) Make sure that diagonal element of pivot is real
            A( K, K ) = REAL( A( K, K ) )
            IF( KSTEP == 2 ) A( K-1, K-1 ) = REAL( A( K-1, K-1 ) )
         END IF
!
!           Update the leading submatrix
!
         IF( KSTEP == 1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
            IF( K > 1 ) THEN
!
!                 Perform a rank-1 update of A(1:k-1,1:k-1) and
!                 store U(k) in column k
!
               IF( ABS( REAL( A( K, K ) ) ) >= SFMIN ) THEN
!
!                    Perform a rank-1 update of A(1:k-1,1:k-1) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*1/D(k)*W(k)**T
!
                  D11 = 1.0E+0 / REAL( A( K, K ) )
                  CALL CHER( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )
!
!                    Store U(k) in column k
!
                  A(1:K-1,K) = D11*A(1:K-1,K)
               ELSE
!
!                    Store L(k) in column K
!
                  D11 = REAL( A( K, K ) )
                  A(1:K-1,K) = A(1:K-1,K) / D11
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                  CALL CHER( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )
               END IF
!
!                 Store the superdiagonal element of D in array E
!
               E( K ) = (0.0E+0,0.0E+0)
!
            END IF
!
         ELSE
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
!                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T
!
!              and store L(k) and L(k+1) in columns k and k+1
!
            IF( K > 2 ) THEN
!                 D = |A12|
               D = SLAPY2( REAL( A( K-1, K ) ), AIMAG( A( K-1, K ) ) )
               D11 = REAL( A( K, K ) / D )
               D22 = REAL( A( K-1, K-1 ) / D )
               D12 = A( K-1, K ) / D
               TT = 1.0E+0 / ( D11*D22-1.0E+0 )
!
               DO J = K - 2, 1, -1
!
!                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
!
                  WKM1 = TT*( D11*A( J, K-1 )-CONJG( D12 )*A( J, K ) )
                  WK = TT*( D22*A( J, K )-D12*A( J, K-1 ) )
!
!                    Perform a rank-2 update of A(1:k-2,1:k-2)
!
                  DO I = J, 1, -1
                     A( I, J ) = A( I, J ) - &
                                 ( A( I, K ) / D )*CONJG( WK ) - &
                                 ( A( I, K-1 ) / D )*CONJG( WKM1 )
                  ENDDO
!
!                    Store U(k) and U(k-1) in cols k and k-1 for row J
!
                  A( J, K ) = WK / D
                  A( J, K-1 ) = WKM1 / D
!                    (*) Make sure that diagonal element of pivot is real
                  A( J, J ) = CMPLX( REAL( A( J, J ) ), 0.0E+0 )
!
               ENDDO
!
            END IF
!
!              Copy superdiagonal elements of D(K) to E(K) and
!              0.0E+0 out superdiagonal entry of A
!
            E( K ) = A( K-1, K )
            E( K-1 ) = (0.0E+0,0.0E+0)
            A( K-1, K ) = (0.0E+0,0.0E+0)
!
         END IF
!
!           End column K is nonsingular
!
      END IF
!
!        Store details of the interchanges in IPIV
!
      IF( KSTEP == 1 ) THEN
         IPIV( K ) = KP
      ELSE
         IPIV( K ) = -P
         IPIV( K-1 ) = -KP
      END IF
!
!        Decrease K and return to the start of the main loop
!
      K = K - KSTEP
      GO TO 10
!
34    CONTINUE
!
   ELSE
!
!        Factorize A as L*D*L**H using the lower triangle of A
!
!        Initialize the unused last entry of the subdiagonal array E.
!
      E( N ) = (0.0E+0,0.0E+0)
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
      K = 1
40    CONTINUE
!
!        If K > N, exit from loop
!
      IF( K > N ) GO TO 64
      KSTEP = 1
      P = K
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
      ABSAKK = ABS( REAL( A( K, K ) ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
      IF( K < N ) THEN
         IMAX = K + ICAMAX( N-K, A( K+1, K ), 1 )
         COLMAX = CABS1( A( IMAX, K ) )
      ELSE
         COLMAX = 0.0E+0
      END IF
!
      IF( MAX( ABSAKK, COLMAX ) == 0.0E+0 ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
         IF( INFO == 0 ) INFO = K
         KP = K
         A( K, K ) = REAL( A( K, K ) )
!
!           Set E( K ) to zero
!
         IF( K < N ) E( K ) = (0.0E+0,0.0E+0)
!
      ELSE
!
!           ============================================================
!
!           BEGIN pivot search
!
!           Case(1)
!           Equivalent to testing for ABSAKK >= ALPHA*COLMAX
!           (used to handle NaN and Inf)
!
         IF( .NOT.( ABSAKK < ALPHA*COLMAX ) ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
            KP = K
!
         ELSE
!
            DONE = .FALSE.
!
!              Loop until pivot found
!
42          CONTINUE
!
!                 BEGIN pivot search loop body
!
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
               IF( IMAX /= K ) THEN
                  JMAX = K - 1 + ICAMAX( IMAX-K, A( IMAX, K ), LDA )
                  ROWMAX = CABS1( A( IMAX, JMAX ) )
               ELSE
                  ROWMAX = 0.0E+0
               END IF
!
               IF( IMAX < N ) THEN
                  ITEMP = IMAX + ICAMAX( N-IMAX, A( IMAX+1, IMAX ), &
                                        1 )
                  STEMP = CABS1( A( ITEMP, IMAX ) )
                  IF( STEMP > ROWMAX ) THEN
                     ROWMAX = STEMP
                     JMAX = ITEMP
                  END IF
               END IF
!
!                 Case(2)
!                 Equivalent to testing for
!                 ABS( REAL( W( IMAX,KW-1 ) ) ) >= ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
               IF( .NOT.( ABS( REAL( A( IMAX, IMAX ) ) ) < ALPHA*ROWMAX ) ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                  KP = IMAX
                  DONE = .TRUE.
!
!                 Case(3)
!                 Equivalent to testing for ROWMAX == COLMAX,
!                 (used to handle NaN and Inf)
!
               ELSE IF( ( P == JMAX ) .OR. ( ROWMAX <= COLMAX ) ) THEN
!
!                    interchange rows and columns K+1 and IMAX,
!                    use 2-by-2 pivot block
!
                  KP = IMAX
                  KSTEP = 2
                  DONE = .TRUE.
!
!                 Case(4)
               ELSE
!
!                    Pivot not found: set params and repeat
!
                  P = IMAX
                  COLMAX = ROWMAX
                  IMAX = JMAX
               END IF
!
!
!                 END pivot search loop body
!
            IF( .NOT.DONE ) GOTO 42
!
         END IF
!
!           END pivot search
!
!           ============================================================
!
!           KK is the column of A where pivoting step stopped
!
         KK = K + KSTEP - 1
!
!           For only a 2x2 pivot, interchange rows and columns K and P
!           in the trailing submatrix A(k:n,k:n)
!
         IF( ( KSTEP == 2 ) .AND. ( P /= K ) ) THEN
!              (1) Swap columnar parts
            IF( P < N ) THEN
              A_TMP(1:N-P) = A(P+1:N,K)
              A(P+1:N,K) = A(P+1:N,P)
              A(P+1:N,P) = A_TMP(1:N-P)
            ENDIF
!              (2) Swap and conjugate middle parts
            A_TMP(1:P-K-1) = CONJG(A(K+1:P-1,K))
            A(K+1:P-1,K) = CONJG(A(P,K+1:P-1))
            A(P,K+1:P-1) = A_TMP(1:P-K-1)
!              (3) Swap and conjugate corner elements at row-col intersection
            A( P, K ) = CONJG( A( P, K ) )
!              (4) Swap diagonal elements at row-col intersection
            R1 = REAL( A( K, K ) )
            A( K, K ) = REAL( A( P, P ) )
            A( P, P ) = R1
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
            IF ( K > 1 ) THEN
               AT_TMP(1:K-1) = A(K,1:K-1)
               A(K,1:K-1) = A(P,1:K-1)
               A(P,1:K-1) = AT_TMP(1:K-1)
            ENDIF
!
         END IF
!
!           For both 1x1 and 2x2 pivots, interchange rows and
!           columns KK and KP in the trailing submatrix A(k:n,k:n)
!
         IF( KP /= KK ) THEN
!              (1) Swap columnar parts
            IF( KP < N ) THEN
               A_TMP(1:N-KP) = A(KP+1:N,KK)
               A(KP+1:N,KK) = A(KP+1:N,KP)
               A(KP+1:N,KP) = A_TMP(1:N-KP)
            ENDIF
!              (2) Swap and conjugate middle parts
            A_TMP(1:KP-KK-1) = CONJG(A(KK+1:KP-1,KK))
            A(KK+1:KP-1,KK) = CONJG(A(KP,KK+1:KP-1))
            A(KP,KK+1:KP-1) = A_TMP(1:KP-KK-1)
!              (3) Swap and conjugate corner elements at row-col intersection
            A( KP, KK ) = CONJG( A( KP, KK ) )
!              (4) Swap diagonal elements at row-col intersection
            R1 = REAL( A( KK, KK ) )
            A( KK, KK ) = REAL( A( KP, KP ) )
            A( KP, KP ) = R1
!
            IF( KSTEP == 2 ) THEN
!                 (*) Make sure that diagonal element of pivot is real
               A( K, K ) = REAL( A( K, K ) )
!                 (5) Swap row elements
               T = A( K+1, K )
               A( K+1, K ) = A( KP, K )
               A( KP, K ) = T
            END IF
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
            IF ( K > 1 ) THEN
               AT_TMP(1:K-1) = A(KK,1:K-1)
               A(KK,1:K-1) = A(KP,1:K-1)
               A(KP,1:K-1) = AT_TMP(1:K-1)
            ENDIF
!
         ELSE
!              (*) Make sure that diagonal element of pivot is real
            A( K, K ) = REAL( A( K, K ) )
            IF( KSTEP == 2 ) &
               A( K+1, K+1 ) = REAL( A( K+1, K+1 ) )
         END IF
!
!           Update the trailing submatrix
!
         IF( KSTEP == 1 ) THEN
!
!              1-by-1 pivot block D(k): column k of A now holds
!
!              W(k) = L(k)*D(k),
!
!              where L(k) is the k-th column of L
!
            IF( K < N ) THEN
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) and
!                 store L(k) in column k
!
!                 Handle division by a small number
!
               IF( ABS( REAL( A( K, K ) ) ) >= SFMIN ) THEN
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!
                  D11 = 1.0E+0 / REAL( A( K, K ) )
                  CALL CHER( UPLO, N-K, -D11, A( K+1, K ), 1, &
                             A( K+1, K+1 ), LDA )
!
!                    Store L(k) in column k
!
                  A(K+1:N,K) = D11*A(K+1:N,K)
               ELSE
!
!                    Store L(k) in column k
!
                  D11 = REAL( A( K, K ) )
                  A(K+1:N,K) = A(K+1:N,K)/D11
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                  CALL CHER( UPLO, N-K, -D11, A( K+1, K ), 1, &
                             A( K+1, K+1 ), LDA )
               END IF
!
!                 Store the subdiagonal element of D in array E
!
               E( K ) = (0.0E+0,0.0E+0)
!
            END IF
!
         ELSE
!
!              2-by-2 pivot block D(k): columns k and k+1 now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
!
!              Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
!                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T
!
!              and store L(k) and L(k+1) in columns k and k+1
!
            IF( K < N-1 ) THEN
!                 D = |A21|
               D = SLAPY2( REAL( A( K+1, K ) ), &
                   AIMAG( A( K+1, K ) ) )
               D11 = REAL( A( K+1, K+1 ) ) / D
               D22 = REAL( A( K, K ) ) / D
               D21 = A( K+1, K ) / D
               TT = 1.0E+0 / ( D11*D22-1.0E+0 )
!
               DO J = K + 2, N
!
!                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
!
                  WK = TT*( D11*A( J, K )-D21*A( J, K+1 ) )
                  WKP1 = TT*( D22*A( J, K+1 )-CONJG( D21 )* &
                         A( J, K ) )
!
!                    Perform a rank-2 update of A(k+2:n,k+2:n)
!
                  A(J:N,J) = A(J:N,J) - (A(J:N,K)/D)*CONJG(WK) - &
                              (A(J:N,K+1)/D )*CONJG(WKP1)
!
!                    Store L(k) and L(k+1) in cols k and k+1 for row J
!
                  A( J, K ) = WK / D
                  A( J, K+1 ) = WKP1 / D
!                    (*) Make sure that diagonal element of pivot is real
                  A( J, J ) = CMPLX( REAL( A( J, J ) ), 0.0E+0 )
!
               ENDDO
!
            END IF
!
!              Copy subdiagonal elements of D(K) to E(K) and
!              0.0E+0 out subdiagonal entry of A
!
            E( K ) = A( K+1, K )
            E( K+1 ) = (0.0E+0,0.0E+0)
            A( K+1, K ) = (0.0E+0,0.0E+0)
!
         END IF
!
!           End column K is nonsingular
!
      END IF
!
!        Store details of the interchanges in IPIV
!
      IF( KSTEP == 1 ) THEN
         IPIV( K ) = KP
      ELSE
         IPIV( K ) = -P
         IPIV( K+1 ) = -KP
      END IF
!
!        Increase K and return to the start of the main loop
!
      K = K + KSTEP
      GO TO 40
!
64    CONTINUE
!
   END IF
!
   RETURN
!
!     End of CHETF2_RK
!
END

