!> \brief \b SSYTF2_RK computes the factorization of a real symmetric indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYTF2_RK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytf2_rk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytf2_rk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytf2_rk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), E ( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SSYTF2_RK computes the factorization of a real symmetric matrix A
!> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
!>
!>    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**T (or L**T) is the transpose of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is symmetric and block
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
!>          symmetric matrix A is stored:
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.
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
!>            a) ONLY diagonal elements of the symmetric block diagonal
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
!>          E is REAL array, dimension (N)
!>          On exit, contains the superdiagonal (or subdiagonal)
!>          elements of the symmetric block diagonal matrix D
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
!>          the symmetric block diagonal matrix D with 1-by-1 or 2-by-2
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
   SUBROUTINE SSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )
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
   REAL               A( LDA, * ), E( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
   REAL               EIGHT, SEVTEN
   PARAMETER          ( EIGHT = 8.0E+0, SEVTEN = 17.0E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER, DONE
   INTEGER            I, IMAX, J, JMAX, ITEMP, K, KK, KP, KSTEP, &
                      P, II
   REAL               ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, &
                      ROWMAX, STEMP, T, WK, WKM1, WKP1, SFMIN
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ISAMAX
   REAL               SLAMCH
   EXTERNAL           LSAME, ISAMAX, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SSCAL, SSWAP, SSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, SQRT
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
      CALL XERBLA( 'SSYTF2_RK', -INFO )
      RETURN
   END IF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
   ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
!     Compute machine safe minimum
!
   SFMIN = SLAMCH( 'S' )
!
   IF( UPPER ) THEN
!
!        Factorize A as U*D*U**T using the upper triangle of A
!
!        Initialize the first entry of array E, where superdiagonal
!        elements of D are stored
!
      E( 1 ) = ZERO
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
      K = N
10    CONTINUE
!
!        If K < 1, exit from loop
!
      IF( K < 1 ) &
         GO TO 34
      KSTEP = 1
      P = K
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
      ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
      IF( K > 1 ) THEN
         IMAX = ISAMAX( K-1, A( 1, K ), 1 )
         COLMAX = ABS( A( IMAX, K ) )
      ELSE
         COLMAX = ZERO
      END IF
!
      IF( (MAX( ABSAKK, COLMAX ) == ZERO) ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
         IF( INFO == 0 ) &
            INFO = K
         KP = K
!
!           Set E( K ) to zero
!
         IF( K > 1 ) &
            E( K ) = ZERO
!
      ELSE
!
!           Test for interchange
!
!           Equivalent to testing for (used to handle NaN and Inf)
!           ABSAKK >= ALPHA*COLMAX
!
         IF( .NOT.( ABSAKK < ALPHA*COLMAX ) ) THEN
!
!              no interchange,
!              use 1-by-1 pivot block
!
            KP = K
         ELSE
!
            DONE = .FALSE.
!
!              Loop until pivot found
!
12          CONTINUE
!
!                 Begin pivot search loop body
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
               IF( IMAX /= K ) THEN
                  JMAX = IMAX + ISAMAX( K-IMAX, A( IMAX, IMAX+1 ), &
                                       LDA )
                  ROWMAX = ABS( A( IMAX, JMAX ) )
               ELSE
                  ROWMAX = ZERO
               END IF
!
               IF( IMAX > 1 ) THEN
                  ITEMP = ISAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  STEMP = ABS( A( ITEMP, IMAX ) )
                  IF( STEMP > ROWMAX ) THEN
                     ROWMAX = STEMP
                     JMAX = ITEMP
                  END IF
               END IF
!
!                 Equivalent to testing for (used to handle NaN and Inf)
!                 ABS( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX
!
               IF( .NOT.( ABS( A( IMAX, IMAX ) ) < ALPHA*ROWMAX ) ) &
               THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                  KP = IMAX
                  DONE = .TRUE.
!
!                 Equivalent to testing for ROWMAX  ==  COLMAX,
!                 used to handle NaN and Inf
!
               ELSE IF( ( P == JMAX ).OR.( ROWMAX <= COLMAX ) ) THEN
!
!                    interchange rows and columns K+1 and IMAX,
!                    use 2-by-2 pivot block
!
                  KP = IMAX
                  KSTEP = 2
                  DONE = .TRUE.
               ELSE
!
!                    Pivot NOT found, set variables and repeat
!
                  P = IMAX
                  COLMAX = ROWMAX
                  IMAX = JMAX
               END IF
!
!                 End pivot search loop body
!
            IF( .NOT. DONE ) GOTO 12
!
         END IF
!
!           Swap TWO rows and TWO columns
!
!           First swap
!
         IF( ( KSTEP == 2 ) .AND. ( P /= K ) ) THEN
!
!              Interchange rows and column K and P in the leading
!              submatrix A(1:k,1:k) if we have a 2-by-2 pivot
!
            IF( P > 1 ) &
               CALL SSWAP( P-1, A( 1, K ), 1, A( 1, P ), 1 )
            IF( P < (K-1) ) &
               CALL SSWAP( K-P-1, A( P+1, K ), 1, A( P, P+1 ), &
                        LDA )
            T = A( K, K )
            A( K, K ) = A( P, P )
            A( P, P ) = T
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
            IF( K < N ) &
               CALL SSWAP( N-K, A( K, K+1 ), LDA, A( P, K+1 ), LDA )
!
         END IF
!
!           Second swap
!
         KK = K - KSTEP + 1
         IF( KP /= KK ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
            IF( KP > 1 ) &
               CALL SSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
            IF( ( KK > 1 ) .AND. ( KP < (KK-1) ) ) &
               CALL SSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), &
                        LDA )
            T = A( KK, KK )
            A( KK, KK ) = A( KP, KP )
            A( KP, KP ) = T
            IF( KSTEP == 2 ) THEN
               T = A( K-1, K )
               A( K-1, K ) = A( KP, K )
               A( KP, K ) = T
            END IF
!
!              Convert upper triangle of A into U form by applying
!              the interchanges in columns k+1:N.
!
            IF( K < N ) &
               CALL SSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), &
                           LDA )
!
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
               IF( ABS( A( K, K ) ) >= SFMIN ) THEN
!
!                    Perform a rank-1 update of A(1:k-1,1:k-1) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*1/D(k)*W(k)**T
!
                  D11 = ONE / A( K, K )
                  CALL SSYR( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )
!
!                    Store U(k) in column k
!
                  CALL SSCAL( K-1, D11, A( 1, K ), 1 )
               ELSE
!
!                    Store L(k) in column K
!
                  D11 = A( K, K )
                  DO II = 1, K - 1
                     A( II, K ) = A( II, K ) / D11
                  ENDDO
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - U(k)*D(k)*U(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                  CALL SSYR( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )
               END IF
!
!                 Store the superdiagonal element of D in array E
!
               E( K ) = ZERO
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
!
               D12 = A( K-1, K )
               D22 = A( K-1, K-1 ) / D12
               D11 = A( K, K ) / D12
               T = ONE / ( D11*D22-ONE )
!
               DO J = K - 2, 1, -1
!
                  WKM1 = T*( D11*A( J, K-1 )-A( J, K ) )
                  WK = T*( D22*A( J, K )-A( J, K-1 ) )
!
                  DO I = J, 1, -1
                     A( I, J ) = A( I, J ) - (A( I, K ) / D12 )*WK - &
                                 ( A( I, K-1 ) / D12 )*WKM1
                  ENDDO
!
!                    Store U(k) and U(k-1) in cols k and k-1 for row J
!
                  A( J, K ) = WK / D12
                  A( J, K-1 ) = WKM1 / D12
!
               ENDDO
!
            END IF
!
!              Copy superdiagonal elements of D(K) to E(K) and
!              ZERO out superdiagonal entry of A
!
            E( K ) = A( K-1, K )
            E( K-1 ) = ZERO
            A( K-1, K ) = ZERO
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
!        Factorize A as L*D*L**T using the lower triangle of A
!
!        Initialize the unused last entry of the subdiagonal array E.
!
      E( N ) = ZERO
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
      K = 1
40    CONTINUE
!
!        If K > N, exit from loop
!
      IF( K > N ) &
         GO TO 64
      KSTEP = 1
      P = K
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
      ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
      IF( K < N ) THEN
         IMAX = K + ISAMAX( N-K, A( K+1, K ), 1 )
         COLMAX = ABS( A( IMAX, K ) )
      ELSE
         COLMAX = ZERO
      END IF
!
      IF( ( MAX( ABSAKK, COLMAX ) == ZERO ) ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
         IF( INFO == 0 ) &
            INFO = K
         KP = K
!
!           Set E( K ) to zero
!
         IF( K < N ) &
            E( K ) = ZERO
!
      ELSE
!
!           Test for interchange
!
!           Equivalent to testing for (used to handle NaN and Inf)
!           ABSAKK >= ALPHA*COLMAX
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
!                 Begin pivot search loop body
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
               IF( IMAX /= K ) THEN
                  JMAX = K - 1 + ISAMAX( IMAX-K, A( IMAX, K ), LDA )
                  ROWMAX = ABS( A( IMAX, JMAX ) )
               ELSE
                  ROWMAX = ZERO
               END IF
!
               IF( IMAX < N ) THEN
                  ITEMP = IMAX + ISAMAX( N-IMAX, A( IMAX+1, IMAX ), &
                                        1 )
                  STEMP = ABS( A( ITEMP, IMAX ) )
                  IF( STEMP > ROWMAX ) THEN
                     ROWMAX = STEMP
                     JMAX = ITEMP
                  END IF
               END IF
!
!                 Equivalent to testing for (used to handle NaN and Inf)
!                 ABS( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX
!
               IF( .NOT.( ABS( A( IMAX, IMAX ) ) < ALPHA*ROWMAX ) ) &
               THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                  KP = IMAX
                  DONE = .TRUE.
!
!                 Equivalent to testing for ROWMAX  ==  COLMAX,
!                 used to handle NaN and Inf
!
               ELSE IF( ( P == JMAX ).OR.( ROWMAX <= COLMAX ) ) THEN
!
!                    interchange rows and columns K+1 and IMAX,
!                    use 2-by-2 pivot block
!
                  KP = IMAX
                  KSTEP = 2
                  DONE = .TRUE.
               ELSE
!
!                    Pivot NOT found, set variables and repeat
!
                  P = IMAX
                  COLMAX = ROWMAX
                  IMAX = JMAX
               END IF
!
!                 End pivot search loop body
!
            IF( .NOT. DONE ) GOTO 42
!
         END IF
!
!           Swap TWO rows and TWO columns
!
!           First swap
!
         IF( ( KSTEP == 2 ) .AND. ( P /= K ) ) THEN
!
!              Interchange rows and column K and P in the trailing
!              submatrix A(k:n,k:n) if we have a 2-by-2 pivot
!
            IF( P < N ) &
               CALL SSWAP( N-P, A( P+1, K ), 1, A( P+1, P ), 1 )
            IF( P > (K+1) ) &
               CALL SSWAP( P-K-1, A( K+1, K ), 1, A( P, K+1 ), LDA )
            T = A( K, K )
            A( K, K ) = A( P, P )
            A( P, P ) = T
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
            IF ( K > 1 ) &
               CALL SSWAP( K-1, A( K, 1 ), LDA, A( P, 1 ), LDA )
!
         END IF
!
!           Second swap
!
         KK = K + KSTEP - 1
         IF( KP /= KK ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
            IF( KP < N ) &
               CALL SSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
            IF( ( KK < N ) .AND. ( KP > (KK+1) ) ) &
               CALL SSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), &
                        LDA )
            T = A( KK, KK )
            A( KK, KK ) = A( KP, KP )
            A( KP, KP ) = T
            IF( KSTEP == 2 ) THEN
               T = A( K+1, K )
               A( K+1, K ) = A( KP, K )
               A( KP, K ) = T
            END IF
!
!              Convert lower triangle of A into L form by applying
!              the interchanges in columns 1:k-1.
!
            IF ( K > 1 ) &
               CALL SSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
!
         END IF
!
!           Update the trailing submatrix
!
         IF( KSTEP == 1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
            IF( K < N ) THEN
!
!              Perform a rank-1 update of A(k+1:n,k+1:n) and
!              store L(k) in column k
!
               IF( ABS( A( K, K ) ) >= SFMIN ) THEN
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!
                  D11 = ONE / A( K, K )
                  CALL SSYR( UPLO, N-K, -D11, A( K+1, K ), 1, &
                             A( K+1, K+1 ), LDA )
!
!                    Store L(k) in column k
!
                  CALL SSCAL( N-K, D11, A( K+1, K ), 1 )
               ELSE
!
!                    Store L(k) in column k
!
                  D11 = A( K, K )
                  DO II = K + 1, N
                     A( II, K ) = A( II, K ) / D11
                  ENDDO
!
!                    Perform a rank-1 update of A(k+1:n,k+1:n) as
!                    A := A - L(k)*D(k)*L(k)**T
!                       = A - W(k)*(1/D(k))*W(k)**T
!                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
!
                  CALL SSYR( UPLO, N-K, -D11, A( K+1, K ), 1, &
                             A( K+1, K+1 ), LDA )
               END IF
!
!                 Store the subdiagonal element of D in array E
!
               E( K ) = ZERO
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
!
               D21 = A( K+1, K )
               D11 = A( K+1, K+1 ) / D21
               D22 = A( K, K ) / D21
               T = ONE / ( D11*D22-ONE )
!
               DO J = K + 2, N
!
!                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
!
                  WK = T*( D11*A( J, K )-A( J, K+1 ) )
                  WKP1 = T*( D22*A( J, K+1 )-A( J, K ) )
!
!                    Perform a rank-2 update of A(k+2:n,k+2:n)
!
                  DO I = J, N
                     A( I, J ) = A( I, J ) - ( A( I, K ) / D21 )*WK - &
                                 ( A( I, K+1 ) / D21 )*WKP1
                  ENDDO
!
!                    Store L(k) and L(k+1) in cols k and k+1 for row J
!
                  A( J, K ) = WK / D21
                  A( J, K+1 ) = WKP1 / D21
!
               ENDDO
!
            END IF
!
!              Copy subdiagonal elements of D(K) to E(K) and
!              ZERO out subdiagonal entry of A
!
            E( K ) = A( K+1, K )
            E( K+1 ) = ZERO
            A( K+1, K ) = ZERO
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
!     End of SSYTF2_RK
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        