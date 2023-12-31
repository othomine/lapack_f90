!> \brief \b CHETF2_ROOK computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETF2_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetf2_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetf2_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetf2_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETF2_ROOK( UPLO, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETF2_ROOK computes the factorization of a complex Hermitian matrix A
!> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:
!>
!>    A = U*D*U**H  or  A = L*D*L**H
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, U**H is the conjugate transpose of U, and D is
!> Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
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
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L (see below for further details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D.
!>
!>          If UPLO = 'U':
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k-1 and -IPIV(k-1) were inerchaged,
!>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k+1 and -IPIV(k+1) were inerchaged,
!>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular, and division by zero will occur if it
!>               is used to solve a system of equations.
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
!> \ingroup hetf2_rook
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', then A = U*D*U**H, where
!>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    v    0   )   k-s
!>     U(k) =  (   0    I    0   )   s
!>             (   0    0    I   )   n-k
!>                k-s   s   n-k
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!>
!>  If UPLO = 'L', then A = L*D*L**H, where
!>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    0     0   )  k-1
!>     L(k) =  (   0    I     0   )  s
!>             (   0    v     I   )  n-k-s+1
!>                k-1   s  n-k-s+1
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!>  01-01-96 - Based on modifications by
!>    J. Lewis, Boeing Computer Services Company
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE CHETF2_ROOK( UPLO, N, A, LDA, IPIV, INFO )
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
   COMPLEX            A( LDA, * )
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
   COMPLEX            D12, D21, T, WK, WKM1, WKP1, Z, A_TMP( LDA )
!     ..
!     .. External Functions ..
!
   LOGICAL            LSAME
   INTEGER            ICAMAX
   REAL               SLAMCH, SLAPY2, CABS1
   EXTERNAL           LSAME, ICAMAX, SLAMCH, SLAPY2, CABS1
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, CSSCAL, CHER, CSWAP
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
      CALL XERBLA( 'CHETF2_ROOK', -INFO )
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
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
      K = N
10    CONTINUE
!
!        If K < 1, exit from loop
!
      IF( K < 1 ) GO TO 70
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
         END IF
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
   ELSE
!
!        Factorize A as L*D*L**H using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
      K = 1
40    CONTINUE
!
!        If K > N, exit from loop
!
      IF( K > N ) GO TO 70
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
         END IF
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
   END IF
!
70 CONTINUE
!
   RETURN
!
!     End of CHETF2_ROOK
!
END

