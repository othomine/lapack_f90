! \brief \b CLAHEF_ROOK computes a partial factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (blocked algorithm, calling Level 3 BLAS).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAHEF_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHEF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KB, LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAHEF_ROOK computes a partial factorization of a complex Hermitian
!> matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting
!> method. The partial factorization has the form:
!>
!> A  =  ( I  U12 ) ( A11  0  ) (  I      0     )  if UPLO = 'U', or:
!>       ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
!>
!> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L'
!>       ( L21  I ) (  0  A22 ) (  0      I     )
!>
!> where the order of D is at most NB. The actual order is returned in
!> the argument KB, and is either NB or NB-1, or N if N <= NB.
!> Note that U**H denotes the conjugate transpose of U.
!>
!> CLAHEF_ROOK is an auxiliary routine called by CHETRF_ROOK. It uses
!> blocked code (calling Level 3 BLAS) to update the submatrix
!> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').
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
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The maximum number of columns of the matrix A that should be
!>          factored.  NB should be at least 2 to allow for 2-by-2 pivot
!>          blocks.
!> \endverbatim
!>
!> \param[out] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of columns of A that were actually factored.
!>          KB is either NB-1 or NB, or N if N <= NB.
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
!>          On exit, A contains details of the partial factorization.
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
!>             Only the last KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k-1 and -IPIV(k-1) were inerchaged,
!>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             Only the first KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k+1 and -IPIV(k+1) were inerchaged,
!>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (LDW,NB)
!> \endverbatim
!>
!> \param[in] LDW
!> \verbatim
!>          LDW is INTEGER
!>          The leading dimension of the array W.  LDW >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular.
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
!> \ingroup lahef_rook
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013, Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE CLAHEF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, &
                           INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, KB, LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   COMPLEX            A( LDA, * ), W( LDW, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            DONE
   INTEGER            IMAX, ITEMP, II, J, JB, JJ, JMAX, JP1, JP2, K, &
                      KK, KKW, KP, KSTEP, KW, P
   REAL               ABSAKK, COLMAX, STEMP, R1, ROWMAX, T, &
                      SFMIN
!
!     Initialize ALPHA for use in choosing pivot block size.
!
   REAL, PARAMETER :: ALPHA = ( 1.0E+0+SQRT( 17.0E+0 ) ) / 8.0E+0
   COMPLEX            D11, D21, D22, Z, A_TMP( N ), W_TMP( NB )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ICAMAX
   REAL               SLAMCH
   EXTERNAL           LSAME, ICAMAX, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CGEMV, CLACGV
!     ..
!     .. Statement Functions ..
   REAL               CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( Z ) = ABS( REAL( Z ) ) + ABS( AIMAG( Z ) )
!     ..
!     .. Executable Statements ..
!
   INFO = 0
!
!     Compute machine safe minimum
!
   SFMIN = SLAMCH( 'S' )
!
   IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11 (note that conjg(W) is actually stored)
!
!        K is the main loop index, decreasing from N in steps of 1 or 2
!
      K = N
10    CONTINUE
!
!        KW is the column of W which corresponds to column K of A
!
      KW = NB + K - N
!
!        Exit from loop
!
      IF( ( K <= N-NB+1 .AND. NB < N ) .OR. K < 1 ) GO TO 30
!
      KSTEP = 1
      P = K
!
!        Copy column K of A to column KW of W and update it
!
      IF( K > 1 ) W(1:K-1,KW) = A(1:K-1,K)
      W( K, KW ) = REAL( A( K, K ) )
      IF( K < N ) THEN
         CALL CGEMV( 'No transpose', K, N-K, -(1.0E+0,0.0E+0), A( 1, K+1 ), LDA, &
                     W( K, KW+1 ), LDW, (1.0E+0,0.0E+0), W( 1, KW ), 1 )
         W( K, KW ) = REAL( W( K, KW ) )
      END IF
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
      ABSAKK = ABS( REAL( W( K, KW ) ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
      IF( K > 1 ) THEN
         IMAX = ICAMAX( K-1, W( 1, KW ), 1 )
         COLMAX = CABS1( W( IMAX, KW ) )
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
         A( K, K ) = REAL( W( K, KW ) )
         IF( K > 1 ) A(1:K-1,K) = W(1:K-1,KW)
      ELSE
!
!           ============================================================
!
!           BEGIN pivot search
!
!           Case(1)
!           Equivalent to testing for ABSAKK >= ALPHA*COLMAX
!           (used to handle NaN and Inf)
         IF( .NOT.( ABSAKK < ALPHA*COLMAX ) ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
            KP = K
!
         ELSE
!
!              Lop until pivot found
!
            DONE = .FALSE.
!
12          CONTINUE
!
!                 BEGIN pivot search loop body
!
!
!                 Copy column IMAX to column KW-1 of W and update it
!
               IF( IMAX > 1 ) W(1:IMAX-1,KW-1) = A(1:IMAX-1,IMAX)
               W( IMAX, KW-1 ) = REAL( A( IMAX, IMAX ) )
!
               W(IMAX+1:K,KW-1) = A(IMAX,IMAX+1:K)
               CALL CLACGV( K-IMAX, W( IMAX+1, KW-1 ), 1 )
!
               IF( K < N ) THEN
                  CALL CGEMV( 'No transpose', K, N-K, -(1.0E+0,0.0E+0), &
                              A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW, &
                              (1.0E+0,0.0E+0), W( 1, KW-1 ), 1 )
                  W( IMAX, KW-1 ) = REAL( W( IMAX, KW-1 ) )
               END IF
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
               IF( IMAX /= K ) THEN
                  JMAX = IMAX + ICAMAX( K-IMAX, W( IMAX+1, KW-1 ), &
                                        1 )
                  ROWMAX = CABS1( W( JMAX, KW-1 ) )
               ELSE
                  ROWMAX = 0.0E+0
               END IF
!
               IF( IMAX > 1 ) THEN
                  ITEMP = ICAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  STEMP = CABS1( W( ITEMP, KW-1 ) )
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
               IF( .NOT.( ABS( REAL( W( IMAX,KW-1 ) ) ) < ALPHA*ROWMAX ) ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                  KP = IMAX
!
!                    copy column KW-1 of W to column KW of W
!
                  W(1:K,KW) = W(1:K,KW-1)
!
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
!
!                    Copy updated JMAXth (next IMAXth) column to Kth of W
!
                  W(1:K,KW) = W(1:K,KW-1)
!
               END IF
!
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
!           KKW is the column of W which corresponds to column KK of A
!
         KKW = NB + KK - N
!
!           Interchange rows and columns P and K.
!           Updated column P is already stored in column KW of W.
!
         IF( ( KSTEP == 2 ) .AND. ( P /= K ) ) THEN
!
!              Copy non-updated column K to column P of submatrix A
!              at step K. No need to copy element into columns
!              K and K-1 of A for 2-by-2 pivot, since these columns
!              will be later overwritten.
!
            A( P, P ) = REAL( A( K, K ) )
            A(P,P+1:K-1) = A(P+1:K-1,K)
            CALL CLACGV( K-1-P, A( P, P+1 ), LDA )
            IF( P > 1 ) A(1:P-1,P) = A(1:P-1,K)
!
!              Interchange rows K and P in the last K+1 to N columns of A
!              (columns K and K-1 of A for 2-by-2 pivot will be
!              later overwritten). Interchange rows K and P
!              in last KKW to NB columns of W.
!
            IF( K < N ) THEN
              A_TMP(1:N-K) = A(K,K+1:N)
              A(K,K+1:N) = A(P,K+1:N)
              A(P,K+1:N) = A_TMP(1:N-K)
            ENDIF
            W_TMP(1:N-KK+1) = W(K,KKW:KKW+N-KK)
            W(K,KKW:KKW+N-KK) = W(P,KKW:KKW+N-KK)
            W(P,KKW:KKW+N-KK) = W_TMP(1:N-KK+1)
         END IF
!
!           Interchange rows and columns KP and KK.
!           Updated column KP is already stored in column KKW of W.
!
         IF( KP /= KK ) THEN
!
!              Copy non-updated column KK to column KP of submatrix A
!              at step K. No need to copy element into column K
!              (or K and K-1 for 2-by-2 pivot) of A, since these columns
!              will be later overwritten.
!
            A( KP, KP ) = REAL( A( KK, KK ) )
            A_TMP(1:KK-1-KP) = A(KP+1:KK-1,KK)
            A(KP+1:KK-1,KK) = A(KP,KP+1:KK-1)
            A(KP,KP+1:KK-1) = A_TMP(1:KK-1-KP)
            CALL CLACGV( KK-1-KP, A( KP, KP+1 ), LDA )
            IF( KP > 1 ) A(1:KP-1,KP) = A(1:KP-1,KK)
!
!              Interchange rows KK and KP in last K+1 to N columns of A
!              (columns K (or K and K-1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in last KKW to NB columns of W.
!
            IF( K < N ) THEN
              A_TMP(1:N-K) = A(KK,K+1:N)
              A(KK,K+1:N) = A(KP,K+1:N)
              A(KP,K+1:N) = A_TMP(1:N-K)
            ENDIF
            W_TMP(1:N-KK+1) = W(KK,KKW:KKW+N-KK)
            W(KK,KKW:KKW+N-KK) = W(KP,KKW:KKW+N-KK)
            W(KP,KKW:KKW+N-KK) = W_TMP(1:N-KK+1)
         END IF
!
         IF( KSTEP == 1 ) THEN
!
!              1-by-1 pivot block D(k): column kw of W now holds
!
!              W(kw) = U(k)*D(k),
!
!              where U(k) is the k-th column of U
!
!              (1) Store subdiag. elements of column U(k)
!              and 1-by-1 block D(k) in column k of A.
!              (NOTE: Diagonal element U(k,k) is a UNIT element
!              and not stored)
!                 A(k,k) := D(k,k) = W(k,kw)
!                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)
!
!              (NOTE: No need to use for Hermitian matrix
!              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal
!              element D(k,k) from W (potentially saves only one load))
            A(1:K,K) = W(1:K,KW)
            IF( K > 1 ) THEN
!
!                 (NOTE: No need to check if A(k,k) is NOT 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  case A(k,k) = 0 falls into 2x2 pivot case(3))
!
!                 Handle division by a small number
!
               A(1:K-1,K) = A(1:K-1,K) / REAL(A(K,K))
!               T = REAL( A( K, K ) )
!               IF( ABS( T ) >= SFMIN ) THEN
!                  R1 = 1.0E+0 / T
!                  CALL CSSCAL( K-1, R1, A( 1, K ), 1 )
!               ELSE
!                  DO II = 1, K-1
!                     A( II, K ) = A( II, K ) / T
!                  ENDDO
!               END IF
!
!                 (2) Conjugate column W(kw)
!
               CALL CLACGV( K-1, W( 1, KW ), 1 )
            END IF
!
         ELSE
!
!              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold
!
!              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              (1) Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
!              block D(k-1:k,k-1:k) in columns k-1 and k of A.
!              (NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
!              block and not stored)
!                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
!                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
!                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )
!
            IF( K > 2 ) THEN
!
!                 Factor out the columns of the inverse of 2-by-2 pivot
!                 block D, so that each column contains 1, to reduce the
!                 number of FLOPS when we multiply panel
!                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).
!
!                 D**(-1) = ( d11 cj(d21) )**(-1) =
!                           ( d21    d22 )
!
!                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
!                                          ( (-d21) (     d11 ) )
!
!                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
!
!                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
!                     (     (      -1 )           ( d11/conj(d21) ) )
!
!                 = 1/(|d21|**2) * 1/(D22*D11-1) *
!
!                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                     (     (  -1 )           ( D22 ) )
!
!                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                                      (     (  -1 )           ( D22 ) )
!
!                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
!                   (               (  -1 )         ( D22 ) )
!
!                 Handle division by a small number. (NOTE: order of
!                 operations is important)
!
!                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) )
!                   (   ((  -1 )          )   (( D22 )     ) ),
!
!                 where D11 = d22/d21,
!                       D22 = d11/conj(d21),
!                       D21 = d21,
!                       T = 1/(D22*D11-1).
!
!                 (NOTE: No need to check for division by 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  (a) d21 != 0 in 2x2 pivot case(4),
!                      since |d21| should be larger than |d11| and |d22|;
!                  (b) (D22*D11 - 1) != 0, since from (a),
!                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)
!
               D21 = W( K-1, KW )
               D11 = W( K, KW ) / CONJG( D21 )
               D22 = W( K-1, KW-1 ) / D21
               T = 1.0E+0 / ( REAL( D11*D22 )-1.0E+0 )
!
!                 Update elements in columns A(k-1) and A(k) as
!                 dot products of rows of ( W(kw-1) W(kw) ) and columns
!                 of D**(-1)
!
               A(1:K-2,K-1) = T*((D11*W(1:K-2,KW-1)-W(1:K-2,KW))/D21 )
               A(1:K-2,K) = T*((D22*W(1:K-2,KW)-W(1:K-2,KW-1))/CONJG( D21 ) )
            END IF
!
!              Copy D(k) to A
!
            A( K-1, K-1 ) = W( K-1, KW-1 )
            A( K-1, K ) = W( K-1, KW )
            A( K, K ) = W( K, KW )
!
!              (2) Conjugate columns W(kw) and W(kw-1)
!
            CALL CLACGV( K-1, W( 1, KW ), 1 )
            CALL CLACGV( K-2, W( 1, KW-1 ), 1 )
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
30    CONTINUE
!
!        Update the upper triangle of A11 (= A(1:k,1:k)) as
!
!        A11 := A11 - U12*D*U12**H = A11 - U12*W**H
!
!        computing blocks of NB columns at a time (note that conjg(W) is
!        actually stored)
!
      DO J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
         JB = MIN( NB, K-J+1 )
!
!           Update the upper triangle of the diagonal block
!
         DO JJ = J, J + JB - 1
            A( JJ, JJ ) = REAL( A( JJ, JJ ) )
            CALL CGEMV( 'No transpose', JJ-J+1, N-K, -(1.0E+0,0.0E+0), &
                        A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, (1.0E+0,0.0E+0), &
                        A( J, JJ ), 1 )
            A( JJ, JJ ) = REAL( A( JJ, JJ ) )
         ENDDO
!
!           Update the rectangular superdiagonal block
!
         IF( J >= 2 ) &
            CALL CGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, &
                        -(1.0E+0,0.0E+0), A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, &
                        (1.0E+0,0.0E+0), A( 1, J ), LDA )
      ENDDO
!
!        Put U12 in standard form by partially undoing the interchanges
!        in of rows in columns k+1:n looping backwards from k+1 to n
!
      J = K + 1
60    CONTINUE
!
!           Undo the interchanges (if any) of rows J and JP2
!           (or J and JP2, and J+1 and JP1) at each step J
!
         KSTEP = 1
         JP1 = 1
!           (Here, J is a diagonal index)
         JJ = J
         JP2 = IPIV( J )
         IF( JP2 < 0 ) THEN
            JP2 = -JP2
!              (Here, J is a diagonal index)
            J = J + 1
            JP1 = -IPIV( J )
            KSTEP = 2
         END IF
!           (NOTE: Here, J is used to determine row length. Length N-J+1
!           of the rows to swap back doesn't include diagonal element)
         J = J + 1
         IF( JP2 /= JJ .AND. J <= N ) THEN
            A_TMP(1:N-J+1) = A(JP2,J:N)
            A(JP2,J:N) = A(JJ,J:N)
            A(JJ,J:N) = A_TMP(1:N-J+1)
         ENDIF
         JJ = JJ + 1
         IF( KSTEP == 2 .AND. JP1 /= JJ .AND. J <= N ) THEN
            A_TMP(1:N-J+1) = A(JP1,J:N)
            A(JP1,J:N) = A(JJ,J:N)
            A(JJ,J:N) = A_TMP(1:N-J+1)
         ENDIF
      IF( J < N ) GO TO 60
!
!        Set KB to the number of columns factorized
!
      KB = N - K
!
   ELSE
!
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22 (note that conjg(W) is actually stored)
!
!        K is the main loop index, increasing from 1 in steps of 1 or 2
!
      K = 1
70    CONTINUE
!
!        Exit from loop
!
      IF( ( K >= NB .AND. NB < N ) .OR. K > N ) GO TO 90
!
      KSTEP = 1
      P = K
!
!        Copy column K of A to column K of W and update column K of W
!
      W( K, K ) = REAL( A( K, K ) )
      IF( K < N ) W(K+1:N,K) = A(K+1:N,K)
      IF( K > 1 ) THEN
         CALL CGEMV( 'No transpose', N-K+1, K-1, -(1.0E+0,0.0E+0), A( K, 1 ), &
                     LDA, W( K, 1 ), LDW, (1.0E+0,0.0E+0), W( K, K ), 1 )
         W( K, K ) = REAL( W( K, K ) )
      END IF
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
      ABSAKK = ABS( REAL( W( K, K ) ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
      IF( K < N ) THEN
         IMAX = K + ICAMAX( N-K, W( K+1, K ), 1 )
         COLMAX = CABS1( W( IMAX, K ) )
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
         A( K, K ) = REAL( W( K, K ) )
         IF( K < N ) A(K+1:N,K) = W(K+1:N,K)
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
72          CONTINUE
!
!                 BEGIN pivot search loop body
!
!
!                 Copy column IMAX to column k+1 of W and update it
!
               W(K:IMAX-1,K+1) = A(IMAX,K:IMAX-1)
               CALL CLACGV( IMAX-K, W( K, K+1 ), 1 )
               W( IMAX, K+1 ) = REAL( A( IMAX, IMAX ) )
!
               IF( IMAX < N ) W(IMAX+1:N,K+1) = A(IMAX+1:N,IMAX)
!
               IF( K > 1 ) THEN
                  CALL CGEMV( 'No transpose', N-K+1, K-1, -(1.0E+0,0.0E+0), &
                               A( K, 1 ), LDA, W( IMAX, 1 ), LDW, &
                               (1.0E+0,0.0E+0), W( K, K+1 ), 1 )
                  W( IMAX, K+1 ) = REAL( W( IMAX, K+1 ) )
               END IF
!
!                 JMAX is the column-index of the largest off-diagonal
!                 element in row IMAX, and ROWMAX is its absolute value.
!                 Determine both ROWMAX and JMAX.
!
               IF( IMAX /= K ) THEN
                  JMAX = K - 1 + ICAMAX( IMAX-K, W( K, K+1 ), 1 )
                  ROWMAX = CABS1( W( JMAX, K+1 ) )
               ELSE
                  ROWMAX = 0.0E+0
               END IF
!
               IF( IMAX < N ) THEN
                  ITEMP = IMAX + ICAMAX( N-IMAX, W( IMAX+1, K+1 ), 1)
                  STEMP = CABS1( W( ITEMP, K+1 ) )
                  IF( STEMP > ROWMAX ) THEN
                     ROWMAX = STEMP
                     JMAX = ITEMP
                  END IF
               END IF
!
!                 Case(2)
!                 Equivalent to testing for
!                 ABS( REAL( W( IMAX,K+1 ) ) ) >= ALPHA*ROWMAX
!                 (used to handle NaN and Inf)
!
               IF( .NOT.( ABS( REAL( W( IMAX,K+1 ) ) ) < ALPHA*ROWMAX ) ) THEN
!
!                    interchange rows and columns K and IMAX,
!                    use 1-by-1 pivot block
!
                  KP = IMAX
!
!                    copy column K+1 of W to column K of W
!
                  W(K:N,K) = W(K:N,K+1)
!
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
!
!                    Copy updated JMAXth (next IMAXth) column to Kth of W
!
                  W(K:N,K) = W(K:N,K+1)
!
               END IF
!
!
!                 End pivot search loop body
!
            IF( .NOT.DONE ) GOTO 72
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
!           Interchange rows and columns P and K (only for 2-by-2 pivot).
!           Updated column P is already stored in column K of W.
!
         IF( ( KSTEP == 2 ) .AND. ( P /= K ) ) THEN
!
!              Copy non-updated column KK-1 to column P of submatrix A
!              at step K. No need to copy element into columns
!              K and K+1 of A for 2-by-2 pivot, since these columns
!              will be later overwritten.
!
            A( P, P ) = REAL( A( K, K ) )
            A(P,K+1:P-1) = A(K+1:P-1,K)
            CALL CLACGV( P-K-1, A( P, K+1 ), LDA )
            IF( P < N ) A(P+1:N,P) = A(P+1:N,K)
!
!              Interchange rows K and P in first K-1 columns of A
!              (columns K and K+1 of A for 2-by-2 pivot will be
!              later overwritten). Interchange rows K and P
!              in first KK columns of W.
!
            IF( K > 1 ) THEN
               A_TMP(1:K-1) = A(K,1:K-1)
               A(K,1:K-1) = A(P,1:K-1)
               A(P,1:K-1) = A_TMP(1:K-1)
            ENDIF
            W_TMP(1:KK) = W(K,1:KK)
            W(K,1:KK) = W(P,1:KK)
            W(P,1:KK) = W_TMP(1:KK)
         END IF
!
!           Interchange rows and columns KP and KK.
!           Updated column KP is already stored in column KK of W.
!
         IF( KP /= KK ) THEN
!
!              Copy non-updated column KK to column KP of submatrix A
!              at step K. No need to copy element into column K
!              (or K and K+1 for 2-by-2 pivot) of A, since these columns
!              will be later overwritten.
!
            A( KP, KP ) = REAL( A( KK, KK ) )
            A(KP,KK+1:KP-1) = A(KK+1:KP-1,KK)
            CALL CLACGV( KP-KK-1, A( KP, KK+1 ), LDA )
            IF( KP < N ) A(KP+1:N,KP) = A(KP+1:N,KK)
!
!              Interchange rows KK and KP in first K-1 columns of A
!              (column K (or K and K+1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in first KK columns of W.
!
            IF( K > 1 ) THEN
               A_TMP(1:K-1) = A(KK,1:K-1)
               A(KK,1:K-1) = A(KP,1:K-1)
               A(KP,1:K-1) = A_TMP(1:K-1)
            ENDIF
            W_TMP(1:KK) = W(KK,1:KK)
            W(KK,1:KK) = W(KP,1:KK)
            W(KP,1:KK) = W_TMP(1:KK)
         END IF
!
         IF( KSTEP == 1 ) THEN
!
!              1-by-1 pivot block D(k): column k of W now holds
!
!              W(k) = L(k)*D(k),
!
!              where L(k) is the k-th column of L
!
!              (1) Store subdiag. elements of column L(k)
!              and 1-by-1 block D(k) in column k of A.
!              (NOTE: Diagonal element L(k,k) is a UNIT element
!              and not stored)
!                 A(k,k) := D(k,k) = W(k,k)
!                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)
!
!              (NOTE: No need to use for Hermitian matrix
!              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal
!              element D(k,k) from W (potentially saves only one load))
            A(K:N,K) = W(K:N,K)
            IF( K < N ) THEN
!
!                 (NOTE: No need to check if A(k,k) is NOT 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  case A(k,k) = 0 falls into 2x2 pivot case(3))
!
!                 Handle division by a small number
!
!                T = REAL( A( K, K ) )
!                IF( ABS( T ) >= SFMIN ) THEN
!                   R1 = 1.0E+0 / T
!                   CALL CSSCAL( N-K, R1, A( K+1, K ), 1 )
!                ELSE
!                   DO II = K + 1, N
!                      A( II, K ) = A( II, K ) / T
!                   ENDDO
!                END IF
               A(K+1:N,K) = A(K+1:N,K) / REAL( A( K, K ) )
!
!                 (2) Conjugate column W(k)
!
               CALL CLACGV( N-K, W( K+1, K ), 1 )
            END IF
!
         ELSE
!
!              2-by-2 pivot block D(k): columns k and k+1 of W now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
!              (1) Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
!              block D(k:k+1,k:k+1) in columns k and k+1 of A.
!              NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
!              block and not stored.
!                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
!                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
!                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )
!
            IF( K < N-1 ) THEN
!
!                 Factor out the columns of the inverse of 2-by-2 pivot
!                 block D, so that each column contains 1, to reduce the
!                 number of FLOPS when we multiply panel
!                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).
!
!                 D**(-1) = ( d11 cj(d21) )**(-1) =
!                           ( d21    d22 )
!
!                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
!                                          ( (-d21) (     d11 ) )
!
!                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
!
!                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
!                     (     (      -1 )           ( d11/conj(d21) ) )
!
!                 = 1/(|d21|**2) * 1/(D22*D11-1) *
!
!                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                     (     (  -1 )           ( D22 ) )
!
!                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
!                                      (     (  -1 )           ( D22 ) )
!
!                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
!                   (               (  -1 )         ( D22 ) )
!
!                 Handle division by a small number. (NOTE: order of
!                 operations is important)
!
!                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) )
!                   (   ((  -1 )          )   (( D22 )     ) ),
!
!                 where D11 = d22/d21,
!                       D22 = d11/conj(d21),
!                       D21 = d21,
!                       T = 1/(D22*D11-1).
!
!                 (NOTE: No need to check for division by 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  (a) d21 != 0 in 2x2 pivot case(4),
!                      since |d21| should be larger than |d11| and |d22|;
!                  (b) (D22*D11 - 1) != 0, since from (a),
!                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)
!
               D21 = W( K+1, K )
               D11 = W( K+1, K+1 ) / D21
               D22 = W( K, K ) / CONJG( D21 )
               T = 1.0E+0 / ( REAL( D11*D22 )-1.0E+0 )
!
!                 Update elements in columns A(k) and A(k+1) as
!                 dot products of rows of ( W(k) W(k+1) ) and columns
!                 of D**(-1)
!
               A(K+2:N,K) = T*((D11*W(K+2:N,K)-W(K+2:N,K+1))/CONJG(D21))
               A(K+2:N,K+1) = T*((D22*W(K+2:N,K+1)-W(K+2:N,K))/D21)
            END IF
!
!              Copy D(k) to A
!
            A( K, K ) = W( K, K )
            A( K+1, K ) = W( K+1, K )
            A( K+1, K+1 ) = W( K+1, K+1 )
!
!              (2) Conjugate columns W(k) and W(k+1)
!
            CALL CLACGV( N-K, W( K+1, K ), 1 )
            CALL CLACGV( N-K-1, W( K+2, K+1 ), 1 )
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
      GO TO 70
!
90    CONTINUE
!
!        Update the lower triangle of A22 (= A(k:n,k:n)) as
!
!        A22 := A22 - L21*D*L21**H = A22 - L21*W**H
!
!        computing blocks of NB columns at a time (note that conjg(W) is
!        actually stored)
!
      DO J = K, N, NB
         JB = MIN( NB, N-J+1 )
!
!           Update the lower triangle of the diagonal block
!
         DO JJ = J, J + JB - 1
            A( JJ, JJ ) = REAL( A( JJ, JJ ) )
            CALL CGEMV( 'No transpose', J+JB-JJ, K-1, -(1.0E+0,0.0E+0), &
                        A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, (1.0E+0,0.0E+0), &
                        A( JJ, JJ ), 1 )
            A( JJ, JJ ) = REAL( A( JJ, JJ ) )
         ENDDO
!
!           Update the rectangular subdiagonal block
!
         IF( J+JB <= N ) &
            CALL CGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB, &
                        K-1, -(1.0E+0,0.0E+0), A( J+JB, 1 ), LDA, W( J, 1 ), &
                        LDW, (1.0E+0,0.0E+0), A( J+JB, J ), LDA )
      ENDDO
!
!        Put L21 in standard form by partially undoing the interchanges
!        of rows in columns 1:k-1 looping backwards from k-1 to 1
!
      J = K - 1
  120    CONTINUE
!
!           Undo the interchanges (if any) of rows J and JP2
!           (or J and JP2, and J-1 and JP1) at each step J
!
         KSTEP = 1
         JP1 = 1
!           (Here, J is a diagonal index)
         JJ = J
         JP2 = IPIV( J )
         IF( JP2 < 0 ) THEN
            JP2 = -JP2
!              (Here, J is a diagonal index)
            J = J - 1
            JP1 = -IPIV( J )
            KSTEP = 2
         END IF
!           (NOTE: Here, J is used to determine row length. Length J
!           of the rows to swap back doesn't include diagonal element)
         J = J - 1
         IF( JP2 /= JJ .AND. J >= 1 ) THEN
            A_TMP(1:J) = A(JP2,1:J)
            A(JP2,1:J) = A(JJ,1:J)
            A(JJ,1:J) = A_TMP(1:J)
         ENDIF
         JJ = JJ -1
         IF( KSTEP == 2 .AND. JP1 /= JJ .AND. J >= 1 ) THEN
            A_TMP(1:J) = A(JP1,1:J)
            A(JP1,1:J) = A(JJ,1:J)
            A(JJ,1:J) = A_TMP(1:J)
         ENDIF
      IF( J > 1 ) GO TO 120
!
!        Set KB to the number of columns factorized
!
      KB = K - 1
!
   END IF
   RETURN
!
!     End of CLAHEF_ROOK
!
END

