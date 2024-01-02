!> \brief \b CLAHEF computes a partial factorization of a complex Hermitian indefinite matrix using the Bunch-Kaufman diagonal pivoting method (blocked algorithm, calling Level 3 BLAS).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAHEF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHEF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
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
!> CLAHEF computes a partial factorization of a complex Hermitian
!> matrix A using the Bunch-Kaufman diagonal pivoting method. The
!> partial factorization has the form:
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
!> CLAHEF is an auxiliary routine called by CHETRF. It uses blocked code
!> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
!> A22 (if UPLO = 'L').
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
!>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
!>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>             is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             Only the first KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
!>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
!>             is a 2-by-2 diagonal block.
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
!> \ingroup lahef
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!> \endverbatim
!
!  =====================================================================
   SUBROUTINE CLAHEF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
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
   INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP, &
                      KSTEP, KW
   REAL               ABSAKK, COLMAX, R1, ROWMAX, T
   COMPLEX            D11, D21, D22, Z, A_TMP( N ), W_TMP( NB )
!
!     Initialize ALPHA for use in choosing pivot block size.
!
   REAL, PARAMETER :: ALPHA = ( 1.0E+0+SQRT( 17.0E+0 ) ) / 8.0E+0
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ICAMAX
   REAL               CABS1
   EXTERNAL           LSAME, ICAMAX, CABS1
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CGEMV
!     ..
!     .. Executable Statements ..
!
   INFO = 0
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
!
!        Copy column K of A to column KW of W and update it
!
      W(1:K-1,KW) = A(1:K-1,K)
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
         A( K, K ) = REAL( A( K, K ) )
      ELSE
!
!           ============================================================
!
!           BEGIN pivot search
!
!           Case(1)
         IF( ABSAKK >= ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
            KP = K
         ELSE
!
!              BEGIN pivot search along IMAX row
!
!
!              Copy column IMAX to column KW-1 of W and update it
!
            W(1:IMAX-1,KW-1) = A(1:IMAX-1,IMAX)
            W( IMAX, KW-1 ) = REAL( A( IMAX, IMAX ) )
            W(IMAX+1:K,KW-1) = CONJG(A(IMAX,IMAX+1:K))
            IF( K < N ) THEN
               CALL CGEMV( 'No transpose', K, N-K, -(1.0E+0,0.0E+0), &
                           A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW, &
                           (1.0E+0,0.0E+0), W( 1, KW-1 ), 1 )
               W( IMAX, KW-1 ) = REAL( W( IMAX, KW-1 ) )
            END IF
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value.
!              Determine only ROWMAX.
!
            JMAX = IMAX + ICAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
            ROWMAX = CABS1( W( JMAX, KW-1 ) )
            IF( IMAX > 1 ) THEN
               JMAX = ICAMAX( IMAX-1, W( 1, KW-1 ), 1 )
               ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, KW-1 ) ) )
            END IF
!
!              Case(2)
            IF( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
               KP = K
!
!              Case(3)
            ELSE IF( ABS( REAL( W( IMAX, KW-1 ) ) ) >= ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
               KP = IMAX
!
!                 copy column KW-1 of W to column KW of W
!
               
               W(1:K,KW) = W(1:K,KW-1)
!
!              Case(4)
            ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
               KP = IMAX
               KSTEP = 2
            END IF
!
!
!              END pivot search along IMAX row
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
            A(KP,KP+1:KK-1) = CONJG(A(KP+1:KK-1,KK))
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
!              A( K, K ) = DBLE( W( K, K) ) to separately copy diagonal
!              element D(k,k) from W (potentially saves only one load))
            A(1:K,K) = W(1:K,KW)
            IF( K > 1 ) THEN
!
!                 (NOTE: No need to check if A(k,k) is NOT 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  case A(k,k) = 0 falls into 2x2 pivot case(4))
!
               A(1:K-1,K) = A(1:K-1,K) / REAL( A( K, K ) )
!
!                 (2) Conjugate column W(kw)
!
               W(1:K-1,KW) = CONJG(W(1:K-1,KW))
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
!                 = ( conj(D21)*( D11 ) D21*(  -1 ) )
!                   (           (  -1 )     ( D22 ) ),
!
!                 where D11 = d22/d21,
!                       D22 = d11/conj(d21),
!                       D21 = T/d21,
!                       T = 1/(D22*D11-1).
!
!                 (NOTE: No need to check for division by 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  (a) d21 != 0, since in 2x2 pivot case(4)
!                      |d21| should be larger than |d11| and |d22|;
!                  (b) (D22*D11 - 1) != 0, since from (a),
!                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)
!
               D21 = W( K-1, KW )
               D11 = W( K, KW ) / CONJG( D21 )
               D22 = W( K-1, KW-1 ) / D21
               T = 1.0E+0 / ( REAL( D11*D22 )-1.0E+0 )
               D21 = T / D21
!
!                 Update elements in columns A(k-1) and A(k) as
!                 dot products of rows of ( W(kw-1) W(kw) ) and columns
!                 of D**(-1)
!
               A(1:K-2,K-1) = D21*(D11*W(1:K-2,KW-1)-W(1:K-2,KW))
               A(1:K-2,K) = CONJG(D21)*(D22*W(1:K-2,KW)-W(1:K-2,KW-1))
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
            W(1:K-1,KW) = CONJG(W(1:K-1,KW))
            W(1:K-2,KW-1) = CONJG(W(1:K-2,KW-1))
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
         IPIV( K ) = -KP
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
!           Undo the interchanges (if any) of rows J and JP
!           at each step J
!
!           (Here, J is a diagonal index)
         JJ = J
         JP = IPIV( J )
         IF( JP < 0 ) THEN
            JP = -JP
!              (Here, J is a diagonal index)
            J = J + 1
         END IF
!           (NOTE: Here, J is used to determine row length. Length N-J+1
!           of the rows to swap back doesn't include diagonal element)
         J = J + 1
         IF( JP /= JJ .AND. J <= N ) THEN
            A_TMP(1:N-J+1) = A(JP,J:N)
            A(JP,J:N) = A(JJ,J:N)
            A(JJ,J:N) = A_TMP(1:N-J+1)
         ENDIF
      IF( J <= N ) GO TO 60
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
!
!        Copy column K of A to column K of W and update it
!
      W( K, K ) = REAL( A( K, K ) )
      IF( K < N ) W(K+1:N,K) = A(K+1:N,K)
      CALL CGEMV( 'No transpose', N-K+1, K-1, -(1.0E+0,0.0E+0), A( K, 1 ), LDA, &
                  W( K, 1 ), LDW, (1.0E+0,0.0E+0), W( K, K ), 1 )
      W( K, K ) = REAL( W( K, K ) )
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
         A( K, K ) = REAL( A( K, K ) )
      ELSE
!
!           ============================================================
!
!           BEGIN pivot search
!
!           Case(1)
         IF( ABSAKK >= ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
            KP = K
         ELSE
!
!              BEGIN pivot search along IMAX row
!
!
!              Copy column IMAX to column K+1 of W and update it
!
            W(K:IMAX-1,K+1) = CONJG(A(IMAX,K:IMAX-1))
            W( IMAX, K+1 ) = REAL( A( IMAX, IMAX ) )
            IF( IMAX < N ) W(IMAX+1:N,K+1) = A(IMAX+1:N,IMAX)
            CALL CGEMV( 'No transpose', N-K+1, K-1, -(1.0E+0,0.0E+0), A( K, 1 ), &
                        LDA, W( IMAX, 1 ), LDW, (1.0E+0,0.0E+0), W( K, K+1 ), 1 )
            W( IMAX, K+1 ) = REAL( W( IMAX, K+1 ) )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value.
!              Determine only ROWMAX.
!
            JMAX = K - 1 + ICAMAX( IMAX-K, W( K, K+1 ), 1 )
            ROWMAX = CABS1( W( JMAX, K+1 ) )
            IF( IMAX < N ) THEN
               JMAX = IMAX + ICAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
               ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, K+1 ) ) )
            END IF
!
!              Case(2)
            IF( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
               KP = K
!
!              Case(3)
            ELSE IF( ABS( REAL( W( IMAX, K+1 ) ) ) >= ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
               KP = IMAX
!
!                 copy column K+1 of W to column K of W
!
               W(K:N,K) = W(K:N,K+1)
!
!              Case(4)
            ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
               KP = IMAX
               KSTEP = 2
            END IF
!
!
!              END pivot search along IMAX row
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
            A(KP,KK+1:KK+KP-KK-1) = CONJG(A(KK+1:KK+KP-KK-1,KK))
            IF( KP < N ) A(KP+1:N,KP) = A(KP+1:N,KK)
!
!              Interchange rows KK and KP in first K-1 columns of A
!              (columns K (or K and K+1 for 2-by-2 pivot) of A will be
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
!              A( K, K ) = DBLE( W( K, K) ) to separately copy diagonal
!              element D(k,k) from W (potentially saves only one load))
            A(K:N,K) = W(K:N,K)
            IF( K < N ) THEN
!
!                 (NOTE: No need to check if A(k,k) is NOT 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  case A(k,k) = 0 falls into 2x2 pivot case(4))
!
               A(K+1:N,K) = A(K+1:N,K) / REAL( A( K, K ) )
!
!                 (2) Conjugate column W(k)
!
               W(K+1:N,K) = CONJG(W(K+1:N,K))
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
!              (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
!              block and not stored)
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
!                 = ( conj(D21)*( D11 ) D21*(  -1 ) )
!                   (           (  -1 )     ( D22 ) )
!
!                 where D11 = d22/d21,
!                       D22 = d11/conj(d21),
!                       D21 = T/d21,
!                       T = 1/(D22*D11-1).
!
!                 (NOTE: No need to check for division by 0.0E+0,
!                  since that was ensured earlier in pivot search:
!                  (a) d21 != 0, since in 2x2 pivot case(4)
!                      |d21| should be larger than |d11| and |d22|;
!                  (b) (D22*D11 - 1) != 0, since from (a),
!                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)
!
               D21 = W( K+1, K )
               D11 = W( K+1, K+1 ) / D21
               D22 = W( K, K ) / CONJG( D21 )
               T = 1.0E+0 / ( REAL( D11*D22 )-1.0E+0 )
               D21 = T / D21
!
!                 Update elements in columns A(k) and A(k+1) as
!                 dot products of rows of ( W(k) W(k+1) ) and columns
!                 of D**(-1)
!
               A(K+2:N,K) = CONJG(D21)*(D11*W(K+2:N,K)-W(K+2:N,K+1))
               A(K+2:N,K+1) = D21*(D22*W(K+2:N,K+1)-W(K+2:N,K))
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
            W(K+1:N,K) = CONJG(W(K+1:N,K))
            W(K+2:N,K+1) = CONJG(W(K+2:N,K+1))
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
         IPIV( K ) = -KP
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
!           Undo the interchanges (if any) of rows J and JP
!           at each step J
!
!           (Here, J is a diagonal index)
         JJ = J
         JP = IPIV( J )
         IF( JP < 0 ) THEN
            JP = -JP
!              (Here, J is a diagonal index)
            J = J - 1
         END IF
!           (NOTE: Here, J is used to determine row length. Length J
!           of the rows to swap back doesn't include diagonal element)
         J = J - 1
         IF( JP /= JJ .AND. J >= 1 ) THEN
            A_TMP(1:J) = A(JP,1:J)
            A(JP,1:J) = A(JJ,1:J)
            A(JJ,1:J) = A_TMP(1:J)
         ENDIF
      IF( J >= 1 ) GO TO 120
!
!        Set KB to the number of columns factorized
!
      KB = K - 1
!
   END IF
   RETURN
!
!     End of CLAHEF
!
END
