!> \brief \b CTPLQT2
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER   INFO, LDA, LDB, LDT, N, M, L
!       ..
!       .. Array Arguments ..
!       COMPLEX   A( LDA, * ), B( LDB, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTPLQT2 computes a LQ a factorization of a complex "triangular-pentagonal"
!> matrix C, which is composed of a triangular block A and pentagonal block B,
!> using the compact WY representation for Q.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The total number of rows of the matrix B.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix B, and the order of
!>          the triangular matrix A.
!>          N >= 0.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of rows of the lower trapezoidal part of B.
!>          MIN(M,N) >= L >= 0.  See Further Details.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,M)
!>          On entry, the lower triangular M-by-M matrix A.
!>          On exit, the elements on and below the diagonal of the array
!>          contain the lower triangular matrix L.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          On entry, the pentagonal M-by-N matrix B.  The first N-L columns
!>          are rectangular, and the last L columns are lower trapezoidal.
!>          On exit, B contains the pentagonal matrix V.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,M)
!>          The N-by-N upper triangular factor T of the block reflector.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max(1,M)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup tplqt2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The input matrix C is a M-by-(M+N) matrix
!>
!>               C = [ A ][ B ]
!>
!>
!>  where A is an lower triangular M-by-M matrix, and B is M-by-N pentagonal
!>  matrix consisting of a M-by-(N-L) rectangular matrix B1 left of a M-by-L
!>  upper trapezoidal matrix B2:
!>
!>               B = [ B1 ][ B2 ]
!>                   [ B1 ]  <-     M-by-(N-L) rectangular
!>                   [ B2 ]  <-     M-by-L lower trapezoidal.
!>
!>  The lower trapezoidal matrix B2 consists of the first L columns of a
!>  N-by-N lower triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
!>  B is rectangular M-by-N; if M=L=N, B is lower triangular.
!>
!>  The matrix W stores the elementary reflectors H(i) in the i-th row
!>  above the diagonal (of A) in the M-by-(M+N) input matrix C
!>
!>               C = [ A ][ B ]
!>                   [ A ]  <- lower triangular M-by-M
!>                   [ B ]  <- M-by-N pentagonal
!>
!>  so that W can be represented as
!>
!>               W = [ I ][ V ]
!>                   [ I ]  <- identity, M-by-M
!>                   [ V ]  <- M-by-N, same form as B.
!>
!>  Thus, all of information needed for W is contained on exit in B, which
!>  we call V above.  Note that V has the same form as B; that is,
!>
!>               W = [ V1 ][ V2 ]
!>                   [ V1 ] <-     M-by-(N-L) rectangular
!>                   [ V2 ] <-     M-by-L lower trapezoidal.
!>
!>  The rows of V represent the vectors which define the H(i)'s.
!>  The (M+N)-by-(M+N) block reflector H is then given by
!>
!>               H = I - W**T * T * W
!>
!>  where W^H is the conjugate transpose of W and T is the upper triangular
!>  factor of the block reflector.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER        INFO, LDA, LDB, LDT, N, M, L
!     ..
!     .. Array Arguments ..
   COMPLEX     A( LDA, * ), B( LDB, * ), T( LDT, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER   I, J, P, MP, NP
   COMPLEX   ALPHA
!     ..
!     .. External Subroutines ..
   EXTERNAL  CLARFG, CGEMV, CGERC, CTRMV, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( L < 0 .OR. L > MIN(M,N) ) THEN
      INFO = -3
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -5
   ELSE IF( LDB < MAX( 1, M ) ) THEN
      INFO = -7
   ELSE IF( LDT < MAX( 1, M ) ) THEN
      INFO = -9
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTPLQT2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 .OR. M == 0 ) RETURN
!
   DO I = 1, M
!
!        Generate elementary reflector H(I) to annihilate B(I,:)
!
      P = N-L+MIN( L, I )
      CALL CLARFG( P+1, A( I, I ), B( I, 1 ), LDB, T( 1, I ) )
      T(1,I)=CONJG(T(1,I))
      IF( I < M ) THEN
         B( I,1:P) = CONJG(B(I,1:P))
!
!           W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]
!
         T( M,1:M-I) = (A( I+1:M, I ))
         CALL CGEMV( 'N', M-I, P, (1.0E+0,0.0E+0), B( I+1, 1 ), LDB, &
                     B( I, 1 ), LDB, (1.0E+0,0.0E+0), T( M, 1 ), LDT )
!
!           C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H
!
         ALPHA = -(T( 1, I ))
         A( I+1:M, I ) = A( I+1:M, I ) + ALPHA*(T( M,1:M-I))
         CALL CGERC( M-I, P, (ALPHA),  T( M, 1 ), LDT, &
             B( I, 1 ), LDB, B( I+1, 1 ), LDB )
         B( I,1:P) = CONJG(B(I,1:P))
      END IF
   END DO
!
   DO I = 2, M
!
!        T(I,1:I-1) := C(I:I-1,1:N)**H * (alpha * C(I,I:N))
!
      ALPHA = -(T( 1, I ))
      T( I,1:I-1) = (0.0E+0,0.0E+0)
      P = MIN( I-1, L )
      NP = MIN( N-L+1, N )
      MP = MIN( P+1, M )
      B(I,1:N-L+P)=CONJG(B(I,1:N-L+P))
!
!        Triangular part of B2
!
      T( I,1:P) = (ALPHA*B( I, N-L+1:N-L+P ))
      CALL CTRMV( 'L', 'N', 'N', P, B( 1, NP ), LDB, T( I, 1 ), LDT )
!
!        Rectangular part of B2
!
      CALL CGEMV( 'N', I-1-P, L,  ALPHA, B( MP, NP ), LDB, &
                  B( I, NP ), LDB, (0.0E+0,0.0E+0), T( I,MP ), LDT )
!
!        B1

!
      CALL CGEMV( 'N', I-1, N-L, ALPHA, B, LDB, B( I, 1 ), LDB, &
                  (1.0E+0,0.0E+0), T( I, 1 ), LDT )
!

!
!        T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)
!
      T(I,1:I-1)=CONJG(T(I,1:I-1))
      CALL CTRMV( 'L', 'C', 'N', I-1, T, LDT, T( I, 1 ), LDT )
      T(I,1:I-1)=CONJG(T(I,1:I-1))
      B(I,1:N-L+P)=CONJG(B(I,1:N-L+P))
!
!        T(I,I) = tau(I)
!
      T( I, I ) = T( 1, I )
      T( 1, I ) = (0.0E+0,0.0E+0)
   END DO
   DO I=1,M
      T(I,I+1:M)=(T(I+1:M,I))
      T(I+1:M,I)=(0.0E+0,0.0E+0)
   END DO

!
!     End of CTPLQT2
!
END
