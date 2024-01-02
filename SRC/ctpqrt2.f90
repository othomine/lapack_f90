!> \brief \b CTPQRT2 computes a QR factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTPQRT2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpqrt2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpqrt2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpqrt2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
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
!> CTPQRT2 computes a QR factorization of a complex "triangular-pentagonal"
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
!>          The number of rows of the upper trapezoidal part of B.
!>          MIN(M,N) >= L >= 0.  See Further Details.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the upper triangular N-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the upper triangular matrix R.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
!>          are rectangular, and the last L rows are upper trapezoidal.
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
!>          T is COMPLEX array, dimension (LDT,N)
!>          The N-by-N upper triangular factor T of the block reflector.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max(1,N)
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
!> \ingroup tpqrt2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The input matrix C is a (N+M)-by-N matrix
!>
!>               C = [ A ]
!>                   [ B ]
!>
!>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
!>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
!>  upper trapezoidal matrix B2:
!>
!>               B = [ B1 ]  <- (M-L)-by-N rectangular
!>                   [ B2 ]  <-     L-by-N upper trapezoidal.
!>
!>  The upper trapezoidal matrix B2 consists of the first L rows of a
!>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
!>  B is rectangular M-by-N; if M=L=N, B is upper triangular.
!>
!>  The matrix W stores the elementary reflectors H(i) in the i-th column
!>  below the diagonal (of A) in the (N+M)-by-N input matrix C
!>
!>               C = [ A ]  <- upper triangular N-by-N
!>                   [ B ]  <- M-by-N pentagonal
!>
!>  so that W can be represented as
!>
!>               W = [ I ]  <- identity, N-by-N
!>                   [ V ]  <- M-by-N, same form as B.
!>
!>  Thus, all of information needed for W is contained on exit in B, which
!>  we call V above.  Note that V has the same form as B; that is,
!>
!>               V = [ V1 ] <- (M-L)-by-N rectangular
!>                   [ V2 ] <-     L-by-N upper trapezoidal.
!>
!>  The columns of V represent the vectors which define the H(i)'s.
!>  The (M+N)-by-(M+N) block reflector H is then given by
!>
!>               H = I - W * T * W**H
!>
!>  where W**H is the conjugate transpose of W and T is the upper triangular
!>  factor of the block reflector.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CTPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER   INFO, LDA, LDB, LDT, N, M, L
!     ..
!     .. Array Arguments ..
   COMPLEX   A( LDA, * ), B( LDB, * ), T( LDT, * )
!     ..
!
!  =====================================================================

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
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -5
   ELSE IF( LDB < MAX( 1, M ) ) THEN
      INFO = -7
   ELSE IF( LDT < MAX( 1, N ) ) THEN
      INFO = -9
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTPQRT2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 .OR. M == 0 ) RETURN
!
   DO I = 1, N
!
!        Generate elementary reflector H(I) to annihilate B(:,I)
!
      P = M-L+MIN( L, I )
      CALL CLARFG( P+1, A( I, I ), B( 1, I ), 1, T( I, 1 ) )
      IF( I < N ) THEN
!
!           W(1:N-I) := C(I:M,I+1:N)**H * C(I:M,I) [use W = T(:,N)]
!
         T(1:N-I, N ) = CONJG(A( I, I+1:N))
         CALL CGEMV( 'C', P, N-I, (1.0E+0,0.0E+0), B( 1, I+1 ), LDB, &
                     B( 1, I ), 1, (1.0E+0,0.0E+0), T( 1, N ), 1 )
!
!           C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)**H
!
         ALPHA = -CONJG(T( I, 1 ))
         A( I, I+1:N ) = A( I, I+1:N ) + ALPHA*CONJG(T(1:N-I, N ))
         CALL CGERC( P, N-I, ALPHA, B( 1, I ), 1, T( 1, N ), 1, B( 1, I+1 ), LDB )
      END IF
   END DO
!
   DO I = 2, N
!
!        T(1:I-1,I) := C(I:M,1:I-1)**H * (alpha * C(I:M,I))
!
      ALPHA = -T( I, 1 )
      T(1:I-1, I ) = (0.0E+0,0.0E+0)
      P = MIN( I-1, L )
      MP = MIN( M-L+1, M )
      NP = MIN( P+1, N )
!
!        Triangular part of B2
!
      T(1:P, I ) = ALPHA*B( M-L+1:M-L+P, I )
      CALL CTRMV( 'U', 'C', 'N', P, B( MP, 1 ), LDB, &
                  T( 1, I ), 1 )
!
!        Rectangular part of B2
!
      CALL CGEMV( 'C', L, I-1-P, ALPHA, B( MP, NP ), LDB, &
                  B( MP, I ), 1, (0.0E+0,0.0E+0), T( NP, I ), 1 )
!
!        B1
!
      CALL CGEMV( 'C', M-L, I-1, ALPHA, B, LDB, B( 1, I ), 1, &
                  (1.0E+0,0.0E+0), T( 1, I ), 1 )
!
!        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
!
      CALL CTRMV( 'U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 )
!
!        T(I,I) = tau(I)
!
      T( I, I ) = T( I, 1 )
      T( I, 1 ) = (0.0E+0,0.0E+0)
   END DO

!
!     End of CTPQRT2
!
END
