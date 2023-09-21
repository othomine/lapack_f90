!> \brief \b ZUNMR2 multiplies a general matrix by the unitary matrix from a RQ factorization determined by cgerqf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNMR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNMR2 overwrites the general complex m-by-n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**H* C  if SIDE = 'L' and TRANS = 'C', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**H if SIDE = 'R' and TRANS = 'C',
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1)**H H(2)**H . . . H(k)**H
!>
!> as returned by ZGERQF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left
!>          = 'R': apply Q or Q**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'C': apply Q**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                               (LDA,M) if SIDE = 'L',
!>                               (LDA,N) if SIDE = 'R'
!>          The i-th row must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZGERQF in the last k rows of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,K).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGERQF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the m-by-n matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
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
!
!> \ingroup unmr2
!
!  =====================================================================
   SUBROUTINE ZUNMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                      WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          SIDE, TRANS
   INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
   COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX*16         ONE
   PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            LEFT, NOTRAN
   INTEGER            I, I1, I2, I3, MI, NI, NQ
   COMPLEX*16         AII, TAUI
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZLACGV, ZLARF
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   LEFT = LSAME( SIDE, 'L' )
   NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
   IF( LEFT ) THEN
      NQ = M
   ELSE
      NQ = N
   END IF
   IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
      INFO = -2
   ELSE IF( M < 0 ) THEN
      INFO = -3
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( K < 0 .OR. K > NQ ) THEN
      INFO = -5
   ELSE IF( LDA < MAX( 1, K ) ) THEN
      INFO = -7
   ELSE IF( LDC < MAX( 1, M ) ) THEN
      INFO = -10
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'ZUNMR2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M == 0 .OR. N == 0 .OR. K == 0 ) &
      RETURN
!
   IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
      I1 = 1
      I2 = K
      I3 = 1
   ELSE
      I1 = K
      I2 = 1
      I3 = -1
   END IF
!
   IF( LEFT ) THEN
      NI = N
   ELSE
      MI = M
   END IF
!
   DO I = I1, I2, I3
      IF( LEFT ) THEN
!
!           H(i) or H(i)**H is applied to C(1:m-k+i,1:n)
!
         MI = M - K + I
      ELSE
!
!           H(i) or H(i)**H is applied to C(1:m,1:n-k+i)
!
         NI = N - K + I
      END IF
!
!        Apply H(i) or H(i)**H
!
      IF( NOTRAN ) THEN
         TAUI = DCONJG( TAU( I ) )
      ELSE
         TAUI = TAU( I )
      END IF
      CALL ZLACGV( NQ-K+I-1, A( I, 1 ), LDA )
      AII = A( I, NQ-K+I )
      A( I, NQ-K+I ) = ONE
      CALL ZLARF( SIDE, MI, NI, A( I, 1 ), LDA, TAUI, C, LDC, WORK )
      A( I, NQ-K+I ) = AII
      CALL ZLACGV( NQ-K+I-1, A( I, 1 ), LDA )
   ENDDO
   RETURN
!
!     End of ZUNMR2
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
