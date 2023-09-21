!> \brief \b CGEMLQT
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT,
!                          C, LDC, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER SIDE, TRANS
!       INTEGER   INFO, K, LDV, LDC, M, N, MB, LDT
!       ..
!       .. Array Arguments ..
!       COMPLEX V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEMLQT overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q C            C Q
!> TRANS = 'C':   Q**H C            C Q**H
!>
!> where Q is a complex unitary matrix defined as the product of K
!> elementary reflectors:
!>
!>       Q = H(1) H(2) . . . H(K) = I - V T V**H
!>
!> generated using the compact WY representation as returned by CGELQT.
!>
!> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'C':  Conjugate transpose, apply Q**H.
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
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The block size used for the storage of T.  K >= MB >= 1.
!>          This must be the same value of MB used to generate T
!>          in CGELQT.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension
!>                               (LDV,M) if SIDE = 'L',
!>                               (LDV,N) if SIDE = 'R'
!>          The i-th row must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          CGELQT in the first K rows of its array argument A.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,K).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,K)
!>          The upper triangular factors of the block reflectors
!>          as returned by CGELQT, stored as a MB-by-K matrix.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= MB.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q C, Q**H C, C Q**H or C Q.
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
!>          WORK is COMPLEX array. The dimension of
!>          WORK is N*MB if SIDE = 'L', or  M*MB if SIDE = 'R'.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup gemlqt
!
!  =====================================================================
   SUBROUTINE CGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT, &
                      C, LDC, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER SIDE, TRANS
   INTEGER   INFO, K, LDV, LDC, M, N, MB, LDT
!     ..
!     .. Array Arguments ..
   COMPLEX   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
   LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
   INTEGER            I, IB, LDWORK, KF, Q
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, CLARFB
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     .. Test the input arguments ..
!
   INFO   = 0
   LEFT   = LSAME( SIDE,  'L' )
   RIGHT  = LSAME( SIDE,  'R' )
   TRAN   = LSAME( TRANS, 'C' )
   NOTRAN = LSAME( TRANS, 'N' )
!
   IF( LEFT ) THEN
      LDWORK = MAX( 1, N )
      Q = M
   ELSE IF ( RIGHT ) THEN
      LDWORK = MAX( 1, M )
      Q = N
   END IF
   IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
      INFO = -1
   ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
      INFO = -2
   ELSE IF( M < 0 ) THEN
      INFO = -3
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( K < 0 .OR. K > Q ) THEN
      INFO = -5
   ELSE IF( MB < 1 .OR. (MB > K .AND. K > 0)) THEN
      INFO = -6
   ELSE IF( LDV < MAX( 1, K ) ) THEN
       INFO = -8
   ELSE IF( LDT < MB ) THEN
      INFO = -10
   ELSE IF( LDC < MAX( 1, M ) ) THEN
      INFO = -12
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CGEMLQT', -INFO )
      RETURN
   END IF
!
!     .. Quick return if possible ..
!
   IF( M == 0 .OR. N == 0 .OR. K == 0 ) RETURN
!
   IF( LEFT .AND. NOTRAN ) THEN
!
      DO I = 1, K, MB
         IB = MIN( MB, K-I+1 )
         CALL CLARFB( 'L', 'C', 'F', 'R', M-I+1, N, IB, &
                      V( I, I ), LDV, T( 1, I ), LDT, &
                      C( I, 1 ), LDC, WORK, LDWORK )
      END DO
!
   ELSE IF( RIGHT .AND. TRAN ) THEN
!
      DO I = 1, K, MB
         IB = MIN( MB, K-I+1 )
         CALL CLARFB( 'R', 'N', 'F', 'R', M, N-I+1, IB, &
                      V( I, I ), LDV, T( 1, I ), LDT, &
                      C( 1, I ), LDC, WORK, LDWORK )
      END DO
!
   ELSE IF( LEFT .AND. TRAN ) THEN
!
      KF = ((K-1)/MB)*MB+1
      DO I = KF, 1, -MB
         IB = MIN( MB, K-I+1 )
         CALL CLARFB( 'L', 'N', 'F', 'R', M-I+1, N, IB, &
                      V( I, I ), LDV, T( 1, I ), LDT, &
                      C( I, 1 ), LDC, WORK, LDWORK )
      END DO
!
   ELSE IF( RIGHT .AND. NOTRAN ) THEN
!
      KF = ((K-1)/MB)*MB+1
      DO I = KF, 1, -MB
         IB = MIN( MB, K-I+1 )
         CALL CLARFB( 'R', 'C', 'F', 'R', M, N-I+1, IB, &
                      V( I, I ), LDV, T( 1, I ), LDT, &
                      C( 1, I ), LDC, WORK, LDWORK )
      END DO
!
   END IF
!
   RETURN
!
!     End of CGEMLQT
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
