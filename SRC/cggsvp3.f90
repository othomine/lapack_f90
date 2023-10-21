!> \brief \b CGGSVP3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGGSVP3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggsvp3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggsvp3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggsvp3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,
!                           TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,
!                           IWORK, RWORK, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBQ, JOBU, JOBV
!       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK
!       REAL               TOLA, TOLB
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGGSVP3 computes unitary matrices U, V and Q such that
!>
!>                    N-K-L  K    L
!>  U**H*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;
!>                 L ( 0     0   A23 )
!>             M-K-L ( 0     0    0  )
!>
!>                  N-K-L  K    L
!>         =     K ( 0    A12  A13 )  if M-K-L < 0;
!>             M-K ( 0     0   A23 )
!>
!>                  N-K-L  K    L
!>  V**H*B*Q =   L ( 0     0   B13 )
!>             P-L ( 0     0    0  )
!>
!> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
!> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
!> otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective
!> numerical rank of the (M+P)-by-N matrix (A**H,B**H)**H.
!>
!> This decomposition is the preprocessing step for computing the
!> Generalized Singular Value Decomposition (GSVD), see subroutine
!> CGGSVD3.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          = 'U':  Unitary matrix U is computed;
!>          = 'N':  U is not computed.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          = 'V':  Unitary matrix V is computed;
!>          = 'N':  V is not computed.
!> \endverbatim
!>
!> \param[in] JOBQ
!> \verbatim
!>          JOBQ is CHARACTER*1
!>          = 'Q':  Unitary matrix Q is computed;
!>          = 'N':  Q is not computed.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows of the matrix B.  P >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, A contains the triangular (or trapezoidal) matrix
!>          described in the Purpose section.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          On entry, the P-by-N matrix B.
!>          On exit, B contains the triangular matrix described in
!>          the Purpose section.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,P).
!> \endverbatim
!>
!> \param[in] TOLA
!> \verbatim
!>          TOLA is REAL
!> \endverbatim
!>
!> \param[in] TOLB
!> \verbatim
!>          TOLB is REAL
!>
!>          TOLA and TOLB are the thresholds to determine the effective
!>          numerical rank of matrix B and a subblock of A. Generally,
!>          they are set to
!>             TOLA = MAX(M,N)*norm(A)*MACHEPS,
!>             TOLB = MAX(P,N)*norm(B)*MACHEPS.
!>          The size of TOLA and TOLB may affect the size of backward
!>          errors of the decomposition.
!> \endverbatim
!>
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is INTEGER
!>
!>          On exit, K and L specify the dimension of the subblocks
!>          described in Purpose section.
!>          K + L = effective numerical rank of (A**H,B**H)**H.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU,M)
!>          If JOBU = 'U', U contains the unitary matrix U.
!>          If JOBU = 'N', U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U. LDU >= max(1,M) if
!>          JOBU = 'U'; LDU >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,P)
!>          If JOBV = 'V', V contains the unitary matrix V.
!>          If JOBV = 'N', V is not referenced.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,P) if
!>          JOBV = 'V'; LDV >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ,N)
!>          If JOBQ = 'Q', Q contains the unitary matrix Q.
!>          If JOBQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= max(1,N) if
!>          JOBQ = 'Q'; LDQ >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup ggsvp3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The subroutine uses LAPACK subroutine CGEQP3 for the QR factorization
!>  with column pivoting to detect the effective numerical rank of the
!>  a matrix. It may be replaced by a better rank determination strategy.
!>
!>  CGGSVP3 replaces the deprecated subroutine CGGSVP.
!>
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, &
                       TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, &
                       IWORK, RWORK, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
   IMPLICIT NONE
!
!     .. Scalar Arguments ..
   CHARACTER          JOBQ, JOBU, JOBV
   INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, &
                      LWORK
   REAL               TOLA, TOLB
!     ..
!     .. Array Arguments ..
   INTEGER            IWORK( * )
   REAL               RWORK( * )
   COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), &
                      TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            FORWRD, WANTQ, WANTU, WANTV, LQUERY
   INTEGER            I, J, LWKOPT
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEQP3, CGEQR2, CGERQ2, CLACPY, CLAPMT, &
                      CLASET, CUNG2R, CUNM2R, CUNMR2, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
   WANTU = LSAME( JOBU, 'U' )
   WANTV = LSAME( JOBV, 'V' )
   WANTQ = LSAME( JOBQ, 'Q' )
   FORWRD = .TRUE.
   LQUERY = ( LWORK == -1 )
   LWKOPT = 1
!
!     Test the input arguments
!
   INFO = 0
   IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
      INFO = -1
   ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
      INFO = -2
   ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
      INFO = -3
   ELSE IF( M < 0 ) THEN
      INFO = -4
   ELSE IF( P < 0 ) THEN
      INFO = -5
   ELSE IF( N < 0 ) THEN
      INFO = -6
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -8
   ELSE IF( LDB < MAX( 1, P ) ) THEN
      INFO = -10
   ELSE IF( LDU < 1 .OR. ( WANTU .AND. LDU < M ) ) THEN
      INFO = -16
   ELSE IF( LDV < 1 .OR. ( WANTV .AND. LDV < P ) ) THEN
      INFO = -18
   ELSE IF( LDQ < 1 .OR. ( WANTQ .AND. LDQ < N ) ) THEN
      INFO = -20
   ELSE IF( LWORK < 1 .AND. .NOT.LQUERY ) THEN
      INFO = -24
   END IF
!
!     Compute workspace
!
   IF( INFO == 0 ) THEN
      CALL CGEQP3( P, N, B, LDB, IWORK, TAU, WORK, -1, RWORK, INFO )
      LWKOPT = INT( WORK ( 1 ) )
      IF( WANTV ) LWKOPT = MAX( LWKOPT, P )
      LWKOPT = MAX( LWKOPT, MIN( N, P ) )
      LWKOPT = MAX( LWKOPT, M )
      IF( WANTQ ) LWKOPT = MAX( LWKOPT, N )
      CALL CGEQP3( M, N, A, LDA, IWORK, TAU, WORK, -1, RWORK, INFO )
      LWKOPT = MAX( LWKOPT, INT( WORK ( 1 ) ) )
      LWKOPT = MAX( 1, LWKOPT )
      WORK( 1 ) = CMPLX( LWKOPT )
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CGGSVP3', -INFO )
      RETURN
   END IF
   IF( LQUERY ) RETURN
!
!     QR with column pivoting of B: B*P = V*( S11 S12 )
!                                           (  0   0  )
!
   IWORK(1:N) = 0
   CALL CGEQP3( P, N, B, LDB, IWORK, TAU, WORK, LWORK, RWORK, INFO )
!
!     Update A := A*P
!
   CALL CLAPMT( FORWRD, M, N, A, LDA, IWORK )
!
!     Determine the effective rank of matrix B.
!
   L = 0
   DO I = 1, MIN( P, N )
      IF( ABS( B( I, I ) ) > TOLB ) L = L + 1
   ENDDO
!
   IF( WANTV ) THEN
!
!        Copy the details of V, and form V.
!
      V(1:P,1:P) = (0.0E+0,0.0E+0)
      IF( P > 1 ) CALL CLACPY( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ), LDV )
      CALL CUNG2R( P, P, MIN( P, N ), V, LDV, TAU, WORK, INFO )
   END IF
!
!     Clean up B
!
   DO J = 1, L - 1
      B(J+1:L,J) = (0.0E+0,0.0E+0)
   ENDDO
   IF( P > L ) B(L+1:P,1:N) = (0.0E+0,0.0E+0)
!
   IF( WANTQ ) THEN
!
!        Set Q = I and Update Q := Q*P
!
      CALL CLASET( 'Full', N, N, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), Q, LDQ )
      CALL CLAPMT( FORWRD, N, N, Q, LDQ, IWORK )
   END IF
!
   IF( P >= L .AND. N /= L ) THEN
!
!        RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z
!
      CALL CGERQ2( L, N, B, LDB, TAU, WORK, INFO )
!
!        Update A := A*Z**H
!
      CALL CUNMR2( 'Right', 'Conjugate transpose', M, N, L, B, LDB, &
                   TAU, A, LDA, WORK, INFO )
      IF( WANTQ ) THEN
!
!           Update Q := Q*Z**H
!
         CALL CUNMR2( 'Right', 'Conjugate transpose', N, N, L, B, &
                      LDB, TAU, Q, LDQ, WORK, INFO )
      END IF
!
!        Clean up B
!
      B(1:L,1:N-L) = (0.0E+0,0.0E+0)
      DO J = N - L + 1, N
         B(J-N+L+1:L,J) = (0.0E+0,0.0E+0)
      ENDDO
!
   END IF
!
!     Let              N-L     L
!                A = ( A11    A12 ) M,
!
!     then the following does the complete QR decomposition of A11:
!
!              A11 = U*(  0  T12 )*P1**H
!                      (  0   0  )
!
   IWORK(1:N-L) = 0
   CALL CGEQP3( M, N-L, A, LDA, IWORK, TAU, WORK, LWORK, RWORK, INFO )
!
!     Determine the effective rank of A11
!
   K = 0
   DO I = 1, MIN( M, N-L )
      IF( ABS( A( I, I ) ) > TOLA ) K = K + 1
   ENDDO
!
!     Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N )
!
   CALL CUNM2R( 'Left', 'Conjugate transpose', M, L, MIN( M, N-L ), &
                A, LDA, TAU, A( 1, N-L+1 ), LDA, WORK, INFO )
!
   IF( WANTU ) THEN
!
!        Copy the details of U, and form U
!
      U(1:M,1:M) = (0.0E+0,0.0E+0)
      IF( M > 1 ) &
         CALL CLACPY( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
      CALL CUNG2R( M, M, MIN( M, N-L ), U, LDU, TAU, WORK, INFO )
   END IF
!
   IF( WANTQ ) THEN
!
!        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1
!
      CALL CLAPMT( FORWRD, N, N-L, Q, LDQ, IWORK )
   END IF
!
!     Clean up A: set the strictly lower triangular part of
!     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.
!
   DO J = 1, K - 1
      A(J+1:K,J) = (0.0E+0,0.0E+0)
   ENDDO
   IF( M > K ) A(K+1:M,1:N-L) = (0.0E+0,0.0E+0)
!
   IF( N-L > K ) THEN
!
!        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1
!
      CALL CGERQ2( K, N-L, A, LDA, TAU, WORK, INFO )
!
      IF( WANTQ ) THEN
!
!           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H
!
         CALL CUNMR2( 'Right', 'Conjugate transpose', N, N-L, K, A, &
                      LDA, TAU, Q, LDQ, WORK, INFO )
      END IF
!
!        Clean up A
!
      A(1:K,1:N-L-K) = (0.0E+0,0.0E+0)
      DO J = N - L - K + 1, N - L
         A(J-N+L+K+1:K,J) = (0.0E+0,0.0E+0)
      ENDDO
!
   END IF
!
   IF( M > K ) THEN
!
!        QR factorization of A( K+1:M,N-L+1:N )
!
      CALL CGEQR2( M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO )
!
      IF( WANTU ) THEN
!
!           Update U(:,K+1:M) := U(:,K+1:M)*U1
!
         CALL CUNM2R( 'Right', 'No transpose', M, M-K, MIN( M-K, L ), &
                      A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU, &
                      WORK, INFO )
      END IF
!
!        Clean up
!
      DO J = N - L + 1, N
         A(J-N+K+L+1:M, J ) = (0.0E+0,0.0E+0)
      ENDDO
!
   END IF
!
   WORK( 1 ) = CMPLX( LWKOPT )
   RETURN
!
!     End of CGGSVP3
!
END
