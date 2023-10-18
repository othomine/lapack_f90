!> \brief \b CGEQP3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEQP3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqp3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqp3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqp3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEQP3 computes a QR factorization with column pivoting of a
!> matrix A:  A*P = Q*R  using Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the upper triangle of the array contains the
!>          min(M,N)-by-N upper trapezoidal matrix R; the elements below
!>          the diagonal, together with the array TAU, represent the
!>          unitary matrix Q as a product of min(M,N) elementary
!>          reflectors.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
!>          to the front of A*P (a leading column); if JPVT(J)=0,
!>          the J-th column of A is a free column.
!>          On exit, if JPVT(J)=K, then the J-th column of A*P was the
!>          the K-th column of A.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= N+1.
!>          For optimal performance LWORK >= ( N+1 )*NB, where NB
!>          is the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit.
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup geqp3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a real/complex vector
!>  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
!>  A(i+1:m,i), and tau in TAU(i).
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!>
!  =====================================================================
   SUBROUTINE CGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, &
                      INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            JPVT( * )
   REAL               RWORK( * )
   COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            INB, INBMIN, IXOVER
   PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY
   INTEGER            FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB, &
                      NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEQRF, CLAQP2, CLAQPS, CSWAP, CUNMQR, XERBLA
!     ..
!     .. External Functions ..
   INTEGER            ILAENV
   REAL               SCNRM2
   EXTERNAL           ILAENV, SCNRM2
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          INT, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!  ====================
!
   INFO = 0
   LQUERY = ( LWORK == -1 )
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -4
   END IF
!
   IF( INFO == 0 ) THEN
      MINMN = MIN( M, N )
      IF( MINMN == 0 ) THEN
         IWS = 1
         LWKOPT = 1
      ELSE
         IWS = N + 1
         NB = ILAENV( INB, 'CGEQRF', ' ', M, N, -1, -1 )
         LWKOPT = ( N + 1 )*NB
      END IF
      WORK( 1 ) = CMPLX( LWKOPT )
!
      IF( ( LWORK < IWS ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
   END IF
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CGEQP3', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Move initial columns up front.
!
   NFXD = 1
   DO J = 1, N
      IF( JPVT( J ) /= 0 ) THEN
         IF( J /= NFXD ) THEN
            CALL CSWAP( M, A( 1, J ), 1, A( 1, NFXD ), 1 )
            JPVT( J ) = JPVT( NFXD )
            JPVT( NFXD ) = J
         ELSE
            JPVT( J ) = J
         END IF
         NFXD = NFXD + 1
      ELSE
         JPVT( J ) = J
      END IF
   ENDDO
   NFXD = NFXD - 1
!
!     Factorize fixed columns
!  =======================
!
!     Compute the QR factorization of fixed columns and update
!     remaining columns.
!
   IF( NFXD > 0 ) THEN
      NA = MIN( M, NFXD )
!CC      CALL CGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
      CALL CGEQRF( M, NA, A, LDA, TAU, WORK, LWORK, INFO )
      IWS = MAX( IWS, INT( WORK( 1 ) ) )
      IF( NA < N ) THEN
!CC         CALL CUNM2R( 'Left', 'Conjugate Transpose', M, N-NA,
!CC  $                   NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK,
!CC  $                   INFO )
         CALL CUNMQR( 'Left', 'Conjugate Transpose', M, N-NA, NA, A, &
                      LDA, TAU, A( 1, NA+1 ), LDA, WORK, LWORK, &
                      INFO )
         IWS = MAX( IWS, INT( WORK( 1 ) ) )
      END IF
   END IF
!
!     Factorize free columns
!  ======================
!
   IF( NFXD < MINMN ) THEN
!
      SM = M - NFXD
      SN = N - NFXD
      SMINMN = MINMN - NFXD
!
!        Determine the block size.
!
      NB = ILAENV( INB, 'CGEQRF', ' ', SM, SN, -1, -1 )
      NBMIN = 2
      NX = 0
!
      IF( ( NB > 1 ) .AND. ( NB < SMINMN ) ) THEN
!
!           Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( IXOVER, 'CGEQRF', ' ', SM, SN, -1, &
              -1 ) )
!
!
         IF( NX < SMINMN ) THEN
!
!              Determine if workspace is large enough for blocked code.
!
            MINWS = ( SN+1 )*NB
            IWS = MAX( IWS, MINWS )
            IF( LWORK < MINWS ) THEN
!
!                 Not enough workspace to use optimal NB: Reduce NB and
!                 determine the minimum value of NB.
!
               NB = LWORK / ( SN+1 )
               NBMIN = MAX( 2, ILAENV( INBMIN, 'CGEQRF', ' ', SM, SN, &
                       -1, -1 ) )
!
!
            END IF
         END IF
      END IF
!
!        Initialize partial column norms. The first N elements of work
!        store the exact column norms.
!
      DO J = NFXD + 1, N
         RWORK( J ) = SCNRM2( SM, A( NFXD+1, J ), 1 )
         RWORK( N+J ) = RWORK( J )
      ENDDO
!
      IF( ( NB >= NBMIN ) .AND. ( NB < SMINMN ) .AND. &
          ( NX < SMINMN ) ) THEN
!
!           Use blocked code initially.
!
         J = NFXD + 1
!
!           Compute factorization: while loop.
!
!
         TOPBMN = MINMN - NX
30       CONTINUE
         IF( J <= TOPBMN ) THEN
            JB = MIN( NB, TOPBMN-J+1 )
!
!              Factorize JB columns among columns J:N.
!
            CALL CLAQPS( M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA, &
                         JPVT( J ), TAU( J ), RWORK( J ), &
                         RWORK( N+J ), WORK( 1 ), WORK( JB+1 ), &
                         N-J+1 )
!
            J = J + FJB
            GO TO 30
         END IF
      ELSE
         J = NFXD + 1
      END IF
!
!        Use unblocked code to factor the last or only block.
!
!
      IF( J <= MINMN ) &
         CALL CLAQP2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ), &
                      TAU( J ), RWORK( J ), RWORK( N+J ), WORK( 1 ) )
!
   END IF
!
   WORK( 1 ) = CMPLX( LWKOPT )
   RETURN
!
!     End of CGEQP3
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
