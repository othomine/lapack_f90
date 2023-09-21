!> \brief \b ZGETSLS
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB,
!     $                     WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGETSLS solves overdetermined or underdetermined complex linear systems
!> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
!> factorization of A.  It is assumed that A has full rank.
!>
!>
!>
!> The following options are provided:
!>
!> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
!>    an overdetermined system, i.e., solve the least squares problem
!>                 minimize || B - A*X ||.
!>
!> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!>    an underdetermined system A * X = B.
!>
!> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
!>    an undetermined system A**T * X = B.
!>
!> 4. If TRANS = 'C' and m < n:  find the least squares solution of
!>    an overdetermined system, i.e., solve the least squares problem
!>                 minimize || B - A**T * X ||.
!>
!> Several right hand side vectors b and solution vectors x can be
!> handled in a single call; they are stored as the columns of the
!> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!> matrix X.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': the linear system involves A;
!>          = 'C': the linear system involves A**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of
!>          columns of the matrices B and X. NRHS >=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          A is overwritten by details of its QR or LQ
!>          factorization as returned by ZGEQR or ZGELQ.
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
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the matrix B of right hand side vectors, stored
!>          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!>          if TRANS = 'C'.
!>          On exit, if INFO = 0, B is overwritten by the solution
!>          vectors, stored columnwise:
!>          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!>          squares solution vectors.
!>          if TRANS = 'N' and m < n, rows 1 to N of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'C' and m >= n, rows 1 to M of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'C' and m < n, rows 1 to M of B contain the
!>          least squares solution vectors.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= MAX(1,M,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
!>          or optimal, if query was assumed) LWORK.
!>          See LWORK for details.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If LWORK = -1 or -2, then a workspace query is assumed.
!>          If LWORK = -1, the routine calculates optimal size of WORK for the
!>          optimal performance and returns this value in WORK(1).
!>          If LWORK = -2, the routine calculates minimal size of WORK and
!>          returns this value in WORK(1).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO =  i, the i-th diagonal element of the
!>                triangular factor of A is zero, so that A does not have
!>                full rank; the least squares solution could not be
!>                computed.
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
!> \ingroup getsls
!
!  =====================================================================
   SUBROUTINE ZGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, &
                       WORK, LWORK, INFO )
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANS
   INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
   COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
   COMPLEX*16         CZERO
   PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY, TRAN
   INTEGER            I, IASCL, IBSCL, J, MAXMN, BROW, &
                      SCLLEN, TSZO, TSZM, LWO, LWM, LW1, LW2, &
                      WSIZEO, WSIZEM, INFO2
   DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM, DUM( 1 )
   COMPLEX*16         TQ( 5 ), WORKQ( 1 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, ZLANGE
   EXTERNAL           LSAME, DLAMCH, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGEQR, ZGEMQR, ZLASCL, ZLASET, &
                      ZTRTRS, XERBLA, ZGELQ, ZGEMLQ
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MAX, MIN, INT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
   INFO = 0
   MAXMN = MAX( M, N )
   TRAN  = LSAME( TRANS, 'C' )
!
   LQUERY = ( LWORK == -1 .OR. LWORK == -2 )
   IF( .NOT.( LSAME( TRANS, 'N' ) .OR. &
       LSAME( TRANS, 'C' ) ) ) THEN
      INFO = -1
   ELSE IF( M < 0 ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -3
   ELSE IF( NRHS < 0 ) THEN
      INFO = -4
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -6
   ELSE IF( LDB < MAX( 1, M, N ) ) THEN
      INFO = -8
   END IF
!
   IF( INFO == 0 ) THEN
!
!     Determine the optimum and minimum LWORK
!
    IF( M >= N ) THEN
      CALL ZGEQR( M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 )
      TSZO = INT( TQ( 1 ) )
      LWO  = INT( WORKQ( 1 ) )
      CALL ZGEMQR( 'L', TRANS, M, NRHS, N, A, LDA, TQ, &
                   TSZO, B, LDB, WORKQ, -1, INFO2 )
      LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
      CALL ZGEQR( M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 )
      TSZM = INT( TQ( 1 ) )
      LWM  = INT( WORKQ( 1 ) )
      CALL ZGEMQR( 'L', TRANS, M, NRHS, N, A, LDA, TQ, &
                   TSZM, B, LDB, WORKQ, -1, INFO2 )
      LWM = MAX( LWM, INT( WORKQ( 1 ) ) )
      WSIZEO = TSZO + LWO
      WSIZEM = TSZM + LWM
    ELSE
      CALL ZGELQ( M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 )
      TSZO = INT( TQ( 1 ) )
      LWO  = INT( WORKQ( 1 ) )
      CALL ZGEMLQ( 'L', TRANS, N, NRHS, M, A, LDA, TQ, &
                   TSZO, B, LDB, WORKQ, -1, INFO2 )
      LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
      CALL ZGELQ( M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 )
      TSZM = INT( TQ( 1 ) )
      LWM  = INT( WORKQ( 1 ) )
      CALL ZGEMLQ( 'L', TRANS, N, NRHS, M, A, LDA, TQ, &
                   TSZM, B, LDB, WORKQ, -1, INFO2 )
      LWM  = MAX( LWM, INT( WORKQ( 1 ) ) )
      WSIZEO = TSZO + LWO
      WSIZEM = TSZM + LWM
    END IF
!
    IF( ( LWORK < WSIZEM ).AND.( .NOT.LQUERY ) ) THEN
       INFO = -10
    END IF
!
    WORK( 1 ) = DBLE( WSIZEO )
!
   END IF
!
   IF( INFO /= 0 ) THEN
     CALL XERBLA( 'ZGETSLS', -INFO )
     RETURN
   END IF
   IF( LQUERY ) THEN
     IF( LWORK == -2 ) WORK( 1 ) = DBLE( WSIZEM )
     RETURN
   END IF
   IF( LWORK < WSIZEO ) THEN
     LW1 = TSZM
     LW2 = LWM
   ELSE
     LW1 = TSZO
     LW2 = LWO
   END IF
!
!     Quick return if possible
!
   IF( MIN( M, N, NRHS ) == 0 ) THEN
        CALL ZLASET( 'FULL', MAX( M, N ), NRHS, CZERO, CZERO, &
                     B, LDB )
        RETURN
   END IF
!
!     Get machine parameters
!
    SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
    BIGNUM = ONE / SMLNUM
!
!     Scale A, B if max element outside range [SMLNUM,BIGNUM]
!
   ANRM = ZLANGE( 'M', M, N, A, LDA, DUM )
   IASCL = 0
   IF( ANRM > ZERO .AND. ANRM < SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
      CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
      IASCL = 1
   ELSE IF( ANRM > BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
      CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
      IASCL = 2
   ELSE IF( ANRM == ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
      CALL ZLASET( 'F', MAXMN, NRHS, CZERO, CZERO, B, LDB )
      GO TO 50
   END IF
!
   BROW = M
   IF ( TRAN ) THEN
     BROW = N
   END IF
   BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, DUM )
   IBSCL = 0
   IF( BNRM > ZERO .AND. BNRM < SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
      CALL ZLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, &
                   INFO )
      IBSCL = 1
   ELSE IF( BNRM > BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
      CALL ZLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, &
                   INFO )
      IBSCL = 2
   END IF
!
   IF ( M >= N ) THEN
!
!        compute QR factorization of A
!
     CALL ZGEQR( M, N, A, LDA, WORK( LW2+1 ), LW1, &
                 WORK( 1 ), LW2, INFO )
     IF ( .NOT.TRAN ) THEN
!
!           Least-Squares Problem min || A * X - B ||
!
!           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
!
       CALL ZGEMQR( 'L' , 'C', M, NRHS, N, A, LDA, &
                    WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, &
                    INFO )
!
!           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
!
       CALL ZTRTRS( 'U', 'N', 'N', N, NRHS, &
                     A, LDA, B, LDB, INFO )
       IF( INFO > 0 ) THEN
         RETURN
       END IF
       SCLLEN = N
     ELSE
!
!           Overdetermined system of equations A**T * X = B
!
!           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
!
         CALL ZTRTRS( 'U', 'C', 'N', N, NRHS, &
                      A, LDA, B, LDB, INFO )
!
         IF( INFO > 0 ) THEN
            RETURN
         END IF
!
!           B(N+1:M,1:NRHS) = CZERO
!
         DO J = 1, NRHS
            DO I = N + 1, M
               B( I, J ) = CZERO
            ENDDO
         ENDDO
!
!           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
!
         CALL ZGEMQR( 'L', 'N', M, NRHS, N, A, LDA, &
                      WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, &
                      INFO )
!
         SCLLEN = M
!
      END IF
!
   ELSE
!
!        Compute LQ factorization of A
!
      CALL ZGELQ( M, N, A, LDA, WORK( LW2+1 ), LW1, &
                  WORK( 1 ), LW2, INFO )
!
!        workspace at least M, optimally M*NB.
!
      IF( .NOT.TRAN ) THEN
!
!           underdetermined system of equations A * X = B
!
!           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
!
         CALL ZTRTRS( 'L', 'N', 'N', M, NRHS, &
                      A, LDA, B, LDB, INFO )
!
         IF( INFO > 0 ) THEN
            RETURN
         END IF
!
!           B(M+1:N,1:NRHS) = 0
!
         DO J = 1, NRHS
            DO I = M + 1, N
               B( I, J ) = CZERO
            ENDDO
         ENDDO
!
!           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
!
         CALL ZGEMLQ( 'L', 'C', N, NRHS, M, A, LDA, &
                      WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, &
                      INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
         SCLLEN = N
!
      ELSE
!
!           overdetermined system min || A**T * X - B ||
!
!           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
!
         CALL ZGEMLQ( 'L', 'N', N, NRHS, M, A, LDA, &
                      WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, &
                      INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
!
         CALL ZTRTRS( 'L', 'C', 'N', M, NRHS, &
                      A, LDA, B, LDB, INFO )
!
         IF( INFO > 0 ) THEN
            RETURN
         END IF
!
         SCLLEN = M
!
      END IF
!
   END IF
!
!     Undo scaling
!
   IF( IASCL == 1 ) THEN
     CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, &
                  INFO )
   ELSE IF( IASCL == 2 ) THEN
     CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, &
                   INFO )
   END IF
   IF( IBSCL == 1 ) THEN
     CALL ZLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, &
                  INFO )
   ELSE IF( IBSCL == 2 ) THEN
     CALL ZLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, &
                  INFO )
   END IF
!
50 CONTINUE
   WORK( 1 ) = DBLE( TSZO + LWO )
   RETURN
!
!     End of ZGETSLS
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        