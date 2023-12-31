!> \brief \b CTREVC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTREVC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrevc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrevc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrevc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!                          LDVR, MM, M, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       REAL               RWORK( * )
!       COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTREVC computes some or all of the right and/or left eigenvectors of
!> a complex upper triangular matrix T.
!> Matrices of this type are produced by the Schur factorization of
!> a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR.
!>
!> The right eigenvector x and the left eigenvector y of T corresponding
!> to an eigenvalue w are defined by:
!>
!>              T*x = w*x,     (y**H)*T = w*(y**H)
!>
!> where y**H denotes the conjugate transpose of the vector y.
!> The eigenvalues are not input to this routine, but are read directly
!> from the diagonal of T.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
!> input matrix.  If Q is the unitary factor that reduces a matrix A to
!> Schur form T, then Q*X and Q*Y are the matrices of right and left
!> eigenvectors of A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  compute right eigenvectors only;
!>          = 'L':  compute left eigenvectors only;
!>          = 'B':  compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A':  compute all right and/or left eigenvectors;
!>          = 'B':  compute all right and/or left eigenvectors,
!>                  backtransformed using the matrices supplied in
!>                  VR and/or VL;
!>          = 'S':  compute selected right and/or left eigenvectors,
!>                  as indicated by the logical array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!>          computed.
!>          The eigenvector corresponding to the j-th eigenvalue is
!>          computed if SELECT(j) = .TRUE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,N)
!>          The upper triangular matrix T.  T is modified, but restored
!>          on exit.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is COMPLEX array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by CHSEQR).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VL, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.  LDVL >= 1, and if
!>          SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is COMPLEX array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by CHSEQR).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*X;
!>          if HOWMNY = 'S', the right eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VR, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1, and if
!>          SIDE = 'R' or 'B'; LDVR >= N.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER
!>          The number of columns in the arrays VL and/or VR. MM >= M.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns in the arrays VL and/or VR actually
!>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
!>          is set to N.  Each selected eigenvector occupies one
!>          column.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup trevc
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The algorithm used in this program is basically backward (forward)
!>  substitution, with scaling to make the the code robust against
!>  possible overflow.
!>
!>  Each eigenvector is normalized so that the element of largest
!>  magnitude has magnitude 1; here the magnitude of a complex number
!>  (x,y) is taken to be |x| + |y|.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                      LDVR, MM, M, WORK, RWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          HOWMNY, SIDE
   INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!     ..
!     .. Array Arguments ..
   LOGICAL            SELECT( * )
   REAL               RWORK( * )
   COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
                      WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV
   INTEGER            I, II, IS, J, K, KI
   REAL               OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
   COMPLEX            CDUM
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLAMCH, CABS1
   EXTERNAL           LSAME, SLAMCH, CABS1
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMV, CLATRS, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
   BOTHV = LSAME( SIDE, 'B' )
   RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
   LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
!
   ALLV = LSAME( HOWMNY, 'A' )
   OVER = LSAME( HOWMNY, 'B' )
   SOMEV = LSAME( HOWMNY, 'S' )
!
!     Set M to the number of columns required to store the selected
!     eigenvectors.
!
   IF( SOMEV ) THEN
      M = COUNT(SELECT(1:N))
   ELSE
      M = N
   END IF
!
   INFO = 0
   IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
      INFO = -1
   ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( LDT < MAX( 1, N ) ) THEN
      INFO = -6
   ELSE IF( LDVL < 1 .OR. ( LEFTV .AND. LDVL < N ) ) THEN
      INFO = -8
   ELSE IF( LDVR < 1 .OR. ( RIGHTV .AND. LDVR < N ) ) THEN
      INFO = -10
   ELSE IF( MM < M ) THEN
      INFO = -11
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTREVC', -INFO )
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF( N == 0 ) RETURN
!
!     Set the constants to control overflow.
!
   UNFL = SLAMCH( 'Safe minimum' )
   OVFL = 1.0E+0 / UNFL
   ULP = SLAMCH( 'Precision' )
   SMLNUM = UNFL*( N / ULP )
!
!     Store the diagonal elements of T in working array WORK.
!
   DO I = 1, N
      WORK( I+N ) = T( I, I )
   ENDDO
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
   RWORK( 1 ) = 0.0E+0
   DO J = 2, N
      RWORK( J ) = sum(ABS(REAL(T(1:J-1,J))) + ABS(AIMAG(T(1:J-1,J))))
   ENDDO
!
   IF( RIGHTV ) THEN
!
!        Compute right eigenvectors.
!
      IS = M
      DO KI = N, 1, -1
!
         IF( SOMEV ) THEN
            IF( .NOT.SELECT( KI ) ) GO TO 80
         END IF
         SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
         WORK( 1 ) = (1.0E+0,0.0E+0)
!
!           Form right-hand side.
!
         WORK(1:KI-1) = -T(1:KI-1,KI)
!
!           Solve the triangular system:
!              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
!
         DO K = 1, KI - 1
            T( K, K ) = T( K, K ) - T( KI, KI )
            IF( CABS1( T( K, K ) ) < SMIN ) T( K, K ) = SMIN
         ENDDO
!
         IF( KI > 1 ) THEN
            CALL CLATRS( 'Upper', 'No transpose', 'Non-unit', 'Y', &
                         KI-1, T, LDT, WORK( 1 ), SCALE, RWORK, INFO )
            WORK( KI ) = SCALE
         END IF
!
!           Copy the vector x or Q*x to VR and normalize.
!
         IF( .NOT.OVER ) THEN
            VR(1:KI,IS) = WORK(1:KI)
!
            II = maxloc(ABS(REAL(VR(1:KI,IS))) + ABS(AIMAG(VR(1:KI,IS))),1)

            REMAX = 1.0E+0 / CABS1( VR( II, IS ) )
            VR(1:KI,IS) = REMAX*VR(1:KI,IS)
!
            VR(KI+1:N,IS) = (0.0E+0,0.0E+0)
         ELSE
            IF( KI > 1 ) &
               CALL CGEMV( 'N', N, KI-1, (1.0E+0,0.0E+0), VR, LDVR, WORK( 1 ), &
                           1, CMPLX( SCALE ), VR( 1, KI ), 1 )
!
            II = maxloc(ABS(REAL(VR(1:N,KI))) + ABS(AIMAG(VR(1:N,KI))),1)
            REMAX = 1.0E+0 / CABS1( VR( II, KI ) )
            VR(1:N,KI) = REMAX*VR(1:N,KI)
         END IF
!
!           Set back the original diagonal elements of T.
!
         DO K = 1, KI - 1
            T( K, K ) = WORK( K+N )
         ENDDO
!
         IS = IS - 1
80    CONTINUE
      ENDDO
   END IF
!
   IF( LEFTV ) THEN
!
!        Compute left eigenvectors.
!
      IS = 1
      DO KI = 1, N
!
         IF( SOMEV ) THEN
            IF( .NOT.SELECT( KI ) ) GO TO 130
         END IF
         SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
         WORK( N ) = (1.0E+0,0.0E+0)
!
!           Form right-hand side.
!
         WORK(KI+1:N) = -CONJG( T( KI, KI+1:N) )
!
!           Solve the triangular system:
!              (T(KI+1:N,KI+1:N) - T(KI,KI))**H*X = SCALE*WORK.
!
         DO K = KI + 1, N
            T( K, K ) = T( K, K ) - T( KI, KI )
            IF( CABS1( T( K, K ) ) < SMIN ) T( K, K ) = SMIN
         ENDDO
!
         IF( KI < N ) THEN
            CALL CLATRS( 'Upper', 'Conjugate transpose', 'Non-unit', &
                         'Y', N-KI, T( KI+1, KI+1 ), LDT, &
                         WORK( KI+1 ), SCALE, RWORK, INFO )
            WORK( KI ) = SCALE
         END IF
!
!           Copy the vector x or Q*x to VL and normalize.
!
         IF( .NOT.OVER ) THEN
            VL(KI:N,IS) = WORK(KI:N)
!
            II = maxloc(ABS(REAL(VL(KI:N,IS))) + ABS(AIMAG(VL(KI:N,IS))),1) + KI - 1

            VL(KI:N,IS) = VL(KI:N,IS) / CABS1( VL( II, IS ) )
!
            VL(1:KI-1, IS ) = (0.0E+0,0.0E+0)
         ELSE
            IF( KI < N ) &
               CALL CGEMV( 'N', N, N-KI, (1.0E+0,0.0E+0), VL( 1, KI+1 ), LDVL, &
                           WORK( KI+1 ), 1, CMPLX( SCALE ), VL( 1, KI ), 1 )
!
            II = maxloc(ABS(REAL(VL(1:N,KI))) + ABS(AIMAG(VL(1:N,KI))),1)
            VL(1:N,KI) = VL(1:N,KI) / CABS1( VL( II, KI ) )
         END IF
!
!           Set back the original diagonal elements of T.
!
         DO K = KI + 1, N
            T( K, K ) = WORK( K+N )
         ENDDO
!
         IS = IS + 1
  130    CONTINUE
      ENDDO
   END IF
!
   RETURN
!
!     End of CTREVC
!
END
