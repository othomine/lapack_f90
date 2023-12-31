!> \brief \b CTREVC3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTREVC3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrevc3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrevc3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrevc3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!                           LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
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
!> CTREVC3 computes some or all of the right and/or left eigenvectors of
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
!> input matrix. If Q is the unitary factor that reduces a matrix A to
!> Schur form T, then Q*X and Q*Y are the matrices of right and left
!> eigenvectors of A.
!>
!> This uses a Level 3 BLAS version of the back transformation.
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
!>          The leading dimension of the array VL.
!>          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N.
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
!>          The leading dimension of the array VR.
!>          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N.
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
!>          used to store the eigenvectors.
!>          If HOWMNY = 'A' or 'B', M is set to N.
!>          Each selected eigenvector occupies one column.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of array WORK. LWORK >= max(1,2*N).
!>          For optimum performance, LWORK >= N + 2*N*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (LRWORK)
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of array RWORK. LRWORK >= max(1,N).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the RWORK array, returns
!>          this value as the first entry of the RWORK array, and no error
!>          message related to LRWORK is issued by XERBLA.
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
!> \ingroup trevc3
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
   SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                       LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
   IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          HOWMNY, SIDE
   INTEGER            INFO, LDT, LDVL, LDVR, LWORK, LRWORK, M, MM, N
!     ..
!     .. Array Arguments ..
   LOGICAL            SELECT( * )
   REAL               RWORK( * )
   COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
                      WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NBMIN, NBMAX
   PARAMETER          ( NBMIN = 8, NBMAX = 128 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV
   INTEGER            I, II, IS, J, K, KI, IV, MAXWRK, NB
   REAL               OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
   COMPLEX            CDUM
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   REAL               SLAMCH, CABS1
   EXTERNAL           LSAME, ILAENV, SLAMCH, CABS1
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, CLASET, CGEMM, CGEMV, &
                      CLATRS, CLACPY
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
   BOTHV  = LSAME( SIDE, 'B' )
   RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
   LEFTV  = LSAME( SIDE, 'L' ) .OR. BOTHV
!
   ALLV  = LSAME( HOWMNY, 'A' )
   OVER  = LSAME( HOWMNY, 'B' )
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
   NB = ILAENV( 1, 'CTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
   MAXWRK = MAX( 1, N + 2*N*NB )
   WORK(1) = MAXWRK
   RWORK(1) = MAX( 1, N )
   LQUERY = ( LWORK == -1 .OR. LRWORK == -1 )
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
   ELSE IF( LWORK < MAX( 1, 2*N ) .AND. .NOT.LQUERY ) THEN
      INFO = -14
   ELSE IF ( LRWORK < MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
      INFO = -16
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTREVC3', -INFO )
      RETURN
   ELSE IF( LQUERY ) THEN
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF( N == 0 ) RETURN
!
!     Use blocked version of back-transformation if sufficient workspace.
!     Zero-out the workspace to avoid potential NaN propagation.
!
   IF( OVER .AND. LWORK  >=  N + 2*N*NBMIN ) THEN
      NB = (LWORK - N) / (2*N)
      NB = MIN( NB, NBMAX )
      CALL CLASET( 'F', N, 1+2*NB, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), WORK, N )
   ELSE
      NB = 1
   END IF
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
      WORK( I ) = T( I, I )
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
!        ============================================================
!        Compute right eigenvectors.
!
!        IV is index of column in current block.
!        Non-blocked version always uses IV=NB=1;
!        blocked     version starts with IV=NB, goes down to 1.
!        (Note the "0-th" column is used to store the original diagonal.)
      IV = NB
      IS = M
      DO KI = N, 1, -1
         IF( SOMEV ) THEN
            IF( .NOT.SELECT( KI ) ) GO TO 80
         END IF
         SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
!           --------------------------------------------------------
!           Complex right eigenvector
!
         WORK( KI + IV*N ) = (1.0E+0,0.0E+0)
!
!           Form right-hand side.
!
         WORK(IV*N+1:IV*N+KI-1 ) = -T(1:KI-1, KI )
!
!           Solve upper triangular system:
!           [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.
!
         DO K = 1, KI - 1
            T( K, K ) = T( K, K ) - T( KI, KI )
            IF( CABS1( T( K, K ) ) < SMIN ) T( K, K ) = SMIN
         ENDDO
!
         IF( KI > 1 ) THEN
            CALL CLATRS( 'Upper', 'No transpose', 'Non-unit', 'Y', &
                         KI-1, T, LDT, WORK( 1 + IV*N ), SCALE, &
                         RWORK, INFO )
            WORK( KI + IV*N ) = SCALE
         END IF
!
!           Copy the vector x or Q*x to VR and normalize.
!
         IF( .NOT.OVER ) THEN
!              ------------------------------
!              no back-transform: copy x to VR and normalize.
            VR(1:KI,IS) = WORK(1+IV*N:IV*N+KI)
!
            II = MAXLOC(ABS(REAL(VR(1:KI,IS))) + ABS(AIMAG(VR(1:KI,IS))),1)

            VR(1:KI,IS) = VR(1:KI,IS) / CABS1( VR( II, IS ) )
!
            VR(KI+1:N, IS ) = (0.0E+0,0.0E+0)
!
         ELSE IF( NB == 1 ) THEN
!              ------------------------------
!              version 1: back-transform each vector with GEMV, Q*x.
            IF( KI > 1 ) &
               CALL CGEMV( 'N', N, KI-1, (1.0E+0,0.0E+0), VR, LDVR, &
                           WORK( 1 + IV*N ), 1, CMPLX( SCALE ), VR( 1, KI ), 1 )
!
            II = MAXLOC(ABS(REAL(VR(1:N,KI))) + ABS(AIMAG(VR(1:N,KI))),1)

            VR(1:N,KI) = VR(1:N,KI) / CABS1( VR( II, KI ) )
!
         ELSE
!              ------------------------------
!              version 2: back-transform block of vectors with GEMM
!              zero out below vector
            WORK(IV*N+KI+1:IV*N+N ) = (0.0E+0,0.0E+0)
!
!              Columns IV:NB of work are valid vectors.
!              When the number of vectors stored reaches NB,
!              or if this was last vector, do the GEMM
            IF( (IV == 1) .OR. (KI == 1) ) THEN
               CALL CGEMM( 'N', 'N', N, NB-IV+1, KI+NB-IV, (1.0E+0,0.0E+0), &
                           VR, LDVR, &
                           WORK( 1 + (IV)*N    ), N, &
                           (0.0E+0,0.0E+0), &
                           WORK( 1 + (NB+IV)*N ), N )
!                 normalize vectors
               DO K = IV, NB
                  II = MAXLOC(ABS(REAL(WORK(1+(NB+K)*N:(NB+K+1)*N))) + &
                    ABS(AIMAG(WORK(1+(NB+K)*N:(NB+K+1)*N))),1)
                  WORK(1+(NB+K)*N:(NB+K+1)*N) = WORK(1+(NB+K)*N:(NB+K+1)*N) / &
                    CABS1( WORK( II + (NB+K)*N ) )
               END DO
               CALL CLACPY( 'F', N, NB-IV+1, &
                            WORK( 1 + (NB+IV)*N ), N, &
                            VR( 1, KI ), LDVR )
               IV = NB
            ELSE
               IV = IV - 1
            END IF
         END IF
!
!           Restore the original diagonal elements of T.
!
         DO K = 1, KI - 1
            T( K, K ) = WORK( K )
         ENDDO
!
         IS = IS - 1
80    CONTINUE
      ENDDO
   END IF
!
   IF( LEFTV ) THEN
!
!        ============================================================
!        Compute left eigenvectors.
!
!        IV is index of column in current block.
!        Non-blocked version always uses IV=1;
!        blocked     version starts with IV=1, goes up to NB.
!        (Note the "0-th" column is used to store the original diagonal.)
      IV = 1
      IS = 1
      DO KI = 1, N
!
         IF( SOMEV ) THEN
            IF( .NOT.SELECT( KI ) ) GO TO 130
         END IF
         SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
!           --------------------------------------------------------
!           Complex left eigenvector
!
         WORK( KI + IV*N ) = (1.0E+0,0.0E+0)
!
!           Form right-hand side.
!
         WORK(IV*N+KI+1:(IV+1)*N) = -CONJG( T( KI,KI+1:N) )
!
!           Solve conjugate-transposed triangular system:
!           [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.
!
         DO K = KI + 1, N
            T( K, K ) = T( K, K ) - T( KI, KI )
            IF( CABS1( T( K, K ) ) < SMIN ) T( K, K ) = SMIN
            ENDDO
!
         IF( KI < N ) THEN
            CALL CLATRS( 'Upper', 'Conjugate transpose', 'Non-unit', &
                         'Y', N-KI, T( KI+1, KI+1 ), LDT, &
                         WORK( KI+1 + IV*N ), SCALE, RWORK, INFO )
            WORK( KI + IV*N ) = SCALE
         END IF
!
!           Copy the vector x or Q*x to VL and normalize.
!
         IF( .NOT.OVER ) THEN
!              ------------------------------
!              no back-transform: copy x to VL and normalize.
            VL(KI:N,IS) = WORK(KI+IV*N:(IV+1)*N)
!
            II = MAXLOC(ABS(REAL(VL(KI:N,IS))) + ABS(AIMAG(VL(KI:N,IS))),1) + KI - 1
            VL(KI:N,IS) = VL(KI:N,IS) / CABS1( VL( II, IS ) )
!
            VL(1:KI-1,IS) = (0.0E+0,0.0E+0)
!
         ELSE IF( NB == 1 ) THEN
!              ------------------------------
!              version 1: back-transform each vector with GEMV, Q*x.
            IF( KI < N ) &
               CALL CGEMV( 'N', N, N-KI, (1.0E+0,0.0E+0), VL( 1, KI+1 ), LDVL, &
                           WORK( KI+1 + IV*N ), 1, CMPLX( SCALE ), &
                           VL( 1, KI ), 1 )
!
            II = MAXLOC(ABS(REAL(VL(1:N,KI))) + ABS(AIMAG(VL(1:N,KI))),1)
            VL(1:N,KI) = VL(1:N,KI) / CABS1( VL( II, KI ) )
!
         ELSE
!              ------------------------------
!              version 2: back-transform block of vectors with GEMM
!              zero out above vector
!              could go from KI-NV+1 to KI-1
            WORK(IV*N+1:IV*N+KI-1) = (0.0E+0,0.0E+0)
!
!              Columns 1:IV of work are valid vectors.
!              When the number of vectors stored reaches NB,
!              or if this was last vector, do the GEMM
            IF( (IV == NB) .OR. (KI == N) ) THEN
               CALL CGEMM( 'N', 'N', N, IV, N-KI+IV, (1.0E+0,0.0E+0), &
                           VL( 1, KI-IV+1 ), LDVL, &
                           WORK( KI-IV+1 + (1)*N ), N, &
                           (0.0E+0,0.0E+0), &
                           WORK( 1 + (NB+1)*N ), N )
!                 normalize vectors
               DO K = 1, IV
                  II = MAXLOC(ABS(REAL(WORK(1+(NB+K)*N:(NB+K+1)*N))) + &
                     ABS(AIMAG(WORK(1+(NB+K)*N:(NB+K+1)*N))),1)
                  WORK(1+(NB+K)*N:(NB+K+1)*N) = WORK(1+(NB+K)*N:(NB+K+1)*N) / &
                     CABS1( WORK( II + (NB+K)*N ) )
               END DO
               CALL CLACPY( 'F', N, IV, &
                            WORK( 1 + (NB+1)*N ), N, &
                            VL( 1, KI-IV+1 ), LDVL )
               IV = 1
            ELSE
               IV = IV + 1
            END IF
         END IF
!
!           Restore the original diagonal elements of T.
!
         DO K = KI + 1, N
            T( K, K ) = WORK( K )
         ENDDO
!
         IS = IS + 1
  130    CONTINUE
      ENDDO
   END IF
!
   RETURN
!
!     End of CTREVC3
!
END
