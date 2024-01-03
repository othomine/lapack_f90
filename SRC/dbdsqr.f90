!> \brief \b DBDSQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DBDSQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,
!                          LDU, C, LDC, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DBDSQR computes the singular values and, optionally, the right and/or
!> left singular vectors from the singular value decomposition (SVD) of
!> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
!> zero-shift QR algorithm.  The SVD of B has the form
!>
!>    B = Q * S * P**T
!>
!> where S is the diagonal matrix of singular values, Q is an orthogonal
!> matrix of left singular vectors, and P is an orthogonal matrix of
!> right singular vectors.  If left singular vectors are requested, this
!> subroutine actually returns U*Q instead of Q, and, if right singular
!> vectors are requested, this subroutine returns P**T*VT instead of
!> P**T, for given real input matrices U and VT.  When U and VT are the
!> orthogonal matrices that reduce a general matrix A to bidiagonal
!> form:  A = U*B*VT, as computed by DGEBRD, then
!>
!>    A = (U*Q) * S * (P**T*VT)
!>
!> is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
!> for a given real input matrix C.
!>
!> See "Computing  Small Singular Values of Bidiagonal Matrices With
!> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
!> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
!> no. 5, pp. 873-912, Sept 1990) and
!> "Accurate singular values and differential qd algorithms," by
!> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
!> Department, University of California at Berkeley, July 1992
!> for a detailed description of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  B is upper bidiagonal;
!>          = 'L':  B is lower bidiagonal.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.  N >= 0.
!> \endverbatim
!>
!> \param[in] NCVT
!> \verbatim
!>          NCVT is INTEGER
!>          The number of columns of the matrix VT. NCVT >= 0.
!> \endverbatim
!>
!> \param[in] NRU
!> \verbatim
!>          NRU is INTEGER
!>          The number of rows of the matrix U. NRU >= 0.
!> \endverbatim
!>
!> \param[in] NCC
!> \verbatim
!>          NCC is INTEGER
!>          The number of columns of the matrix C. NCC >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the n diagonal elements of the bidiagonal matrix B.
!>          On exit, if INFO=0, the singular values of B in decreasing
!>          order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the N-1 offdiagonal elements of the bidiagonal
!>          matrix B.
!>          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
!>          will contain the diagonal and superdiagonal elements of a
!>          bidiagonal matrix orthogonally equivalent to the one given
!>          as input.
!> \endverbatim
!>
!> \param[in,out] VT
!> \verbatim
!>          VT is DOUBLE PRECISION array, dimension (LDVT, NCVT)
!>          On entry, an N-by-NCVT matrix VT.
!>          On exit, VT is overwritten by P**T * VT.
!>          Not referenced if NCVT = 0.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.
!>          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
!> \endverbatim
!>
!> \param[in,out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU, N)
!>          On entry, an NRU-by-N matrix U.
!>          On exit, U is overwritten by U * Q.
!>          Not referenced if NRU = 0.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,NRU).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC, NCC)
!>          On entry, an N-by-NCC matrix C.
!>          On exit, C is overwritten by Q**T * C.
!>          Not referenced if NCC = 0.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.
!>          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (4*(N-1))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  If INFO = -i, the i-th argument had an illegal value
!>          > 0:
!>             if NCVT = NRU = NCC = 0,
!>                = 1, a split was marked by a positive value in E
!>                = 2, current block of Z not diagonalized after 30*N
!>                     iterations (in inner while loop)
!>                = 3, termination criterion of outer while loop not met
!>                     (program created more than N unreduced blocks)
!>             else NCVT = NRU = NCC = 0,
!>                   the algorithm did not converge; D and E contain the
!>                   elements of a bidiagonal matrix which is orthogonally
!>                   similar to the input matrix B;  if INFO = i, i
!>                   elements of E have not converged to zero.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
!>          TOLMUL controls the convergence criterion of the QR loop.
!>          If it is positive, TOLMUL*EPS is the desired relative
!>             precision in the computed singular values.
!>          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
!>             desired absolute accuracy in the computed singular
!>             values (corresponds to relative accuracy
!>             abs(TOLMUL*EPS) in the largest singular value.
!>          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
!>             between 10 (for fast convergence) and .1/EPS
!>             (for there to be some accuracy in the results).
!>          Default is to lose at either one eighth or 2 of the
!>             available decimal digits in each computed singular value
!>             (whichever is smaller).
!>
!>  MAXITR  INTEGER, default = 6
!>          MAXITR controls the maximum number of passes of the
!>          algorithm through its inner loop. The algorithms stops
!>          (and so fails to converge) if the number of passes
!>          through the inner loop exceeds MAXITR*N**2.
!>
!> \endverbatim
!
!> \par Note:
!  ===========
!>
!> \verbatim
!>  Bug report from Cezary Dendek.
!>  On March 23rd 2017, the INTEGER variable MAXIT = MAXITR*N**2 is
!>  removed since it can overflow pretty easily (for N larger or equal
!>  than 18,919). We instead use MAXITDIVN = MAXITR*N.
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
!> \ingroup bdsqr
!
!  =====================================================================
   SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
                      LDU, C, LDC, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ), &
                      VT( LDVT, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXITR
   PARAMETER          ( MAXITR = 6 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LOWER, ROTATE
   INTEGER            I, IDIR, ISUB, ITER, ITERDIVN, J, LL, LLL, M, &
                      MAXITDIVN, NM1, NM12, NM13, OLDLL, OLDM
   DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU, &
                      OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL, &
                      SINR, SLL, SMAX, SMIN, SMINOA, &
                      SN, THRESH, TOL, TOLMUL, UNFL
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           LSAME, DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLARTG, DLAS2, DLASQ1, DLASR, DLASV2, DROT, &
                      DSCAL, DSWAP, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   LOWER = LSAME( UPLO, 'L' )
   IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( NCVT < 0 ) THEN
      INFO = -3
   ELSE IF( NRU < 0 ) THEN
      INFO = -4
   ELSE IF( NCC < 0 ) THEN
      INFO = -5
   ELSE IF( ( NCVT == 0 .AND. LDVT < 1 ) .OR. &
            ( NCVT > 0 .AND. LDVT < MAX( 1, N ) ) ) THEN
      INFO = -9
   ELSE IF( LDU < MAX( 1, NRU ) ) THEN
      INFO = -11
   ELSE IF( ( NCC == 0 .AND. LDC < 1 ) .OR. &
            ( NCC > 0 .AND. LDC < MAX( 1, N ) ) ) THEN
      INFO = -13
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DBDSQR', -INFO )
      RETURN
   END IF
   IF( N == 0 ) RETURN
   IF( N == 1 ) GO TO 160
!
!     ROTATE is true if any singular vectors desired, false otherwise
!
   ROTATE = ( NCVT > 0 ) .OR. ( NRU > 0 ) .OR. ( NCC > 0 )
!
!     If no singular vectors desired, use qd algorithm
!
   IF( .NOT.ROTATE ) THEN
      CALL DLASQ1( N, D, E, WORK, INFO )
!
!     If INFO equals 2, dqds didn't finish, try to finish
!
      IF( INFO  /=  2 ) RETURN
      INFO = 0
   END IF
!
   NM1 = N - 1
   NM12 = NM1 + NM1
   NM13 = NM12 + NM1
   IDIR = 0
!
!     Get machine constants
!
   EPS = DLAMCH( 'Epsilon' )
   UNFL = DLAMCH( 'Safe minimum' )
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
   IF( LOWER ) THEN
      DO I = 1, N - 1
         CALL DLARTG( D( I ), E( I ), CS, SN, R )
         D( I ) = R
         E( I ) = SN*D( I+1 )
         D( I+1 ) = CS*D( I+1 )
         WORK( I ) = CS
         WORK( NM1+I ) = SN
      ENDDO
!
!        Update singular vectors if desired
!
      IF( NRU > 0 ) &
         CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U, &
                     LDU )
      IF( NCC > 0 ) &
         CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C, &
                     LDC )
   END IF
!
!     Compute singular values to relative accuracy TOL
!     (By setting TOL to be negative, algorithm will compute
!     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
!
   TOLMUL = MAX( 10.0D0, MIN( 100.0D0, EPS**-0.125D0 ) )
   TOL = TOLMUL*EPS
!
!     Compute approximate maximum, minimum singular values
!
   SMAX = MAXVAL(ABS( D( 1:N ) ))
   SMAX = MAX(SMAX,MAXVAL(ABS( E( 1:N-1 ) )))
   SMIN = 0.0D0
   IF( TOL >= 0.0D0 ) THEN
!
!        Relative accuracy desired
!
      SMINOA = ABS( D( 1 ) )
      IF( SMINOA == 0.0D0 ) GO TO 50
      MU = SMINOA
      DO I = 2, N
         MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
         SMINOA = MIN( SMINOA, MU )
         IF( SMINOA == 0.0D0 ) GO TO 50
      ENDDO
50    CONTINUE
      SMINOA = SMINOA / SQRT( DBLE( N ) )
      THRESH = MAX( TOL*SMINOA, MAXITR*(N*(N*UNFL)) )
   ELSE
!
!        Absolute accuracy desired
!
      THRESH = MAX( ABS( TOL )*SMAX, MAXITR*(N*(N*UNFL)) )
   END IF
!
!     Prepare for main iteration loop for the singular values
!     (MAXIT is the maximum number of passes through the inner
!     loop permitted before nonconvergence signalled.)
!
   MAXITDIVN = MAXITR*N
   ITERDIVN = 0
   ITER = -1
   OLDLL = -1
   OLDM = -1
!
!     M points to last element of unconverged part of matrix
!
   M = N
!
!     Begin main iteration loop
!
60 CONTINUE
!
!     Check for convergence or exceeding iteration count
!
   IF( M <= 1 ) GO TO 160
!
   IF( ITER >= N ) THEN
      ITER = ITER - N
      ITERDIVN = ITERDIVN + 1
      IF( ITERDIVN >= MAXITDIVN ) THEN
!
!     Maximum number of iterations exceeded, failure to converge
!
         INFO = COUNT(E(1:N-1) /= 0.0D0)
         RETURN
      ENDIF
   END IF
!
!     Find diagonal block of matrix to work on
!
   IF( TOL < 0.0D0 .AND. ABS( D( M ) ) <= THRESH ) D( M ) = 0.0D0
   SMAX = ABS( D( M ) )
   DO LLL = 1, M - 1
      LL = M - LLL
      ABSS = ABS( D( LL ) )
      ABSE = ABS( E( LL ) )
      IF( TOL < 0.0D0 .AND. ABSS <= THRESH ) D( LL ) = 0.0D0
      IF( ABSE <= THRESH ) GO TO 80
      SMAX = MAX( SMAX, ABSS, ABSE )
   ENDDO
   LL = 0
   GO TO 90
80 CONTINUE
   E( LL ) = 0.0D0
!
!     Matrix splits since E(LL) = 0
!
   IF( LL == M-1 ) THEN
!
!        Convergence of bottom singular value, return to top of loop
!
      M = M - 1
      GO TO 60
   END IF
90 CONTINUE
   LL = LL + 1
!
!     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
!
   IF( LL == M-1 ) THEN
!
!        2 by 2 block, handle separately
!
      CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR, &
                   COSR, SINL, COSL )
      D( M-1 ) = SIGMX
      E( M-1 ) = 0.0D0
      D( M ) = SIGMN
!
!        Compute singular vectors, if desired
!
      IF( NCVT > 0 ) CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR, SINR )
      IF( NRU > 0 ) CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
      IF( NCC > 0 ) CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL, SINL )
      M = M - 2
      GO TO 60
   END IF
!
!     If working on new submatrix, choose shift direction
!     (from larger end diagonal element towards smaller)
!
   IF( LL > OLDM .OR. M < OLDLL ) THEN
      IF( ABS( D( LL ) ) >= ABS( D( M ) ) ) THEN
!
!           Chase bulge from top (big end) to bottom (small end)
!
         IDIR = 1
      ELSE
!
!           Chase bulge from bottom (big end) to top (small end)
!
         IDIR = 2
      END IF
   END IF
!
!     Apply convergence tests
!
   IF( IDIR == 1 ) THEN
!
!        Run convergence test in forward direction
!        First apply standard test to bottom of matrix
!
      IF( ABS( E( M-1 ) ) <= ABS( TOL )*ABS( D( M ) ) .OR. &
          ( TOL < 0.0D0 .AND. ABS( E( M-1 ) ) <= THRESH ) ) THEN
         E( M-1 ) = 0.0D0
         GO TO 60
      END IF
!
      IF( TOL >= 0.0D0 ) THEN
!
!           If relative accuracy desired,
!           apply convergence criterion forward
!
         MU = ABS( D( LL ) )
         SMIN = MU
         DO LLL = LL, M - 1
            IF( ABS( E( LLL ) ) <= TOL*MU ) THEN
               E( LLL ) = 0.0D0
               GO TO 60
            END IF
            MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
            SMIN = MIN( SMIN, MU )
            ENDDO
      END IF
!
   ELSE
!
!        Run convergence test in backward direction
!        First apply standard test to top of matrix
!
      IF( ABS( E( LL ) ) <= ABS( TOL )*ABS( D( LL ) ) .OR. &
          ( TOL < 0.0D0 .AND. ABS( E( LL ) ) <= THRESH ) ) THEN
         E( LL ) = 0.0D0
         GO TO 60
      END IF
!
      IF( TOL >= 0.0D0 ) THEN
!
!           If relative accuracy desired,
!           apply convergence criterion backward
!
         MU = ABS( D( M ) )
         SMIN = MU
         DO LLL = M - 1, LL, -1
            IF( ABS( E( LLL ) ) <= TOL*MU ) THEN
               E( LLL ) = 0.0D0
               GO TO 60
            END IF
            MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
            SMIN = MIN( SMIN, MU )
            ENDDO
      END IF
   END IF
   OLDLL = LL
   OLDM = M
!
!     Compute shift.  First, test if shifting would ruin relative
!     accuracy, and if so set the shift to zero.
!
   IF( TOL >= 0.0D0 .AND. N*TOL*( SMIN / SMAX ) <= MAX( EPS, 0.01D0*TOL ) ) THEN
!
!        Use a zero shift to avoid loss of relative accuracy
!
      SHIFT = 0.0D0
   ELSE
!
!        Compute the shift from 2-by-2 block at end of matrix
!
      IF( IDIR == 1 ) THEN
         SLL = ABS( D( LL ) )
         CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
      ELSE
         SLL = ABS( D( M ) )
         CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
      END IF
!
!        Test if shift negligible, and if so set to zero
!
      IF( SLL > 0.0D0 ) THEN
         IF( ( SHIFT / SLL )**2 < EPS ) SHIFT = 0.0D0
      END IF
   END IF
!
!     Increment iteration count
!
   ITER = ITER + M - LL
!
!     If SHIFT = 0, do simplified QR iteration
!
   IF( SHIFT == 0.0D0 ) THEN
      IF( IDIR == 1 ) THEN
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
         CS = 1.0D0
         OLDCS = 1.0D0
         DO I = LL, M - 1
            CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
            IF( I > LL ) E( I-1 ) = OLDSN*R
            CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
            WORK( I-LL+1 ) = CS
            WORK( I-LL+1+NM1 ) = SN
            WORK( I-LL+1+NM12 ) = OLDCS
            WORK( I-LL+1+NM13 ) = OLDSN
            ENDDO
         H = D( M )*CS
         D( M ) = H*OLDCS
         E( M-1 ) = H*OLDSN
!
!           Update singular vectors
!
         IF( NCVT > 0 ) &
            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), WORK( N ), VT( LL, 1 ), LDVT )
         IF( NRU > 0 ) &
            CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), WORK( NM13+1 ), U( 1, LL ), LDU )
         IF( NCC > 0 ) &
            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), WORK( NM13+1 ), C( LL, 1 ), LDC )
!
!           Test convergence
!
         IF( ABS( E( M-1 ) ) <= THRESH ) E( M-1 ) = 0.0D0
!
      ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
         CS = 1.0D0
         OLDCS = 1.0D0
         DO I = M, LL + 1, -1
            CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
            IF( I < M ) E( I ) = OLDSN*R
            CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
            WORK( I-LL ) = CS
            WORK( I-LL+NM1 ) = -SN
            WORK( I-LL+NM12 ) = OLDCS
            WORK( I-LL+NM13 ) = -OLDSN
         ENDDO
         H = D( LL )*CS
         D( LL ) = H*OLDCS
         E( LL ) = H*OLDSN
!
!           Update singular vectors
!
         IF( NCVT > 0 ) CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                        WORK( NM13+1 ), VT( LL, 1 ), LDVT )
         IF( NRU > 0 ) CALL DLASR( 'R', 'V', 'B',NRU, M-LL+1, WORK( 1 ), &
                        WORK( N ), U( 1, LL ), LDU )
         IF( NCC > 0 ) CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                        WORK( N ), C( LL, 1 ), LDC )
!
!           Test convergence
!
         IF( ABS( E( LL ) ) <= THRESH ) E( LL ) = 0.0D0
      END IF
   ELSE
!
!        Use nonzero shift
!
      IF( IDIR == 1 ) THEN
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
         F = ( ABS( D( LL ) )-SHIFT )* ( SIGN( 1.0D0, D( LL ) )+SHIFT / D( LL ) )
         G = E( LL )
         DO I = LL, M - 1
            CALL DLARTG( F, G, COSR, SINR, R )
            IF( I > LL ) E( I-1 ) = R
            F = COSR*D( I ) + SINR*E( I )
            E( I ) = COSR*E( I ) - SINR*D( I )
            G = SINR*D( I+1 )
            D( I+1 ) = COSR*D( I+1 )
            CALL DLARTG( F, G, COSL, SINL, R )
            D( I ) = R
            F = COSL*E( I ) + SINL*D( I+1 )
            D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
            IF( I < M-1 ) THEN
               G = SINL*E( I+1 )
               E( I+1 ) = COSL*E( I+1 )
            END IF
            WORK( I-LL+1 ) = COSR
            WORK( I-LL+1+NM1 ) = SINR
            WORK( I-LL+1+NM12 ) = COSL
            WORK( I-LL+1+NM13 ) = SINL
         ENDDO
         E( M-1 ) = F
!
!           Update singular vectors
!
         IF( NCVT > 0 ) CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), &
                        WORK( N ), VT( LL, 1 ), LDVT )
         IF( NRU > 0 ) CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), &
                        WORK( NM13+1 ), U( 1, LL ), LDU )
         IF( NCC > 0 ) CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), &
                        WORK( NM13+1 ), C( LL, 1 ), LDC )
!
!           Test convergence
!
         IF( ABS( E( M-1 ) ) <= THRESH ) E( M-1 ) = 0.0D0
!
      ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
         F = ( ABS( D( M ) )-SHIFT )*( SIGN( 1.0D0, D( M ) )+SHIFT / D( M ) )
         G = E( M-1 )
         DO I = M, LL + 1, -1
            CALL DLARTG( F, G, COSR, SINR, R )
            IF( I < M ) E( I ) = R
            F = COSR*D( I ) + SINR*E( I-1 )
            E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
            G = SINR*D( I-1 )
            D( I-1 ) = COSR*D( I-1 )
            CALL DLARTG( F, G, COSL, SINL, R )
            D( I ) = R
            F = COSL*E( I-1 ) + SINL*D( I-1 )
            D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
            IF( I > LL+1 ) THEN
               G = SINL*E( I-2 )
               E( I-2 ) = COSL*E( I-2 )
            END IF
            WORK( I-LL ) = COSR
            WORK( I-LL+NM1 ) = -SINR
            WORK( I-LL+NM12 ) = COSL
            WORK( I-LL+NM13 ) = -SINL
         ENDDO
         E( LL ) = F
!
!           Test convergence
!
         IF( ABS( E( LL ) ) <= THRESH ) E( LL ) = 0.0D0
!
!           Update singular vectors if desired
!
         IF( NCVT > 0 ) CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                        WORK( NM13+1 ), VT( LL, 1 ), LDVT )
         IF( NRU > 0 ) CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), &
                        WORK( N ), U( 1, LL ), LDU )
         IF( NCC > 0 ) CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                        WORK( N ), C( LL, 1 ), LDC )
      END IF
   END IF
!
!     QR iteration finished, go back and check convergence
!
   GO TO 60
!
!     All singular values converged, so make them positive
!
  160 CONTINUE
   DO I = 1, N
      IF( D( I ) < 0.0D0 ) THEN
         D( I ) = -D( I )
!
!           Change sign of singular vectors, if desired
!
         IF( NCVT > 0 ) VT(I,1:NCVT) = -VT(I,1:NCVT)
      END IF
      ENDDO
!
!     Sort the singular values into decreasing order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
   DO I = 1, N - 1
!
!        Scan for smallest D(I)
!
      ISUB = 1
      SMIN = D( 1 )
      DO J = 2, N + 1 - I
         IF( D( J ) <= SMIN ) THEN
            ISUB = J
            SMIN = D( J )
         END IF
         ENDDO
      IF( ISUB /= N+1-I ) THEN
!
!           Swap singular values and vectors
!
         D( ISUB ) = D( N+1-I )
         D( N+1-I ) = SMIN
         IF( NCVT > 0 ) CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ), LDVT )
         IF( NRU > 0 ) CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
         IF( NCC > 0 ) CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
      END IF
      ENDDO
   RETURN
!
!     End of DBDSQR
!
END
