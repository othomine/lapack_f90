!> \brief \b STGEVC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STGEVC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgevc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgevc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgevc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
!                          LDVL, VR, LDVR, MM, M, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       REAL               P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
!      $                   VR( LDVR, * ), WORK( * )
!       ..
!
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STGEVC computes some or all of the right and/or left eigenvectors of
!> a pair of real matrices (S,P), where S is a quasi-triangular matrix
!> and P is upper triangular.  Matrix pairs of this type are produced by
!> the generalized Schur factorization of a matrix pair (A,B):
!>
!>    A = Q*S*Z**T,  B = Q*P*Z**T
!>
!> as computed by SGGHRD + SHGEQZ.
!>
!> The right eigenvector x and the left eigenvector y of (S,P)
!> corresponding to an eigenvalue w are defined by:
!>
!>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
!>
!> where y**H denotes the conjugate transpose of y.
!> The eigenvalues are not input to this routine, but are computed
!> directly from the diagonal blocks of S and P.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
!> where Z and Q are input matrices.
!> If Q and Z are the orthogonal factors from the generalized Schur
!> factorization of a matrix pair (A,B), then Z*X and Q*Y
!> are the matrices of right and left eigenvectors of (A,B).
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R': compute right eigenvectors only;
!>          = 'L': compute left eigenvectors only;
!>          = 'B': compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A': compute all right and/or left eigenvectors;
!>          = 'B': compute all right and/or left eigenvectors,
!>                 backtransformed by the matrices in VR and/or VL;
!>          = 'S': compute selected right and/or left eigenvectors,
!>                 specified by the logical array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY='S', SELECT specifies the eigenvectors to be
!>          computed.  If w(j) is a real eigenvalue, the corresponding
!>          real eigenvector is computed if SELECT(j) is .TRUE..
!>          If w(j) and w(j+1) are the real and imaginary parts of a
!>          complex eigenvalue, the corresponding complex eigenvector
!>          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,
!>          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is
!>          set to .FALSE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices S and P.  N >= 0.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (LDS,N)
!>          The upper quasi-triangular matrix S from a generalized Schur
!>          factorization, as computed by SHGEQZ.
!> \endverbatim
!>
!> \param[in] LDS
!> \verbatim
!>          LDS is INTEGER
!>          The leading dimension of array S.  LDS >= max(1,N).
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is REAL array, dimension (LDP,N)
!>          The upper triangular matrix P from a generalized Schur
!>          factorization, as computed by SHGEQZ.
!>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks
!>          of S must be in positive diagonal form.
!> \endverbatim
!>
!> \param[in] LDP
!> \verbatim
!>          LDP is INTEGER
!>          The leading dimension of array P.  LDP >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is REAL array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!>          of left Schur vectors returned by SHGEQZ).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
!>                      SELECT, stored consecutively in the columns of
!>                      VL, in the same order as their eigenvalues.
!>
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part, and the second the imaginary part.
!>
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of array VL.  LDVL >= 1, and if
!>          SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is REAL array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Z (usually the orthogonal matrix Z
!>          of right Schur vectors returned by SHGEQZ).
!>
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);
!>          if HOWMNY = 'B' or 'b', the matrix Z*X;
!>          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)
!>                      specified by SELECT, stored consecutively in the
!>                      columns of VR, in the same order as their
!>                      eigenvalues.
!>
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part and the second the imaginary part.
!>
!>          Not referenced if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1, and if
!>          SIDE = 'R' or 'B', LDVR >= N.
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
!>          is set to N.  Each selected real eigenvector occupies one
!>          column and each selected complex eigenvector occupies two
!>          columns.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (6*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex
!>                eigenvalue.
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
!> \ingroup tgevc
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Allocation of workspace:
!>  ---------- -- ---------
!>
!>     WORK( j ) = 1-norm of j-th column of A, above the diagonal
!>     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal
!>     WORK( 2*N+1:3*N ) = real part of eigenvector
!>     WORK( 3*N+1:4*N ) = imaginary part of eigenvector
!>     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector
!>     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector
!>
!>  Rowwise vs. columnwise solution methods:
!>  ------- --  ---------- -------- -------
!>
!>  Finding a generalized eigenvector consists basically of solving the
!>  singular triangular system
!>
!>   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)
!>
!>  Consider finding the i-th right eigenvector (assume all eigenvalues
!>  are real). The equation to be solved is:
!>       n                   i
!>  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1
!>      k=j                 k=j
!>
!>  where  C = (A - w B)  (The components v(i+1:n) are 0.)
!>
!>  The "rowwise" method is:
!>
!>  (1)  v(i) := 1
!>  for j = i-1,. . .,1:
!>                          i
!>      (2) compute  s = - sum C(j,k) v(k)   and
!>                        k=j+1
!>
!>      (3) v(j) := s / C(j,j)
!>
!>  Step 2 is sometimes called the "dot product" step, since it is an
!>  inner product between the j-th row and the portion of the eigenvector
!>  that has been computed so far.
!>
!>  The "columnwise" method consists basically in doing the sums
!>  for all the rows in parallel.  As each v(j) is computed, the
!>  contribution of v(j) times the j-th column of C is added to the
!>  partial sums.  Since FORTRAN arrays are stored columnwise, this has
!>  the advantage that at each step, the elements of C that are accessed
!>  are adjacent to one another, whereas with the rowwise method, the
!>  elements accessed at a step are spaced LDS (and LDP) words apart.
!>
!>  When finding left eigenvectors, the matrix in question is the
!>  transpose of the one in storage, so the rowwise method then
!>  actually accesses columns of A and B at each step, and so is the
!>  preferred method.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE STGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, &
                      LDVL, VR, LDVR, MM, M, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          HOWMNY, SIDE
   INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
!     ..
!     .. Array Arguments ..
   LOGICAL            SELECT( * )
   REAL               P( LDP, * ), S( LDS, * ), VL( LDVL, * ), &
                      VR( LDVR, * ), WORK( * )
!     ..
!
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE, SAFETY
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, &
                      SAFETY = 1.0E+2 )
!     ..
!     .. Local Scalars ..
   LOGICAL            COMPL, COMPR, IL2BY2, ILABAD, ILALL, ILBACK, &
                      ILBBAD, ILCOMP, ILCPLX, LSA, LSB
   INTEGER            I, IBEG, IEIG, IEND, IHWMNY, IINFO, IM, ISIDE, &
                      J, JA, JC, JE, JR, JW, NA, NW
   REAL               ACOEF, ACOEFA, ANORM, ASCALE, BCOEFA, BCOEFI, &
                      BCOEFR, BIG, BIGNUM, BNORM, BSCALE, CIM2A, &
                      CIM2B, CIMAGA, CIMAGB, CRE2A, CRE2B, CREALA, &
                      CREALB, DMIN, SAFMIN, SALFAR, SBETA, SCALE, &
                      SMALL, TEMP, TEMP2, TEMP2I, TEMP2R, ULP, XMAX, &
                      XSCALE
!     ..
!     .. Local Arrays ..
   REAL               BDIAG( 2 ), SUM( 2, 2 ), SUMS( 2, 2 ), &
                      SUMP( 2, 2 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLAMCH
   EXTERNAL           LSAME, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMV, SLACPY, SLAG2, SLALN2, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
   IF( LSAME( HOWMNY, 'A' ) ) THEN
      IHWMNY = 1
      ILALL = .TRUE.
      ILBACK = .FALSE.
   ELSE IF( LSAME( HOWMNY, 'S' ) ) THEN
      IHWMNY = 2
      ILALL = .FALSE.
      ILBACK = .FALSE.
   ELSE IF( LSAME( HOWMNY, 'B' ) ) THEN
      IHWMNY = 3
      ILALL = .TRUE.
      ILBACK = .TRUE.
   ELSE
      IHWMNY = -1
      ILALL = .TRUE.
   END IF
!
   IF( LSAME( SIDE, 'R' ) ) THEN
      ISIDE = 1
      COMPL = .FALSE.
      COMPR = .TRUE.
   ELSE IF( LSAME( SIDE, 'L' ) ) THEN
      ISIDE = 2
      COMPL = .TRUE.
      COMPR = .FALSE.
   ELSE IF( LSAME( SIDE, 'B' ) ) THEN
      ISIDE = 3
      COMPL = .TRUE.
      COMPR = .TRUE.
   ELSE
      ISIDE = -1
   END IF
!
   INFO = 0
   IF( ISIDE < 0 ) THEN
      INFO = -1
   ELSE IF( IHWMNY < 0 ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -4
   ELSE IF( LDS < MAX( 1, N ) ) THEN
      INFO = -6
   ELSE IF( LDP < MAX( 1, N ) ) THEN
      INFO = -8
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'STGEVC', -INFO )
      RETURN
   END IF
!
!     Count the number of eigenvectors to be computed
!
   IF( .NOT.ILALL ) THEN
      IM = 0
      ILCPLX = .FALSE.
      DO J = 1, N
         IF( ILCPLX ) THEN
            ILCPLX = .FALSE.
            GO TO 10
         END IF
         IF( J < N ) THEN
            IF( S( J+1, J ) /= ZERO ) &
               ILCPLX = .TRUE.
         END IF
         IF( ILCPLX ) THEN
            IF( SELECT( J ) .OR. SELECT( J+1 ) ) &
               IM = IM + 2
         ELSE
            IF( SELECT( J ) ) &
               IM = IM + 1
         END IF
10    CONTINUE
      ENDDO
   ELSE
      IM = N
   END IF
!
!     Check 2-by-2 diagonal blocks of A, B
!
   ILABAD = .FALSE.
   ILBBAD = .FALSE.
   DO J = 1, N - 1
      IF( S( J+1, J ) /= ZERO ) THEN
         IF( P( J, J ) == ZERO .OR. P( J+1, J+1 ) == ZERO .OR. &
             P( J, J+1 ) /= ZERO )ILBBAD = .TRUE.
         IF( J < N-1 ) THEN
            IF( S( J+2, J+1 ) /= ZERO ) &
               ILABAD = .TRUE.
         END IF
      END IF
   ENDDO
!
   IF( ILABAD ) THEN
      INFO = -5
   ELSE IF( ILBBAD ) THEN
      INFO = -7
   ELSE IF( COMPL .AND. LDVL < N .OR. LDVL < 1 ) THEN
      INFO = -10
   ELSE IF( COMPR .AND. LDVR < N .OR. LDVR < 1 ) THEN
      INFO = -12
   ELSE IF( MM < IM ) THEN
      INFO = -13
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'STGEVC', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   M = IM
   IF( N == 0 ) &
      RETURN
!
!     Machine Constants
!
   SAFMIN = SLAMCH( 'Safe minimum' )
   BIG = ONE / SAFMIN
   ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
   SMALL = SAFMIN*N / ULP
   BIG = ONE / SMALL
   BIGNUM = ONE / ( SAFMIN*N )
!
!     Compute the 1-norm of each column of the strictly upper triangular
!     part (i.e., excluding all elements belonging to the diagonal
!     blocks) of A and B to check for possible overflow in the
!     triangular solver.
!
   ANORM = ABS( S( 1, 1 ) )
   IF( N > 1 ) &
      ANORM = ANORM + ABS( S( 2, 1 ) )
   BNORM = ABS( P( 1, 1 ) )
   WORK( 1 ) = ZERO
   WORK( N+1 ) = ZERO
!
   DO J = 2, N
      TEMP = ZERO
      TEMP2 = ZERO
      IF( S( J, J-1 ) == ZERO ) THEN
         IEND = J - 1
      ELSE
         IEND = J - 2
      END IF
      DO I = 1, IEND
         TEMP = TEMP + ABS( S( I, J ) )
         TEMP2 = TEMP2 + ABS( P( I, J ) )
      ENDDO
      WORK( J ) = TEMP
      WORK( N+J ) = TEMP2
      DO I = IEND + 1, MIN( J+1, N )
         TEMP = TEMP + ABS( S( I, J ) )
         TEMP2 = TEMP2 + ABS( P( I, J ) )
      ENDDO
      ANORM = MAX( ANORM, TEMP )
      BNORM = MAX( BNORM, TEMP2 )
   ENDDO
!
   ASCALE = ONE / MAX( ANORM, SAFMIN )
   BSCALE = ONE / MAX( BNORM, SAFMIN )
!
!     Left eigenvectors
!
   IF( COMPL ) THEN
      IEIG = 0
!
!        Main loop over eigenvalues
!
      ILCPLX = .FALSE.
      DO JE = 1, N
!
!           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
!           (b) this would be the second of a complex pair.
!           Check for complex eigenvalue, so as to be sure of which
!           entry(-ies) of SELECT to look at.
!
         IF( ILCPLX ) THEN
            ILCPLX = .FALSE.
            GO TO 220
         END IF
         NW = 1
         IF( JE < N ) THEN
            IF( S( JE+1, JE ) /= ZERO ) THEN
               ILCPLX = .TRUE.
               NW = 2
            END IF
         END IF
         IF( ILALL ) THEN
            ILCOMP = .TRUE.
         ELSE IF( ILCPLX ) THEN
            ILCOMP = SELECT( JE ) .OR. SELECT( JE+1 )
         ELSE
            ILCOMP = SELECT( JE )
         END IF
         IF( .NOT.ILCOMP ) &
            GO TO 220
!
!           Decide if (a) singular pencil, (b) real eigenvalue, or
!           (c) complex eigenvalue.
!
         IF( .NOT.ILCPLX ) THEN
            IF( ABS( S( JE, JE ) ) <= SAFMIN .AND. &
                ABS( P( JE, JE ) ) <= SAFMIN ) THEN
!
!                 Singular matrix pencil -- return unit eigenvector
!
               IEIG = IEIG + 1
               DO JR = 1, N
                  VL( JR, IEIG ) = ZERO
               ENDDO
               VL( IEIG, IEIG ) = ONE
               GO TO 220
            END IF
         END IF
!
!           Clear vector
!
         DO JR = 1, NW*N
            WORK( 2*N+JR ) = ZERO
         ENDDO
!                                                 T
!           Compute coefficients in  ( a A - b B )  y = 0
!              a  is  ACOEF
!              b  is  BCOEFR + i*BCOEFI
!
         IF( .NOT.ILCPLX ) THEN
!
!              Real eigenvalue
!
            TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE, &
                   ABS( P( JE, JE ) )*BSCALE, SAFMIN )
            SALFAR = ( TEMP*S( JE, JE ) )*ASCALE
            SBETA = ( TEMP*P( JE, JE ) )*BSCALE
            ACOEF = SBETA*ASCALE
            BCOEFR = SALFAR*BSCALE
            BCOEFI = ZERO
!
!              Scale to avoid underflow
!
            SCALE = ONE
            LSA = ABS( SBETA ) >= SAFMIN .AND. ABS( ACOEF ) < SMALL
            LSB = ABS( SALFAR ) >= SAFMIN .AND. ABS( BCOEFR ) < &
                  SMALL
            IF( LSA ) &
               SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
            IF( LSB ) &
               SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )* &
                       MIN( BNORM, BIG ) )
            IF( LSA .OR. LSB ) THEN
               SCALE = MIN( SCALE, ONE / &
                       ( SAFMIN*MAX( ONE, ABS( ACOEF ), &
                       ABS( BCOEFR ) ) ) )
               IF( LSA ) THEN
                  ACOEF = ASCALE*( SCALE*SBETA )
               ELSE
                  ACOEF = SCALE*ACOEF
               END IF
               IF( LSB ) THEN
                  BCOEFR = BSCALE*( SCALE*SALFAR )
               ELSE
                  BCOEFR = SCALE*BCOEFR
               END IF
            END IF
            ACOEFA = ABS( ACOEF )
            BCOEFA = ABS( BCOEFR )
!
!              First component is 1
!
            WORK( 2*N+JE ) = ONE
            XMAX = ONE
         ELSE
!
!              Complex eigenvalue
!
            CALL SLAG2( S( JE, JE ), LDS, P( JE, JE ), LDP, &
                        SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2, &
                        BCOEFI )
            BCOEFI = -BCOEFI
            IF( BCOEFI == ZERO ) THEN
               INFO = JE
               RETURN
            END IF
!
!              Scale to avoid over/underflow
!
            ACOEFA = ABS( ACOEF )
            BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
            SCALE = ONE
            IF( ACOEFA*ULP < SAFMIN .AND. ACOEFA >= SAFMIN ) &
               SCALE = ( SAFMIN / ULP ) / ACOEFA
            IF( BCOEFA*ULP < SAFMIN .AND. BCOEFA >= SAFMIN ) &
               SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA )
            IF( SAFMIN*ACOEFA > ASCALE ) &
               SCALE = ASCALE / ( SAFMIN*ACOEFA )
            IF( SAFMIN*BCOEFA > BSCALE ) &
               SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
            IF( SCALE /= ONE ) THEN
               ACOEF = SCALE*ACOEF
               ACOEFA = ABS( ACOEF )
               BCOEFR = SCALE*BCOEFR
               BCOEFI = SCALE*BCOEFI
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
            END IF
!
!              Compute first two components of eigenvector
!
            TEMP = ACOEF*S( JE+1, JE )
            TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE )
            TEMP2I = -BCOEFI*P( JE, JE )
            IF( ABS( TEMP ) > ABS( TEMP2R )+ABS( TEMP2I ) ) THEN
               WORK( 2*N+JE ) = ONE
               WORK( 3*N+JE ) = ZERO
               WORK( 2*N+JE+1 ) = -TEMP2R / TEMP
               WORK( 3*N+JE+1 ) = -TEMP2I / TEMP
            ELSE
               WORK( 2*N+JE+1 ) = ONE
               WORK( 3*N+JE+1 ) = ZERO
               TEMP = ACOEF*S( JE, JE+1 )
               WORK( 2*N+JE ) = ( BCOEFR*P( JE+1, JE+1 )-ACOEF* &
                                S( JE+1, JE+1 ) ) / TEMP
               WORK( 3*N+JE ) = BCOEFI*P( JE+1, JE+1 ) / TEMP
            END IF
            XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ), &
                   ABS( WORK( 2*N+JE+1 ) )+ABS( WORK( 3*N+JE+1 ) ) )
         END IF
!
         DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
!
!                                           T
!           Triangular solve of  (a A - b B)  y = 0
!
!                                   T
!           (rowwise in  (a A - b B) , or columnwise in (a A - b B) )
!
         IL2BY2 = .FALSE.
!
         DO J = JE + NW, N
            IF( IL2BY2 ) THEN
               IL2BY2 = .FALSE.
               GO TO 160
            END IF
!
            NA = 1
            BDIAG( 1 ) = P( J, J )
            IF( J < N ) THEN
               IF( S( J+1, J ) /= ZERO ) THEN
                  IL2BY2 = .TRUE.
                  BDIAG( 2 ) = P( J+1, J+1 )
                  NA = 2
               END IF
            END IF
!
!              Check whether scaling is necessary for dot products
!
            XSCALE = ONE / MAX( ONE, XMAX )
            TEMP = MAX( WORK( J ), WORK( N+J ), &
                   ACOEFA*WORK( J )+BCOEFA*WORK( N+J ) )
            IF( IL2BY2 ) &
               TEMP = MAX( TEMP, WORK( J+1 ), WORK( N+J+1 ), &
                      ACOEFA*WORK( J+1 )+BCOEFA*WORK( N+J+1 ) )
            IF( TEMP > BIGNUM*XSCALE ) THEN
               DO JW = 0, NW - 1
                  DO JR = JE, J - 1
                     WORK( ( JW+2 )*N+JR ) = XSCALE* &
                        WORK( ( JW+2 )*N+JR )
                  ENDDO
               ENDDO
               XMAX = XMAX*XSCALE
            END IF
!
!              Compute dot products
!
!                    j-1
!              SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
!                    k=je
!
!              To reduce the op count, this is done as
!
!              _        j-1                  _        j-1
!              a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) )
!                       k=je                          k=je
!
!              which may cause underflow problems if A or B are close
!              to underflow.  (E.g., less than SMALL.)
!
!
            DO JW = 1, NW
               DO JA = 1, NA
                  SUMS( JA, JW ) = ZERO
                  SUMP( JA, JW ) = ZERO
!
                  DO JR = JE, J - 1
                     SUMS( JA, JW ) = SUMS( JA, JW ) + &
                                      S( JR, J+JA-1 )* &
                                      WORK( ( JW+1 )*N+JR )
                     SUMP( JA, JW ) = SUMP( JA, JW ) + &
                                      P( JR, J+JA-1 )* &
                                      WORK( ( JW+1 )*N+JR )
                     ENDDO
                  ENDDO
               ENDDO
!
            DO JA = 1, NA
               IF( ILCPLX ) THEN
                  SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) + &
                                 BCOEFR*SUMP( JA, 1 ) - &
                                 BCOEFI*SUMP( JA, 2 )
                  SUM( JA, 2 ) = -ACOEF*SUMS( JA, 2 ) + &
                                 BCOEFR*SUMP( JA, 2 ) + &
                                 BCOEFI*SUMP( JA, 1 )
               ELSE
                  SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) + &
                                 BCOEFR*SUMP( JA, 1 )
               END IF
               ENDDO
!
!                                  T
!              Solve  ( a A - b B )  y = SUM(,)
!              with scaling and perturbation of the denominator
!
            CALL SLALN2( .TRUE., NA, NW, DMIN, ACOEF, S( J, J ), LDS, &
                         BDIAG( 1 ), BDIAG( 2 ), SUM, 2, BCOEFR, &
                         BCOEFI, WORK( 2*N+J ), N, SCALE, TEMP, &
                         IINFO )
            IF( SCALE < ONE ) THEN
               DO JW = 0, NW - 1
                  DO JR = JE, J - 1
                     WORK( ( JW+2 )*N+JR ) = SCALE* &
                        WORK( ( JW+2 )*N+JR )
                     ENDDO
                  ENDDO
               XMAX = SCALE*XMAX
            END IF
            XMAX = MAX( XMAX, TEMP )
  160       CONTINUE
            ENDDO
!
!           Copy eigenvector to VL, back transforming if
!           HOWMNY='B'.
!
         IEIG = IEIG + 1
         IF( ILBACK ) THEN
            DO JW = 0, NW - 1
               CALL SGEMV( 'N', N, N+1-JE, ONE, VL( 1, JE ), LDVL, &
                           WORK( ( JW+2 )*N+JE ), 1, ZERO, &
                           WORK( ( JW+4 )*N+1 ), 1 )
               ENDDO
            CALL SLACPY( ' ', N, NW, WORK( 4*N+1 ), N, VL( 1, JE ), &
                         LDVL )
            IBEG = 1
         ELSE
            CALL SLACPY( ' ', N, NW, WORK( 2*N+1 ), N, VL( 1, IEIG ), &
                         LDVL )
            IBEG = JE
         END IF
!
!           Scale eigenvector
!
         XMAX = ZERO
         IF( ILCPLX ) THEN
            DO J = IBEG, N
               XMAX = MAX( XMAX, ABS( VL( J, IEIG ) )+ &
                      ABS( VL( J, IEIG+1 ) ) )
               ENDDO
         ELSE
            DO J = IBEG, N
               XMAX = MAX( XMAX, ABS( VL( J, IEIG ) ) )
               ENDDO
         END IF
!
         IF( XMAX > SAFMIN ) THEN
            XSCALE = ONE / XMAX
!
            DO JW = 0, NW - 1
               DO JR = IBEG, N
                  VL( JR, IEIG+JW ) = XSCALE*VL( JR, IEIG+JW )
                  ENDDO
               ENDDO
         END IF
         IEIG = IEIG + NW - 1
!
  220    CONTINUE
         ENDDO
   END IF
!
!     Right eigenvectors
!
   IF( COMPR ) THEN
      IEIG = IM + 1
!
!        Main loop over eigenvalues
!
      ILCPLX = .FALSE.
      DO JE = N, 1, -1
!
!           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
!           (b) this would be the second of a complex pair.
!           Check for complex eigenvalue, so as to be sure of which
!           entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
!           or SELECT(JE-1).
!           If this is a complex pair, the 2-by-2 diagonal block
!           corresponding to the eigenvalue is in rows/columns JE-1:JE
!
         IF( ILCPLX ) THEN
            ILCPLX = .FALSE.
            GO TO 500
         END IF
         NW = 1
         IF( JE > 1 ) THEN
            IF( S( JE, JE-1 ) /= ZERO ) THEN
               ILCPLX = .TRUE.
               NW = 2
            END IF
         END IF
         IF( ILALL ) THEN
            ILCOMP = .TRUE.
         ELSE IF( ILCPLX ) THEN
            ILCOMP = SELECT( JE ) .OR. SELECT( JE-1 )
         ELSE
            ILCOMP = SELECT( JE )
         END IF
         IF( .NOT.ILCOMP ) &
            GO TO 500
!
!           Decide if (a) singular pencil, (b) real eigenvalue, or
!           (c) complex eigenvalue.
!
         IF( .NOT.ILCPLX ) THEN
            IF( ABS( S( JE, JE ) ) <= SAFMIN .AND. &
                ABS( P( JE, JE ) ) <= SAFMIN ) THEN
!
!                 Singular matrix pencil -- unit eigenvector
!
               IEIG = IEIG - 1
               DO JR = 1, N
                  VR( JR, IEIG ) = ZERO
                  ENDDO
               VR( IEIG, IEIG ) = ONE
               GO TO 500
            END IF
         END IF
!
!           Clear vector
!
         DO JW = 0, NW - 1
            DO JR = 1, N
               WORK( ( JW+2 )*N+JR ) = ZERO
               ENDDO
            ENDDO
!
!           Compute coefficients in  ( a A - b B ) x = 0
!              a  is  ACOEF
!              b  is  BCOEFR + i*BCOEFI
!
         IF( .NOT.ILCPLX ) THEN
!
!              Real eigenvalue
!
            TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE, &
                   ABS( P( JE, JE ) )*BSCALE, SAFMIN )
            SALFAR = ( TEMP*S( JE, JE ) )*ASCALE
            SBETA = ( TEMP*P( JE, JE ) )*BSCALE
            ACOEF = SBETA*ASCALE
            BCOEFR = SALFAR*BSCALE
            BCOEFI = ZERO
!
!              Scale to avoid underflow
!
            SCALE = ONE
            LSA = ABS( SBETA ) >= SAFMIN .AND. ABS( ACOEF ) < SMALL
            LSB = ABS( SALFAR ) >= SAFMIN .AND. ABS( BCOEFR ) < &
                  SMALL
            IF( LSA ) &
               SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
            IF( LSB ) &
               SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )* &
                       MIN( BNORM, BIG ) )
            IF( LSA .OR. LSB ) THEN
               SCALE = MIN( SCALE, ONE / &
                       ( SAFMIN*MAX( ONE, ABS( ACOEF ), &
                       ABS( BCOEFR ) ) ) )
               IF( LSA ) THEN
                  ACOEF = ASCALE*( SCALE*SBETA )
               ELSE
                  ACOEF = SCALE*ACOEF
               END IF
               IF( LSB ) THEN
                  BCOEFR = BSCALE*( SCALE*SALFAR )
               ELSE
                  BCOEFR = SCALE*BCOEFR
               END IF
            END IF
            ACOEFA = ABS( ACOEF )
            BCOEFA = ABS( BCOEFR )
!
!              First component is 1
!
            WORK( 2*N+JE ) = ONE
            XMAX = ONE
!
!              Compute contribution from column JE of A and B to sum
!              (See "Further Details", above.)
!
            DO JR = 1, JE - 1
               WORK( 2*N+JR ) = BCOEFR*P( JR, JE ) - &
                                ACOEF*S( JR, JE )
               ENDDO
         ELSE
!
!              Complex eigenvalue
!
            CALL SLAG2( S( JE-1, JE-1 ), LDS, P( JE-1, JE-1 ), LDP, &
                        SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2, &
                        BCOEFI )
            IF( BCOEFI == ZERO ) THEN
               INFO = JE - 1
               RETURN
            END IF
!
!              Scale to avoid over/underflow
!
            ACOEFA = ABS( ACOEF )
            BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
            SCALE = ONE
            IF( ACOEFA*ULP < SAFMIN .AND. ACOEFA >= SAFMIN ) &
               SCALE = ( SAFMIN / ULP ) / ACOEFA
            IF( BCOEFA*ULP < SAFMIN .AND. BCOEFA >= SAFMIN ) &
               SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA )
            IF( SAFMIN*ACOEFA > ASCALE ) &
               SCALE = ASCALE / ( SAFMIN*ACOEFA )
            IF( SAFMIN*BCOEFA > BSCALE ) &
               SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
            IF( SCALE /= ONE ) THEN
               ACOEF = SCALE*ACOEF
               ACOEFA = ABS( ACOEF )
               BCOEFR = SCALE*BCOEFR
               BCOEFI = SCALE*BCOEFI
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
            END IF
!
!              Compute first two components of eigenvector
!              and contribution to sums
!
            TEMP = ACOEF*S( JE, JE-1 )
            TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE )
            TEMP2I = -BCOEFI*P( JE, JE )
            IF( ABS( TEMP ) >= ABS( TEMP2R )+ABS( TEMP2I ) ) THEN
               WORK( 2*N+JE ) = ONE
               WORK( 3*N+JE ) = ZERO
               WORK( 2*N+JE-1 ) = -TEMP2R / TEMP
               WORK( 3*N+JE-1 ) = -TEMP2I / TEMP
            ELSE
               WORK( 2*N+JE-1 ) = ONE
               WORK( 3*N+JE-1 ) = ZERO
               TEMP = ACOEF*S( JE-1, JE )
               WORK( 2*N+JE ) = ( BCOEFR*P( JE-1, JE-1 )-ACOEF* &
                                S( JE-1, JE-1 ) ) / TEMP
               WORK( 3*N+JE ) = BCOEFI*P( JE-1, JE-1 ) / TEMP
            END IF
!
            XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ), &
                   ABS( WORK( 2*N+JE-1 ) )+ABS( WORK( 3*N+JE-1 ) ) )
!
!              Compute contribution from columns JE and JE-1
!              of A and B to the sums.
!
            CREALA = ACOEF*WORK( 2*N+JE-1 )
            CIMAGA = ACOEF*WORK( 3*N+JE-1 )
            CREALB = BCOEFR*WORK( 2*N+JE-1 ) - &
                     BCOEFI*WORK( 3*N+JE-1 )
            CIMAGB = BCOEFI*WORK( 2*N+JE-1 ) + &
                     BCOEFR*WORK( 3*N+JE-1 )
            CRE2A = ACOEF*WORK( 2*N+JE )
            CIM2A = ACOEF*WORK( 3*N+JE )
            CRE2B = BCOEFR*WORK( 2*N+JE ) - BCOEFI*WORK( 3*N+JE )
            CIM2B = BCOEFI*WORK( 2*N+JE ) + BCOEFR*WORK( 3*N+JE )
            DO JR = 1, JE - 2
               WORK( 2*N+JR ) = -CREALA*S( JR, JE-1 ) + &
                                CREALB*P( JR, JE-1 ) - &
                                CRE2A*S( JR, JE ) + CRE2B*P( JR, JE )
               WORK( 3*N+JR ) = -CIMAGA*S( JR, JE-1 ) + &
                                CIMAGB*P( JR, JE-1 ) - &
                                CIM2A*S( JR, JE ) + CIM2B*P( JR, JE )
               ENDDO
         END IF
!
         DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
!
!           Columnwise triangular solve of  (a A - b B)  x = 0
!
         IL2BY2 = .FALSE.
         DO J = JE - NW, 1, -1
!
!              If a 2-by-2 block, is in position j-1:j, wait until
!              next iteration to process it (when it will be j:j+1)
!
            IF( .NOT.IL2BY2 .AND. J > 1 ) THEN
               IF( S( J, J-1 ) /= ZERO ) THEN
                  IL2BY2 = .TRUE.
                  GO TO 370
               END IF
            END IF
            BDIAG( 1 ) = P( J, J )
            IF( IL2BY2 ) THEN
               NA = 2
               BDIAG( 2 ) = P( J+1, J+1 )
            ELSE
               NA = 1
            END IF
!
!              Compute x(j) (and x(j+1), if 2-by-2 block)
!
            CALL SLALN2( .FALSE., NA, NW, DMIN, ACOEF, S( J, J ), &
                         LDS, BDIAG( 1 ), BDIAG( 2 ), WORK( 2*N+J ), &
                         N, BCOEFR, BCOEFI, SUM, 2, SCALE, TEMP, &
                         IINFO )
            IF( SCALE < ONE ) THEN
!
               DO JW = 0, NW - 1
                  DO JR = 1, JE
                     WORK( ( JW+2 )*N+JR ) = SCALE* &
                        WORK( ( JW+2 )*N+JR )
                     ENDDO
                  ENDDO
            END IF
            XMAX = MAX( SCALE*XMAX, TEMP )
!
            DO JW = 1, NW
               DO JA = 1, NA
                  WORK( ( JW+1 )*N+J+JA-1 ) = SUM( JA, JW )
                  ENDDO
               ENDDO
!
!              w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
!
            IF( J > 1 ) THEN
!
!                 Check whether scaling is necessary for sum.
!
               XSCALE = ONE / MAX( ONE, XMAX )
               TEMP = ACOEFA*WORK( J ) + BCOEFA*WORK( N+J )
               IF( IL2BY2 ) &
                  TEMP = MAX( TEMP, ACOEFA*WORK( J+1 )+BCOEFA* &
                         WORK( N+J+1 ) )
               TEMP = MAX( TEMP, ACOEFA, BCOEFA )
               IF( TEMP > BIGNUM*XSCALE ) THEN
!
                  DO JW = 0, NW - 1
                     DO JR = 1, JE
                        WORK( ( JW+2 )*N+JR ) = XSCALE* &
                           WORK( ( JW+2 )*N+JR )
                        ENDDO
                     ENDDO
                  XMAX = XMAX*XSCALE
               END IF
!
!                 Compute the contributions of the off-diagonals of
!                 column j (and j+1, if 2-by-2 block) of A and B to the
!                 sums.
!
!
               DO JA = 1, NA
                  IF( ILCPLX ) THEN
                     CREALA = ACOEF*WORK( 2*N+J+JA-1 )
                     CIMAGA = ACOEF*WORK( 3*N+J+JA-1 )
                     CREALB = BCOEFR*WORK( 2*N+J+JA-1 ) - &
                              BCOEFI*WORK( 3*N+J+JA-1 )
                     CIMAGB = BCOEFI*WORK( 2*N+J+JA-1 ) + &
                              BCOEFR*WORK( 3*N+J+JA-1 )
                     DO JR = 1, J - 1
                        WORK( 2*N+JR ) = WORK( 2*N+JR ) - &
                                         CREALA*S( JR, J+JA-1 ) + &
                                         CREALB*P( JR, J+JA-1 )
                        WORK( 3*N+JR ) = WORK( 3*N+JR ) - &
                                         CIMAGA*S( JR, J+JA-1 ) + &
                                         CIMAGB*P( JR, J+JA-1 )
                        ENDDO
                  ELSE
                     CREALA = ACOEF*WORK( 2*N+J+JA-1 )
                     CREALB = BCOEFR*WORK( 2*N+J+JA-1 )
                     DO JR = 1, J - 1
                        WORK( 2*N+JR ) = WORK( 2*N+JR ) - &
                                         CREALA*S( JR, J+JA-1 ) + &
                                         CREALB*P( JR, J+JA-1 )
                        ENDDO
                  END IF
                  ENDDO
            END IF
!
            IL2BY2 = .FALSE.
  370       CONTINUE
            ENDDO
!
!           Copy eigenvector to VR, back transforming if
!           HOWMNY='B'.
!
         IEIG = IEIG - NW
         IF( ILBACK ) THEN
!
            DO JW = 0, NW - 1
               DO JR = 1, N
                  WORK( ( JW+4 )*N+JR ) = WORK( ( JW+2 )*N+1 )* &
                                          VR( JR, 1 )
                  ENDDO
!
!                 A series of compiler directives to defeat
!                 vectorization for the next loop
!
!
               DO JC = 2, JE
                  DO JR = 1, N
                     WORK( ( JW+4 )*N+JR ) = WORK( ( JW+4 )*N+JR ) + &
                        WORK( ( JW+2 )*N+JC )*VR( JR, JC )
                     ENDDO
                  ENDDO
               ENDDO
!
            DO JW = 0, NW - 1
               DO JR = 1, N
                  VR( JR, IEIG+JW ) = WORK( ( JW+4 )*N+JR )
                  ENDDO
               ENDDO
!
            IEND = N
         ELSE
            DO JW = 0, NW - 1
               DO JR = 1, N
                  VR( JR, IEIG+JW ) = WORK( ( JW+2 )*N+JR )
                  ENDDO
               ENDDO
!
            IEND = JE
         END IF
!
!           Scale eigenvector
!
         XMAX = ZERO
         IF( ILCPLX ) THEN
            DO J = 1, IEND
               XMAX = MAX( XMAX, ABS( VR( J, IEIG ) )+ &
                      ABS( VR( J, IEIG+1 ) ) )
               ENDDO
         ELSE
            DO J = 1, IEND
               XMAX = MAX( XMAX, ABS( VR( J, IEIG ) ) )
               ENDDO
         END IF
!
         IF( XMAX > SAFMIN ) THEN
            XSCALE = ONE / XMAX
            DO JW = 0, NW - 1
               DO JR = 1, IEND
                  VR( JR, IEIG+JW ) = XSCALE*VR( JR, IEIG+JW )
                  ENDDO
               ENDDO
         END IF
  500    CONTINUE
         ENDDO
   END IF
!
   RETURN
!
!     End of STGEVC
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
