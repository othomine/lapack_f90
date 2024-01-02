!> \brief \b CTGEVC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTGEVC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgevc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgevc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgevc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
!                          LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       REAL               RWORK( * )
!       COMPLEX            P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
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
!> CTGEVC computes some or all of the right and/or left eigenvectors of
!> a pair of complex matrices (S,P), where S and P are upper triangular.
!> Matrix pairs of this type are produced by the generalized Schur
!> factorization of a complex matrix pair (A,B):
!>
!>    A = Q*S*Z**H,  B = Q*P*Z**H
!>
!> as computed by CGGHRD + CHGEQZ.
!>
!> The right eigenvector x and the left eigenvector y of (S,P)
!> corresponding to an eigenvalue w are defined by:
!>
!>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
!>
!> where y**H denotes the conjugate transpose of y.
!> The eigenvalues are not input to this routine, but are computed
!> directly from the diagonal elements of S and P.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
!> where Z and Q are input matrices.
!> If Q and Z are the unitary factors from the generalized Schur
!> factorization of a matrix pair (A,B), then Z*X and Q*Y
!> are the matrices of right and left eigenvectors of (A,B).
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
!>          computed.  The eigenvector corresponding to the j-th
!>          eigenvalue is computed if SELECT(j) = .TRUE..
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
!>          S is COMPLEX array, dimension (LDS,N)
!>          The upper triangular matrix S from a generalized Schur
!>          factorization, as computed by CHGEQZ.
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
!>          P is COMPLEX array, dimension (LDP,N)
!>          The upper triangular matrix P from a generalized Schur
!>          factorization, as computed by CHGEQZ.  P must have real
!>          diagonal elements.
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
!>          VL is COMPLEX array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q
!>          of left Schur vectors returned by CHGEQZ).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
!>                      SELECT, stored consecutively in the columns of
!>                      VL, in the same order as their eigenvalues.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of array VL.  LDVL >= 1, and if
!>          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is COMPLEX array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Z (usually the unitary matrix Z
!>          of right Schur vectors returned by CHGEQZ).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);
!>          if HOWMNY = 'B', the matrix Z*X;
!>          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by
!>                      SELECT, stored consecutively in the columns of
!>                      VR, in the same order as their eigenvalues.
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
!>          is set to N.  Each selected eigenvector occupies one column.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
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
!>          = 0:  successful exit.
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
!> \ingroup tgevc
!
!  =====================================================================
   SUBROUTINE CTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, &
                      LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
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
   REAL               RWORK( * )
   COMPLEX            P( LDP, * ), S( LDS, * ), VL( LDVL, * ), &
                      VR( LDVR, * ), WORK( * )
!     ..
!
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            COMPL, COMPR, ILALL, ILBACK, ILBBAD, ILCOMP, &
                      LSA, LSB
   INTEGER            I, IBEG, IEIG, IEND, IHWMNY, IM, ISIDE, ISRC, &
                      J, JE, JR
   REAL               ACOEFA, ACOEFF, ANORM, ASCALE, BCOEFA, BIG, &
                      BIGNUM, BNORM, BSCALE, DMIN, SAFMIN, SBETA, &
                      SCALE, SMALL, TEMP, ULP, XMAX
   COMPLEX            BCOEFF, CA, CB, D, SALPHA, SOMME, SUMA, SUMB, X
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLAMCH, CABS1
   COMPLEX            CLADIV
   EXTERNAL           LSAME, SLAMCH, CLADIV, CABS1
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMV, XERBLA

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
      CALL XERBLA( 'CTGEVC', -INFO )
      RETURN
   END IF
!
!     Count the number of eigenvectors
!
   IF( .NOT.ILALL ) THEN
      IM = COUNT(SELECT(1:N))
   ELSE
      IM = N
   END IF
!
!     Check diagonal of B
!
   ILBBAD = .FALSE.
   DO J = 1, N
      IF( AIMAG( P( J, J ) ) /= 0.0E+0 ) ILBBAD = .TRUE.
   ENDDO
!
   IF( ILBBAD ) THEN
      INFO = -7
   ELSE IF( COMPL .AND. LDVL < N .OR. LDVL < 1 ) THEN
      INFO = -10
   ELSE IF( COMPR .AND. LDVR < N .OR. LDVR < 1 ) THEN
      INFO = -12
   ELSE IF( MM < IM ) THEN
      INFO = -13
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTGEVC', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   M = IM
   IF( N == 0 ) RETURN
!
!     Machine Constants
!
   SAFMIN = SLAMCH( 'Safe minimum' )
   BIG = 1.0E+0 / SAFMIN
   ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
   SMALL = SAFMIN*N / ULP
   BIG = 1.0E+0 / SMALL
   BIGNUM = 1.0E+0 / ( SAFMIN*N )
!
!     Compute the 1-norm of each column of the strictly upper triangular
!     part of A and B to check for possible overflow in the triangular
!     solver.
!
   ANORM = CABS1( S( 1, 1 ) )
   BNORM = CABS1( P( 1, 1 ) )
   RWORK( 1 ) = 0.0E+0
   RWORK( N+1 ) = 0.0E+0
   DO J = 2, N
      RWORK( J ) = 0.0E+0
      RWORK( N+J ) = 0.0E+0
      DO I = 1, J - 1
         RWORK( J ) = RWORK( J ) + CABS1( S( I, J ) )
         RWORK( N+J ) = RWORK( N+J ) + CABS1( P( I, J ) )
      ENDDO
      ANORM = MAX( ANORM, RWORK( J )+CABS1( S( J, J ) ) )
      BNORM = MAX( BNORM, RWORK( N+J )+CABS1( P( J, J ) ) )
   ENDDO
!
   ASCALE = 1.0E+0 / MAX( ANORM, SAFMIN )
   BSCALE = 1.0E+0 / MAX( BNORM, SAFMIN )
!
!     Left eigenvectors
!
   IF( COMPL ) THEN
      IEIG = 0
!
!        Main loop over eigenvalues
!
      DO JE = 1, N
         IF( ILALL ) THEN
            ILCOMP = .TRUE.
         ELSE
            ILCOMP = SELECT( JE )
         END IF
         IF( ILCOMP ) THEN
            IEIG = IEIG + 1
!
            IF( CABS1( S( JE, JE ) ) <= SAFMIN .AND. &
                ABS( REAL( P( JE, JE ) ) ) <= SAFMIN ) THEN
!
!                 Singular matrix pencil -- return unit eigenvector
!
               VL(1:N, IEIG ) = (0.0E+0,0.0E+0)
               VL( IEIG, IEIG ) = (1.0E+0,0.0E+0)
               GO TO 140
            END IF
!
!              Non-singular eigenvalue:
!              Compute coefficients  a  and  b  in
!                   H
!                 y  ( a A - b B ) = 0
!
            TEMP = 1.0E+0 / MAX( CABS1( S( JE, JE ) )*ASCALE, &
                   ABS( REAL( P( JE, JE ) ) )*BSCALE, SAFMIN )
            SALPHA = ( TEMP*S( JE, JE ) )*ASCALE
            SBETA = ( TEMP*REAL( P( JE, JE ) ) )*BSCALE
            ACOEFF = SBETA*ASCALE
            BCOEFF = SALPHA*BSCALE
!
!              Scale to avoid underflow
!
            LSA = ABS( SBETA ) >= SAFMIN .AND. ABS( ACOEFF ) < SMALL
            LSB = CABS1( SALPHA ) >= SAFMIN .AND. CABS1( BCOEFF ) < &
                  SMALL
!
            SCALE = 1.0E+0
            IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
            IF( LSB ) &
               SCALE = MAX( SCALE, ( SMALL / CABS1( SALPHA ) )* MIN( BNORM, BIG ) )
            IF( LSA .OR. LSB ) THEN
               SCALE = MIN( SCALE, 1.0E+0 / &
                       ( SAFMIN*MAX( 1.0E+0, ABS( ACOEFF ), CABS1( BCOEFF ) ) ) )
               IF( LSA ) THEN
                  ACOEFF = ASCALE*( SCALE*SBETA )
               ELSE
                  ACOEFF = SCALE*ACOEFF
               END IF
               IF( LSB ) THEN
                  BCOEFF = BSCALE*( SCALE*SALPHA )
               ELSE
                  BCOEFF = SCALE*BCOEFF
               END IF
            END IF
!
            ACOEFA = ABS( ACOEFF )
            BCOEFA = CABS1( BCOEFF )
            XMAX = 1.0E+0
            WORK(1:N) = (0.0E+0,0.0E+0)
            WORK( JE ) = (1.0E+0,0.0E+0)
            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
!
!                                              H
!              Triangular solve of  (a A - b B)  y = 0
!
!                                      H
!              (rowwise in  (a A - b B) , or columnwise in a A - b B)
!
            DO J = JE + 1, N
!
!                 Compute
!                       j-1
!                 SOMME = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
!                       k=je
!                 (Scale if necessary)
!
               TEMP = 1.0E+0 / XMAX
               IF( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ) > BIGNUM* TEMP ) THEN
                  WORK(JE:J-1) = TEMP*WORK(JE:J-1)
                  XMAX = 1.0E+0
               END IF
               SUMA = SUM(CONJG( S( JE:J-1, J ) )*WORK( JE:J-1 ))
               SUMB = SUM(CONJG( P( JE:J-1, J ) )*WORK( JE:J-1 ))
               SOMME = ACOEFF*SUMA - CONJG( BCOEFF )*SUMB
!
!                 Form x(j) = - SOMME / conjg( a*S(j,j) - b*P(j,j) )
!
!                 with scaling and perturbation of the denominator
!
               D = CONJG( ACOEFF*S( J, J )-BCOEFF*P( J, J ) )
               IF( CABS1( D ) <= DMIN ) D = CMPLX( DMIN )
!
               IF( CABS1( D ) < 1.0E+0 ) THEN
                  IF( CABS1( SOMME ) >= BIGNUM*CABS1( D ) ) THEN
                     TEMP = 1.0E+0 / CABS1( SOMME )
                     WORK(JE:J-1) = TEMP*WORK(JE:J-1)
                     XMAX = TEMP*XMAX
                     SOMME = TEMP*SOMME
                  END IF
               END IF
               WORK( J ) = CLADIV( -SOMME, D )
               XMAX = MAX( XMAX, CABS1( WORK( J ) ) )
            ENDDO
!
!              Back transform eigenvector if HOWMNY='B'.
!
            IF( ILBACK ) THEN
               CALL CGEMV( 'N', N, N+1-JE, (1.0E+0,0.0E+0), VL( 1, JE ), LDVL, &
                           WORK( JE ), 1, (0.0E+0,0.0E+0), WORK( N+1 ), 1 )
               ISRC = 2
               IBEG = 1
            ELSE
               ISRC = 1
               IBEG = JE
            END IF
!
!              Copy and scale eigenvector into column of VL
!
            XMAX = 0.0E+0
            DO JR = IBEG, N
               XMAX = MAX( XMAX, CABS1( WORK( ( ISRC-1 )*N+JR ) ) )
            ENDDO
!
            IF( XMAX > SAFMIN ) THEN
               VL(IBEG:N, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+IBEG:ISRC*N) / XMAX
            ELSE
               IBEG = N + 1
            END IF
!
            VL(1:IBEG-1, IEIG ) = (0.0E+0,0.0E+0)
!
         END IF
  140    CONTINUE
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
      DO JE = N, 1, -1
         IF( ILALL ) THEN
            ILCOMP = .TRUE.
         ELSE
            ILCOMP = SELECT( JE )
         END IF
         IF( ILCOMP ) THEN
            IEIG = IEIG - 1
!
            IF( CABS1( S( JE, JE ) ) <= SAFMIN .AND. &
                ABS( REAL( P( JE, JE ) ) ) <= SAFMIN ) THEN
!
!                 Singular matrix pencil -- return unit eigenvector
!
               VR(1:N, IEIG ) = (0.0E+0,0.0E+0)
               VR( IEIG, IEIG ) = (1.0E+0,0.0E+0)
               GO TO 250
            END IF
!
!              Non-singular eigenvalue:
!              Compute coefficients  a  and  b  in
!
!              ( a A - b B ) x  = 0
!
            TEMP = 1.0E+0 / MAX( CABS1( S( JE, JE ) )*ASCALE, &
                   ABS( REAL( P( JE, JE ) ) )*BSCALE, SAFMIN )
            SALPHA = ( TEMP*S( JE, JE ) )*ASCALE
            SBETA = ( TEMP*REAL( P( JE, JE ) ) )*BSCALE
            ACOEFF = SBETA*ASCALE
            BCOEFF = SALPHA*BSCALE
!
!              Scale to avoid underflow
!
            LSA = ABS( SBETA ) >= SAFMIN .AND. ABS( ACOEFF ) < SMALL
            LSB = CABS1( SALPHA ) >= SAFMIN .AND. CABS1( BCOEFF ) < SMALL
!
            SCALE = 1.0E+0
            IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
            IF( LSB ) SCALE = MAX( SCALE, ( SMALL / CABS1( SALPHA ) )* MIN( BNORM, BIG ) )
            IF( LSA .OR. LSB ) THEN
               SCALE = MIN( SCALE, 1.0E+0 / ( SAFMIN*MAX( 1.0E+0, ABS( ACOEFF ), &
                       CABS1( BCOEFF ) ) ) )
               IF( LSA ) THEN
                  ACOEFF = ASCALE*( SCALE*SBETA )
               ELSE
                  ACOEFF = SCALE*ACOEFF
               END IF
               IF( LSB ) THEN
                  BCOEFF = BSCALE*( SCALE*SALPHA )
               ELSE
                  BCOEFF = SCALE*BCOEFF
               END IF
            END IF
!
            ACOEFA = ABS( ACOEFF )
            BCOEFA = CABS1( BCOEFF )
            XMAX = 1.0E+0
            WORK(1:N) = (0.0E+0,0.0E+0)
            WORK( JE ) = (1.0E+0,0.0E+0)
            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
!
!              Triangular solve of  (a A - b B) x = 0  (columnwise)
!
!              WORK(1:j-1) contains sums w,
!              WORK(j+1:JE) contains x
!
            WORK(1:JE-1) = ACOEFF*S(1:JE-1, JE ) - BCOEFF*P(1:JE-1, JE )
            WORK( JE ) = (1.0E+0,0.0E+0)
!
            DO J = JE - 1, 1, -1
!
!                 Form x(j) := - w(j) / d
!                 with scaling and perturbation of the denominator
!
               D = ACOEFF*S( J, J ) - BCOEFF*P( J, J )
               IF( CABS1( D ) <= DMIN ) D = CMPLX( DMIN )
!
               IF( CABS1( D ) < 1.0E+0 ) THEN
                  IF( CABS1( WORK( J ) ) >= BIGNUM*CABS1( D ) ) THEN
                     TEMP = 1.0E+0 / CABS1( WORK( J ) )
                     WORK(1:JE) = TEMP*WORK(1:JE)
                  END IF
               END IF
!
               WORK( J ) = CLADIV( -WORK( J ), D )
!
               IF( J > 1 ) THEN
!
!                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
!
                  IF( CABS1( WORK( J ) ) > 1.0E+0 ) THEN
                     TEMP = 1.0E+0 / CABS1( WORK( J ) )
                     IF( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ) >= BIGNUM*TEMP ) THEN
                        WORK(1:JE) = TEMP*WORK(1:JE)
                     END IF
                  END IF
!
                  CA = ACOEFF*WORK( J )
                  CB = BCOEFF*WORK( J )
                  WORK(1:J-1) = WORK(1:J-1) + CA*S(1:J-1, J ) - CB*P(1:J-1, J )
               END IF
            ENDDO
!
!              Back transform eigenvector if HOWMNY='B'.
!
            IF( ILBACK ) THEN
               CALL CGEMV( 'N', N, JE, (1.0E+0,0.0E+0), VR, LDVR, WORK, 1, &
                           (0.0E+0,0.0E+0), WORK( N+1 ), 1 )
               ISRC = 2
               IEND = N
            ELSE
               ISRC = 1
               IEND = JE
            END IF
!
!              Copy and scale eigenvector into column of VR
!
            XMAX = 0.0E+0
            DO JR = 1, IEND
               XMAX = MAX( XMAX, CABS1( WORK( ( ISRC-1 )*N+JR ) ) )
            ENDDO
!
            IF( XMAX > SAFMIN ) THEN
               TEMP = 1.0E+0 / XMAX
               VR(1:IEND,IEIG ) = TEMP*WORK((ISRC-1)*N+1:(ISRC-1)*N+IEND )
            ELSE
               IEND = 0
            END IF
!
            VR(IEND+1:N, IEIG ) = (0.0E+0,0.0E+0)
!
         END IF
  250    CONTINUE
      ENDDO
   END IF
!
   RETURN
!
!     End of CTGEVC
!
END
