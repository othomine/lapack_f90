!> \brief \b CLATPS solves a triangular system of equations with the matrix held in packed storage.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLATPS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatps.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatps.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatps.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE,
!                          CNORM, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORMIN, TRANS, UPLO
!       INTEGER            INFO, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       REAL               CNORM( * )
!       COMPLEX            AP( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLATPS solves one of the triangular systems
!>
!>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
!>
!> with scaling to prevent overflow, where A is an upper or lower
!> triangular matrix stored in packed form.  Here A**T denotes the
!> transpose of A, A**H denotes the conjugate transpose of A, x and b
!> are n-element vectors, and s is a scaling factor, usually less than
!> or equal to 1, chosen so that the components of x will be less than
!> the overflow threshold.  If the unscaled problem will not cause
!> overflow, the Level 2 BLAS routine CTPSV is called. If the matrix A
!> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
!> non-trivial solution to A*x = 0 is returned.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to A.
!>          = 'N':  Solve A * x = s*b     (No transpose)
!>          = 'T':  Solve A**T * x = s*b  (Transpose)
!>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] NORMIN
!> \verbatim
!>          NORMIN is CHARACTER*1
!>          Specifies whether CNORM has been set or not.
!>          = 'Y':  CNORM contains the column norms on entry
!>          = 'N':  CNORM is not set on entry.  On exit, the norms will
!>                  be computed and stored in CNORM.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          The upper or lower triangular matrix A, packed columnwise in
!>          a linear array.  The j-th column of A is stored in the array
!>          AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (N)
!>          On entry, the right hand side b of the triangular system.
!>          On exit, X is overwritten by the solution vector x.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          The scaling factor s for the triangular system
!>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
!>          If SCALE = 0, the matrix A is singular or badly scaled, and
!>          the vector x is an exact or approximate solution to A*x = 0.
!> \endverbatim
!>
!> \param[in,out] CNORM
!> \verbatim
!>          CNORM is REAL array, dimension (N)
!>
!>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
!>          contains the norm of the off-diagonal part of the j-th column
!>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
!>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
!>          must be greater than or equal to the 1-norm.
!>
!>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
!>          returns the 1-norm of the offdiagonal part of the j-th column
!>          of A.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -k, the k-th argument had an illegal value
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
!> \ingroup latps
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  A rough bound on x is computed; if that is less than overflow, CTPSV
!>  is called, otherwise, specific code is used which checks for possible
!>  overflow or divide-by-zero at every operation.
!>
!>  A columnwise scheme is used for solving A*x = b.  The basic algorithm
!>  if A is lower triangular is
!>
!>       x[1:n] := b[1:n]
!>       for j = 1, ..., n
!>            x(j) := x(j) / A(j,j)
!>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
!>       end
!>
!>  Define bounds on the components of x after j iterations of the loop:
!>     M(j) = bound on x[1:j]
!>     G(j) = bound on x[j+1:n]
!>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
!>
!>  Then for iteration j+1 we have
!>     M(j+1) <= G(j) / | A(j+1,j+1) |
!>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
!>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
!>
!>  where CNORM(j+1) is greater than or equal to the infinity-norm of
!>  column j+1 of A, not counting the diagonal.  Hence
!>
!>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
!>                  1<=i<=j
!>  and
!>
!>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
!>                                   1<=i< j
!>
!>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTPSV if the
!>  reciprocal of the largest M(j), j=1,..,n, is larger than
!>  max(underflow, 1/overflow).
!>
!>  The bound on x(j) is also used to determine when a step in the
!>  columnwise method can be performed without fear of overflow.  If
!>  the computed bound is greater than a large constant, x is scaled to
!>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
!>
!>  Similarly, a row-wise scheme is used to solve A**T *x = b  or
!>  A**H *x = b.  The basic algorithm for A upper triangular is
!>
!>       for j = 1, ..., n
!>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
!>       end
!>
!>  We simultaneously compute two bounds
!>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
!>       M(j) = bound on x(i), 1<=i<=j
!>
!>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
!>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
!>  Then the bound on x(j) is
!>
!>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
!>
!>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
!>                      1<=i<=j
!>
!>  and we can safely call CTPSV if 1/M(n) and 1/G(n) are both greater
!>  than max(underflow, 1/overflow).
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, &
                      CNORM, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, NORMIN, TRANS, UPLO
   INTEGER            INFO, N
   REAL               SCALE
!     ..
!     .. Array Arguments ..
   REAL               CNORM( * )
   COMPLEX            AP( * ), X( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            NOTRAN, NOUNIT, UPPER
   INTEGER            I, IMAX, IP, J, JFIRST, JINC, JLAST, JLEN
   REAL               BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL, &
                      XBND, XJ, XMAX
   COMPLEX            CSUMJ, TJJS, USCAL, ZDUM
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ICAMAX, ISAMAX
   REAL               SCASUM, SLAMCH
   COMPLEX            CDOTC, CDOTU, CLADIV
   EXTERNAL           LSAME, ICAMAX, ISAMAX, SCASUM, SLAMCH, CDOTC, &
                      CDOTU, CLADIV
!     ..
!     .. External Subroutines ..
   EXTERNAL           CTPSV, XERBLA
!     ..
!     .. Statement Functions ..
   REAL               CABS1, CABS2
!     ..
!     .. Statement Function definitions ..
   CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
   CABS2( ZDUM ) = ABS( REAL( ZDUM ) / 2. ) + &
                   ABS( AIMAG( ZDUM ) / 2. )
!     ..
!     .. Executable Statements ..
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   NOTRAN = LSAME( TRANS, 'N' )
   NOUNIT = LSAME( DIAG, 'N' )
!
!     Test the input parameters.
!
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
            LSAME( TRANS, 'C' ) ) THEN
      INFO = -2
   ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
      INFO = -3
   ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. &
            LSAME( NORMIN, 'N' ) ) THEN
      INFO = -4
   ELSE IF( N < 0 ) THEN
      INFO = -5
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CLATPS', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Determine machine dependent parameters to control overflow.
!
   SMLNUM = SLAMCH( 'Safe minimum' )
   BIGNUM = 1.0E+0 / SMLNUM
   SMLNUM = SMLNUM / SLAMCH( 'Precision' )
   BIGNUM = 1.0E+0 / SMLNUM
   SCALE = 1.0E+0
!
   IF( LSAME( NORMIN, 'N' ) ) THEN
!
!        Compute the 1-norm of each column, not including the diagonal.
!
      IF( UPPER ) THEN
!
!           A is upper triangular.
!
         IP = 1
         DO J = 1, N
            CNORM( J ) = SCASUM( J-1, AP( IP ), 1 )
            IP = IP + J
         ENDDO
      ELSE
!
!           A is lower triangular.
!
         IP = 1
         DO J = 1, N - 1
            CNORM( J ) = SCASUM( N-J, AP( IP+1 ), 1 )
            IP = IP + N - J + 1
         ENDDO
         CNORM( N ) = 0.0E+0
      END IF
   END IF
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM/2.
!
   IMAX = ISAMAX( N, CNORM, 1 )
   TMAX = CNORM( IMAX )
   IF( TMAX <= BIGNUM*0.5E+0 ) THEN
      TSCAL = 1.0E+0
   ELSE
      TSCAL = 0.5E+0 / ( SMLNUM*TMAX )
      CNORM(1:N) = TSCAL*CNORM(1:N)
   END IF
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine CTPSV can be used.
!
   XMAX = 0.0E+0
   DO J = 1, N
      XMAX = MAX( XMAX, CABS2( X( J ) ) )
   ENDDO
   XBND = XMAX
   IF( NOTRAN ) THEN
!
!        Compute the growth in A * x = b.
!
      IF( UPPER ) THEN
         JFIRST = N
         JLAST = 1
         JINC = -1
      ELSE
         JFIRST = 1
         JLAST = N
         JINC = 1
      END IF
!
      IF( TSCAL /= 1.0E+0 ) THEN
         GROW = 0.0E+0
         GO TO 60
      END IF
!
      IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, G(0) = max{x(i), i=1,...,n}.
!
         GROW = 0.5E+0 / MAX( XBND, SMLNUM )
         XBND = GROW
         IP = JFIRST*( JFIRST+1 ) / 2
         JLEN = N
         DO J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
            IF( GROW <= SMLNUM ) GO TO 60
!
            TJJS = AP( IP )
            TJJ = CABS1( TJJS )
!
            IF( TJJ >= SMLNUM ) THEN
!
!                 M(j) = G(j-1) / abs(A(j,j))
!
               XBND = MIN( XBND, MIN( 1.0E+0, TJJ )*GROW )
            ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
               XBND = 0.0E+0
            END IF
!
            IF( TJJ+CNORM( J ) >= SMLNUM ) THEN
!
!                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
!
               GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
            ELSE
!
!                 G(j) could overflow, set GROW to 0.
!
               GROW = 0.0E+0
            END IF
            IP = IP + JINC*JLEN
            JLEN = JLEN - 1
         ENDDO
         GROW = XBND
      ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
         GROW = MIN( 1.0E+0, 0.5E+0 / MAX( XBND, SMLNUM ) )
         DO J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
            IF( GROW <= SMLNUM ) GO TO 60
!
!              G(j) = G(j-1)*( 1 + CNORM(j) )
!
            GROW = GROW*( 1.0E+0 / ( 1.0E+0+CNORM( J ) ) )
         ENDDO
      END IF
60    CONTINUE
!
   ELSE
!
!        Compute the growth in A**T * x = b  or  A**H * x = b.
!
      IF( UPPER ) THEN
         JFIRST = 1
         JLAST = N
         JINC = 1
      ELSE
         JFIRST = N
         JLAST = 1
         JINC = -1
      END IF
!
      IF( TSCAL /= 1.0E+0 ) THEN
         GROW = 0.0E+0
         GO TO 90
      END IF
!
      IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, M(0) = max{x(i), i=1,...,n}.
!
         GROW = 0.5E+0 / MAX( XBND, SMLNUM )
         XBND = GROW
         IP = JFIRST*( JFIRST+1 ) / 2
         JLEN = 1
         DO J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
            IF( GROW <= SMLNUM ) GO TO 90
!
!              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
!
            XJ = 1.0E+0 + CNORM( J )
            GROW = MIN( GROW, XBND / XJ )
!
            TJJS = AP( IP )
            TJJ = CABS1( TJJS )
!
            IF( TJJ >= SMLNUM ) THEN
!
!                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
!
               IF( XJ > TJJ ) XBND = XBND*( TJJ / XJ )
            ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
               XBND = 0.0E+0
            END IF
            JLEN = JLEN + 1
            IP = IP + JINC*JLEN
         ENDDO
         GROW = MIN( GROW, XBND )
      ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
         GROW = MIN( 1.0E+0, 0.5E+0 / MAX( XBND, SMLNUM ) )
         DO J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
            IF( GROW <= SMLNUM ) GO TO 90
!
!              G(j) = ( 1 + CNORM(j) )*G(j-1)
!
            XJ = 1.0E+0 + CNORM( J )
            GROW = GROW / XJ
         ENDDO
      END IF
90    CONTINUE
   END IF
!
   IF( ( GROW*TSCAL ) > SMLNUM ) THEN
!
!        Use the Level 2 BLAS solve if the reciprocal of the bound on
!        elements of X is not too small.
!
      CALL CTPSV( UPLO, TRANS, DIAG, N, AP, X, 1 )
   ELSE
!
!        Use a Level 1 BLAS solve, scaling intermediate results.
!
      IF( XMAX > BIGNUM*0.5E+0 ) THEN
!
!           Scale X so that its components are less than or equal to
!           BIGNUM in absolute value.
!
         SCALE = ( BIGNUM*0.5E+0 ) / XMAX
         X(1:N) = SCALE*X(1:N)
         XMAX = BIGNUM
      ELSE
         XMAX = XMAX*2.0E+0
      END IF
!
      IF( NOTRAN ) THEN
!
!           Solve A * x = b
!
         IP = JFIRST*( JFIRST+1 ) / 2
         DO J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
!
            XJ = CABS1( X( J ) )
            IF( NOUNIT ) THEN
               TJJS = AP( IP )*TSCAL
            ELSE
               TJJS = TSCAL
               IF( TSCAL == 1.0E+0 ) GO TO 105
            END IF
               TJJ = CABS1( TJJS )
               IF( TJJ > SMLNUM ) THEN
!
!                    abs(A(j,j)) > SMLNUM:
!
                  IF( TJJ < 1.0E+0 ) THEN
                     IF( XJ > TJJ*BIGNUM ) THEN
!
!                          Scale x by 1/b(j).
!
                        REC = 1.0E+0 / XJ
                        X(1:N) = REC*X(1:N)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = CLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE IF( TJJ > 0.0E+0 ) THEN
!
!                    0 < abs(A(j,j)) <= SMLNUM:
!
                  IF( XJ > TJJ*BIGNUM ) THEN
!
!                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
!                       to avoid overflow when dividing by A(j,j).
!
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ) > 1.0E+0 ) THEN
!
!                          Scale by 1/CNORM(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
                        REC = REC / CNORM( J )
                     END IF
                     X(1:N) = REC*X(1:N)
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = CLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE
!
!                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to A*x = 0.
!
                  X(1:N) = 0.0E+0
                  X( J ) = 1.0E+0
                  XJ = 1.0E+0
                  SCALE = 0.0E+0
                  XMAX = 0.0E+0
               END IF
  105          CONTINUE
!
!              Scale x if necessary to avoid overflow when adding a
!              multiple of column j of A.
!
            IF( XJ > 1.0E+0 ) THEN
               REC = 1.0E+0 / XJ
               IF( CNORM( J ) > ( BIGNUM-XMAX )*REC ) THEN
!
!                    Scale x by 1/(2*abs(x(j))).
!
                  REC = REC*0.5E+0
                  X(1:N) = REC*X(1:N)
                  SCALE = SCALE*REC
               END IF
            ELSE IF( XJ*CNORM( J ) > ( BIGNUM-XMAX ) ) THEN
!
!                 Scale x by 1/2.
!
                X(1:N) = 0.5E+0*X(1:N)
               SCALE = SCALE*0.5E+0
            END IF
!
            IF( UPPER ) THEN
               IF( J > 1 ) THEN
!
!                    Compute the update
!                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
!
                  X(1:J-1) = X(1:J-1)-X( J )*TSCAL*AP(IP-J+1:IP-1)
                  I = ICAMAX( J-1, X, 1 )
                  XMAX = CABS1( X( I ) )
               END IF
               IP = IP - J
            ELSE
               IF( J < N ) THEN
!
!                    Compute the update
!                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
!
                  X(J+1:N) = X(J+1:N)-X( J )*TSCAL*AP(IP+1:IP+N-J)
                  I = J + ICAMAX( N-J, X( J+1 ), 1 )
                  XMAX = CABS1( X( I ) )
               END IF
               IP = IP + N - J + 1
            END IF
            ENDDO
!
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
!
!           Solve A**T * x = b
!
         IP = JFIRST*( JFIRST+1 ) / 2
         JLEN = 1
         DO J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
            XJ = CABS1( X( J ) )
            USCAL = TSCAL
            REC = 1.0E+0 / MAX( XMAX, 1.0E+0 )
            IF( CNORM( J ) > ( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
               REC = REC*0.5E+0
               IF( NOUNIT ) THEN
                  TJJS = AP( IP )*TSCAL
               ELSE
                  TJJS = TSCAL
               END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ > 1.0E+0 ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( 1.0E+0, REC*TJJ )
                     USCAL = CLADIV( USCAL, TJJS )
                  END IF
               IF( REC < 1.0E+0 ) THEN
                  X(1:N) = REC*X(1:N)
                  SCALE = SCALE*REC
                  XMAX = XMAX*REC
               END IF
            END IF
!
            CSUMJ = 0.0E+0
            IF( USCAL == CMPLX( 1.0E+0 ) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call CDOTU to perform the dot product.
!
               IF( UPPER ) THEN
                  CSUMJ = CDOTU( J-1, AP( IP-J+1 ), 1, X, 1 )
               ELSE IF( J < N ) THEN
                  CSUMJ = CDOTU( N-J, AP( IP+1 ), 1, X( J+1 ), 1 )
               END IF
            ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
               IF( UPPER ) THEN
                  CSUMJ = CSUMJ + USCAL*SUM(AP(IP-J+1:IP-1)*X(1:J-1))
               ELSE IF( J < N ) THEN
                  CSUMJ = CSUMJ + USCAL*SUM(AP(IP+1:IP+N-J)*X(J+1:N))
               END IF
            END IF
!
            IF( USCAL == CMPLX( TSCAL ) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
               X( J ) = X( J ) - CSUMJ
               XJ = CABS1( X( J ) )
               IF( NOUNIT ) THEN
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                  TJJS = AP( IP )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL == 1.0E+0 ) GO TO 145
               END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ > SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ < 1.0E+0 ) THEN
                        IF( XJ > TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = 1.0E+0 / XJ
                           X(1:N) = REC*X(1:N)
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = CLADIV( X( J ), TJJS )
                  ELSE IF( TJJ > 0.0E+0 ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ > TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        X(1:N) = REC*X(1:N)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = CLADIV( X( J ), TJJS )
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**T *x = 0.
!
                     X(1:N) = 0.0E+0
                     X( J ) = 1.0E+0
                     SCALE = 0.0E+0
                     XMAX = 0.0E+0
                  END IF
  145             CONTINUE
            ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
               X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ
            END IF
            XMAX = MAX( XMAX, CABS1( X( J ) ) )
            JLEN = JLEN + 1
            IP = IP + JINC*JLEN
            ENDDO
!
      ELSE
!
!           Solve A**H * x = b
!
         IP = JFIRST*( JFIRST+1 ) / 2
         JLEN = 1
         DO J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
            XJ = CABS1( X( J ) )
            USCAL = TSCAL
            REC = 1.0E+0 / MAX( XMAX, 1.0E+0 )
            IF( CNORM( J ) > ( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
               REC = REC*0.5E+0
               IF( NOUNIT ) THEN
                  TJJS = CONJG( AP( IP ) )*TSCAL
               ELSE
                  TJJS = TSCAL
               END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ > 1.0E+0 ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( 1.0E+0, REC*TJJ )
                     USCAL = CLADIV( USCAL, TJJS )
                  END IF
               IF( REC < 1.0E+0 ) THEN
                  X(1:N) = REC*X(1:N)
                  SCALE = SCALE*REC
                  XMAX = XMAX*REC
               END IF
            END IF
!
            CSUMJ = 0.0E+0
            IF( USCAL == CMPLX( 1.0E+0 ) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call CDOTC to perform the dot product.
!
               IF( UPPER ) THEN
                  CSUMJ = CDOTC( J-1, AP( IP-J+1 ), 1, X, 1 )
               ELSE IF( J < N ) THEN
                  CSUMJ = CDOTC( N-J, AP( IP+1 ), 1, X( J+1 ), 1 )
               END IF
            ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
               IF( UPPER ) THEN
                  CSUMJ = CSUMJ + USCAL*SUM(CONJG(AP(IP-J+1:IP-1))*X(1:J-1))
               ELSE IF( J < N ) THEN
                  CSUMJ = CSUMJ + USCAL*SUM(CONJG(AP(IP+1:IP+N-J))*X(J+1:N))
               END IF
            END IF
!
            IF( USCAL == CMPLX( TSCAL ) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
               X( J ) = X( J ) - CSUMJ
               XJ = CABS1( X( J ) )
               IF( NOUNIT ) THEN
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                  TJJS = CONJG( AP( IP ) )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL == 1.0E+0 ) GO TO 185
               END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ > SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ < 1.0E+0 ) THEN
                        IF( XJ > TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = 1.0E+0 / XJ
                           X(1:N) = REC*X(1:N)
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = CLADIV( X( J ), TJJS )
                  ELSE IF( TJJ > 0.0E+0 ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ > TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        X(1:N) = REC*X(1:N)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = CLADIV( X( J ), TJJS )
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**H *x = 0.
!
                     X(1:N) = 0.0E+0
                     X( J ) = 1.0E+0
                     SCALE = 0.0E+0
                     XMAX = 0.0E+0
                  END IF
  185             CONTINUE
            ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
               X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ
            END IF
            XMAX = MAX( XMAX, CABS1( X( J ) ) )
            JLEN = JLEN + 1
            IP = IP + JINC*JLEN
         ENDDO
      END IF
      SCALE = SCALE / TSCAL
   END IF
!
!     Scale the column norms by 1/TSCAL for return.
!
   IF( TSCAL /= 1.0E+0 ) CNORM(1:N) = CNORM(1:N)/TSCAL
!
   RETURN
!
!     End of CLATPS
!
END
