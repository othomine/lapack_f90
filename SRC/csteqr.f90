!> \brief \b CSTEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSTEQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csteqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csteqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csteqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER            INFO, LDZ, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), WORK( * )
!       COMPLEX            Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the implicit QL or QR method.
!> The eigenvectors of a full or band complex Hermitian matrix can also
!> be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this
!> matrix to tridiagonal form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'V':  Compute eigenvalues and eigenvectors of the original
!>                  Hermitian matrix.  On entry, Z must contain the
!>                  unitary matrix used to reduce the original matrix
!>                  to tridiagonal form.
!>          = 'I':  Compute eigenvalues and eigenvectors of the
!>                  tridiagonal matrix.  Z is initialized to the identity
!>                  matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
!>          On entry, if  COMPZ = 'V', then Z contains the unitary
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original Hermitian matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          eigenvectors are desired, then  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (max(1,2*N-2))
!>          If COMPZ = 'N', then WORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  the algorithm has failed to find all the eigenvalues in
!>                a total of 30*N iterations; if INFO = i, then i
!>                elements of E have not converged to zero; on exit, D
!>                and E contain the elements of a symmetric tridiagonal
!>                matrix which is unitarily similar to the original
!>                matrix.
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
!> \ingroup steqr
!
!  =====================================================================
   SUBROUTINE CSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
   IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          COMPZ
   INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
   REAL               D( * ), E( * ), WORK( * )
   COMPLEX            Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXIT
   PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Array ..
   COMPLEX            Z_tmp( LDZ)

!     ..
!     .. Local Scalars ..
   INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                      LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, &
                      NM1, NMAXIT
   REAL               ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                      S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               SLAMCH, SLANST, SLAPY2
   EXTERNAL           LSAME, SLAMCH, SLANST, SLAPY2
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLASET, CLASR, SLAE2, SLAEV2, SLARTG, &
                      SLASCL, SLASRT, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
!
   IF( LSAME( COMPZ, 'N' ) ) THEN
      ICOMPZ = 0
   ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
      ICOMPZ = 1
   ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
      ICOMPZ = 2
   ELSE
      ICOMPZ = -1
   END IF
   IF( ICOMPZ < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( ( LDZ < 1 ) .OR. ( ICOMPZ > 0 .AND. LDZ < MAX( 1, &
            N ) ) ) THEN
      INFO = -6
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSTEQR', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
   IF( N == 1 ) THEN
      IF( ICOMPZ == 2 ) Z( 1, 1 ) = (1.0E+0,0.0E+0 )
      RETURN
   END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
   EPS = SLAMCH( 'E' )
   EPS2 = EPS**2
   SAFMIN = SLAMCH( 'S' )
   SAFMAX = 1.0E+0 / SAFMIN
   SSFMAX = SQRT( SAFMAX ) / 3.0E+0
   SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
   IF( ICOMPZ == 2 ) &
      CALL CLASET( 'Full', N, N, (0.0E+0,0.0E+0 ), (1.0E+0,0.0E+0 ), Z, LDZ )
!
   NMAXIT = N*MAXIT
   JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
   L1 = 1
   NM1 = N - 1
!
10 CONTINUE
   IF( L1 > N ) GO TO 160
   IF( L1 > 1 ) E( L1-1 ) = 0.0E+0
   IF( L1 <= NM1 ) THEN
      DO M = L1, NM1
         TST = ABS( E( M ) )
         IF( TST == 0.0E+0 ) GO TO 30
         IF( TST <= ( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ 1 ) ) ) )*EPS ) THEN
            E( M ) = 0.0E+0
            GO TO 30
         END IF
      ENDDO
   END IF
   M = N
!
30 CONTINUE
   L = L1
   LSV = L
   LEND = M
   LENDSV = LEND
   L1 = M + 1
   IF( LEND == L ) GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
   ANORM = SLANST( 'I', LEND-L+1, D( L ), E( L ) )
   ISCALE = 0
   IF( ANORM == 0.0E+0 ) GO TO 10
   IF( ANORM > SSFMAX ) THEN
      ISCALE = 1
      CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO )
      CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO )
   ELSE IF( ANORM < SSFMIN ) THEN
      ISCALE = 2
      CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO )
      CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO )
   END IF
!
!     Choose between QL and QR iteration
!
   IF( ABS( D( LEND ) ) < ABS( D( L ) ) ) THEN
      LEND = LSV
      L = LENDSV
   END IF
!
   IF( LEND > L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
40    CONTINUE
      IF( L /= LEND ) THEN
         LENDM1 = LEND - 1
         DO M = L, LENDM1
            TST = ABS( E( M ) )**2
            IF( TST <= ( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ SAFMIN )GO TO 60
         ENDDO
      END IF
!
      M = LEND
!
60    CONTINUE
      IF( M < LEND ) E( M ) = 0.0E+0
      P = D( L )
      IF( M == L ) GO TO 80
!
!        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
!        to compute its eigensystem.
!
      IF( M == L+1 ) THEN
         IF( ICOMPZ > 0 ) THEN
            CALL SLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
            WORK( L ) = C
            WORK( N-1+L ) = S
            CALL CLASR( 'R', 'V', 'B', N, 2, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ )
         ELSE
            CALL SLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
         END IF
         D( L ) = RT1
         D( L+1 ) = RT2
         E( L ) = 0.0E+0
         L = L + 2
         IF( L <= LEND ) GO TO 40
         GO TO 140
      END IF
!
      IF( JTOT == NMAXIT ) GO TO 140
      JTOT = JTOT + 1
!
!        Form shift.
!
      G = ( D( L+1 )-P ) / ( 2.0E+0*E( L ) )
      R = SLAPY2( G, 1.0E+0 )
      G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
      S = 1.0E+0
      C = 1.0E+0
      P = 0.0E+0
!
!        Inner loop
!
      MM1 = M - 1
      DO I = MM1, L, -1
         F = S*E( I )
         B = C*E( I )
         CALL SLARTG( G, F, C, S, R )
         IF( I /= M-1 ) E( I+1 ) = R
         G = D( I+1 ) - P
         R = ( D( I )-G )*S + 2.0E+0*C*B
         P = S*R
         D( I+1 ) = G + P
         G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
         IF( ICOMPZ > 0 ) THEN
            WORK( I ) = C
            WORK( N-1+I ) = -S
         END IF
!
      ENDDO
!
!        If eigenvectors are desired, then apply saved rotations.
!
      IF( ICOMPZ > 0 ) THEN
         MM = M - L + 1
         CALL CLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ )
      END IF
!
      D( L ) = D( L ) - P
      E( L ) = G
      GO TO 40
!
!        Eigenvalue found.
!
80    CONTINUE
      D( L ) = P
!
      L = L + 1
      IF( L <= LEND ) GO TO 40
      GO TO 140
!
   ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
90    CONTINUE
      IF( L /= LEND ) THEN
         LENDP1 = LEND + 1
         DO M = L, LENDP1, -1
            TST = ABS( E( M-1 ) )**2
            IF( TST <= ( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ SAFMIN )GO TO 110
            ENDDO
      END IF
!
      M = LEND
!
  110    CONTINUE
      IF( M > LEND ) E( M-1 ) = 0.0E+0
      P = D( L )
      IF( M == L ) GO TO 130
!
!        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
!        to compute its eigensystem.
!
      IF( M == L-1 ) THEN
         IF( ICOMPZ > 0 ) THEN
            CALL SLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
            WORK( M ) = C
            WORK( N-1+M ) = S
            CALL CLASR( 'R', 'V', 'F', N, 2, WORK( M ), WORK( N-1+M ), Z( 1, L-1 ), LDZ )
         ELSE
            CALL SLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
         END IF
         D( L-1 ) = RT1
         D( L ) = RT2
         E( L-1 ) = 0.0E+0
         L = L - 2
         IF( L >= LEND ) GO TO 90
         GO TO 140
      END IF
!
      IF( JTOT == NMAXIT ) GO TO 140
      JTOT = JTOT + 1
!
!        Form shift.
!
      G = ( D( L-1 )-P ) / ( 2.0E+0*E( L-1 ) )
      R = SLAPY2( G, 1.0E+0 )
      G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
      S = 1.0E+0
      C = 1.0E+0
      P = 0.0E+0
!
!        Inner loop
!
      LM1 = L - 1
      DO I = M, LM1
         F = S*E( I )
         B = C*E( I )
         CALL SLARTG( G, F, C, S, R )
         IF( I /= M ) E( I-1 ) = R
         G = D( I ) - P
         R = ( D( I+1 )-G )*S + 2.0E+0*C*B
         P = S*R
         D( I ) = G + P
         G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
         IF( ICOMPZ > 0 ) THEN
            WORK( I ) = C
            WORK( N-1+I ) = S
         END IF
!
         ENDDO
!
!        If eigenvectors are desired, then apply saved rotations.
!
      IF( ICOMPZ > 0 ) THEN
         MM = L - M + 1
         CALL CLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), Z( 1, M ), LDZ )
      END IF
!
      D( L ) = D( L ) - P
      E( LM1 ) = G
      GO TO 90
!
!        Eigenvalue found.
!
  130    CONTINUE
      D( L ) = P
!
      L = L - 1
      IF( L >= LEND ) GO TO 90
      GO TO 140
!
   END IF
!
!     Undo scaling if necessary
!
  140 CONTINUE
   IF( ISCALE == 1 ) THEN
      CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )
      CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO )
   ELSE IF( ISCALE == 2 ) THEN
      CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )
      CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO )
   END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
   IF( JTOT == NMAXIT ) THEN
      INFO = INFO + COUNT(E(1:N-1)/= 0.0E+0)
      RETURN
   END IF
   GO TO 10
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
   IF( ICOMPZ == 0 ) THEN
!
!        Use Quick Sort
!
      CALL SLASRT( 'I', N, D, INFO )
!
   ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
      DO II = 2, N
         I = II - 1
         K = I
         P = D( I )
         DO J = II, N
            IF( D( J ) < P ) THEN
               K = J
               P = D( J )
            END IF
         ENDDO
         IF( K /= I ) THEN
            D( K ) = D( I )
            D( I ) = P
            Z_tmp(1:N) = Z(1:N,I)
            Z(1:N,I) = Z(1:N,K)
            Z(1:N,K) = Z_tmp(1:N)
         END IF
      ENDDO
   END IF
   RETURN
!
!     End of CSTEQR
!
END
