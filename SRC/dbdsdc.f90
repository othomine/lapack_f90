!> \brief \b DBDSDC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DBDSDC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsdc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsdc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsdc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ,
!                          WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, UPLO
!       INTEGER            INFO, LDU, LDVT, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IQ( * ), IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), Q( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DBDSDC computes the singular value decomposition (SVD) of a real
!> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
!> using a divide and conquer method, where S is a diagonal matrix
!> with non-negative diagonal elements (the singular values of B), and
!> U and VT are orthogonal matrices of left and right singular vectors,
!> respectively. DBDSDC can be used to compute all singular values,
!> and optionally, singular vectors or singular vectors in compact form.
!>
!> The code currently calls DLASDQ if singular values only are desired.
!> However, it can be slightly modified to compute singular values
!> using the divide and conquer method.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  B is upper bidiagonal.
!>          = 'L':  B is lower bidiagonal.
!> \endverbatim
!>
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          Specifies whether singular vectors are to be computed
!>          as follows:
!>          = 'N':  Compute singular values only;
!>          = 'P':  Compute singular values and compute singular
!>                  vectors in compact form;
!>          = 'I':  Compute singular values and singular vectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the n diagonal elements of the bidiagonal matrix B.
!>          On exit, if INFO=0, the singular values of B.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the elements of E contain the offdiagonal
!>          elements of the bidiagonal matrix whose SVD is desired.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU,N)
!>          If  COMPQ = 'I', then:
!>             On exit, if INFO = 0, U contains the left singular vectors
!>             of the bidiagonal matrix.
!>          For other values of COMPQ, U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= 1.
!>          If singular vectors are desired, then LDU >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is DOUBLE PRECISION array, dimension (LDVT,N)
!>          If  COMPQ = 'I', then:
!>             On exit, if INFO = 0, VT**T contains the right singular
!>             vectors of the bidiagonal matrix.
!>          For other values of COMPQ, VT is not referenced.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= 1.
!>          If singular vectors are desired, then LDVT >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ)
!>          If  COMPQ = 'P', then:
!>             On exit, if INFO = 0, Q and IQ contain the left
!>             and right singular vectors in a compact form,
!>             requiring O(N log N) space instead of 2*N**2.
!>             In particular, Q contains all the DOUBLE PRECISION data in
!>             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
!>             words of memory, where SMLSIZ is returned by ILAENV and
!>             is equal to the maximum size of the subproblems at the
!>             bottom of the computation tree (usually about 25).
!>          For other values of COMPQ, Q is not referenced.
!> \endverbatim
!>
!> \param[out] IQ
!> \verbatim
!>          IQ is INTEGER array, dimension (LDIQ)
!>          If  COMPQ = 'P', then:
!>             On exit, if INFO = 0, Q and IQ contain the left
!>             and right singular vectors in a compact form,
!>             requiring O(N log N) space instead of 2*N**2.
!>             In particular, IQ contains all INTEGER data in
!>             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
!>             words of memory, where SMLSIZ is returned by ILAENV and
!>             is equal to the maximum size of the subproblems at the
!>             bottom of the computation tree (usually about 25).
!>          For other values of COMPQ, IQ is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          If COMPQ = 'N' then LWORK >= (4 * N).
!>          If COMPQ = 'P' then LWORK >= (6 * N).
!>          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (8*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute a singular value.
!>                The update process of divide and conquer failed.
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
!> \ingroup bdsdc
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
   SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, WORK, IWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          COMPQ, UPLO
   INTEGER            INFO, LDU, LDVT, N
!     ..
!     .. Array Arguments ..
   INTEGER            IQ( * ), IWORK( * )
   DOUBLE PRECISION   D( * ), E( * ), Q( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
!     ..
!
!  =====================================================================
!  Changed dimension statement in comment describing E from (N) to
!  (N-1).  Sven, 17 Feb 05.
!  =====================================================================
!     ..
!     .. Local Array ..
   DOUBLE PRECISION   U_tmp( N ), VT_tmp( N )
!     .. Local Scalars ..
   INTEGER            DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, IC, &
                      ICOMPQ, IERR, II, IS, IU, IUPLO, IVT, J, K, KK, &
                      MLVL, NM1, NSIZE, PERM, POLES, QSTART, SMLSIZ, &
                      SMLSZP, SQRE, START, WSTART, Z
   DOUBLE PRECISION   CS, EPS, ORGNRM, P, R, SN
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   INTEGER            ILAENV
   DOUBLE PRECISION   DLAMCH, DLANST
   EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DLARTG, DLASCL, DLASD0, DLASDA, DLASDQ, &
                      DLASET, DLASR, DSWAP, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
!
   IUPLO = 0
   IF( LSAME( UPLO, 'U' ) ) IUPLO = 1
   IF( LSAME( UPLO, 'L' ) ) IUPLO = 2
   IF( LSAME( COMPQ, 'N' ) ) THEN
      ICOMPQ = 0
   ELSE IF( LSAME( COMPQ, 'P' ) ) THEN
      ICOMPQ = 1
   ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
      ICOMPQ = 2
   ELSE
      ICOMPQ = -1
   END IF
   IF( IUPLO == 0 ) THEN
      INFO = -1
   ELSE IF( ICOMPQ < 0 ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -3
   ELSE IF( ( LDU < 1 ) .OR. ( ( ICOMPQ == 2 ) .AND. ( LDU < N ) ) ) THEN
      INFO = -7
   ELSE IF( ( LDVT < 1 ) .OR. ( ( ICOMPQ == 2 ) .AND. ( LDVT < N ) ) ) THEN
      INFO = -9
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DBDSDC', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
   SMLSIZ = ILAENV( 9, 'DBDSDC', ' ', 0, 0, 0, 0 )
   IF( N == 1 ) THEN
      IF( ICOMPQ == 1 ) THEN
         Q( 1 ) = SIGN( 1.0D0, D( 1 ) )
         Q( 1+SMLSIZ*N ) = 1.0D0
      ELSE IF( ICOMPQ == 2 ) THEN
         U( 1, 1 ) = SIGN( 1.0D0, D( 1 ) )
         VT( 1, 1 ) = 1.0D0
      END IF
      D( 1 ) = ABS( D( 1 ) )
      RETURN
   END IF
   NM1 = N - 1
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
   WSTART = 1
   QSTART = 3
   IF( ICOMPQ == 1 ) THEN
      Q(1:N) = D(1:N)
      Q(N+1:2*N) = E(1:N-1)
   END IF
   IF( IUPLO == 2 ) THEN
      QSTART = 5
      IF( ICOMPQ  ==  2 ) WSTART = 2*N - 1
      DO I = 1, N - 1
         CALL DLARTG( D( I ), E( I ), CS, SN, R )
         D( I ) = R
         E( I ) = SN*D( I+1 )
         D( I+1 ) = CS*D( I+1 )
         IF( ICOMPQ == 1 ) THEN
            Q( I+2*N ) = CS
            Q( I+3*N ) = SN
         ELSE IF( ICOMPQ == 2 ) THEN
            WORK( I ) = CS
            WORK( NM1+I ) = -SN
         END IF
      ENDDO
   END IF
!
!     If ICOMPQ = 0, use DLASDQ to compute the singular values.
!
   IF( ICOMPQ == 0 ) THEN
!        Ignore WSTART, instead using WORK( 1 ), since the two vectors
!        for CS and -SN above are added only if ICOMPQ == 2,
!        and adding them exceeds documented WORK size of 4*n.
      CALL DLASDQ( 'U', 0, N, 0, 0, 0, D, E, VT, LDVT, U, LDU, U, &
                   LDU, WORK( 1 ), INFO )
      GO TO 40
   END IF
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
   IF( N <= SMLSIZ ) THEN
      IF( ICOMPQ == 2 ) THEN
         U(1:N,1:N) = 0.0D0
         VT(1:N,1:N) = 0.0D0
         DO I=1, N
            U(I,I) = 1.0D0
            VT(I,I) = 1.0D0
         ENDDO
         CALL DLASDQ( 'U', 0, N, N, N, 0, D, E, VT, LDVT, U, LDU, U, &
                      LDU, WORK( WSTART ), INFO )
      ELSE IF( ICOMPQ == 1 ) THEN
         IU = 1
         IVT = IU + N
         CALL DLASET( 'A', N, N, 0.0D0, 1.0D0, Q( IU+( QSTART-1 )*N ), N )
         CALL DLASET( 'A', N, N, 0.0D0, 1.0D0, Q( IVT+( QSTART-1 )*N ), N )
         CALL DLASDQ( 'U', 0, N, N, N, 0, D, E, &
                      Q( IVT+( QSTART-1 )*N ), N, &
                      Q( IU+( QSTART-1 )*N ), N, &
                      Q( IU+( QSTART-1 )*N ), N, WORK( WSTART ), &
                      INFO )
      END IF
      GO TO 40
   END IF
!
   IF( ICOMPQ == 2 ) THEN
      U(1:N,1:N) = 0.0D0
      VT(1:N,1:N) = 0.0D0
      DO I=1, N
         U(I,I) = 1.0D0
         VT(I,I) = 1.0D0
      ENDDO
   END IF
!
!     Scale.
!
   ORGNRM = DLANST( 'M', N, D, E )
   IF( ORGNRM == 0.0D0 ) RETURN
   CALL DLASCL( 'G', 0, 0, ORGNRM, 1.0D0, N, 1, D, N, IERR )
   CALL DLASCL( 'G', 0, 0, ORGNRM, 1.0D0, NM1, 1, E, NM1, IERR )
!
   EPS = (0.9D+0)*DLAMCH( 'Epsilon' )
!
   MLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( 2.0D0 ) ) + 1
   SMLSZP = SMLSIZ + 1
!
   IF( ICOMPQ == 1 ) THEN
      IU = 1
      IVT = 1 + SMLSIZ
      DIFL = IVT + SMLSZP
      DIFR = DIFL + MLVL
      Z = DIFR + MLVL*2
      IC = Z + MLVL
      IS = IC + 1
      POLES = IS + 1
      GIVNUM = POLES + 2*MLVL
!
      K = 1
      GIVPTR = 2
      PERM = 3
      GIVCOL = PERM + MLVL
   END IF
!
   DO I = 1, N
      IF( ABS( D( I ) ) < EPS ) THEN
         D( I ) = SIGN( EPS, D( I ) )
      END IF
   ENDDO
!
   START = 1
   SQRE = 0
!
   DO I = 1, NM1
      IF( ( ABS( E( I ) ) < EPS ) .OR. ( I == NM1 ) ) THEN
!
!           Subproblem found. First determine its size and then
!           apply divide and conquer on it.
!
         IF( I < NM1 ) THEN
!
!              A subproblem with E(I) small for I < NM1.
!
            NSIZE = I - START + 1
         ELSE IF( ABS( E( I ) ) >= EPS ) THEN
!
!              A subproblem with E(NM1) not too small but I = NM1.
!
            NSIZE = N - START + 1
         ELSE
!
!              A subproblem with E(NM1) small. This implies an
!              1-by-1 subproblem at D(N). Solve this 1-by-1 problem
!              first.
!
            NSIZE = I - START + 1
            IF( ICOMPQ == 2 ) THEN
               U( N, N ) = SIGN( 1.0D0, D( N ) )
               VT( N, N ) = 1.0D0
            ELSE IF( ICOMPQ == 1 ) THEN
               Q( N+( QSTART-1 )*N ) = SIGN( 1.0D0, D( N ) )
               Q( N+( SMLSIZ+QSTART-1 )*N ) = 1.0D0
            END IF
            D( N ) = ABS( D( N ) )
         END IF
         IF( ICOMPQ == 2 ) THEN
            CALL DLASD0( NSIZE, SQRE, D( START ), E( START ), &
                         U( START, START ), LDU, VT( START, START ), &
                         LDVT, SMLSIZ, IWORK, WORK( WSTART ), INFO )
         ELSE
            CALL DLASDA( ICOMPQ, SMLSIZ, NSIZE, SQRE, D( START ), &
                         E( START ), Q( START+( IU+QSTART-2 )*N ), N, &
                         Q( START+( IVT+QSTART-2 )*N ), &
                         IQ( START+K*N ), Q( START+( DIFL+QSTART-2 )* &
                         N ), Q( START+( DIFR+QSTART-2 )*N ), &
                         Q( START+( Z+QSTART-2 )*N ), &
                         Q( START+( POLES+QSTART-2 )*N ), &
                         IQ( START+GIVPTR*N ), IQ( START+GIVCOL*N ), &
                         N, IQ( START+PERM*N ), &
                         Q( START+( GIVNUM+QSTART-2 )*N ), &
                         Q( START+( IC+QSTART-2 )*N ), &
                         Q( START+( IS+QSTART-2 )*N ), &
                         WORK( WSTART ), IWORK, INFO )
         END IF
         IF( INFO /= 0 ) RETURN
         START = I + 1
      END IF
   ENDDO
!
!     Unscale
!
   CALL DLASCL( 'G', 0, 0, 1.0D0, ORGNRM, N, 1, D, N, IERR )
40 CONTINUE
!
!     Use Selection Sort to minimize swaps of singular vectors
!
   DO II = 2, N
      I = II - 1
      KK = I
      P = D( I )
      DO J = II, N
         IF( D( J ) > P ) THEN
            KK = J
            P = D( J )
         END IF
      ENDDO
      IF( KK /= I ) THEN
         D( KK ) = D( I )
         D( I ) = P
         IF( ICOMPQ == 1 ) THEN
            IQ( I ) = KK
         ELSE IF( ICOMPQ == 2 ) THEN
            U_tmp(1:N) = U(1:N,I)
            U(1:N,I) = U(1:N,KK)
            U(1:N,KK) = U_tmp(1:N)
            VT_tmp(1:N) = VT(I,1:N)
            VT(I,1:N) = VT(KK,1:N)
            VT(KK,1:N) = VT_tmp(1:N)
         END IF
      ELSE IF( ICOMPQ == 1 ) THEN
         IQ( I ) = I
      END IF
   ENDDO
!
!     If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO
!
   IF( ICOMPQ == 1 ) THEN
      IF( IUPLO == 1 ) THEN
         IQ( N ) = 1
      ELSE
         IQ( N ) = 0
      END IF
   END IF
!
!     If B is lower bidiagonal, update U by those Givens rotations
!     which rotated B to be upper bidiagonal
!
   IF( ( IUPLO == 2 ) .AND. ( ICOMPQ == 2 ) ) &
      CALL DLASR( 'L', 'V', 'B', N, N, WORK( 1 ), WORK( N ), U, LDU )
!
   RETURN
!
!     End of DBDSDC
!
END
