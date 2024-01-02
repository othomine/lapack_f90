!> \brief \b CTRSYL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTRSYL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrsyl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrsyl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrsyl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
!                          LDC, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANA, TRANB
!       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRSYL solves the complex Sylvester matrix equation:
!>
!>    op(A)*X + X*op(B) = scale*C or
!>    op(A)*X - X*op(B) = scale*C,
!>
!> where op(A) = A or A**H, and A and B are both upper triangular. A is
!> M-by-M and B is N-by-N; the right hand side C and the solution X are
!> M-by-N; and scale is an output scale factor, set <= 1 to avoid
!> overflow in X.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANA
!> \verbatim
!>          TRANA is CHARACTER*1
!>          Specifies the option op(A):
!>          = 'N': op(A) = A    (No transpose)
!>          = 'C': op(A) = A**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] TRANB
!> \verbatim
!>          TRANB is CHARACTER*1
!>          Specifies the option op(B):
!>          = 'N': op(B) = B    (No transpose)
!>          = 'C': op(B) = B**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] ISGN
!> \verbatim
!>          ISGN is INTEGER
!>          Specifies the sign in the equation:
!>          = +1: solve op(A)*X + X*op(B) = scale*C
!>          = -1: solve op(A)*X - X*op(B) = scale*C
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrix A, and the number of rows in the
!>          matrices X and C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B, and the number of columns in the
!>          matrices X and C. N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,M)
!>          The upper triangular matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          The upper triangular matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          On entry, the M-by-N right hand side matrix C.
!>          On exit, C is overwritten by the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M)
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          The scale factor, scale, set <= 1 to avoid overflow in X.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          = 1: A and B have common or very close eigenvalues; perturbed
!>               values were used to solve the equation (but the matrices
!>               A and B are unchanged).
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
!> \ingroup trsyl
!
!  =====================================================================
   SUBROUTINE CTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANA, TRANB
   INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
   REAL               SCALE
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            NOTRNA, NOTRNB
   INTEGER            J, K, L
   REAL               BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, &
                      SMLNUM
   COMPLEX            A11, SUML, SUMR, VEC, X11
!     ..
!     .. Local Arrays ..
   REAL               DUM( 1 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   REAL               CLANGE, SLAMCH
   COMPLEX            CLADIV
   EXTERNAL           LSAME, CLANGE, SLAMCH, CLADIV
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Decode and Test input parameters
!
   NOTRNA = LSAME( TRANA, 'N' )
   NOTRNB = LSAME( TRANB, 'N' )
!
   INFO = 0
   IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'C' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'C' ) ) THEN
      INFO = -2
   ELSE IF( ISGN /= 1 .AND. ISGN /= -1 ) THEN
      INFO = -3
   ELSE IF( M < 0 ) THEN
      INFO = -4
   ELSE IF( N < 0 ) THEN
      INFO = -5
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -7
   ELSE IF( LDB < MAX( 1, N ) ) THEN
      INFO = -9
   ELSE IF( LDC < MAX( 1, M ) ) THEN
      INFO = -11
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CTRSYL', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   SCALE = 1.0E+0
   IF( M == 0 .OR. N == 0 ) RETURN
!
!     Set constants to control overflow
!
   EPS = SLAMCH( 'P' )
   SMLNUM = SLAMCH( 'S' )
   BIGNUM = 1.0E+0 / SMLNUM
   SMLNUM = SMLNUM*REAL( M*N ) / EPS
   BIGNUM = 1.0E+0 / SMLNUM
   SMIN = MAX( SMLNUM, EPS*CLANGE( 'M', M, M, A, LDA, DUM ), &
          EPS*CLANGE( 'M', N, N, B, LDB, DUM ) )
   SGN = ISGN
!
   IF( NOTRNA .AND. NOTRNB ) THEN
!
!        Solve    A*X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    M                        L-1
!          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
!                  I=K+1                      J=1
!
      DO L = 1, N
         DO K = M, 1, -1
!
            IF (K == M) THEN
               SUML = (0.0E+0,0.0E+0)
            ELSE
               SUML = SUM(A(K,K+1:M)*C(K+1:M,L))
            ENDIF
            SUMR = SUM(C(K,1:L-1)*B(1:L-1,L))
            VEC = C( K, L ) - ( SUML+SGN*SUMR )
!
            SCALOC = 1.0E+0
            A11 = A( K, K ) + SGN*B( L, L )
            DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
            IF( DA11 <= SMIN ) THEN
               A11 = SMIN
               DA11 = SMIN
               INFO = 1
            END IF
            DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
            IF( DA11 < 1.0E+0 .AND. DB > 1.0E+0 ) THEN
               IF( DB > BIGNUM*DA11 ) SCALOC = 1.0E+0 / DB
            END IF
            X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
            IF( SCALOC /= 1.0E+0 ) THEN
               C(1:M,1:N) = SCALOC*C(1:M,1:N)
               SCALE = SCALE*SCALOC
            END IF
            C( K, L ) = X11
!
         ENDDO
      ENDDO
!
   ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
!
!        Solve    A**H *X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        upper-left corner column by column by
!
!            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                   K-1                           L-1
!          R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
!                   I=1                           J=1
!
      DO L = 1, N
         DO K = 1, M
!
            SUML = SUM(CONJG(A(1:K-1,K))*C(1:K-1,L))
            SUMR = SUM(C(K,1:L-1)*B(1:L-1,L))
            VEC = C( K, L ) - ( SUML+SGN*SUMR )
!
            SCALOC = 1.0E+0
            A11 = CONJG( A( K, K ) ) + SGN*B( L, L )
            DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
            IF( DA11 <= SMIN ) THEN
               A11 = SMIN
               DA11 = SMIN
               INFO = 1
            END IF
            DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
            IF( DA11 < 1.0E+0 .AND. DB > 1.0E+0 ) THEN
               IF( DB > BIGNUM*DA11 ) SCALOC = 1.0E+0 / DB
            END IF
!
            X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
            IF( SCALOC /= 1.0E+0 ) THEN
               C(1:M,1:N) = SCALOC*C(1:M,1:N)
               SCALE = SCALE*SCALOC
            END IF
            C( K, L ) = X11
!
         ENDDO
      ENDDO
!
   ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
!
!        Solve    A**H*X + ISGN*X*B**H = C.
!
!        The (K,L)th block of X is determined starting from
!        upper-right corner column by column by
!
!            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    K-1
!           R(K,L) = SUM [A**H(I,K)*X(I,L)] +
!                    I=1
!                           N
!                     ISGN*SUM [X(K,J)*B**H(L,J)].
!                          J=L+1
!
      DO L = N, 1, -1
         DO K = 1, M
!
            SUML = SUM(CONJG(A(1:K-1,K))*C(1:K-1,L))
            IF (L == N) THEN
               SUMR = (0.0E+0,0.0E+0)
            ELSE
               SUMR = SUM((C(K,L+1:N))*B(L,L+1:N))
            ENDIF
            VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
!
            SCALOC = 1.0E+0
            A11 = CONJG( A( K, K )+SGN*B( L, L ) )
            DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
            IF( DA11 <= SMIN ) THEN
               A11 = SMIN
               DA11 = SMIN
               INFO = 1
            END IF
            DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
            IF( DA11 < 1.0E+0 .AND. DB > 1.0E+0 ) THEN
               IF( DB > BIGNUM*DA11 ) SCALOC = 1.0E+0 / DB
            END IF
!
            X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
            IF( SCALOC /= 1.0E+0 ) THEN
               C(1:M,1:N) = SCALOC*C(1:M,1:N)
               SCALE = SCALE*SCALOC
            END IF
            C( K, L ) = X11
!
         ENDDO
      ENDDO
!
   ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
!
!        Solve    A*X + ISGN*X*B**H = C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!           A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    M                          N
!          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
!                  I=K+1                      J=L+1
!
      DO L = N, 1, -1
         DO K = M, 1, -1
!
            IF (K == M) THEN
               SUML = (0.0E+0,0.0E+0)
            ELSE
               SUML = SUM(A(K,K+1:M)*C(K+1:M,L))
            ENDIF
            IF (L == N) THEN
               SUMR = (0.0E+0,0.0E+0)
            ELSE
               SUMR = SUM(CONJG(C(K,L+1:N))*B(L,L+1:N))
            ENDIF
            VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
!
            SCALOC = 1.0E+0
            A11 = A( K, K ) + SGN*CONJG( B( L, L ) )
            DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
            IF( DA11 <= SMIN ) THEN
               A11 = SMIN
               DA11 = SMIN
               INFO = 1
            END IF
            DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
            IF( DA11 < 1.0E+0 .AND. DB > 1.0E+0 ) THEN
               IF( DB > BIGNUM*DA11 ) SCALOC = 1.0E+0 / DB
            END IF
!
            X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
            IF( SCALOC /= 1.0E+0 ) THEN
               C(1:M,1:N) = SCALOC*C(1:M,1:N)
               SCALE = SCALE*SCALOC
            END IF
            C( K, L ) = X11
!
            ENDDO
         ENDDO
!
   END IF
!
   RETURN
!
!     End of CTRSYL
!
END
