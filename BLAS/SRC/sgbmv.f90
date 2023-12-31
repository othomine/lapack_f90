!> \brief \b SGBMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
!       INTEGER INCX,INCY,KL,KU,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGBMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n band matrix, with kl sub-diagonals and ku super-diagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!>
!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!>
!>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           On entry, KL specifies the number of sub-diagonals of the
!>           matrix A. KL must satisfy  0 .le. KL.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           On entry, KU specifies the number of super-diagonals of the
!>           matrix A. KU must satisfy  0 .le. KU.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension ( LDA, N )
!>           Before entry, the leading ( kl + ku + 1 ) by n part of the
!>           array A must contain the matrix of coefficients, supplied
!>           column by column, with the leading diagonal of the matrix in
!>           row ( ku + 1 ) of the array, the first super-diagonal
!>           starting at position 2 in row ku, the first sub-diagonal
!>           starting at position 1 in row ( ku + 2 ), and so on.
!>           Elements in the array A that do not correspond to elements
!>           in the band matrix (such as the top left ku by ku triangle)
!>           are not referenced.
!>           The following program segment will transfer a band matrix
!>           from conventional full matrix storage to band storage:
!>
!>                 DO  J = 1, N
!>                    K = KU + 1 - J
!>                    DO  I = MAX( 1, J - KU ), MIN( M, J + KL )
!>                       A( K + I, J ) = matrix( I, J )
!>                10    ENDDO
!>                20 ENDDO
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           ( kl + ku + 1 ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!>           Before entry, the incremented array X must contain the
!>           vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry, the incremented array Y must contain the
!>           vector y. On exit, Y is overwritten by the updated vector y.
!>           If either m or n is zero, then Y not referenced and the function
!>           performs a quick return.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Olivier Thomine
!
!> \ingroup gbmv
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   REAL ALPHA,BETA
   INTEGER INCX,INCY,KL,KU,LDA,M,N,MINI,MAXI
   CHARACTER TRANS
!     ..
!     .. Array Arguments ..
   REAL A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   REAL TEMP
   INTEGER I,INFO,IX,IY,J,JX,JY,K,KUP1,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
   LOGICAL LSAME
   EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL XERBLA
!     ..
!
!     Test the input parameters.
!
   INFO = 0
   IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
       .NOT.LSAME(TRANS,'C')) THEN
       INFO = 1
   ELSE IF (M < 0) THEN
       INFO = 2
   ELSE IF (N < 0) THEN
       INFO = 3
   ELSE IF (KL < 0) THEN
       INFO = 4
   ELSE IF (KU < 0) THEN
       INFO = 5
   ELSE IF (LDA <  (KL+KU+1)) THEN
       INFO = 8
   ELSE IF (INCX == 0) THEN
       INFO = 10
   ELSE IF (INCY == 0) THEN
       INFO = 13
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('SGBMV ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((M == 0) .OR. (N == 0) .OR. &
       ((ALPHA == 0.0E+0).AND. (BETA == 1.0E+0))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
   IF (LSAME(TRANS,'N')) THEN
       LENX = N
       LENY = M
   ELSE
       LENX = M
       LENY = N
   END IF
   IF (INCX > 0) THEN
       KX = 1
   ELSE
       KX = 1 - (LENX-1)*INCX
   END IF
   IF (INCY > 0) THEN
       KY = 1
   ELSE
       KY = 1 - (LENY-1)*INCY
   END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the band part of A.
!
!     First form  y := beta*y.
!
   IF (BETA == 0.0E+0) THEN
       IF (INCY == 1) THEN
           Y(1:LENY) = 0.0E+0
       ELSE
           IY = KY
           DO I = 1,LENY
               Y(IY) = 0.0E+0
               IY = IY + INCY
           ENDDO
       END IF
   ELSEIF (BETA /= 1.0E+0) THEN
       IF (INCY == 1) THEN
           Y(1:LENY) = 0.0E+0
       ELSE
           IY = KY
           DO I = 1,LENY
               Y(IY) = BETA*Y(IY)
               IY = IY + INCY
           ENDDO
       END IF
   END IF
   IF (ALPHA == 0.0E+0) RETURN
   KUP1 = KU + 1
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
       JX = KX
       IF (INCY == 1) THEN
           DO J = 1,N
               K = KUP1 - J
               MINI = MAX(1,J-KU)
               MAXI = MIN(M,J+KL)
               Y(MINI:MAXI) = Y(MINI:MAXI) + ALPHA*X(JX)*A(K+MINI:K+MAXI,J)
               JX = JX + INCX
           ENDDO
       ELSE
           DO J = 1,N
               TEMP = ALPHA*X(JX)
               IY = KY
               K = KUP1 - J
               MINI = MAX(1,J-KU)
               MAXI = MIN(M,J+KL)
               DO I = MINI, MAXI
                   Y(IY) = Y(IY) + TEMP*A(K+I,J)
                   IY = IY + INCY
               ENDDO
               JX = JX + INCX
               IF (J > KU) KY = KY + INCY
           ENDDO
       END IF
   ELSE
!
!        Form  y := alpha*A**T*x + y.
!
       JY = KY
       IF (INCX == 1) THEN
           DO J = 1,N
               K = KUP1 - J
               MINI = MAX(1,J-KU)
               MAXI = MIN(M,J+KL)
               Y(JY) = Y(JY) + ALPHA*sum(A(K+MINI:K+MAXI,J)*X(MINI:MAXI))
               JY = JY + INCY
           ENDDO
       ELSE
           DO J = 1,N
               TEMP = 0.0E+0
               IX = KX
               K = KUP1 - J
               MINI = MAX(1,J-KU)
               MAXI = MIN(M,J+KL)
               DO I = MINI, MAXI
                   TEMP = TEMP + A(K+I,J)*X(IX)
                   IX = IX + INCX
               ENDDO
               Y(JY) = Y(JY) + ALPHA*TEMP
               JY = JY + INCY
               IF (J > KU) KX = KX + INCX
          ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of SGBMV
!
END
