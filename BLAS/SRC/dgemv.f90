!> \brief \b DGEMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
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
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension at least
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
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry with BETA non-zero, the incremented array Y
!>           must contain the vector y. On exit, Y is overwritten by the
!>           updated vector y.
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
!> \ingroup gemv
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
   SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION ALPHA,BETA
   INTEGER INCX,INCY,LDA,M,N
   CHARACTER TRANS
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION TEMP
   INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
   LOGICAL LSAME
   EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC MAX
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
   ELSE IF (LDA < MAX(1,M)) THEN
       INFO = 6
   ELSE IF (INCX == 0) THEN
       INFO = 8
   ELSE IF (INCY == 0) THEN
       INFO = 11
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('DGEMV ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((M == 0) .OR. (N == 0) .OR. &
       ((ALPHA == 0.0D+0).AND. (BETA == 1.0D+0))) RETURN
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
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
   IF (BETA == 0.0D+0) THEN
       IF (INCY == 1) THEN
           Y(1:LENY) = 0.0D+0
       ELSE
           IY = KY
           DO I = 1,LENY
               Y(IY) = BETA*Y(IY)
               IY = IY + INCY
           ENDDO
       END IF
   ELSEIF (BETA /= 1.0D+0) THEN
       IF (INCY == 1) THEN
           Y(1:LENY) = BETA*Y(1:LENY)
       ELSE
           IY = KY
           DO I = 1,LENY
               Y(IY) = BETA*Y(IY)
               IY = IY + INCY
           ENDDO
       END IF
   ENDIF
   IF (ALPHA == 0.0D+0) RETURN
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
       JX = KX
       IF (INCY == 1) THEN
           DO J = 1,N
               Y(1:M) = Y(1:M) + ALPHA*X(JX)*A(1:M,J)
               JX = JX + INCX
           ENDDO
       ELSE
           DO J = 1,N
               TEMP = ALPHA*X(JX)
               IY = KY
               DO I = 1,M
                   Y(IY) = Y(IY) + TEMP*A(I,J)
                   IY = IY + INCY
               ENDDO
               JX = JX + INCX
           ENDDO
       END IF
   ELSE
!
!        Form  y := alpha*A**T*x + y.
!
       JY = KY
       IF (INCX == 1) THEN
           DO J = 1,N
               Y(JY) = Y(JY) + ALPHA*sum(A(1:M,J)*X(1:M))
               JY = JY + INCY
           ENDDO
       ELSE
           DO J = 1,N
               TEMP = 0.0D+0
               IX = KX
               DO I = 1,M
                   TEMP = TEMP + A(I,J)*X(IX)
                   IX = IX + INCX
               ENDDO
               Y(JY) = Y(JY) + ALPHA*TEMP
               JY = JY + INCY
           ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of DGEMV
!
END