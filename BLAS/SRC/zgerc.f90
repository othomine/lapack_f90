!> \brief \b ZGERC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGERC  performs the rank 1 operation
!>
!>    A := alpha*x*y**H + A,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and A is an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>          ALPHA is COMPLEX*16
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the m
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients. On exit, A is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
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
!> \ingroup ger
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
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
   SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16 ALPHA
   INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
   COMPLEX*16 A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   COMPLEX*16 TEMP
   INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
   EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC DCONJG,MAX
!     ..
!
!     Test the input parameters.
!
   INFO = 0
   IF (M < 0) THEN
       INFO = 1
   ELSE IF (N < 0) THEN
       INFO = 2
   ELSE IF (INCX == 0) THEN
       INFO = 5
   ELSE IF (INCY == 0) THEN
       INFO = 7
   ELSE IF (LDA < MAX(1,M)) THEN
       INFO = 9
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('ZGERC ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((M == 0) .OR. (N == 0) .OR. (ALPHA == (0.0D+0,0.0D+0))) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
   IF (INCY > 0) THEN
       JY = 1
   ELSE
       JY = 1 - (N-1)*INCY
   END IF
   IF (INCX == 1) THEN
       DO J = 1,N
           IF (Y(JY) /= (0.0D+0,0.0D+0)) THEN
               A(1:M,J) = A(1:M,J) + X(1:M)*ALPHA*DCONJG(Y(JY))
           END IF
           JY = JY + INCY
       ENDDO
   ELSE
       IF (INCX > 0) THEN
           KX = 1
       ELSE
           KX = 1 - (M-1)*INCX
       END IF
       DO J = 1,N
           IF (Y(JY) /= (0.0D+0,0.0D+0)) THEN
               TEMP = ALPHA*DCONJG(Y(JY))
               IX = KX
               DO I = 1,M
                   A(I,J) = A(I,J) + X(IX)*TEMP
                   IX = IX + INCX
               ENDDO
           END IF
           JY = JY + INCY
       ENDDO
   END IF
!
   RETURN
!
!     End of ZGERC
!
END
