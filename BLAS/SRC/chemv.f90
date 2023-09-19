!> \brief \b CHEMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA,BETA
!       INTEGER INCX,INCY,LDA,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHEMV  performs the matrix-vector  operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n hermitian matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the array A is to be referenced as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of A
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of A
!>                                  is to be referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular part of the hermitian matrix and the strictly
!>           lower triangular part of A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular part of the hermitian matrix and the strictly
!>           upper triangular part of A is not referenced.
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set and are assumed to be zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
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
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y. On exit, Y is overwritten by the updated
!>           vector y.
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
!> \ingroup hemv
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
!>     converted to F90 and optimized 2023, olivier thomine
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX ALPHA,BETA
   INTEGER INCX,INCY,LDA,N
   CHARACTER UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   COMPLEX TEMP1,TEMP2
   INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
!     ..
!     .. External Functions ..
   LOGICAL LSAME
   EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC CONJG,MAX,REAL
!     ..
!
!     Test the input parameters.
!
   INFO = 0
   IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
       INFO = 1
   ELSE IF (N < 0) THEN
       INFO = 2
   ELSE IF (LDA < MAX(1,N)) THEN
       INFO = 5
   ELSE IF (INCX == 0) THEN
       INFO = 7
   ELSE IF (INCY == 0) THEN
       INFO = 10
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('CHEMV ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N == 0) .OR. ((ALPHA == (0.0E+0,0.0E+0)).AND. (BETA == (1.0E+0,0.0E+0)))) RETURN
!
!     Set up the start points in  X  and  Y.
!
   IF (INCX > 0) THEN
       KX = 1
   ELSE
       KX = 1 - (N-1)*INCX
   END IF
   IF (INCY > 0) THEN
       KY = 1
   ELSE
       KY = 1 - (N-1)*INCY
   END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
!     First form  y := beta*y.
!
   IF (BETA == (0.0E+0,0.0E+0)) THEN
       IF (INCY == 1) THEN
           Y(1:N) = (0.0E+0,0.0E+0)
       ELSE
           IY = KY
           DO I = 1,N
               Y(IY) = (0.0E+0,0.0E+0)
               IY = IY + INCY
           ENDDO
       END IF
   ELSEIF (BETA /= (1.0E+0,0.0E+0)) THEN
       IF (INCY == 1) THEN
           Y(1:N) = BETA*Y(1:N)
       ELSE
           IY = KY
           DO I = 1,N
               Y(IY) = BETA*Y(IY)
               IY = IY + INCY
           ENDDO
       END IF
   ENDIF
   IF (ALPHA == (0.0E+0,0.0E+0)) RETURN
   IF (LSAME(UPLO,'U')) THEN
!
!        Form  y  when A is stored in upper triangle.
!
       IF ((INCX == 1) .AND. (INCY == 1)) THEN
           DO J = 1,N
               TEMP1 = ALPHA*X(J)
               Y(1:J - 1) = Y(1:J - 1) + TEMP1*A(1:J - 1,J)
               Y(J) = Y(J) + TEMP1*REAL(A(J,J)) + ALPHA*sum(CONJG(A(1:J - 1,J))*X(1:J - 1))
           ENDDO
       ELSE
           JX = KX
           JY = KY
           DO J = 1,N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = (0.0E+0,0.0E+0)
               IX = KX
               IY = KY
               DO I = 1,J - 1
                   Y(IY) = Y(IY) + TEMP1*A(I,J)
                   TEMP2 = TEMP2 + CONJG(A(I,J))*X(IX)
                   IX = IX + INCX
                   IY = IY + INCY
               ENDDO
               Y(JY) = Y(JY) + TEMP1*REAL(A(J,J)) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
           ENDDO
       END IF
   ELSE
!
!        Form  y  when A is stored in lower triangle.
!
       IF ((INCX == 1) .AND. (INCY == 1)) THEN
           DO J = 1,N
               TEMP1 = ALPHA*X(J)
               Y(J) = Y(J) + TEMP1*REAL(A(J,J))
               Y(J + 1:N) = Y(J + 1:N) + TEMP1*A(J + 1:N,J)
               Y(J) = Y(J) + ALPHA*sum(CONJG(A(J + 1:N,J))*X(J + 1:N))
           ENDDO
       ELSE
           JX = KX
           JY = KY
           DO J = 1,N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = (0.0E+0,0.0E+0)
               Y(JY) = Y(JY) + TEMP1*REAL(A(J,J))
               IX = JX
               IY = JY
               DO I = J + 1,N
                   IX = IX + INCX
                   IY = IY + INCY
                   Y(IY) = Y(IY) + TEMP1*A(I,J)
                   TEMP2 = TEMP2 + CONJG(A(I,J))*X(IX)
               ENDDO
               Y(JY) = Y(JY) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
           ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of CHEMV
!
END
