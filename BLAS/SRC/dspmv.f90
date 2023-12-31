!> \brief \b DSPMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER INCX,INCY,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION AP(*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPMV  performs the matrix-vector operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix, supplied in packed form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the matrix A is supplied in the packed
!>           array AP as follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  supplied in AP.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  supplied in AP.
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
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension at least
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
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension at least
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
!> \ingroup hpmv
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
   SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION ALPHA,BETA
   INTEGER INCX,INCY,N
   CHARACTER UPLO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION AP(*),X(*),Y(*)
!     ..
!
!  =====================================================================
!     .. Local Scalars ..
   DOUBLE PRECISION TEMP1,TEMP2
   INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
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
   IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
       INFO = 1
   ELSE IF (N < 0) THEN
       INFO = 2
   ELSE IF (INCX == 0) THEN
       INFO = 6
   ELSE IF (INCY == 0) THEN
       INFO = 9
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('DSPMV ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N == 0) .OR. ((ALPHA == 0.0D+0).AND. (BETA == 1.0D+0))) RETURN
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
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
!     First form  y := beta*y.
!
   IF (BETA /= 1.0D+0) THEN
       IF (INCY == 1) THEN
           Y(1:N) = BETA*Y(1:N)
       ELSE
           IY = KY
           IF (BETA == 0.0D+0) THEN
               DO I = 1,N
                   Y(IY) = 0.0D+0
                   IY = IY + INCY
               ENDDO
           ELSE
               DO I = 1,N
                   Y(IY) = BETA*Y(IY)
                   IY = IY + INCY
               ENDDO
           END IF
       END IF
   END IF
   IF (ALPHA == 0.0D+0) RETURN
   KK = 1
   IF (LSAME(UPLO,'U')) THEN
!
!        Form  y  when AP contains the upper triangle.
!
       IF ((INCX == 1) .AND. (INCY == 1)) THEN
           DO J = 1,N
               TEMP1 = ALPHA*X(J)
               Y(1:J-1) = Y(1:J-1) + TEMP1*AP(KK:KK+J-2)
               Y(J) = Y(J) + TEMP1*AP(KK+J-1) + ALPHA*sum(AP(KK:KK+J-2)*X(1:J-1))
               KK = KK + J
           ENDDO
       ELSE
           JX = KX
           JY = KY
           DO J = 1,N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = 0.0D+0
               IX = KX
               IY = KY
               DO K = KK,KK + J - 2
                   Y(IY) = Y(IY) + TEMP1*AP(K)
                   TEMP2 = TEMP2 + AP(K)*X(IX)
                   IX = IX + INCX
                   IY = IY + INCY
               ENDDO
               Y(JY) = Y(JY) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + J
           ENDDO
       END IF
   ELSE
!
!        Form  y  when AP contains the lower triangle.
!
       IF ((INCX == 1) .AND. (INCY == 1)) THEN
           DO J = 1,N
               TEMP1 = ALPHA*X(J)
               Y(J) = Y(J) + TEMP1*AP(KK) + ALPHA*sum(AP(KK+1:KK-J+N)*X(J+1:N))
               Y(J+1:N) = Y(J+1:N) + TEMP1*AP(KK+1:KK-J+N)
               KK = KK + (N-J+1)
           ENDDO
       ELSE
           JX = KX
           JY = KY
           DO J = 1,N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = 0.0D+0
               Y(JY) = Y(JY) + TEMP1*AP(KK)
               IX = JX
               IY = JY
               DO K = KK + 1,KK + N - J
                   IX = IX + INCX
                   IY = IY + INCY
                   Y(IY) = Y(IY) + TEMP1*AP(K)
                   TEMP2 = TEMP2 + AP(K)*X(IX)
               ENDDO
               Y(JY) = Y(JY) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + (N-J+1)
           ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of DSPMV
!
END
