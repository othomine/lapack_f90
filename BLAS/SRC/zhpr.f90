!> \brief \b ZHPR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER INCX,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 AP(*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHPR    performs the hermitian rank 1 operation
!>
!>    A := alpha*x*x**H + A,
!>
!> where alpha is a real scalar, x is an n element vector and A is an
!> n by n hermitian matrix, supplied in packed form.
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
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with  UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the hermitian matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the upper triangular part of the
!>           updated matrix.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the hermitian matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the lower triangular part of the
!>           updated matrix.
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set, they are assumed to be zero, and on exit they
!>           are set to zero.
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
!> \ingroup hpr
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
   SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION ALPHA
   INTEGER INCX,N
   CHARACTER UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX*16 AP(*),X(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   COMPLEX*16 TEMP
   INTEGER I,INFO,IX,J,JX,K,KK,KX
!     ..
!     .. External Functions ..
   LOGICAL LSAME
   EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC DBLE,DCONJG
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
       INFO = 5
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('ZHPR  ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N == 0) .OR. (ALPHA == DBLE((0.0D+0,0.0D+0)))) RETURN
!
!     Set the start point in X if the increment is not unity.
!
   IF (INCX <= 0) THEN
       KX = 1 - (N-1)*INCX
   ELSE IF (INCX /= 1) THEN
       KX = 1
   END IF
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
   KK = 1
   IF (LSAME(UPLO,'U')) THEN
!
!        Form  A  when upper triangle is stored in AP.
!
       IF (INCX == 1) THEN
           DO J = 1,N
               IF (X(J) /= (0.0D+0,0.0D+0)) THEN
                   TEMP = ALPHA*DCONJG(X(J))
                   K = KK
                   AP(KK:KK+J-2) = AP(KK:KK+J-2) + X(1:J-1)*TEMP
                   AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(J)*TEMP)
               ELSE
                   AP(KK+J-1) = DBLE(AP(KK+J-1))
               END IF
               KK = KK + J
           ENDDO
       ELSE
           JX = KX
           DO J = 1,N
               IF (X(JX) /= (0.0D+0,0.0D+0)) THEN
                   TEMP = ALPHA*DCONJG(X(JX))
                   IX = KX
                   DO K = KK,KK + J - 2
                       AP(K) = AP(K) + X(IX)*TEMP
                       IX = IX + INCX
                   ENDDO
                   AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(JX)*TEMP)
               ELSE
                   AP(KK+J-1) = DBLE(AP(KK+J-1))
               END IF
               JX = JX + INCX
               KK = KK + J
           ENDDO
       END IF
   ELSE
!
!        Form  A  when lower triangle is stored in AP.
!
       IF (INCX == 1) THEN
           DO J = 1,N
               IF (X(J) /= (0.0D+0,0.0D+0)) THEN
                   TEMP = ALPHA*DCONJG(X(J))
                   AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(J))
                   AP(KK+1:KK-J+N) = AP(KK+1:KK-J+N) + X(J+1:N)*TEMP
               ELSE
                   AP(KK) = DBLE(AP(KK))
               END IF
               KK = KK + N - J + 1
           ENDDO
       ELSE
           JX = KX
           DO J = 1,N
               IF (X(JX) /= (0.0D+0,0.0D+0)) THEN
                   TEMP = ALPHA*DCONJG(X(JX))
                   AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(JX))
                   IX = JX
                   DO K = KK + 1,KK + N - J
                       IX = IX + INCX
                       AP(K) = AP(K) + X(IX)*TEMP
                   ENDDO
               ELSE
                   AP(KK) = DBLE(AP(KK))
               END IF
               JX = JX + INCX
               KK = KK + N - J + 1
           ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of ZHPR
!
END
