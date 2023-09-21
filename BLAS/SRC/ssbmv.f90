!> \brief \b SSBMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
!       INTEGER INCX,INCY,K,LDA,N
!       CHARACTER UPLO
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
!> SSBMV  performs the matrix-vector  operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric band matrix, with k super-diagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the band matrix A is being supplied as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  being supplied.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  being supplied.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry, K specifies the number of super-diagonals of the
!>           matrix A. K must satisfy  0 .le. K.
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
!>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!>           by n part of the array A must contain the upper triangular
!>           band part of the symmetric matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row
!>           ( k + 1 ) of the array, the first super-diagonal starting at
!>           position 2 in row k, and so on. The top left k by k triangle
!>           of the array A is not referenced.
!>           The following program segment will transfer the upper
!>           triangular part of a symmetric band matrix from conventional
!>           full matrix storage to band storage:
!>
!>                 DO  J = 1, N
!>                    M = K + 1 - J
!>                    DO  I = MAX( 1, J - K ), J
!>                       A( M + I, J ) = matrix( I, J )
!>              10    ENDDO
!>              20 ENDDO
!>
!>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!>           by n part of the array A must contain the lower triangular
!>           band part of the symmetric matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row 1 of
!>           the array, the first sub-diagonal starting at position 1 in
!>           row 2, and so on. The bottom right k by k triangle of the
!>           array A is not referenced.
!>           The following program segment will transfer the lower
!>           triangular part of a symmetric band matrix from conventional
!>           full matrix storage to band storage:
!>
!>                 DO  J = 1, N
!>                    M = 1 - J
!>                    DO  I = J, MIN( N, J + K )
!>                       A( M + I, J ) = matrix( I, J )
!>              10    ENDDO
!>              20 ENDDO
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           ( k + 1 ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
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
!>           On entry, BETA specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the
!>           vector y. On exit, Y is overwritten by the updated vector y.
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
!> \ingroup hbmv
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
   SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   REAL ALPHA,BETA
   INTEGER INCX,INCY,K,LDA,N
   CHARACTER UPLO
!     ..
!     .. Array Arguments ..
   REAL A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   REAL TEMP1,TEMP2
   INTEGER I,INFO,IX,IY,J,JX,JY,KPLUS1,KX,KY,L,MINI
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
   ELSE IF (K < 0) THEN
       INFO = 3
   ELSE IF (LDA <  (K+1)) THEN
       INFO = 6
   ELSE IF (INCX == 0) THEN
       INFO = 8
   ELSE IF (INCY == 0) THEN
       INFO = 11
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('SSBMV ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N == 0) .OR. ((ALPHA == 0.0E+0).AND. (BETA == 1.0E+0))) RETURN
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
!     Start the operations. In this version the elements of the array A
!     are accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
   IF (BETA /= 1.0E+0) THEN
       IF (INCY == 1) THEN
           IF (BETA == 0.0E+0) THEN
               Y(1:N) = 0.0E+0
           ELSE
               Y(1:N) = BETA*Y(1:N)
           END IF
       ELSE
           IY = KY
           IF (BETA == 0.0E+0) THEN
               DO I = 1,N
                   Y(IY) = 0.0E+0
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
   IF (ALPHA == 0.0E+0) RETURN
   IF (LSAME(UPLO,'U')) THEN
!
!        Form  y  when upper triangle of A is stored.
!
       KPLUS1 = K + 1
       IF ((INCX == 1) .AND. (INCY == 1)) THEN
           DO J = 1,N
               MINI = MAX(1,J-K)
               TEMP1 = ALPHA*X(J)
               TEMP2 = sum(A(KPLUS1 - J+MINI:KPLUS1 - 1,J)*X(MINI:J - 1))
               Y(MINI:J - 1) = Y(MINI:J - 1) + TEMP1*A(KPLUS1 - J+MINI:KPLUS1 - 1,J)
               Y(J) = Y(J) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
           ENDDO
       ELSE
           JX = KX
           JY = KY
           DO J = 1,N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = 0.0E+0
               IX = KX
               IY = KY
               L = KPLUS1 - J
               DO I = MAX(1,J-K),J - 1
                   Y(IY) = Y(IY) + TEMP1*A(L+I,J)
                   TEMP2 = TEMP2 + A(L+I,J)*X(IX)
                   IX = IX + INCX
                   IY = IY + INCY
               ENDDO
               Y(JY) = Y(JY) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
               IF (J > K) THEN
                   KX = KX + INCX
                   KY = KY + INCY
               END IF
           ENDDO
       END IF
   ELSE
!
!        Form  y  when lower triangle of A is stored.
!
       IF ((INCX == 1) .AND. (INCY == 1)) THEN
           DO J = 1,N
               MINI = MIN(N,J+K)
               TEMP1 = ALPHA*X(J)
               TEMP2 = sum(A(2:1-J+MINI,J)*X(J+1:MINI))
               Y(J) = Y(J) + TEMP1*A(1,J) + ALPHA*TEMP2
               Y(J + 1:MINI) = Y(J + 1:MINI) + TEMP1*A(2:1-J+MINI,J)
           ENDDO
       ELSE
           JX = KX
           JY = KY
           DO J = 1,N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = 0.0E+0
               Y(JY) = Y(JY) + TEMP1*A(1,J)
               L = 1 - J
               IX = JX
               IY = JY
               DO I = J + 1,MIN(N,J+K)
                   IX = IX + INCX
                   IY = IY + INCY
                   Y(IY) = Y(IY) + TEMP1*A(L+I,J)
                   TEMP2 = TEMP2 + A(L+I,J)*X(IX)
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
!     End of SSBMV
!
END
