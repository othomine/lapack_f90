!> \brief \b CTBSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,K,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTBSV  solves one of the systems of equations
!>
!>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
!>
!> where b and x are n element vectors and A is an n by n unit, or
!> non-unit, upper or lower triangular band matrix, with ( k + 1 )
!> diagonals.
!>
!> No test for singularity or near-singularity is included in this
!> routine. Such tests must be performed before calling this routine.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the equations to be solved as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   A*x = b.
!>
!>              TRANS = 'T' or 't'   A**T*x = b.
!>
!>              TRANS = 'C' or 'c'   A**H*x = b.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit
!>           triangular as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
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
!>           On entry with UPLO = 'U' or 'u', K specifies the number of
!>           super-diagonals of the matrix A.
!>           On entry with UPLO = 'L' or 'l', K specifies the number of
!>           sub-diagonals of the matrix A.
!>           K must satisfy  0 .le. K.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!>           by n part of the array A must contain the upper triangular
!>           band part of the matrix of coefficients, supplied column by
!>           column, with the leading diagonal of the matrix in row
!>           ( k + 1 ) of the array, the first super-diagonal starting at
!>           position 2 in row k, and so on. The top left k by k triangle
!>           of the array A is not referenced.
!>           The following program segment will transfer an upper
!>           triangular band matrix from conventional full matrix storage
!>           to band storage:
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
!>           band part of the matrix of coefficients, supplied column by
!>           column, with the leading diagonal of the matrix in row 1 of
!>           the array, the first sub-diagonal starting at position 1 in
!>           row 2, and so on. The bottom right k by k triangle of the
!>           array A is not referenced.
!>           The following program segment will transfer a lower
!>           triangular band matrix from conventional full matrix storage
!>           to band storage:
!>
!>                 DO  J = 1, N
!>                    M = 1 - J
!>                    DO  I = J, MIN( N, J + K )
!>                       A( M + I, J ) = matrix( I, J )
!>              10    ENDDO
!>              20 ENDDO
!>
!>           Note that when DIAG = 'U' or 'u' the elements of the array A
!>           corresponding to the diagonal elements of the matrix are not
!>           referenced, but are assumed to be unity.
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
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element right-hand side vector b. On exit, X is overwritten
!>           with the solution vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
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
!> \ingroup tbsv
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
   SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER INCX,K,LDA,N
   CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX ZERO
   PARAMETER (ZERO= (0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
   COMPLEX TEMP
   INTEGER I,INFO,IX,J,JX,KX,L
   LOGICAL NOCONJ,NOUNIT
!     ..
!     .. External Functions ..
   LOGICAL LSAME
   EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC CONJG,MAX,MIN
!     ..
!
!     Test the input parameters.
!
   INFO = 0
   IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
       INFO = 1
   ELSE IF (.NOT.LSAME(TRANS,'N') .AND. &
            .NOT.LSAME(TRANS,'T') .AND. &
            .NOT.LSAME(TRANS,'C')) THEN
       INFO = 2
   ELSE IF (.NOT.LSAME(DIAG,'U') .AND. &
            .NOT.LSAME(DIAG,'N')) THEN
       INFO = 3
   ELSE IF (N < 0) THEN
       INFO = 4
   ELSE IF (K < 0) THEN
       INFO = 5
   ELSE IF (LDA <  (K+1)) THEN
       INFO = 7
   ELSE IF (INCX == 0) THEN
       INFO = 9
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('CTBSV ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF (N == 0) RETURN
!
   NOCONJ = LSAME(TRANS,'T')
   NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
   IF (INCX <= 0) THEN
       KX = 1 - (N-1)*INCX
   ELSE IF (INCX /= 1) THEN
       KX = 1
   END IF
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := inv( A )*x.
!
       IF (LSAME(UPLO,'U')) THEN
           IF (INCX == 1) THEN
               DO J = N,1,-1
                   IF (X(J) /= (0.0E+0,0.0E+0)) THEN
                       IF (NOUNIT) X(J) = X(J)/A(K + 1,J)
                       MINI = MAX(1,J-K)
                       X(MINI:J-1) = X(MINI:J-1) - X(J)*A(K+1-J+MINI:K,J)
                   END IF
               ENDDO
           ELSE
               KX = KX + (N-1)*INCX
               JX = KX
               DO J = N,1,-1
                   KX = KX - INCX
                   IF (X(JX) /= (0.0E+0,0.0E+0)) THEN
                       IX = KX
                       L = K + 1 - J
                       IF (NOUNIT) X(JX) = X(JX)/A(K + 1,J)
                       TEMP = X(JX)
                       MINI = MAX(1,J-K)
                       DO I = J - 1,MINI,-1
                           X(IX) = X(IX) - TEMP*A(L+I,J)
                           IX = IX - INCX
                       ENDDO
                   END IF
                   JX = JX - INCX
               ENDDO
           END IF
       ELSE
           IF (INCX == 1) THEN
               DO J = 1,N
                   IF (X(J) /= (0.0E+0,0.0E+0)) THEN
                       IF (NOUNIT) X(J) = X(J)/A(1,J)
                       MAXI = MIN(N,J+K)
                       X(J+1:MAXI) = X(J+1:MAXI) - X(J)*A(2:1-J+MAXI,J)
                   END IF
               ENDDO
           ELSE
               JX = KX
               DO J = 1,N
                   KX = KX + INCX
                   IF (X(JX) /= (0.0E+0,0.0E+0)) THEN
                       IX = KX
                       L = 1 - J
                       IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                       TEMP = X(JX)
                       MAXI = MIN(N,J+K)
                       DO I = J + 1,MAXI
                           X(IX) = X(IX) - TEMP*A(L+I,J)
                           IX = IX + INCX
                       ENDDO
                   END IF
                   JX = JX + INCX
               ENDDO
           END IF
       END IF
   ELSE
!
!        Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
!
       IF (LSAME(UPLO,'U')) THEN
           IF (INCX == 1) THEN
               DO J = 1,N
                   MINI = MAX(1,J-K)
                   IF (NOCONJ) THEN
                       TEMP = X(J) - sum(A(K+1-J+MINI:K,J)*X(MINI:J-1))
                       IF (NOUNIT) TEMP = TEMP/A(K+1,J)
                   ELSE
                       TEMP = X(J) - sum(CONJG(A(K+1-J+MINI:K,J))*X(MINI:J-1))
                       IF (NOUNIT) TEMP = TEMP/CONJG(A(K+1,J))
                   END IF
                   X(J) = TEMP
               ENDDO
           ELSE
               JX = KX
               DO J = 1,N
                   TEMP = X(JX)
                   IX = KX
                   L = K+1-J
                   MINI = MAX(1,J-K)
                   IF (NOCONJ) THEN
                       DO I = MINI,J - 1
                           TEMP = TEMP - A(L+I,J)*X(IX)
                           IX = IX + INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/A(K+1,J)
                   ELSE
                       DO I = MINI,J - 1
                           TEMP = TEMP - CONJG(A(L+I,J))*X(IX)
                           IX = IX + INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/CONJG(A(K+1,J))
                   END IF
                   X(JX) = TEMP
                   JX = JX + INCX
                   IF (J > K) KX = KX + INCX
               ENDDO
           END IF
       ELSE
           IF (INCX == 1) THEN
               DO J = N,1,-1
                   MAXI = MIN(N,J+K)
                   IF (NOCONJ) THEN
                       TEMP = X(J) - sum(A(2:1-J+MAXI,J)*X(J+1:MAXI))
                       IF (NOUNIT) TEMP = TEMP/A(1,J)
                   ELSE
                       TEMP = X(J) - sum(CONJG(A(2:1-J+MAXI,J))*X(J+1:MAXI))
                       IF (NOUNIT) TEMP = TEMP/CONJG(A(1,J))
                   END IF
                   X(J) = TEMP
               ENDDO
           ELSE
               KX = KX + (N-1)*INCX
               JX = KX
               DO J = N,1,-1
                   TEMP = X(JX)
                   IX = KX
                   L = 1 - J
                   MAXI = MIN(N,J+K)
                   IF (NOCONJ) THEN
                       DO I = MAXI, J + 1, -1
                           TEMP = TEMP - A(L+I,J)*X(IX)
                           IX = IX - INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/A(1,J)
                   ELSE
                       DO I = MAXI, J + 1, -1
                           TEMP = TEMP - CONJG(A(L+I,J))*X(IX)
                           IX = IX - INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/CONJG(A(1,J))
                   END IF
                   X(JX) = TEMP
                   JX = JX - INCX
                   IF ((N-J) >= K) KX = KX - INCX
               ENDDO
           END IF
       END IF
   END IF
!
   RETURN
!
!     End of CTBSV
!
END
