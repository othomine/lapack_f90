!> \brief \b ZTRSM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),B(LDB,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRSM  solves one of the matrix equations
!>
!>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!>
!> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
!>
!> The matrix X is overwritten on B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry, SIDE specifies whether op( A ) appears on the left
!>           or right of X as follows:
!>
!>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!>
!>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix A is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n'   op( A ) = A.
!>
!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c'   op( A ) = A**H.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit triangular
!>           as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, k ),
!>           where k is m when SIDE = 'L' or 'l'
!>             and k is n when SIDE = 'R' or 'r'.
!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!>           upper triangular part of the array  A must contain the upper
!>           triangular matrix  and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!>           lower triangular part of the array  A must contain the lower
!>           triangular matrix  and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!>           A  are not referenced either,  but are assumed to be  unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension ( LDB, N )
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain  the  right-hand  side  matrix  B,  and  on exit  is
!>           overwritten by the solution matrix  X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
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
!> \ingroup trsm
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16 ALPHA
   INTEGER LDA,LDB,M,N
   CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX*16 A(LDA,*),B(LDB,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
   LOGICAL LSAME
   EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL XERBLA
!     ..
!     .. Local Scalars ..
   COMPLEX*16 TEMP
   INTEGER I,INFO,J,K,NROWA
   LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
!     ..
!
!     Test the input parameters.
!
   LSIDE = LSAME(SIDE,'L')
   IF (LSIDE) THEN
       NROWA = M
   ELSE
       NROWA = N
   END IF
   NOCONJ = LSAME(TRANSA,'T')
   NOUNIT = LSAME(DIAG,'N')
   UPPER = LSAME(UPLO,'U')
!
   INFO = 0
   IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
       INFO = 1
   ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
       INFO = 2
   ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
            (.NOT.LSAME(TRANSA,'T')) .AND. &
            (.NOT.LSAME(TRANSA,'C'))) THEN
       INFO = 3
   ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. &
            (.NOT.LSAME(DIAG,'N'))) THEN
       INFO = 4
   ELSE IF (M < 0) THEN
       INFO = 5
   ELSE IF (N < 0) THEN
       INFO = 6
   ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 9
   ELSE IF (LDB < MAX(1,M)) THEN
       INFO = 11
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('ZTRSM ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF (M == 0 .OR. N == 0) RETURN
!
!     And when  alpha.eq.zero.
!
   IF (ALPHA == (0.0D+0,0.0D+0)) THEN
       B(1:M,1:N) = (0.0D+0,0.0D+0)
       RETURN
   END IF
!
!     Start the operations.
!
   IF (LSIDE) THEN
       IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*inv( A )*B.
!
           IF (UPPER) THEN
               DO J = 1,N
                   IF (ALPHA /= (1.0D+0,0.0D+0)) THEN
                       B(1:M,J) = ALPHA*B(1:M,J)
                   END IF
                   DO K = M,1,-1
                       IF (B(K,J) /= (0.0D+0,0.0D+0)) THEN
                           IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                           B(1:K-1,J) = B(1:K-1,J) - B(K,J)*A(1:K-1,K)
                       END IF
                   ENDDO
               ENDDO
           ELSE
               DO J = 1,N
                   IF (ALPHA /= (1.0D+0,0.0D+0)) THEN
                       B(1:M,J) = ALPHA*B(1:M,J)
                   END IF
                   DO K = 1,M
                       IF (B(K,J) /= (0.0D+0,0.0D+0)) THEN
                           IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                           B(K+1:M,J) = B(K+1:M,J) - B(K,J)*A(K+1:M,K)
                       END IF
                   ENDDO
               ENDDO
           END IF
       ELSE
!
!           Form  B := alpha*inv( A**T )*B
!           or    B := alpha*inv( A**H )*B.
!
           IF (UPPER) THEN
               DO J = 1,N
                   DO I = 1,M
                       TEMP = ALPHA*B(I,J)
                       IF (NOCONJ) THEN
                           TEMP = TEMP - sum(A(1:I-1,I)*B(1:I-1,J))
                           IF (NOUNIT) TEMP = TEMP/A(I,I)
                       ELSE
                           TEMP = TEMP - sum(DCONJG(A(1:I-1,I))*B(1:I-1,J))
                           IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                       END IF
                       B(I,J) = TEMP
                   ENDDO
               ENDDO
           ELSE
               DO J = 1,N
                   DO I = M,1,-1
                       TEMP = ALPHA*B(I,J)
                       IF (NOCONJ) THEN
                           TEMP = TEMP - sum(A(I+1:M,I)*B(I+1:M,J))
                           IF (NOUNIT) TEMP = TEMP/A(I,I)
                       ELSE
                           TEMP = TEMP - sum(DCONJG(A(I+1:M,I))*B(I+1:M,J))
                           IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                       END IF
                       B(I,J) = TEMP
                   ENDDO
               ENDDO
           END IF
       END IF
   ELSE
       IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*inv( A ).
!
           IF (UPPER) THEN
               DO J = 1,N
                   IF (ALPHA /= (1.0D+0,0.0D+0)) THEN
                       B(1:M,J) = ALPHA*B(1:M,J)
                   END IF
                   DO K = 1,J - 1
                       IF (A(K,J) /= (0.0D+0,0.0D+0)) THEN
                           B(1:M,J) = B(1:M,J) - A(K,J)*B(1:M,K)
                       END IF
                   ENDDO
                   IF (NOUNIT) THEN
                       B(1:M,J) = B(1:M,J)/A(J,J)
                   END IF
               ENDDO
           ELSE
               DO J = N,1,-1
                   IF (ALPHA /= (1.0D+0,0.0D+0)) THEN
                       B(1:M,J) = ALPHA*B(1:M,J)
                   END IF
                   DO K = J + 1,N
                       IF (A(K,J) /= (0.0D+0,0.0D+0)) THEN
                           B(1:M,J) = B(1:M,J) - A(K,J)*B(1:M,K)
                       END IF
                   ENDDO
                   IF (NOUNIT) THEN
                       TEMP = (1.0D+0,0.0D+0)/A(J,J)
                       B(1:M,J) = TEMP*B(1:M,J)
                   END IF
               ENDDO
           END IF
       ELSE
!
!           Form  B := alpha*B*inv( A**T )
!           or    B := alpha*B*inv( A**H ).
!
           IF (UPPER) THEN
               DO K = N,1,-1
                   IF (NOUNIT) THEN
                       IF (NOCONJ) THEN
                           TEMP = (1.0D+0,0.0D+0)/A(K,K)
                       ELSE
                           TEMP = (1.0D+0,0.0D+0)/DCONJG(A(K,K))
                       END IF
                       B(1:M,K) = TEMP*B(1:M,K)
                   END IF
                   DO J = 1,K - 1
                       IF (A(J,K) /= (0.0D+0,0.0D+0)) THEN
                           IF (NOCONJ) THEN
                               TEMP = A(J,K)
                           ELSE
                               TEMP = DCONJG(A(J,K))
                           END IF
                           B(1:M,J) = B(1:M,J) - TEMP*B(1:M,K)
                       END IF
                   ENDDO
                   IF (ALPHA /= (1.0D+0,0.0D+0)) THEN
                       B(1:M,K) = ALPHA*B(1:M,K)
                   END IF
               ENDDO
           ELSE
               DO K = 1,N
                   IF (NOUNIT) THEN
                       IF (NOCONJ) THEN
                           TEMP = (1.0D+0,0.0D+0)/A(K,K)
                       ELSE
                           TEMP = (1.0D+0,0.0D+0)/DCONJG(A(K,K))
                       END IF
                       B(1:M,K) = TEMP*B(1:M,K)
                   END IF
                   DO J = K + 1,N
                       IF (A(J,K) /= (0.0D+0,0.0D+0)) THEN
                           IF (NOCONJ) THEN
                               TEMP = A(J,K)
                           ELSE
                               TEMP = DCONJG(A(J,K))
                           END IF
                           B(1:M,J) = B(1:M,J) - TEMP*B(1:M,K)
                       END IF
                   ENDDO
                   IF (ALPHA /= (1.0D+0,0.0D+0)) THEN
                       B(1:M,K) = ALPHA*B(1:M,K)
                   END IF
               ENDDO
           END IF
       END IF
   END IF
!
   RETURN
!
!     End of ZTRSM
!
END
