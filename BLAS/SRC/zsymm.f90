!> \brief \b ZSYMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA,BETA
!       INTEGER LDA,LDB,LDC,M,N
!       CHARACTER SIDE,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSYMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*A*B + beta*C,
!>
!> or
!>
!>    C := alpha*B*A + beta*C,
!>
!> where  alpha and beta are scalars, A is a symmetric matrix and  B and
!> C are m by n matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry,  SIDE  specifies whether  the  symmetric matrix  A
!>           appears on the  left or right  in the  operation as follows:
!>
!>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
!>
!>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!>           triangular  part  of  the  symmetric  matrix   A  is  to  be
!>           referenced as follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of the
!>                                  symmetric matrix is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of the
!>                                  symmetric matrix is to be referenced.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies the number of rows of the matrix  C.
!>           M  must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix C.
!>           N  must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, ka ), where ka is
!>           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
!>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
!>           the array  A  must contain the  symmetric matrix,  such that
!>           when  UPLO = 'U' or 'u', the leading m by m upper triangular
!>           part of the array  A  must contain the upper triangular part
!>           of the  symmetric matrix and the  strictly  lower triangular
!>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!>           the leading  m by m  lower triangular part  of the  array  A
!>           must  contain  the  lower triangular part  of the  symmetric
!>           matrix and the  strictly upper triangular part of  A  is not
!>           referenced.
!>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
!>           the array  A  must contain the  symmetric matrix,  such that
!>           when  UPLO = 'U' or 'u', the leading n by n upper triangular
!>           part of the array  A  must contain the upper triangular part
!>           of the  symmetric matrix and the  strictly  lower triangular
!>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!>           the leading  n by n  lower triangular part  of the  array  A
!>           must  contain  the  lower triangular part  of the  symmetric
!>           matrix and the  strictly upper triangular part of  A  is not
!>           referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least max( 1, n ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension ( LDB, N )
!>           Before entry, the leading  m by n part of the array  B  must
!>           contain the matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16
!>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!>           supplied as zero then C need not be set on input.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension ( LDC, N )
!>           Before entry, the leading  m by n  part of the array  C must
!>           contain the matrix  C,  except when  beta  is zero, in which
!>           case C need not be set on entry.
!>           On exit, the array  C  is overwritten by the  m by n updated
!>           matrix.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
!> \ingroup hemm
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
   SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16 ALPHA,BETA
   INTEGER LDA,LDB,LDC,M,N
   CHARACTER SIDE,UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
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
   COMPLEX*16 TEMP1
   INTEGER I,INFO,J,K,NROWA
   LOGICAL UPPER
!     ..
!
!     Set NROWA as the number of rows of A.
!
   IF (LSAME(SIDE,'L')) THEN
       NROWA = M
   ELSE
       NROWA = N
   END IF
   UPPER = LSAME(UPLO,'U')
!
!     Test the input parameters.
!
   INFO = 0
   IF ((.NOT.LSAME(SIDE,'L')) .AND. &
       (.NOT.LSAME(SIDE,'R'))) THEN
       INFO = 1
   ELSE IF ((.NOT.UPPER) .AND. &
            (.NOT.LSAME(UPLO,'L'))) THEN
       INFO = 2
   ELSE IF (M < 0) THEN
       INFO = 3
   ELSE IF (N < 0) THEN
       INFO = 4
   ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 7
   ELSE IF (LDB < MAX(1,M)) THEN
       INFO = 9
   ELSE IF (LDC < MAX(1,M)) THEN
       INFO = 12
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('ZSYMM ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((M == 0) .OR. (N == 0) .OR. &
       ((ALPHA == (0.0D+0,0.0D+0)).AND. (BETA == (1.0D+0,0.0D+0)))) RETURN
!
!     And when  alpha.eq.zero.
!
   IF (BETA == (0.0D+0,0.0D+0)) THEN
       C(1:M,1:N) = (0.0D+0,0.0D+0)
   ELSE
       C(1:M,1:N) = BETA*C(1:M,1:N)
   END IF
   IF (ALPHA == (0.0D+0,0.0D+0)) RETURN
!
!     Start the operations.
!
   IF (LSAME(SIDE,'L')) THEN
!
!        Form  C := alpha*A*B + beta*C.
!
       IF (UPPER) THEN
           DO J = 1,N
               DO I = 1,M
                   TEMP1 = ALPHA*B(I,J)
                   C(1:I-1,J) = C(1:I-1,J) + TEMP1*A(1:I-1,I)
                   C(I,J) = C(I,J) + TEMP1*A(I,I) + ALPHA*sum(B(1:I-1,J)*A(1:I-1,I))
               ENDDO
           ENDDO
       ELSE
           DO J = 1,N
               DO I = M,1,-1
                   TEMP1 = ALPHA*B(I,J)
                   C(I+1:M,J) = C(I+1:M,J) + TEMP1*A(I+1:M,I)
                   C(I,J) = C(I,J) + TEMP1*A(I,I) + ALPHA*sum(B(I+1:M,J)*A(I+1:M,I))
               ENDDO
           ENDDO
       END IF
   ELSE
!
!        Form  C := alpha*B*A + beta*C.
!
       IF (UPPER) THEN
           DO J = 1,N
               C(1:M,J) = C(1:M,J) + ALPHA*A(J,J)*B(1:M,J)
               DO K = 1,J - 1
                   C(1:M,J) = C(1:M,J) + ALPHA*A(K,J)*B(1:M,K)
               ENDDO
               DO K = J + 1,N
                   C(1:M,J) = C(1:M,J) + ALPHA*A(J,K)*B(1:M,K)
               ENDDO
           ENDDO
       ELSE
           DO J = 1,N
               C(1:M,J) = C(1:M,J) + ALPHA*A(J,J)*B(1:M,J)
               DO K = 1,J - 1
                   C(1:M,J) = C(1:M,J) + ALPHA*A(J,K)*B(1:M,K)
               ENDDO
               DO K = J + 1,N
                   C(1:M,J) = C(1:M,J) + ALPHA*A(K,J)*B(1:M,K)
               ENDDO
           ENDDO
       ENDIF
   END IF
!
   RETURN
!
!     End of ZSYMM
!
END
