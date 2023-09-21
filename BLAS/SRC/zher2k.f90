!> \brief \b ZHER2K
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA
!       DOUBLE PRECISION BETA
!       INTEGER K,LDA,LDB,LDC,N
!       CHARACTER TRANS,UPLO
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
!> ZHER2K  performs one of the hermitian rank 2k operations
!>
!>    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
!>
!> or
!>
!>    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
!>
!> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
!> hermitian matrix and  A and B  are  n by k matrices in the first case
!> and  k by n  matrices in the second case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!>           triangular  part  of the  array  C  is to be  referenced  as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!>                                  is to be referenced.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry,  TRANS  specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'    C := alpha*A*B**H          +
!>                                         conjg( alpha )*B*A**H +
!>                                         beta*C.
!>
!>              TRANS = 'C' or 'c'    C := alpha*A**H*B          +
!>                                         conjg( alpha )*B**H*A +
!>                                         beta*C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N specifies the order of the matrix C.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!>           of  columns  of the  matrices  A and B,  and on  entry  with
!>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
!>           matrices  A and B.  K must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16 .
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, ka ), where ka is
!>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by n  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!>           be at least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension ( LDB, kb ), where kb is
!>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  k by n  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!>           then  LDB must be at least  max( 1, n ), otherwise  LDB must
!>           be at least  max( 1, k ).
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION .
!>           On entry, BETA specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension ( LDC, N )
!>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!>           upper triangular part of the array C must contain the upper
!>           triangular part  of the  hermitian matrix  and the strictly
!>           lower triangular part of C is not referenced.  On exit, the
!>           upper triangular part of the array  C is overwritten by the
!>           upper triangular part of the updated matrix.
!>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!>           lower triangular part of the array C must contain the lower
!>           triangular part  of the  hermitian matrix  and the strictly
!>           upper triangular part of C is not referenced.  On exit, the
!>           lower triangular part of the array  C is overwritten by the
!>           lower triangular part of the updated matrix.
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set,  they are assumed to be zero,  and on exit they
!>           are set to zero.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!>           max( 1, n ).
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
!> \ingroup her2k
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
!>  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.
!>     Ed Anderson, Cray Research Inc.
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16 ALPHA
   DOUBLE PRECISION BETA
   INTEGER K,LDA,LDB,LDC,N,I
   CHARACTER TRANS,UPLO
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
   COMPLEX*16 TEMP1,TEMP2
   INTEGER INFO,J,L,NROWA
   LOGICAL UPPER
!     ..
!
!     Test the input parameters.
!
   IF (LSAME(TRANS,'N')) THEN
       NROWA = N
   ELSE
       NROWA = K
   END IF
   UPPER = LSAME(UPLO,'U')
!
   INFO = 0
   IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
       INFO = 1
   ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND. &
            (.NOT.LSAME(TRANS,'C'))) THEN
       INFO = 2
   ELSE IF (N < 0) THEN
       INFO = 3
   ELSE IF (K < 0) THEN
       INFO = 4
   ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 7
   ELSE IF (LDB < MAX(1,NROWA)) THEN
       INFO = 9
   ELSE IF (LDC < MAX(1,N)) THEN
       INFO = 12
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('ZHER2K',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N == 0) .OR. (((ALPHA == (0.0D+0,0.0D+0)).OR. &
       (K == 0)).AND. (BETA == 1.0D+0))) RETURN
!
!     And when  alpha.eq.zero.
!
   IF (UPPER) THEN
       IF (BETA == DBLE((0.0D+0,0.0D+0))) THEN
           DO J = 1,N
               C(1:J,J) = (0.0D+0,0.0D+0)
           ENDDO
       ELSE
           DO J = 1,N
               C(1:J-1,J) = BETA*C(1:J-1,J)
               C(J,J) = BETA*DBLE(C(J,J))
           ENDDO
       END IF
   ELSE
       IF (BETA == DBLE((0.0D+0,0.0D+0))) THEN
           DO J = 1,N
               C(J:N,J) = (0.0D+0,0.0D+0)
           ENDDO
       ELSE
           DO J = 1,N
               C(J,J) = BETA*DBLE(C(J,J))
               C(J+1:N,J) = BETA*C(J+1:N,J)
           ENDDO
       END IF
   END IF
   IF (ALPHA == (0.0D+0,0.0D+0)) RETURN
!
!     Start the operations.
!
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
!                   C.
!
       IF (UPPER) THEN
           DO J = 1,N
               DO L = 1,K
                   IF ((A(J,L) /= (0.0D+0,0.0D+0)) .OR. (B(J,L) /= (0.0D+0,0.0D+0))) THEN
                       TEMP1 = ALPHA*DCONJG(B(J,L))
                       TEMP2 = DCONJG(ALPHA*A(J,L))
                       C(1:J-1,J) = C(1:J-1,J) + A(1:J-1,L)*TEMP1 + B(1:J-1,L)*TEMP2
                       C(J,J) = DBLE(C(J,J)) + DBLE(A(J,L)*TEMP1+B(J,L)*TEMP2)
                   END IF
               ENDDO
           ENDDO
       ELSE
           DO J = 1,N
               DO L = 1,K
                   IF ((A(J,L) /= (0.0D+0,0.0D+0)) .OR. (B(J,L) /= (0.0D+0,0.0D+0))) THEN
                       TEMP1 = ALPHA*DCONJG(B(J,L))
                       TEMP2 = DCONJG(ALPHA*A(J,L))
                       C(J+1:N,J) = C(J+1:N,J) + A(J+1:N,L)*TEMP1 + B(J+1:N,L)*TEMP2
                       C(J,J) = DBLE(C(J,J)) + DBLE(A(J,L)*TEMP1+B(J,L)*TEMP2)
                   END IF
               ENDDO
           ENDDO
       END IF
   ELSE
!
!        Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
!                   C.
!
       IF (UPPER) THEN
           DO J = 1,N
               DO I = 1,J - 1
                   C(I,J) = C(I,J) + ALPHA*sum(DCONJG(A(1:K,I))*B(1:K,J)) + DCONJG(ALPHA)*sum(DCONJG(B(1:K,I))*A(1:K,J))
               ENDDO
               C(J,J) = DBLE(C(J,J)) + DBLE(ALPHA*sum(DCONJG(A(1:K,J))*B(1:K,J)) + DCONJG(ALPHA)*sum(DCONJG(B(1:K,J))*A(1:K,J)))
           ENDDO
       ELSE
           DO J = 1,N
               C(J,J) = DBLE(C(J,J)) + DBLE(ALPHA*sum(DCONJG(A(1:K,J))*B(1:K,J))+DCONJG(ALPHA)*sum(DCONJG(B(1:K,J))*A(1:K,J)))
               DO I = J + 1,N
                   C(I,J) = C(I,J) + ALPHA*sum(DCONJG(A(1:K,I))*B(1:K,J)) + DCONJG(ALPHA)*sum(DCONJG(B(1:K,I))*A(1:K,J))
               ENDDO
           ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of ZHER2K
!
END
