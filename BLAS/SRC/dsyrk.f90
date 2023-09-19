!> \brief \b DSYRK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER K,LDA,LDC,N
!       CHARACTER TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYRK  performs one of the symmetric rank k operations
!>
!>    C := alpha*A*A**T + beta*C,
!>
!> or
!>
!>    C := alpha*A**T*A + beta*C,
!>
!> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
!> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
!> in the second case.
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
!>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
!>
!>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
!>
!>              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.
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
!>           of  columns   of  the   matrix   A,   and  on   entry   with
!>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!>           of rows of the matrix  A.  K must be at least zero.
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
!>          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
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
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension ( LDC, N )
!>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!>           upper triangular part of the array C must contain the upper
!>           triangular part  of the  symmetric matrix  and the strictly
!>           lower triangular part of C is not referenced.  On exit, the
!>           upper triangular part of the array  C is overwritten by the
!>           upper triangular part of the updated matrix.
!>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!>           lower triangular part of the array C must contain the lower
!>           triangular part  of the  symmetric matrix  and the strictly
!>           upper triangular part of C is not referenced.  On exit, the
!>           lower triangular part of the array  C is overwritten by the
!>           lower triangular part of the updated matrix.
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
!> \ingroup herk
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
   SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION ALPHA,BETA
   INTEGER K,LDA,LDC,N
   CHARACTER TRANS,UPLO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION A(LDA,*),C(LDC,*)
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
!     .. Intrinsic Functions ..
   INTRINSIC MAX
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION TEMP
   INTEGER I,INFO,J,L,NROWA
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
            (.NOT.LSAME(TRANS,'T')) .AND. &
            (.NOT.LSAME(TRANS,'C'))) THEN
       INFO = 2
   ELSE IF (N < 0) THEN
       INFO = 3
   ELSE IF (K < 0) THEN
       INFO = 4
   ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 7
   ELSE IF (LDC < MAX(1,N)) THEN
       INFO = 10
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('DSYRK ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N == 0) .OR. (((ALPHA == 0.0D+0).OR. &
       (K == 0)).AND. (BETA == 1.0D+0))) RETURN
!
!     And when  alpha.eq.zero.
!
   IF (ALPHA == 0.0D+0) THEN
       IF (UPPER) THEN
           DO J = 1,N
               C(1:J,J) = BETA*C(1:J,J)
           ENDDO
       ELSE
           DO J = 1,N
               C(J:N,J) = BETA*C(J:N,J)
           ENDDO
       END IF
       RETURN
   END IF
!
!     Start the operations.
!
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  C := alpha*A*A**T + beta*C.
!
       IF (UPPER) THEN
           DO J = 1,N
               C(1:J,J) = BETA*C(1:J,J)
               DO L = 1,K
                   C(1:J,J) = C(1:J,J) + ALPHA*A(J,L)*A(1:J,L)
               ENDDO
           ENDDO
       ELSE
           DO J = 1,N
               C(J:N,J) = BETA*C(J:N,J)
               DO L = 1,K
                   C(J:N,J) = C(J:N,J) + ALPHA*A(J,L)*A(J:N,L)
               ENDDO
           ENDDO
       END IF
   ELSE
!
!        Form  C := alpha*A**T*A + beta*C.
!
       IF (UPPER) THEN
           DO J = 1,N
               DO I = 1,J
                   C(I,J) = ALPHA*sum(A(1:K,I)*A(1:K,J)) + BETA*C(I,J)
               ENDDO
           ENDDO
       ELSE
           DO J = 1,N
               DO I = J,N
                   C(I,J) = ALPHA*sum(A(1:K,I)*A(1:KL,J)) + BETA*C(I,J)
               ENDDO
           ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of DSYRK
!
END
