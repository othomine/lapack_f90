!> \brief \b CHER2K
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA
!       REAL BETA
!       INTEGER K,LDA,LDB,LDC,N
!       CHARACTER TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHER2K  performs one of the hermitian rank 2k operations
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
!>          ALPHA is COMPLEX
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, ka ), where ka is
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
!>          B is COMPLEX array, dimension ( LDB, kb ), where kb is
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
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL
!>           On entry, BETA specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension ( LDC, N )
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
!>  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.
!>     Ed Anderson, Cray Research Inc.
!>
!>     converted to F90 and optimized 2023, Olivier Thomine
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX ALPHA
   REAL BETA
   INTEGER K,LDA,LDB,LDC,N
   CHARACTER TRANS,UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
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
   INTRINSIC CONJG,MAX,REAL
!     ..
!     .. Local Scalars ..
   COMPLEX TEMP1,TEMP2
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
       CALL XERBLA('CHER2K',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N == 0) .OR. (((ALPHA == (0.0E+0,0.0E+0)).OR. &
       (K == 0)).AND. (BETA == 1.0E+0))) RETURN
!
!     And when  alpha.eq.zero.
!
   IF (ALPHA == (0.0E+0,0.0E+0)) THEN
       IF (UPPER) THEN
           IF (BETA == 0.0E+0) THEN
               DO J = 1,N
                   C(1:J,J) = (0.0E+0,0.0E+0)
               ENDDO
           ELSE
               DO J = 1,N
                   C(1:J - 1,J) = BETA*C(1:J - 1,J)
                   C(J,J) = BETA*REAL(C(J,J))
               ENDDO
           END IF
       ELSE
           IF (BETA == 0.0E+0) THEN
               DO J = 1,N
                   C(J:N,J) = (0.0E+0,0.0E+0)
               ENDDO
           ELSE
               DO J = 1,N
                   C(J + 1:N,J) = BETA*C(J + 1:N,J)
                   C(J,J) = BETA*REAL(C(J,J))
               ENDDO
           END IF
       END IF
       RETURN
   END IF
!
!     Start the operations.
!
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
!                   C.
!
       IF (UPPER) THEN
           IF (BETA == 0.0E+0) THEN
               DO J = 1,N
                   C(1:J,J) = (0.0E+0,0.0E+0)
               ENDDO
           ELSE IF (BETA /= 1.0E+0) THEN
               DO J = 1,N
                   C(1:J - 1,J) = BETA*C(1:J - 1,J)
                   C(J,J) = BETA*REAL(C(J,J))
               ENDDO
           ELSE
               DO J = 1,N
                   C(J,J) = REAL(C(J,J))
               ENDDO
           ENDIF
           DO J = 1,N
               DO L = 1,K
                   IF ((A(J,L) /= (0.0E+0,0.0E+0)) .OR. (B(J,L) /= (0.0E+0,0.0E+0))) THEN
                       TEMP1 = ALPHA*CONJG(B(J,L))
                       TEMP2 = CONJG(ALPHA*A(J,L))
                       C(1:J - 1,J) = C(1:J - 1,J) + A(1:J - 1,L)*TEMP1 + B(1:J - 1,L)*TEMP2
                       C(J,J) = REAL(C(J,J)) + REAL(A(J,L)*TEMP1+B(J,L)*TEMP2)
                   END IF
               ENDDO
           ENDDO
       ELSE
           IF (BETA == 0.0E+0) THEN
               DO J = 1,N
                   C(J:N,J) = (0.0E+0,0.0E+0)
               ENDDO
           ELSE IF (BETA /= 1.0E+0) THEN
               DO J = 1,N
                   C(J + 1:N,J) = BETA*C(J + 1:N,J)
                   C(J,J) = BETA*REAL(C(J,J))
               ENDDO
           ELSE
               DO J = 1,N
                   C(J,J) = REAL(C(J,J))
               ENDDO
           ENDIF
           DO J = 1,N
               DO L = 1,K
                   IF ((A(J,L) /= (0.0E+0,0.0E+0)) .OR. (B(J,L) /= (0.0E+0,0.0E+0))) THEN
                       TEMP1 = ALPHA*CONJG(B(J,L))
                       TEMP2 = CONJG(ALPHA*A(J,L))
                       C(J + 1:N,J) = C(J + 1:N,J) + A(J + 1:N,L)*TEMP1 + B(J + 1:N,L)*TEMP2
                       C(J,J) = REAL(C(J,J)) + REAL(A(J,L)*TEMP1+B(J,L)*TEMP2)
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
           IF (BETA == 0.0E+0) THEN
               DO J = 1,N
                   DO I = 1, J - 1
                       TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,J))
                       TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,J))
                       C(I,J) = ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                   ENDDO
                   TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,I))
                   TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,I))
                   C(J,J) = REAL(ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2)
               ENDDO
           ELSEIF (BETA /= 1.0E+0) THEN
               DO J = 1,N
                   DO I = 1, J - 1
                       TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,J))
                       TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,J))
                       C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                   ENDDO
                   TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,I))
                   TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,I))
                   C(J,J) = BETA*REAL(C(J,J)) + REAL(ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2)
               ENDDO
           ELSE
               DO J = 1,N
                   DO I = 1, J - 1
                       TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,J))
                       TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,J))
                       C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                   ENDDO
                   TEMP1 = sum(CONJG(A(1:K,J))*B(1:K,J))
                   TEMP2 = sum(CONJG(B(1:K,J))*A(1:K,J))
                   C(J,J) = REAL(C(J,J)) + REAL(ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2)
               ENDDO
           ENDIF
       ELSE
           IF (BETA == 0.0E+0) THEN
               DO J = 1, N
                   TEMP1 = sum(CONJG(A(1:K,J))*B(1:K,J))
                   TEMP2 = sum(CONJG(B(1:K,J))*A(1:K,J))
                   C(J,J) = REAL(ALPHA*TEMP1+CONJG(ALPHA)*TEMP2)
                   DO I = J + 1,N
                       TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,J))
                       TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,J))
                       C(I,J) = ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                   ENDDO
               ENDDO
           ELSEIF (BETA /= 1.0E+0) THEN
               DO J = 1, N
                   TEMP1 = sum(CONJG(A(1:K,J))*B(1:K,J))
                   TEMP2 = sum(CONJG(B(1:K,J))*A(1:K,J))
                   C(J,J) = BETA*REAL(C(J,J)) + REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2)
                   DO I = J + 1, N
                       TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,J))
                       TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,J))
                       C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                   ENDDO
               ENDDO
           ELSE
               DO J = 1, N
                   TEMP1 = sum(CONJG(A(1:K,J))*B(1:K,J))
                   TEMP2 = sum(CONJG(B(1:K,J))*A(1:K,J))
                   C(J,J) = REAL(C(J,J)) + REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2)
                   DO I = J + 1, N
                       TEMP1 = sum(CONJG(A(1:K,I))*B(1:K,J))
                       TEMP2 = sum(CONJG(B(1:K,I))*A(1:K,J))
                       C(I,J) = C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                   ENDDO
               ENDDO
           ENDIF
       END IF
   END IF
!
   RETURN
!
!     End of CHER2K
!
END
