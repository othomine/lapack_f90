!> \brief \b ZGEMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
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
!> ZGEMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*op( A )*op( B ) + beta*C,
!>
!> where  op( X ) is one of
!>
!>    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
!>
!> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n',  op( A ) = A.
!>
!>              TRANSA = 'T' or 't',  op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c',  op( A ) = A**H.
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER*1
!>           On entry, TRANSB specifies the form of op( B ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSB = 'N' or 'n',  op( B ) = B.
!>
!>              TRANSB = 'T' or 't',  op( B ) = B**T.
!>
!>              TRANSB = 'C' or 'c',  op( B ) = B**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies  the number  of rows  of the  matrix
!>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N  specifies the number  of columns of the matrix
!>           op( B ) and the number of columns of the matrix C. N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry,  K  specifies  the number of columns of the matrix
!>           op( A ) and the number of rows of the matrix op( B ). K must
!>           be at least  zero.
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
!>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by m  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension ( LDB, kb ), where kb is
!>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  n by k  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!>           least  max( 1, n ).
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
!>           On exit, the array  C  is overwritten by the  m by n  matrix
!>           ( alpha*op( A )*op( B ) + beta*C ).
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
!> \ingroup gemm
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
   SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16 ALPHA,BETA
   INTEGER K,LDA,LDB,LDC,M,N
   CHARACTER TRANSA,TRANSB
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
   INTEGER INFO,J,L,NROWA,NROWB
   LOGICAL CONJA,CONJB,NOTA,NOTB
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!     B  respectively are to be  transposed but  not conjugated  and set
!     NROWA and NROWB  as the number of rows  of  A  and  B  respectively.
!
   NOTA = LSAME(TRANSA,'N')
   NOTB = LSAME(TRANSB,'N')
   CONJA = LSAME(TRANSA,'C')
   CONJB = LSAME(TRANSB,'C')
   IF (NOTA) THEN
       NROWA = M
   ELSE
       NROWA = K
   END IF
   IF (NOTB) THEN
       NROWB = K
   ELSE
       NROWB = N
   END IF
!
!     Test the input parameters.
!
   INFO = 0
   IF ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND. &
       (.NOT.LSAME(TRANSA,'T'))) THEN
       INFO = 1
   ELSE IF ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND. &
            (.NOT.LSAME(TRANSB,'T'))) THEN
       INFO = 2
   ELSE IF (M < 0) THEN
       INFO = 3
   ELSE IF (N < 0) THEN
       INFO = 4
   ELSE IF (K < 0) THEN
       INFO = 5
   ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 8
   ELSE IF (LDB < MAX(1,NROWB)) THEN
       INFO = 10
   ELSE IF (LDC < MAX(1,M)) THEN
       INFO = 13
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('ZGEMM ',INFO)
       RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((M == 0) .OR. (N == 0) .OR. &
       (((ALPHA == (0.0D+0,0.0D+0)).OR. (K == 0)).AND. (BETA == (1.0D+0,0.0D+0)))) RETURN
!
!     And when  alpha.eq.zero.
!
!
!     Start the operations.
!
   IF (BETA == (0.0D+0,0.0D+0)) THEN
       C(1:M,1:N) = (0.0D+0,0.0D+0)
   ELSE IF (BETA /= (1.0D+0,0.0D+0)) THEN
       C(1:M,1:N) = BETA*C(1:M,1:N)
   END IF
   IF (ALPHA == (0.0D+0,0.0D+0)) RETURN
   IF (NOTB) THEN
       IF (NOTA) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*B(L,J)*A(1:M,L)
               ENDDO
           ENDDO
       ELSE IF (CONJA) THEN
!
!           Form  C := alpha*A**H*B + beta*C.
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*DCONJG(A(L,1:M))*B(L,J)
               ENDDO
           ENDDO
       ELSE
!
!           Form  C := alpha*A**T*B + beta*C
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*A(L,1:M)*B(L,J)
               ENDDO
           ENDDO
       END IF
   ELSE IF (NOTA) THEN
       IF (CONJB) THEN
!
!           Form  C := alpha*A*B**H + beta*C.
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*DCONJG(B(J,L))*A(1:M,L)
               ENDDO
           ENDDO
       ELSE
!
!           Form  C := alpha*A*B**T + beta*C
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*B(J,L)*A(1:M,L)
               ENDDO
           ENDDO
       END IF
   ELSE IF (CONJA) THEN
       IF (CONJB) THEN
!
!           Form  C := alpha*A**H*B**H + beta*C.
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*DCONJG(A(L,1:M))*DCONJG(B(J,L))
               ENDDO
           ENDDO
       ELSE
!
!           Form  C := alpha*A**H*B**T + beta*C
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*DCONJG(A(L,1:M))*B(J,L)
               ENDDO
           ENDDO
       END IF
   ELSE
       IF (CONJB) THEN
!
!           Form  C := alpha*A**T*B**H + beta*C
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*A(L,1:M)*DCONJG(B(J,L))
               ENDDO
           ENDDO
       ELSE
!
!           Form  C := alpha*A**T*B**T + beta*C
!
           DO J = 1,N
               DO L = 1,K
                   C(1:M,J) = C(1:M,J) + ALPHA*A(L,1:M)*B(J,L)
               ENDDO
           ENDDO
       END IF
   END IF
!
   RETURN
!
!     End of ZGEMM
!
END
