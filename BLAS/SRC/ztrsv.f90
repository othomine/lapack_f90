!> \brief \b ZTRSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRSV  solves one of the systems of equations
!>
!>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
!>
!> where b and x are n element vectors and A is an n by n unit, or
!> non-unit, upper or lower triangular matrix.
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular matrix and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular matrix and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!>           A are not referenced either, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
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
!> \ingroup trsv
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
   SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER INCX,LDA,N
   CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
   COMPLEX*16 A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   COMPLEX*16 TEMP
   INTEGER I,INFO,IX,J,JX,KX
   LOGICAL NOCONJ,NOUNIT
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
   ELSE IF (.NOT.LSAME(TRANS,'N') .AND. &
            .NOT.LSAME(TRANS,'T') .AND. &
            .NOT.LSAME(TRANS,'C')) THEN
       INFO = 2
   ELSE IF (.NOT.LSAME(DIAG,'U') .AND. &
            .NOT.LSAME(DIAG,'N')) THEN
       INFO = 3
   ELSE IF (N < 0) THEN
       INFO = 4
   ELSE IF (LDA < MAX(1,N)) THEN
       INFO = 6
   ELSE IF (INCX == 0) THEN
       INFO = 8
   END IF
   IF (INFO /= 0) THEN
       CALL XERBLA('ZTRSV ',INFO)
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
!     accessed sequentially with one pass through A.
!
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := inv( A )*x.
!
       IF (LSAME(UPLO,'U')) THEN
           IF (INCX == 1) THEN
               DO J = N,1,-1
                   IF (X(J) /= (0.0D+0,0.0D+0)) THEN
                       IF (NOUNIT) X(J) = X(J)/A(J,J)
                       X(1:J-1) = X(1:J-1) - X(J)*A(1:J-1,J)
                   END IF
               ENDDO
           ELSE
               JX = KX + (N-1)*INCX
               DO J = N,1,-1
                   IF (X(JX) /= (0.0D+0,0.0D+0)) THEN
                       IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                       TEMP = X(JX)
                       IX = JX
                       DO I = J - 1,1,-1
                           IX = IX - INCX
                           X(IX) = X(IX) - TEMP*A(I,J)
                       ENDDO
                   END IF
                   JX = JX - INCX
               ENDDO
           END IF
       ELSE
           IF (INCX == 1) THEN
               DO J = 1,N
                   IF (X(J) /= (0.0D+0,0.0D+0)) THEN
                       IF (NOUNIT) X(J) = X(J)/A(J,J)
                       X(J+1:N) = X(J+1:N) - X(J)*A(J+1:N,J)
                   END IF
               ENDDO
           ELSE
               JX = KX
               DO J = 1,N
                   IF (X(JX) /= (0.0D+0,0.0D+0)) THEN
                       IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                       TEMP = X(JX)
                       IX = JX
                       DO I = J + 1,N
                           IX = IX + INCX
                           X(IX) = X(IX) - TEMP*A(I,J)
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
                   TEMP = X(J)
                   IF (NOCONJ) THEN
                       TEMP = TEMP - sum(A(1:J-1,J)*X(1:J-1))
                       IF (NOUNIT) TEMP = TEMP/A(J,J)
                   ELSE
                       TEMP = TEMP - sum(DCONJG(A(1:J-1,J))*X(1:J-1))
                       IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                   END IF
                   X(J) = TEMP
               ENDDO
           ELSE
               JX = KX
               DO J = 1,N
                   IX = KX
                   TEMP = X(JX)
                   IF (NOCONJ) THEN
                       DO I = 1,J - 1
                           TEMP = TEMP - A(I,J)*X(IX)
                           IX = IX + INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/A(J,J)
                   ELSE
                       DO I = 1,J - 1
                           TEMP = TEMP - DCONJG(A(I,J))*X(IX)
                           IX = IX + INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                   END IF
                   X(JX) = TEMP
                   JX = JX + INCX
               ENDDO
           END IF
       ELSE
           IF (INCX == 1) THEN
               DO J = N,1,-1
                   TEMP = X(J)
                   IF (NOCONJ) THEN
                       TEMP = TEMP - sum(A(J+1:N,J)*X(J+1:N))
                       IF (NOUNIT) TEMP = TEMP/A(J,J)
                   ELSE
                       TEMP = TEMP - sum(DCONJG(A(J+1:N,J))*X(J+1:N))
                       IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                   END IF
                   X(J) = TEMP
               ENDDO
           ELSE
               KX = KX + (N-1)*INCX
               JX = KX
               DO J = N,1,-1
                   IX = KX
                   TEMP = X(JX)
                   IF (NOCONJ) THEN
                       DO I = N,J + 1,-1
                           TEMP = TEMP - A(I,J)*X(IX)
                           IX = IX - INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/A(J,J)
                   ELSE
                       DO I = N,J + 1,-1
                           TEMP = TEMP - DCONJG(A(I,J))*X(IX)
                           IX = IX - INCX
                       ENDDO
                       IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                   END IF
                   X(JX) = TEMP
                   JX = JX - INCX
               ENDDO
           END IF
       END IF
   END IF
!
   RETURN
!
!     End of ZTRSV
!
END
