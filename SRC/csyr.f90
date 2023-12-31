!> \brief \b CSYR performs the symmetric rank-1 update of a complex symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCX, LDA, N
!       COMPLEX            ALPHA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYR   performs the symmetric rank 1 operation
!>
!>    A := alpha*x*x**H + A,
!>
!> where alpha is a complex scalar, x is an n element vector and A is an
!> n by n symmetric matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the array A is to be referenced as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of A
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of A
!>                                  is to be referenced.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the N-
!>           element vector x.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           Before entry, with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular part of the symmetric matrix and the strictly
!>           lower triangular part of A is not referenced. On exit, the
!>           upper triangular part of the array A is overwritten by the
!>           upper triangular part of the updated matrix.
!>           Before entry, with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular part of the symmetric matrix and the strictly
!>           upper triangular part of A is not referenced. On exit, the
!>           lower triangular part of the array A is overwritten by the
!>           lower triangular part of the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, N ).
!>           Unchanged on exit.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup her
!
!  =====================================================================
   SUBROUTINE CSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INCX, LDA, N
   COMPLEX            ALPHA
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), X( * )
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, IX, J, JX, KX
   COMPLEX            TEMP
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = 1
   ELSE IF( N < 0 ) THEN
      INFO = 2
   ELSE IF( INCX == 0 ) THEN
      INFO = 5
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = 7
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSYR  ', INFO )
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF( ( N == 0 ) .OR. ( ALPHA == (0.0E+0,0.0E+0) ) ) RETURN
!
!     Set the start point in X if the increment is not unity.
!
   IF( INCX <= 0 ) THEN
      KX = 1 - ( N-1 )*INCX
   ELSE IF( INCX /= 1 ) THEN
      KX = 1
   END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
   IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Form  A  when A is stored in upper triangle.
!
      IF( INCX == 1 ) THEN
         DO J = 1, N
            IF( X( J ) /= (0.0E+0,0.0E+0) ) THEN
               TEMP = ALPHA*X( J )
               A(1:J,J) = A(1:J,J) + X(1:J)*TEMP
            END IF
         ENDDO
      ELSE
         JX = KX
         DO J = 1, N
            IF( X( JX ) /= (0.0E+0,0.0E+0) ) THEN
               TEMP = ALPHA*X( JX )
               IX = KX
               DO I = 1, J
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX = IX + INCX
               ENDDO
            END IF
            JX = JX + INCX
         ENDDO
      END IF
   ELSE
!
!        Form  A  when A is stored in lower triangle.
!
      IF( INCX == 1 ) THEN
         DO J = 1, N
            IF( X( J ) /= (0.0E+0,0.0E+0) ) THEN
               TEMP = ALPHA*X( J )
               A(J:N,J) = A(J:N,J) + X(J:N)*TEMP
            END IF
         ENDDO
      ELSE
         JX = KX
         DO J = 1, N
            IF( X( JX ) /= (0.0E+0,0.0E+0) ) THEN
               TEMP = ALPHA*X( JX )
               IX = JX
               DO I = J, N
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX = IX + INCX
               ENDDO
            END IF
            JX = JX + INCX
         ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of CSYR
!
END
