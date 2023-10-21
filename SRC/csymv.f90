!> \brief \b CSYMV computes a matrix-vector product for a complex symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYMV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csymv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csymv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csymv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCX, INCY, LDA, N
!       COMPLEX            ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYMV  performs the matrix-vector  operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix.
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           Before entry, with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular part of the symmetric matrix and the strictly
!>           lower triangular part of A is not referenced.
!>           Before entry, with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular part of the symmetric matrix and the strictly
!>           upper triangular part of A is not referenced.
!>           Unchanged on exit.
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
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y. On exit, Y is overwritten by the updated
!>           vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
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
!> \ingroup hemv
!
!  =====================================================================
   SUBROUTINE CSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INCX, INCY, LDA, N
   COMPLEX            ALPHA, BETA
!     ..
!     .. Array Arguments ..
   COMPLEX            A( LDA, * ), X( * ), Y( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
   COMPLEX            ONE
   PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
   COMPLEX            ZERO
   PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
   COMPLEX            TEMP1, TEMP2
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX
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
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = 5
   ELSE IF( INCX == 0 ) THEN
      INFO = 7
   ELSE IF( INCY == 0 ) THEN
      INFO = 10
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSYMV ', INFO )
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF( ( N == 0 ) .OR. ( ( ALPHA == ZERO ) .AND. ( BETA == ONE ) ) ) &
      RETURN
!
!     Set up the start points in  X  and  Y.
!
   IF( INCX > 0 ) THEN
      KX = 1
   ELSE
      KX = 1 - ( N-1 )*INCX
   END IF
   IF( INCY > 0 ) THEN
      KY = 1
   ELSE
      KY = 1 - ( N-1 )*INCY
   END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
!     First form  y := beta*y.
!
   IF( BETA /= ONE ) THEN
      IF( INCY == 1 ) THEN
         IF( BETA == ZERO ) THEN
            DO I = 1, N
               Y( I ) = ZERO
            ENDDO
         ELSE
            DO I = 1, N
               Y( I ) = BETA*Y( I )
            ENDDO
         END IF
      ELSE
         IY = KY
         IF( BETA == ZERO ) THEN
            DO I = 1, N
               Y( IY ) = ZERO
               IY = IY + INCY
            ENDDO
         ELSE
            DO I = 1, N
               Y( IY ) = BETA*Y( IY )
               IY = IY + INCY
            ENDDO
         END IF
      END IF
   END IF
   IF( ALPHA == ZERO ) &
      RETURN
   IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Form  y  when A is stored in upper triangle.
!
      IF( ( INCX == 1 ) .AND. ( INCY == 1 ) ) THEN
         DO J = 1, N
            TEMP1 = ALPHA*X( J )
            TEMP2 = ZERO
            DO I = 1, J - 1
               Y( I ) = Y( I ) + TEMP1*A( I, J )
               TEMP2 = TEMP2 + A( I, J )*X( I )
            ENDDO
            Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
         ENDDO
      ELSE
         JX = KX
         JY = KY
         DO J = 1, N
            TEMP1 = ALPHA*X( JX )
            TEMP2 = ZERO
            IX = KX
            IY = KY
            DO I = 1, J - 1
               Y( IY ) = Y( IY ) + TEMP1*A( I, J )
               TEMP2 = TEMP2 + A( I, J )*X( IX )
               IX = IX + INCX
               IY = IY + INCY
            ENDDO
            Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
            JX = JX + INCX
            JY = JY + INCY
         ENDDO
      END IF
   ELSE
!
!        Form  y  when A is stored in lower triangle.
!
      IF( ( INCX == 1 ) .AND. ( INCY == 1 ) ) THEN
         DO J = 1, N
            TEMP1 = ALPHA*X( J )
            TEMP2 = ZERO
            Y( J ) = Y( J ) + TEMP1*A( J, J )
            DO I = J + 1, N
               Y( I ) = Y( I ) + TEMP1*A( I, J )
               TEMP2 = TEMP2 + A( I, J )*X( I )
            ENDDO
            Y( J ) = Y( J ) + ALPHA*TEMP2
            ENDDO
      ELSE
         JX = KX
         JY = KY
         DO J = 1, N
            TEMP1 = ALPHA*X( JX )
            TEMP2 = ZERO
            Y( JY ) = Y( JY ) + TEMP1*A( J, J )
            IX = JX
            IY = JY
            DO I = J + 1, N
               IX = IX + INCX
               IY = IY + INCY
               Y( IY ) = Y( IY ) + TEMP1*A( I, J )
               TEMP2 = TEMP2 + A( I, J )*X( IX )
               ENDDO
            Y( JY ) = Y( JY ) + ALPHA*TEMP2
            JX = JX + INCX
            JY = JY + INCY
            ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of CSYMV
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

