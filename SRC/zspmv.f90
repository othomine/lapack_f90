!> \brief \b ZSPMV computes a matrix-vector product for complex vectors using a complex symmetric packed matrix
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSPMV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspmv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspmv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspmv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCX, INCY, N
!       COMPLEX*16         ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         AP( * ), X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSPMV  performs the matrix-vector operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix, supplied in packed form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the matrix A is supplied in the packed
!>           array AP as follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  supplied in AP.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  supplied in AP.
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
!>          ALPHA is COMPLEX*16
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension at least
!>           ( ( N*( N + 1 ) )/2 ).
!>           Before entry, with UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on.
!>           Before entry, with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
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
!>          BETA is COMPLEX*16
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension at least
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
!
!> \ingroup hpmv
!
!  =====================================================================
   SUBROUTINE ZSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INCX, INCY, N
   COMPLEX*16         ALPHA, BETA
!     ..
!     .. Array Arguments ..
   COMPLEX*16         AP( * ), X( * ), Y( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
   COMPLEX*16         ONE
   PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
   COMPLEX*16         ZERO
   PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
   COMPLEX*16         TEMP1, TEMP2
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
      INFO = 6
   ELSE IF( INCY == 0 ) THEN
      INFO = 9
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'ZSPMV ', INFO )
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
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
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
   KK = 1
   IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Form  y  when AP contains the upper triangle.
!
      IF( ( INCX == 1 ) .AND. ( INCY == 1 ) ) THEN
         DO J = 1, N
            TEMP1 = ALPHA*X( J )
            TEMP2 = ZERO
            K = KK
            DO I = 1, J - 1
               Y( I ) = Y( I ) + TEMP1*AP( K )
               TEMP2 = TEMP2 + AP( K )*X( I )
               K = K + 1
            ENDDO
            Y( J ) = Y( J ) + TEMP1*AP( KK+J-1 ) + ALPHA*TEMP2
            KK = KK + J
         ENDDO
      ELSE
         JX = KX
         JY = KY
         DO J = 1, N
            TEMP1 = ALPHA*X( JX )
            TEMP2 = ZERO
            IX = KX
            IY = KY
            DO K = KK, KK + J - 2
               Y( IY ) = Y( IY ) + TEMP1*AP( K )
               TEMP2 = TEMP2 + AP( K )*X( IX )
               IX = IX + INCX
               IY = IY + INCY
            ENDDO
            Y( JY ) = Y( JY ) + TEMP1*AP( KK+J-1 ) + ALPHA*TEMP2
            JX = JX + INCX
            JY = JY + INCY
            KK = KK + J
         ENDDO
      END IF
   ELSE
!
!        Form  y  when AP contains the lower triangle.
!
      IF( ( INCX == 1 ) .AND. ( INCY == 1 ) ) THEN
         DO J = 1, N
            TEMP1 = ALPHA*X( J )
            TEMP2 = ZERO
            Y( J ) = Y( J ) + TEMP1*AP( KK )
            K = KK + 1
            DO I = J + 1, N
               Y( I ) = Y( I ) + TEMP1*AP( K )
               TEMP2 = TEMP2 + AP( K )*X( I )
               K = K + 1
            ENDDO
            Y( J ) = Y( J ) + ALPHA*TEMP2
            KK = KK + ( N-J+1 )
            ENDDO
      ELSE
         JX = KX
         JY = KY
         DO J = 1, N
            TEMP1 = ALPHA*X( JX )
            TEMP2 = ZERO
            Y( JY ) = Y( JY ) + TEMP1*AP( KK )
            IX = JX
            IY = JY
            DO K = KK + 1, KK + N - J
               IX = IX + INCX
               IY = IY + INCY
               Y( IY ) = Y( IY ) + TEMP1*AP( K )
               TEMP2 = TEMP2 + AP( K )*X( IX )
               ENDDO
            Y( JY ) = Y( JY ) + ALPHA*TEMP2
            JX = JX + INCX
            JY = JY + INCY
            KK = KK + ( N-J+1 )
            ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of ZSPMV
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
