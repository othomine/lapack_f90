!> \brief \b ZLA_GEAMV computes a matrix-vector product using a general matrix to calculate error bounds.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLA_GEAMV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_geamv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_geamv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_geamv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLA_GEAMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA,
!                             Y, INCY )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   ALPHA, BETA
!       INTEGER            INCX, INCY, LDA, M, N
!       INTEGER            TRANS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), X( * )
!       DOUBLE PRECISION   Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLA_GEAMV  performs one of the matrix-vector operations
!>
!>         y := alpha*abs(A)*abs(x) + beta*abs(y),
!>    or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
!>
!> This function is primarily used in calculating error bounds.
!> To protect against underflow during evaluation, components in
!> the resulting vector are perturbed away from zero by (N+1)
!> times the underflow threshold.  To prevent unnecessarily large
!> errors for block-structure embedded in general matrices,
!> "symbolically" zero components are not perturbed.  A zero
!> entry is considered "symbolic" if all multiplications involved
!> in computing that entry have at least one zero multiplicand.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is INTEGER
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)
!>             BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y)
!>             BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y)
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, n )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!>           Before entry, the incremented array X must contain the
!>           vector x.
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
!>          BETA is DOUBLE PRECISION
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry with BETA non-zero, the incremented array Y
!>           must contain the vector y. On exit, Y is overwritten by the
!>           updated vector y.
!>           If either m or n is zero, then Y not referenced and the function
!>           performs a quick return.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!>           Unchanged on exit.
!>
!>  Level 2 Blas routine.
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
!> \ingroup la_geamv
!
!  =====================================================================
   SUBROUTINE ZLA_GEAMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, &
                         Y, INCY )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION   ALPHA, BETA
   INTEGER            INCX, INCY, LDA, M, N
   INTEGER            TRANS
!     ..
!     .. Array Arguments ..
   COMPLEX*16         A( LDA, * ), X( * )
   DOUBLE PRECISION   Y( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   COMPLEX*16         ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            SYMB_ZERO
   DOUBLE PRECISION   TEMP, SAFE1
   INTEGER            I, INFO, IY, J, JX, KX, KY, LENX, LENY
   COMPLEX*16         CDUM
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, DLAMCH
   DOUBLE PRECISION   DLAMCH
!     ..
!     .. External Functions ..
   EXTERNAL           ILATRANS
   INTEGER            ILATRANS
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, ABS, REAL, DIMAG, SIGN
!     ..
!     .. Statement Functions ..
   DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function Definitions ..
   CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   IF     ( .NOT.( ( TRANS == ILATRANS( 'N' ) ) &
              .OR. ( TRANS == ILATRANS( 'T' ) ) &
              .OR. ( TRANS == ILATRANS( 'C' ) ) ) ) THEN
      INFO = 1
   ELSE IF( M < 0 )THEN
      INFO = 2
   ELSE IF( N < 0 )THEN
      INFO = 3
   ELSE IF( LDA < MAX( 1, M ) )THEN
      INFO = 6
   ELSE IF( INCX == 0 )THEN
      INFO = 8
   ELSE IF( INCY == 0 )THEN
      INFO = 11
   END IF
   IF( INFO /= 0 )THEN
      CALL XERBLA( 'ZLA_GEAMV ', INFO )
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF( ( M == 0 ).OR.( N == 0 ).OR. &
       ( ( ALPHA == ZERO ).AND.( BETA == ONE ) ) ) &
      RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
   IF( TRANS == ILATRANS( 'N' ) )THEN
      LENX = N
      LENY = M
   ELSE
      LENX = M
      LENY = N
   END IF
   IF( INCX > 0 )THEN
      KX = 1
   ELSE
      KX = 1 - ( LENX - 1 )*INCX
   END IF
   IF( INCY > 0 )THEN
      KY = 1
   ELSE
      KY = 1 - ( LENY - 1 )*INCY
   END IF
!
!     Set SAFE1 essentially to be the underflow threshold times the
!     number of additions in each row.
!
   SAFE1 = DLAMCH( 'Safe minimum' )
   SAFE1 = (N+1)*SAFE1
!
!     Form  y := alpha*abs(A)*abs(x) + beta*abs(y).
!
!     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
!     the inexact flag.  Still doesn't help change the iteration order
!     to per-column.
!
   IY = KY
   IF ( INCX == 1 ) THEN
      IF( TRANS == ILATRANS( 'N' ) )THEN
         DO I = 1, LENY
            IF ( BETA  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
               Y( IY ) = 0.0D+0
            ELSE IF ( Y( IY )  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
            ELSE
               SYMB_ZERO = .FALSE.
               Y( IY ) = BETA * ABS( Y( IY ) )
            END IF
            IF ( ALPHA  /=  0.0D+0 ) THEN
               DO J = 1, LENX
                  TEMP = CABS1( A( I, J ) )
                  SYMB_ZERO = SYMB_ZERO .AND. &
                       ( X( J )  ==  ZERO .OR. TEMP  ==  ZERO )

                  Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
               END DO
            END IF

            IF ( .NOT.SYMB_ZERO ) Y( IY ) = &
                 Y( IY ) + SIGN( SAFE1, Y( IY ) )

            IY = IY + INCY
         END DO
      ELSE
         DO I = 1, LENY
            IF ( BETA  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
               Y( IY ) = 0.0D+0
            ELSE IF ( Y( IY )  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
            ELSE
               SYMB_ZERO = .FALSE.
               Y( IY ) = BETA * ABS( Y( IY ) )
            END IF
            IF ( ALPHA  /=  0.0D+0 ) THEN
               DO J = 1, LENX
                  TEMP = CABS1( A( J, I ) )
                  SYMB_ZERO = SYMB_ZERO .AND. &
                       ( X( J )  ==  ZERO .OR. TEMP  ==  ZERO )

                  Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
               END DO
            END IF

            IF ( .NOT.SYMB_ZERO ) Y( IY ) = &
                 Y( IY ) + SIGN( SAFE1, Y( IY ) )

            IY = IY + INCY
         END DO
      END IF
   ELSE
      IF( TRANS == ILATRANS( 'N' ) )THEN
         DO I = 1, LENY
            IF ( BETA  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
               Y( IY ) = 0.0D+0
            ELSE IF ( Y( IY )  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
            ELSE
               SYMB_ZERO = .FALSE.
               Y( IY ) = BETA * ABS( Y( IY ) )
            END IF
            IF ( ALPHA  /=  0.0D+0 ) THEN
               JX = KX
               DO J = 1, LENX
                  TEMP = CABS1( A( I, J ) )
                  SYMB_ZERO = SYMB_ZERO .AND. &
                       ( X( JX )  ==  ZERO .OR. TEMP  ==  ZERO )

                  Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                  JX = JX + INCX
               END DO
            END IF

            IF ( .NOT.SYMB_ZERO ) Y( IY ) = &
                 Y( IY ) + SIGN( SAFE1, Y( IY ) )

            IY = IY + INCY
         END DO
      ELSE
         DO I = 1, LENY
            IF ( BETA  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
               Y( IY ) = 0.0D+0
            ELSE IF ( Y( IY )  ==  0.0D+0 ) THEN
               SYMB_ZERO = .TRUE.
            ELSE
               SYMB_ZERO = .FALSE.
               Y( IY ) = BETA * ABS( Y( IY ) )
            END IF
            IF ( ALPHA  /=  0.0D+0 ) THEN
               JX = KX
               DO J = 1, LENX
                  TEMP = CABS1( A( J, I ) )
                  SYMB_ZERO = SYMB_ZERO .AND. &
                       ( X( JX )  ==  ZERO .OR. TEMP  ==  ZERO )

                  Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                  JX = JX + INCX
               END DO
            END IF

            IF ( .NOT.SYMB_ZERO ) Y( IY ) = &
                 Y( IY ) + SIGN( SAFE1, Y( IY ) )

            IY = IY + INCY
         END DO
      END IF

   END IF
!
   RETURN
!
!     End of ZLA_GEAMV
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

