!> \brief \b DORBDB6
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORBDB6 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb6.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb6.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb6.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,
!                           LDQ2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,
!      $                   N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!>\verbatim
!>
!> DORBDB6 orthogonalizes the column vector
!>      X = [ X1 ]
!>          [ X2 ]
!> with respect to the columns of
!>      Q = [ Q1 ] .
!>          [ Q2 ]
!> The Euclidean norm of X must be one and the columns of Q must be
!> orthonormal. The orthogonalized vector will be zero if and only if it
!> lies entirely in the range of Q.
!>
!> The projection is computed with at most two iterations of the
!> classical Gram-Schmidt algorithm, see
!> * L. Giraud, J. Langou, M. Rozložník. "On the round-off error
!>   analysis of the Gram-Schmidt algorithm with reorthogonalization."
!>   2002. CERFACS Technical Report No. TR/PA/02/33. URL:
!>   https://www.cerfacs.fr/algor/reports/2002/TR_PA_02_33.pdf
!>
!>\endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M1
!> \verbatim
!>          M1 is INTEGER
!>           The dimension of X1 and the number of rows in Q1. 0 <= M1.
!> \endverbatim
!>
!> \param[in] M2
!> \verbatim
!>          M2 is INTEGER
!>           The dimension of X2 and the number of rows in Q2. 0 <= M2.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The number of columns in Q1 and Q2. 0 <= N.
!> \endverbatim
!>
!> \param[in,out] X1
!> \verbatim
!>          X1 is DOUBLE PRECISION array, dimension (M1)
!>           On entry, the top part of the vector to be orthogonalized.
!>           On exit, the top part of the projected vector.
!> \endverbatim
!>
!> \param[in] INCX1
!> \verbatim
!>          INCX1 is INTEGER
!>           Increment for entries of X1.
!> \endverbatim
!>
!> \param[in,out] X2
!> \verbatim
!>          X2 is DOUBLE PRECISION array, dimension (M2)
!>           On entry, the bottom part of the vector to be
!>           orthogonalized. On exit, the bottom part of the projected
!>           vector.
!> \endverbatim
!>
!> \param[in] INCX2
!> \verbatim
!>          INCX2 is INTEGER
!>           Increment for entries of X2.
!> \endverbatim
!>
!> \param[in] Q1
!> \verbatim
!>          Q1 is DOUBLE PRECISION array, dimension (LDQ1, N)
!>           The top part of the orthonormal basis matrix.
!> \endverbatim
!>
!> \param[in] LDQ1
!> \verbatim
!>          LDQ1 is INTEGER
!>           The leading dimension of Q1. LDQ1 >= M1.
!> \endverbatim
!>
!> \param[in] Q2
!> \verbatim
!>          Q2 is DOUBLE PRECISION array, dimension (LDQ2, N)
!>           The bottom part of the orthonormal basis matrix.
!> \endverbatim
!>
!> \param[in] LDQ2
!> \verbatim
!>          LDQ2 is INTEGER
!>           The leading dimension of Q2. LDQ2 >= M2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK. LWORK >= N.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit.
!>           < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup unbdb6
!
!  =====================================================================
   SUBROUTINE DORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, &
                       LDQ2, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, &
                      N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ALPHA, REALONE, REALZERO
   PARAMETER          ( ALPHA = 0.01D0, REALONE = 1.0D0, &
                        REALZERO = 0.0D0 )
   DOUBLE PRECISION   NEGONE, ONE, ZERO
   PARAMETER          ( NEGONE = -1.0D0, ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IX
   DOUBLE PRECISION   EPS, NORM, NORM_NEW, SCL, SSQ
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMV, DLASSQ, XERBLA
!     ..
!     .. Intrinsic Function ..
   INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!
   INFO = 0
   IF( M1  <  0 ) THEN
      INFO = -1
   ELSE IF( M2  <  0 ) THEN
      INFO = -2
   ELSE IF( N  <  0 ) THEN
      INFO = -3
   ELSE IF( INCX1  <  1 ) THEN
      INFO = -5
   ELSE IF( INCX2  <  1 ) THEN
      INFO = -7
   ELSE IF( LDQ1  <  MAX( 1, M1 ) ) THEN
      INFO = -9
   ELSE IF( LDQ2  <  MAX( 1, M2 ) ) THEN
      INFO = -11
   ELSE IF( LWORK  <  N ) THEN
      INFO = -13
   END IF
!
   IF( INFO  /=  0 ) THEN
      CALL XERBLA( 'DORBDB6', -INFO )
      RETURN
   END IF
!
   EPS = DLAMCH( 'Precision' )
!
!     First, project X onto the orthogonal complement of Q's column
!     space
!
!     Christoph Conrads: In debugging mode the norm should be computed
!     and an assertion added comparing the norm with one. Alas, Fortran
!     never made it into 1989 when assert() was introduced into the C
!     programming language.
   NORM = REALONE
!
   IF( M1  ==  0 ) THEN
      DO I = 1, N
         WORK(I) = ZERO
      END DO
   ELSE
      CALL DGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, &
                  1 )
   END IF
!
   CALL DGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )
!
   CALL DGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, &
               INCX1 )
   CALL DGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, &
               INCX2 )
!
   SCL = REALZERO
   SSQ = REALZERO
   CALL DLASSQ( M1, X1, INCX1, SCL, SSQ )
   CALL DLASSQ( M2, X2, INCX2, SCL, SSQ )
   NORM_NEW = SCL * SQRT(SSQ)
!
!     If projection is sufficiently large in norm, then stop.
!     If projection is zero, then stop.
!     Otherwise, project again.
!
   IF( NORM_NEW  >=  ALPHA * NORM ) THEN
      RETURN
   END IF
!
   IF( NORM_NEW  <=  N * EPS * NORM ) THEN
      DO IX = 1, 1 + (M1-1)*INCX1, INCX1
        X1( IX ) = ZERO
      END DO
      DO IX = 1, 1 + (M2-1)*INCX2, INCX2
        X2( IX ) = ZERO
      END DO
      RETURN
   END IF
!
   NORM = NORM_NEW
!
   DO I = 1, N
      WORK(I) = ZERO
   END DO
!
   IF( M1  ==  0 ) THEN
      DO I = 1, N
         WORK(I) = ZERO
      END DO
   ELSE
      CALL DGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, &
                  1 )
   END IF
!
   CALL DGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )
!
   CALL DGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, &
               INCX1 )
   CALL DGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, &
               INCX2 )
!
   SCL = REALZERO
   SSQ = REALZERO
   CALL DLASSQ( M1, X1, INCX1, SCL, SSQ )
   CALL DLASSQ( M2, X2, INCX2, SCL, SSQ )
   NORM_NEW = SCL * SQRT(SSQ)
!
!     If second projection is sufficiently large in norm, then do
!     nothing more. Alternatively, if it shrunk significantly, then
!     truncate it to zero.
!
   IF( NORM_NEW  <  ALPHA * NORM ) THEN
      DO IX = 1, 1 + (M1-1)*INCX1, INCX1
         X1(IX) = ZERO
      END DO
      DO IX = 1, 1 + (M2-1)*INCX2, INCX2
         X2(IX) = ZERO
      END DO
   END IF
!
   RETURN
!
!     End of DORBDB6
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
