!> \brief \b SGEEQU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEEQU + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeequ.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeequ.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeequ.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       REAL               AMAX, COLCND, ROWCND
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), C( * ), R( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEEQU computes row and column scalings intended to equilibrate an
!> M-by-N matrix A and reduce its condition number.  R returns the row
!> scale factors and C the column scale factors, chosen to try to make
!> the largest element in each row and column of the matrix B with
!> elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
!>
!> R(i) and C(j) are restricted to be between SMLNUM = smallest safe
!> number and BIGNUM = largest safe number.  Use of these scaling
!> factors is not guaranteed to reduce the condition number of A but
!> works well in practice.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The M-by-N matrix whose equilibration factors are
!>          to be computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is REAL array, dimension (M)
!>          If INFO = 0 or INFO > M, R contains the row scale factors
!>          for A.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL array, dimension (N)
!>          If INFO = 0,  C contains the column scale factors for A.
!> \endverbatim
!>
!> \param[out] ROWCND
!> \verbatim
!>          ROWCND is REAL
!>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
!>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
!>          AMAX is neither too large nor too small, it is not worth
!>          scaling by R.
!> \endverbatim
!>
!> \param[out] COLCND
!> \verbatim
!>          COLCND is REAL
!>          If INFO = 0, COLCND contains the ratio of the smallest
!>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
!>          worth scaling by C.
!> \endverbatim
!>
!> \param[out] AMAX
!> \verbatim
!>          AMAX is REAL
!>          Absolute value of largest matrix element.  If AMAX is very
!>          close to overflow or very close to underflow, the matrix
!>          should be scaled.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i,  and i is
!>                <= M:  the i-th row of A is exactly zero
!>                >  M:  the (i-M)-th column of A is exactly zero
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
!> \ingroup geequ
!
!  =====================================================================
   SUBROUTINE SGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
                      INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, M, N
   REAL               AMAX, COLCND, ROWCND
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), C( * ), R( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               BIGNUM, RCMAX, RCMIN, SMLNUM
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   IF( M < 0 ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, M ) ) THEN
      INFO = -4
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SGEEQU', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M == 0 .OR. N == 0 ) THEN
      ROWCND = ONE
      COLCND = ONE
      AMAX = ZERO
      RETURN
   END IF
!
!     Get machine constants.
!
   SMLNUM = SLAMCH( 'S' )
   BIGNUM = ONE / SMLNUM
!
!     Compute row scale factors.
!
   DO I = 1, M
      R( I ) = ZERO
   ENDDO
!
!     Find the maximum element in each row.
!
   DO J = 1, N
      DO I = 1, M
         R( I ) = MAX( R( I ), ABS( A( I, J ) ) )
      ENDDO
   ENDDO
!
!     Find the maximum and minimum scale factors.
!
   RCMIN = BIGNUM
   RCMAX = ZERO
   DO I = 1, M
      RCMAX = MAX( RCMAX, R( I ) )
      RCMIN = MIN( RCMIN, R( I ) )
   ENDDO
   AMAX = RCMAX
!
   IF( RCMIN == ZERO ) THEN
!
!        Find the first zero scale factor and return an error code.
!
      DO I = 1, M
         IF( R( I ) == ZERO ) THEN
            INFO = I
            RETURN
         END IF
      ENDDO
   ELSE
!
!        Invert the scale factors.
!
      DO I = 1, M
         R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
      ENDDO
!
!        Compute ROWCND = min(R(I)) / max(R(I))
!
      ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
   END IF
!
!     Compute column scale factors
!
   DO J = 1, N
      C( J ) = ZERO
   ENDDO
!
!     Find the maximum element in each column,
!     assuming the row scaling computed above.
!
   DO J = 1, N
      DO I = 1, M
         C( J ) = MAX( C( J ), ABS( A( I, J ) )*R( I ) )
      ENDDO
   ENDDO
!
!     Find the maximum and minimum scale factors.
!
   RCMIN = BIGNUM
   RCMAX = ZERO
   DO J = 1, N
      RCMIN = MIN( RCMIN, C( J ) )
      RCMAX = MAX( RCMAX, C( J ) )
      ENDDO
!
   IF( RCMIN == ZERO ) THEN
!
!        Find the first zero scale factor and return an error code.
!
      DO J = 1, N
         IF( C( J ) == ZERO ) THEN
            INFO = M + J
            RETURN
         END IF
         ENDDO
   ELSE
!
!        Invert the scale factors.
!
      DO J = 1, N
         C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
         ENDDO
!
!        Compute COLCND = min(C(J)) / max(C(J))
!
      COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
   END IF
!
   RETURN
!
!     End of SGEEQU
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

