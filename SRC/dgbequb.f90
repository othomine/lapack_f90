!> \brief \b DGBEQUB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBEQUB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbequb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbequb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbequb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,
!                           AMAX, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, M, N
!       DOUBLE PRECISION   AMAX, COLCND, ROWCND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AB( LDAB, * ), C( * ), R( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBEQUB computes row and column scalings intended to equilibrate an
!> M-by-N matrix A and reduce its condition number.  R returns the row
!> scale factors and C the column scale factors, chosen to try to make
!> the largest element in each row and column of the matrix B with
!> elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most
!> the radix.
!>
!> R(i) and C(j) are restricted to be a power of the radix between
!> SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use
!> of these scaling factors is not guaranteed to reduce the condition
!> number of A but works well in practice.
!>
!> This routine differs from DGEEQU by restricting the scaling factors
!> to a power of the radix.  Barring over- and underflow, scaling by
!> these factors introduces no additional rounding errors.  However, the
!> scaled entries' magnitudes are no longer approximately 1 but lie
!> between sqrt(radix) and 1/sqrt(radix).
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array A.  LDAB >= max(1,M).
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (M)
!>          If INFO = 0 or INFO > M, R contains the row scale factors
!>          for A.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0,  C contains the column scale factors for A.
!> \endverbatim
!>
!> \param[out] ROWCND
!> \verbatim
!>          ROWCND is DOUBLE PRECISION
!>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
!>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
!>          AMAX is neither too large nor too small, it is not worth
!>          scaling by R.
!> \endverbatim
!>
!> \param[out] COLCND
!> \verbatim
!>          COLCND is DOUBLE PRECISION
!>          If INFO = 0, COLCND contains the ratio of the smallest
!>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
!>          worth scaling by C.
!> \endverbatim
!>
!> \param[out] AMAX
!> \verbatim
!>          AMAX is DOUBLE PRECISION
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
!> \ingroup gbequb
!
!  =====================================================================
   SUBROUTINE DGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, &
                       AMAX, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, KL, KU, LDAB, M, N
   DOUBLE PRECISION   AMAX, COLCND, ROWCND
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AB( LDAB, * ), C( * ), R( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, KD
   DOUBLE PRECISION   BIGNUM, RCMAX, RCMIN, SMLNUM, RADIX, LOGRDX
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
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
   ELSE IF( KL < 0 ) THEN
      INFO = -3
   ELSE IF( KU < 0 ) THEN
      INFO = -4
   ELSE IF( LDAB < KL+KU+1 ) THEN
      INFO = -6
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'DGBEQUB', -INFO )
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF( M == 0 .OR. N == 0 ) THEN
      ROWCND = 1.0D0
      COLCND = 1.0D0
      AMAX = 0.0D0
      RETURN
   END IF
!
!     Get machine constants.  Assume SMLNUM is a power of the radix.
!
   SMLNUM = DLAMCH( 'S' )
   BIGNUM = 1.0D0 / SMLNUM
   RADIX = DLAMCH( 'B' )
   LOGRDX = LOG(RADIX)
!
!     Compute row scale factors.
!
   R(1:M) = 0.0D0
!
!     Find the maximum element in each row.
!
   KD = KU + 1
   DO J = 1, N
      DO I = MAX( J-KU, 1 ), MIN( J+KL, M )
         R( I ) = MAX( R( I ), ABS( AB( KD+I-J, J ) ) )
      ENDDO
   ENDDO
   DO I = 1, M
      IF( R( I ) > 0.0D0 ) THEN
         R( I ) = RADIX**INT( LOG( R( I ) ) / LOGRDX )
      END IF
   END DO
!
!     Find the maximum and minimum scale factors.
!
   RCMIN = MINVAL(R(1:M))
   RCMAX = MAXVAL(R(1:M))
   AMAX = RCMAX
!
   IF( RCMIN == 0.0D0 ) THEN
!
!        Find the first zero scale factor and return an error code.
!
      DO I = 1, M
         IF( R( I ) == 0.0D0 ) THEN
            INFO = I
            RETURN
         END IF
      ENDDO
   ELSE
!
!        Invert the scale factors.
!
      DO I = 1, M
         R( I ) = 1.0D0 / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
      ENDDO
!
!        Compute ROWCND = min(R(I)) / max(R(I)).
!
      ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
   END IF
!
!     Compute column scale factors.
!
   C(1:N) = 0.0D0
!
!     Find the maximum element in each column,
!     assuming the row scaling computed above.
!
   DO J = 1, N
      DO I = MAX( J-KU, 1 ), MIN( J+KL, M )
         C( J ) = MAX( C( J ), ABS( AB( KD+I-J, J ) )*R( I ) )
      ENDDO
      IF( C( J ) > 0.0D0 ) THEN
         C( J ) = RADIX**INT( LOG( C( J ) ) / LOGRDX )
      END IF
   ENDDO
!
!     Find the maximum and minimum scale factors.
!
   RCMIN = MINVAL(C(1:N))
   RCMAX = MAXVAL(C(1:N))
!
   IF( RCMIN == 0.0D0 ) THEN
!
!        Find the first zero scale factor and return an error code.
!
      DO J = 1, N
         IF( C( J ) == 0.0D0 ) THEN
            INFO = M + J
            RETURN
         END IF
      ENDDO
   ELSE
!
!        Invert the scale factors.
!
      DO J = 1, N
         C( J ) = 1.0D0 / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
      ENDDO
!
!        Compute COLCND = min(C(J)) / max(C(J)).
!
      COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
   END IF
!
   RETURN
!
!     End of DGBEQUB
!
END
