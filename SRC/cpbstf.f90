!> \brief \b CPBSTF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPBSTF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbstf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbstf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbstf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPBSTF( UPLO, N, KD, AB, LDAB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPBSTF computes a split Cholesky factorization of a complex
!> Hermitian positive definite band matrix A.
!>
!> This routine is designed to be used in conjunction with CHBGST.
!>
!> The factorization has the form  A = S**H*S  where S is a band matrix
!> of the same bandwidth as A and the following structure:
!>
!>   S = ( U    )
!>       ( M  L )
!>
!> where U is upper triangular of order m = (n+kd)/2, and L is lower
!> triangular of order n-m.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the Hermitian band
!>          matrix A, stored in the first kd+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, if INFO = 0, the factor S from the split Cholesky
!>          factorization A = S**H*S. See Further Details.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, the factorization could not be completed,
!>               because the updated element a(i,i) was negative; the
!>               matrix A is not positive definite.
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
!> \ingroup pbstf
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  N = 7, KD = 2:
!>
!>  S = ( s11  s12  s13                     )
!>      (      s22  s23  s24                )
!>      (           s33  s34                )
!>      (                s44                )
!>      (           s53  s54  s55           )
!>      (                s64  s65  s66      )
!>      (                     s75  s76  s77 )
!>
!>  If UPLO = 'U', the array AB holds:
!>
!>  on entry:                          on exit:
!>
!>   *    *   a13  a24  a35  a46  a57   *    *   s13  s24  s53**H s64**H s75**H
!>   *   a12  a23  a34  a45  a56  a67   *   s12  s23  s34  s54**H s65**H s76**H
!>  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55    s66    s77
!>
!>  If UPLO = 'L', the array AB holds:
!>
!>  on entry:                          on exit:
!>
!>  a11  a22  a33  a44  a55  a66  a77  s11    s22    s33    s44  s55  s66  s77
!>  a21  a32  a43  a54  a65  a76   *   s12**H s23**H s34**H s54  s65  s76   *
!>  a31  a42  a53  a64  a64   *    *   s13**H s24**H s53    s64  s75   *    *
!>
!>  Array elements marked * are not used by the routine; s12**H denotes
!>  conjg(s12); the diagonal elements of S are real.
!> \endverbatim
!>
!  =====================================================================
   SUBROUTINE CPBSTF( UPLO, N, KD, AB, LDAB, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
   COMPLEX            AB( LDAB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            J, KLD, KM, M
   REAL               AJJ
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CHER, CLACGV, CSSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN, REAL, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( KD < 0 ) THEN
      INFO = -3
   ELSE IF( LDAB < KD+1 ) THEN
      INFO = -5
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CPBSTF', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) &
      RETURN
!
   KLD = MAX( 1, LDAB-1 )
!
!     Set the splitting point m.
!
   M = ( N+KD ) / 2
!
   IF( UPPER ) THEN
!
!        Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m).
!
      DO J = N, M + 1, -1
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
         AJJ = REAL( AB( KD+1, J ) )
         IF( AJJ <= ZERO ) THEN
            AB( KD+1, J ) = AJJ
            GO TO 50
         END IF
         AJJ = SQRT( AJJ )
         AB( KD+1, J ) = AJJ
         KM = MIN( J-1, KD )
!
!           Compute elements j-km:j-1 of the j-th column and update the
!           the leading submatrix within the band.
!
         CALL CSSCAL( KM, ONE / AJJ, AB( KD+1-KM, J ), 1 )
         CALL CHER( 'Upper', KM, -ONE, AB( KD+1-KM, J ), 1, &
                    AB( KD+1, J-KM ), KLD )
      ENDDO
!
!        Factorize the updated submatrix A(1:m,1:m) as U**H*U.
!
      DO J = 1, M
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
         AJJ = REAL( AB( KD+1, J ) )
         IF( AJJ <= ZERO ) THEN
            AB( KD+1, J ) = AJJ
            GO TO 50
         END IF
         AJJ = SQRT( AJJ )
         AB( KD+1, J ) = AJJ
         KM = MIN( KD, M-J )
!
!           Compute elements j+1:j+km of the j-th row and update the
!           trailing submatrix within the band.
!
         IF( KM > 0 ) THEN
            CALL CSSCAL( KM, ONE / AJJ, AB( KD, J+1 ), KLD )
            CALL CLACGV( KM, AB( KD, J+1 ), KLD )
            CALL CHER( 'Upper', KM, -ONE, AB( KD, J+1 ), KLD, &
                       AB( KD+1, J+1 ), KLD )
            CALL CLACGV( KM, AB( KD, J+1 ), KLD )
         END IF
      ENDDO
   ELSE
!
!        Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m).
!
      DO J = N, M + 1, -1
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
         AJJ = REAL( AB( 1, J ) )
         IF( AJJ <= ZERO ) THEN
            AB( 1, J ) = AJJ
            GO TO 50
         END IF
         AJJ = SQRT( AJJ )
         AB( 1, J ) = AJJ
         KM = MIN( J-1, KD )
!
!           Compute elements j-km:j-1 of the j-th row and update the
!           trailing submatrix within the band.
!
         CALL CSSCAL( KM, ONE / AJJ, AB( KM+1, J-KM ), KLD )
         CALL CLACGV( KM, AB( KM+1, J-KM ), KLD )
         CALL CHER( 'Lower', KM, -ONE, AB( KM+1, J-KM ), KLD, &
                    AB( 1, J-KM ), KLD )
         CALL CLACGV( KM, AB( KM+1, J-KM ), KLD )
      ENDDO
!
!        Factorize the updated submatrix A(1:m,1:m) as U**H*U.
!
      DO J = 1, M
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
         AJJ = REAL( AB( 1, J ) )
         IF( AJJ <= ZERO ) THEN
            AB( 1, J ) = AJJ
            GO TO 50
         END IF
         AJJ = SQRT( AJJ )
         AB( 1, J ) = AJJ
         KM = MIN( KD, M-J )
!
!           Compute elements j+1:j+km of the j-th column and update the
!           trailing submatrix within the band.
!
         IF( KM > 0 ) THEN
            CALL CSSCAL( KM, ONE / AJJ, AB( 2, J ), 1 )
            CALL CHER( 'Lower', KM, -ONE, AB( 2, J ), 1, &
                       AB( 1, J+1 ), KLD )
         END IF
      ENDDO
   END IF
   RETURN
!
50 CONTINUE
   INFO = J
   RETURN
!
!     End of CPBSTF
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        