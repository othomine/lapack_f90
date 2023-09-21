!> \brief \b SLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular matrix supplied in packed form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLANTP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLANTP( NORM, UPLO, DIAG, N, AP, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       REAL               AP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLANTP  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> triangular matrix A, supplied in packed form.
!> \endverbatim
!>
!> \return SLANTP
!> \verbatim
!>
!>    SLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in SLANTP as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, SLANTP is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          The upper or lower triangular matrix A, packed columnwise in
!>          a linear array.  The j-th column of A is stored in the array
!>          AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>          Note that when DIAG = 'U', the elements of the array AP
!>          corresponding to the diagonal elements of the matrix A are
!>          not referenced, but are assumed to be one.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!>          referenced.
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
!> \ingroup lantp
!
!  =====================================================================
   REAL             FUNCTION SLANTP( NORM, UPLO, DIAG, N, AP, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, NORM, UPLO
   INTEGER            N
!     ..
!     .. Array Arguments ..
   REAL               AP( * ), WORK( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UDIAG
   INTEGER            I, J, K
   REAL               SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SLASSQ
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, SISNAN
   EXTERNAL           LSAME, SISNAN
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
   IF( N == 0 ) THEN
      VALUE = ZERO
   ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
      K = 1
      IF( LSAME( DIAG, 'U' ) ) THEN
         VALUE = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = K, K + J - 2
                  SUM = ABS( AP( I ) )
                  IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               ENDDO
               K = K + J
            ENDDO
         ELSE
            DO J = 1, N
               DO I = K + 1, K + N - J
                  SUM = ABS( AP( I ) )
                  IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               ENDDO
               K = K + N - J + 1
            ENDDO
         END IF
      ELSE
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 1, N
               DO I = K, K + J - 1
                  SUM = ABS( AP( I ) )
                  IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               ENDDO
               K = K + J
            ENDDO
         ELSE
            DO J = 1, N
               DO I = K, K + N - J
                  SUM = ABS( AP( I ) )
                  IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               ENDDO
               K = K + N - J + 1
            ENDDO
         END IF
      END IF
   ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN
!
!        Find norm1(A).
!
      VALUE = ZERO
      K = 1
      UDIAG = LSAME( DIAG, 'U' )
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            IF( UDIAG ) THEN
               SUM = ONE
               DO I = K, K + J - 2
                  SUM = SUM + ABS( AP( I ) )
               ENDDO
            ELSE
               SUM = ZERO
               DO I = K, K + J - 1
                  SUM = SUM + ABS( AP( I ) )
                  ENDDO
            END IF
            K = K + J
            IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            ENDDO
      ELSE
         DO J = 1, N
            IF( UDIAG ) THEN
               SUM = ONE
               DO I = K + 1, K + N - J
                  SUM = SUM + ABS( AP( I ) )
                  ENDDO
            ELSE
               SUM = ZERO
               DO I = K, K + N - J
                  SUM = SUM + ABS( AP( I ) )
                  ENDDO
            END IF
            K = K + N - J + 1
            IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            ENDDO
      END IF
   ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
      K = 1
      IF( LSAME( UPLO, 'U' ) ) THEN
         IF( LSAME( DIAG, 'U' ) ) THEN
            DO I = 1, N
               WORK( I ) = ONE
               ENDDO
            DO J = 1, N
               DO I = 1, J - 1
                  WORK( I ) = WORK( I ) + ABS( AP( K ) )
                  K = K + 1
                  ENDDO
               K = K + 1
               ENDDO
         ELSE
            DO I = 1, N
               WORK( I ) = ZERO
               ENDDO
            DO J = 1, N
               DO I = 1, J
                  WORK( I ) = WORK( I ) + ABS( AP( K ) )
                  K = K + 1
                  ENDDO
               ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            DO I = 1, N
               WORK( I ) = ONE
               ENDDO
            DO J = 1, N
               K = K + 1
               DO I = J + 1, N
                  WORK( I ) = WORK( I ) + ABS( AP( K ) )
                  K = K + 1
                  ENDDO
               ENDDO
         ELSE
            DO I = 1, N
               WORK( I ) = ZERO
               ENDDO
            DO J = 1, N
               DO I = J, N
                  WORK( I ) = WORK( I ) + ABS( AP( K ) )
                  K = K + 1
                  ENDDO
               ENDDO
         END IF
      END IF
      VALUE = ZERO
      DO I = 1, N
         SUM = WORK( I )
         IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
         ENDDO
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         IF( LSAME( DIAG, 'U' ) ) THEN
            SCALE = ONE
            SUM = N
            K = 2
            DO J = 2, N
               CALL SLASSQ( J-1, AP( K ), 1, SCALE, SUM )
               K = K + J
               ENDDO
         ELSE
            SCALE = ZERO
            SUM = ONE
            K = 1
            DO J = 1, N
               CALL SLASSQ( J, AP( K ), 1, SCALE, SUM )
               K = K + J
               ENDDO
         END IF
      ELSE
         IF( LSAME( DIAG, 'U' ) ) THEN
            SCALE = ONE
            SUM = N
            K = 2
            DO J = 1, N - 1
               CALL SLASSQ( N-J, AP( K ), 1, SCALE, SUM )
               K = K + N - J + 1
               ENDDO
         ELSE
            SCALE = ZERO
            SUM = ONE
            K = 1
            DO J = 1, N
               CALL SLASSQ( N-J+1, AP( K ), 1, SCALE, SUM )
               K = K + N - J + 1
               ENDDO
         END IF
      END IF
      VALUE = SCALE*SQRT( SUM )
   END IF
!
   SLANTP = VALUE
   RETURN
!
!     End of SLANTP
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
