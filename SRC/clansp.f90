!> \brief \b CLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric matrix supplied in packed form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANSP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANSP( NORM, UPLO, N, AP, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANSP  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex symmetric matrix A,  supplied in packed form.
!> \endverbatim
!>
!> \return CLANSP
!> \verbatim
!>
!>    CLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANSP as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is supplied.
!>          = 'U':  Upper triangular part of A is supplied
!>          = 'L':  Lower triangular part of A is supplied
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANSP is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          The upper or lower triangle of the symmetric matrix A, packed
!>          columnwise in a linear array.  The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!>          WORK is not referenced.
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
!> \ingroup lanhp
!
!  =====================================================================
   REAL             FUNCTION CLANSP( NORM, UPLO, N, AP, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          NORM, UPLO
   INTEGER            N
!     ..
!     .. Array Arguments ..
   REAL               WORK( * )
   COMPLEX            AP( * )
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, K
   REAL               ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Functions ..
   LOGICAL            LSAME, SISNAN
   EXTERNAL           LSAME, SISNAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLASSQ
!     ..
!     .. Executable Statements ..
!
   IF( N == 0 ) THEN
      VALUE = 0.0E+0
   ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
      VALUE = 0.0E+0
      IF( LSAME( UPLO, 'U' ) ) THEN
         K = 1
         DO J = 1, N
            DO I = K, K + J - 1
               SUM = ABS( AP( I ) )
               IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            ENDDO
            K = K + J
         ENDDO
      ELSE
         K = 1
         DO J = 1, N
            DO I = K, K + N - J
               SUM = ABS( AP( I ) )
               IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            ENDDO
            K = K + N - J + 1
         ENDDO
      END IF
   ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
      VALUE = 0.0E+0
      K = 1
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            SUM = 0.0E+0
            DO I = 1, J - 1
               ABSA = ABS( AP( K ) )
               SUM = SUM + ABSA
               WORK( I ) = WORK( I ) + ABSA
               K = K + 1
            ENDDO
            WORK( J ) = SUM + ABS( AP( K ) )
            K = K + 1
         ENDDO
         DO I = 1, N
            SUM = WORK( I )
            IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
         ENDDO
      ELSE
         WORK(1:N) = 0.0E+0
         DO J = 1, N
            SUM = WORK( J ) + ABS( AP( K ) )
            K = K + 1
            DO I = J + 1, N
               ABSA = ABS( AP( K ) )
               SUM = SUM + ABSA
               WORK( I ) = WORK( I ) + ABSA
               K = K + 1
            ENDDO
            IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
         ENDDO
      END IF
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      SCALE = 0.0E+0
      SUM = 1.0E+0
      K = 2
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 2, N
            CALL CLASSQ( J-1, AP( K ), 1, SCALE, SUM )
            K = K + J
         ENDDO
      ELSE
         DO J = 1, N - 1
            CALL CLASSQ( N-J, AP( K ), 1, SCALE, SUM )
            K = K + N - J + 1
         ENDDO
      END IF
      SUM = 2*SUM
      K = 1
      DO I = 1, N
         IF( REAL( AP( K ) ) /= 0.0E+0 ) THEN
            ABSA = ABS( REAL( AP( K ) ) )
            IF( SCALE < ABSA ) THEN
               SUM = 1.0E+0 + SUM*( SCALE / ABSA )**2
               SCALE = ABSA
            ELSE
               SUM = SUM + ( ABSA / SCALE )**2
            END IF
         END IF
         IF( AIMAG( AP( K ) ) /= 0.0E+0 ) THEN
            ABSA = ABS( AIMAG( AP( K ) ) )
            IF( SCALE < ABSA ) THEN
               SUM = 1.0E+0 + SUM*( SCALE / ABSA )**2
               SCALE = ABSA
            ELSE
               SUM = SUM + ( ABSA / SCALE )**2
            END IF
         END IF
         IF( LSAME( UPLO, 'U' ) ) THEN
            K = K + I + 1
         ELSE
            K = K + N - I + 1
         END IF
      ENDDO
      VALUE = SCALE*SQRT( SUM )
   END IF
!
   CLANSP = VALUE
   RETURN
!
!     End of CLANSP
!
END
