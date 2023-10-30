!> \brief \b CLANHB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a Hermitian band matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANHB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANHB( NORM, UPLO, N, K, AB, LDAB,
!                        WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            K, LDAB, N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANHB  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the element of  largest absolute value  of an
!> n by n hermitian band matrix A,  with k super-diagonals.
!> \endverbatim
!>
!> \return CLANHB
!> \verbatim
!>
!>    CLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANHB as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          band matrix A is supplied.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANHB is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of super-diagonals or sub-diagonals of the
!>          band matrix A.  K >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          The upper or lower triangle of the hermitian band matrix A,
!>          stored in the first K+1 rows of AB.  The j-th column of A is
!>          stored in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
!>          Note that the imaginary parts of the diagonal elements need
!>          not be set and are assumed to be zero.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= K+1.
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
!> \ingroup lanhb
!
!  =====================================================================
   REAL             FUNCTION CLANHB( NORM, UPLO, N, K, AB, LDAB, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          NORM, UPLO
   INTEGER            K, LDAB, N
!     ..
!     .. Array Arguments ..
   REAL               WORK( * )
   COMPLEX            AB( LDAB, * )
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J, L
   REAL               ABSA, SCALE, SOMME, VALUE
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
         DO J = 1, N
            DO I = MAX( K+2-J, 1 ), K
               SOMME = ABS( AB( I, J ) )
               IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
            ENDDO
            SOMME = ABS( REAL( AB( K+1, J ) ) )
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
         ENDDO
      ELSE
         DO J = 1, N
            SOMME = ABS( REAL( AB( 1, J ) ) )
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
            DO I = 2, MIN( N+1-J, K+1 )
               SOMME = ABS( AB( I, J ) )
               IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
            ENDDO
         ENDDO
      END IF
   ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
            ( NORM == '1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is hermitian).
!
      VALUE = 0.0E+0
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            SOMME = 0.0E+0
            L = K + 1 - J
            DO I = MAX( 1, J-K ), J - 1
               ABSA = ABS( AB( L+I, J ) )
               SOMME = SOMME + ABSA
               WORK( I ) = WORK( I ) + ABSA
            ENDDO
            WORK( J ) = SOMME + ABS( REAL( AB( K+1, J ) ) )
         ENDDO
         DO I = 1, N
            IF( VALUE  <  WORK( I ) .OR. SISNAN( WORK( I ) ) ) VALUE = WORK( I )
         ENDDO
      ELSE
         DO I = 1, N
            WORK( I ) = 0.0E+0
         ENDDO
         DO J = 1, N
            SOMME = WORK( J ) + ABS( REAL( AB( 1, J ) ) )
            L = 1 - J
            DO I = J + 1, MIN( N, J+K )
               ABSA = ABS( AB( L+I, J ) )
               SOMME = SOMME + ABSA
               WORK( I ) = WORK( I ) + ABSA
            ENDDO
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
            ENDDO
      END IF
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      SCALE = 0.0E+0
      SOMME = 1.0E+0
      IF( K > 0 ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO J = 2, N
               CALL CLASSQ( MIN( J-1, K ), AB( MAX( K+2-J, 1 ), J ), 1, SCALE, SOMME )
               ENDDO
            L = K + 1
         ELSE
            DO J = 1, N - 1
               CALL CLASSQ( MIN( N-J, K ), AB( 2, J ), 1, SCALE, SOMME )
               ENDDO
            L = 1
         END IF
         SOMME = 2*SOMME
      ELSE
         L = 1
      END IF
      DO J = 1, N
         IF( REAL( AB( L, J ) ) /= 0.0E+0 ) THEN
            ABSA = ABS( REAL( AB( L, J ) ) )
            IF( SCALE < ABSA ) THEN
               SOMME = 1.0E+0 + SOMME*( SCALE / ABSA )**2
               SCALE = ABSA
            ELSE
               SOMME = SOMME + ( ABSA / SCALE )**2
            END IF
         END IF
         ENDDO
      VALUE = SCALE*SQRT( SOMME )
   END IF
!
   CLANHB = VALUE
   RETURN
!
!     End of CLANHB
!
END
