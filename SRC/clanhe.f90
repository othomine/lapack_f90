!> \brief \b CLANHE returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a complex Hermitian matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANHE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhe.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhe.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhe.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANHE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex hermitian matrix A.
!> \endverbatim
!>
!> \return CLANHE
!> \verbatim
!>
!>    CLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANHE as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          hermitian matrix A is to be referenced.
!>          = 'U':  Upper triangular part of A is referenced
!>          = 'L':  Lower triangular part of A is referenced
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANHE is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The hermitian matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of A contains the lower triangular part of
!>          the matrix A, and the strictly upper triangular part of A is
!>          not referenced. Note that the imaginary parts of the diagonal
!>          elements need not be set and are assumed to be zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(N,1).
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
!> \ingroup lanhe
!
!  =====================================================================
   REAL             FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          NORM, UPLO
   INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
   REAL               WORK( * )
   COMPLEX            A( LDA, * )
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
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
            DO I = 1, J - 1
               IF( VALUE  <  ABS( A( I, J ) ) .OR. SISNAN( ABS( A( I, J ) ) ) ) VALUE = ABS( A( I, J ) )
            ENDDO
            IF( VALUE<ABS(REAL(A(J,J))) .OR. SISNAN(ABS(REAL(A(J,J))))) VALUE = ABS(REAL(A(J,J)))
         ENDDO
      ELSE
         DO J = 1, N
            IF( VALUE  <  ABS(REAL(A(J,J))) .OR. SISNAN( ABS(REAL(A(J,J)))) ) VALUE = ABS(REAL(A(J,J)))
            DO I = J + 1, N
               IF( VALUE < ABS(A(I,J)) .OR. SISNAN( ABS(A(I,J))) ) VALUE = ABS(A(I,J))
            ENDDO
         ENDDO
      END IF
   ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is hermitian).
!
      VALUE = 0.0E+0
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            WORK(1:J-1) = WORK(1:J-1) + ABS(A(1:J-1,J))
            WORK( J ) = SUM(ABS(A(1:J-1,J))) + ABS(REAL(A(J,J)))
         ENDDO
         DO I = 1, N
            IF( VALUE  <  WORK(I) .OR. SISNAN( WORK(I) ) ) VALUE = WORK(I)
         ENDDO
      ELSE
         WORK(1:N) = 0.0E+0
         DO J = 1, N
            SOMME = WORK( J ) + ABS( REAL( A( J, J ) ) )
            WORK(J+1:N) = WORK(J+1:N) + ABS(A(J+1:N,J))
            SOMME = SOMME + SUM(ABS(A(J+1:N,J)))
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
         ENDDO
      END IF
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      SCALE = 0.0E+0
      SOMME = 1.0E+0
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 2, N
            CALL CLASSQ( J-1, A( 1, J ), 1, SCALE, SOMME )
         ENDDO
      ELSE
         DO J = 1, N - 1
            CALL CLASSQ( N-J, A( J+1, J ), 1, SCALE, SOMME )
         ENDDO
      END IF
      SOMME = 2*SOMME
      DO I = 1, N
         IF( REAL( A( I, I ) ) /= 0.0E+0 ) THEN
            ABSA = ABS( REAL( A( I, I ) ) )
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
   CLANHE = VALUE
   RETURN
!
!     End of CLANHE
!
END
