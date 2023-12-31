!> \brief \b CLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a complex symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANSY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANSY( NORM, UPLO, N, A, LDA, WORK )
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
!> CLANSY  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex symmetric matrix A.
!> \endverbatim
!>
!> \return CLANSY
!> \verbatim
!>
!>    CLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANSY as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is to be referenced.
!>          = 'U':  Upper triangular part of A is referenced
!>          = 'L':  Lower triangular part of A is referenced
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANSY is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of A contains the lower triangular part of
!>          the matrix A, and the strictly upper triangular part of A is
!>          not referenced.
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
   REAL             FUNCTION CLANSY( NORM, UPLO, N, A, LDA, WORK )
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
            DO I = 1, J
               SOMME = ABS( A( I, J ) )
               IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
            ENDDO
         ENDDO
      ELSE
         DO J = 1, N
            DO I = J, N
               SOMME = ABS( A( I, J ) )
               IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
            ENDDO
         ENDDO
      END IF
   ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
            ( NORM == '1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
      VALUE = 0.0E+0
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            WORK(1:J-1) = WORK(1:J-1) + ABS(A(1:J-1,J))
            WORK(J) = SUM(ABS(A(1:J,J)))
         ENDDO
         DO I = 1, N
            SOMME = WORK( I )
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
         ENDDO
      ELSE
         WORK(1:N) = 0.0E+0
         DO J = 1, N
            SOMME = WORK( J ) + SUM(ABS(A(J:N,J)))
            WORK(J+1:N) = WORK(J+1:N) + ABS(A(J+1:N,J))
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
      CALL CLASSQ( N, A, LDA+1, SCALE, SOMME )
      VALUE = SCALE*SQRT( SOMME )
   END IF
!
   CLANSY = VALUE
   RETURN
!
!     End of CLANSY
!
END
