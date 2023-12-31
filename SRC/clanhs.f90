!> \brief \b CLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLANHS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANHS( NORM, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
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
!> CLANHS  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> Hessenberg matrix A.
!> \endverbatim
!>
!> \return CLANHS
!> \verbatim
!>
!>    CLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANHS as described
!>          above.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, CLANHS is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The n by n upper Hessenberg matrix A; the part of A below the
!>          first sub-diagonal is not referenced.
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup lanhs
!
!  =====================================================================
   REAL             FUNCTION CLANHS( NORM, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          NORM
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
   REAL               SCALE, SOMME, VALUE
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
      DO J = 1, N
         DO I = 1, MIN( N, J+1 )
            SOMME = ABS( A( I, J ) )
            IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
         ENDDO
      ENDDO
   ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN
!
!        Find norm1(A).
!
      VALUE = 0.0E+0
      DO J = 1, N
         SOMME = SUM(ABS(A(1:MIN(N,J+1),J)))
         IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
      ENDDO
   ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
      WORK(1:N) = 0.0E+0
      DO J = 1, N
         WORK(1:MIN(N,J+1)) = WORK(1:MIN(N,J+1)) + ABS(A(1:MIN(N,J+1),J))
      ENDDO
      VALUE = 0.0E+0
      DO I = 1, N
         SOMME = WORK( I )
         IF( VALUE  <  SOMME .OR. SISNAN( SOMME ) ) VALUE = SOMME
      ENDDO
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      SCALE = 0.0E+0
      SOMME = 1.0E+0
      DO J = 1, N
         CALL CLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SOMME )
      ENDDO
      VALUE = SCALE*SQRT( SOMME )
   END IF
!
   CLANHS = VALUE
   RETURN
!
!     End of CLANHS
!
END
