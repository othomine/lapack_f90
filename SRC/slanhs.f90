!> \brief \b SLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLANHS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slanhs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slanhs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slanhs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLANHS( NORM, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLANHS  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> Hessenberg matrix A.
!> \endverbatim
!>
!> \return SLANHS
!> \verbatim
!>
!>    SLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in SLANHS as described
!>          above.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, SLANHS is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
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
   REAL             FUNCTION SLANHS( NORM, N, A, LDA, WORK )
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
   REAL               A( LDA, * ), WORK( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
   REAL               ONE, ZERO
   PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
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
   INTRINSIC          ABS, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
   IF( N == 0 ) THEN
      VALUE = ZERO
   ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
      VALUE = ZERO
      DO J = 1, N
         DO I = 1, MIN( N, J+1 )
            SUM = ABS( A( I, J ) )
            IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
         ENDDO
      ENDDO
   ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN
!
!        Find norm1(A).
!
      VALUE = ZERO
      DO J = 1, N
         SUM = ZERO
         DO I = 1, MIN( N, J+1 )
            SUM = SUM + ABS( A( I, J ) )
         ENDDO
         IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
      ENDDO
   ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
      DO I = 1, N
         WORK( I ) = ZERO
      ENDDO
      DO J = 1, N
         DO I = 1, MIN( N, J+1 )
            WORK( I ) = WORK( I ) + ABS( A( I, J ) )
         ENDDO
      ENDDO
      VALUE = ZERO
      DO I = 1, N
         SUM = WORK( I )
         IF( VALUE  <  SUM .OR. SISNAN( SUM ) ) VALUE = SUM
      ENDDO
   ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
      SCALE = ZERO
      SUM = ONE
      DO J = 1, N
         CALL SLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
      ENDDO
      VALUE = SCALE*SQRT( SUM )
   END IF
!
   SLANHS = VALUE
   RETURN
!
!     End of SLANHS
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

