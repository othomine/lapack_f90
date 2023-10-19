!> \brief \b SSTECT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSTECT( N, A, B, SHIFT, NUM )
!
!       .. Scalar Arguments ..
!       INTEGER            N, NUM
!       REAL               SHIFT
!       ..
!       .. Array Arguments ..
!       REAL               A( * ), B( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SSTECT counts the number NUM of eigenvalues of a tridiagonal
!>    matrix T which are less than or equal to SHIFT. T has
!>    diagonal entries A(1), ... , A(N), and offdiagonal entries
!>    B(1), ..., B(N-1).
!>    See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!>    Matrix", Report CS41, Computer Science Dept., Stanford
!>    University, July 21, 1966
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (N)
!>          The diagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (N-1)
!>          The offdiagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] SHIFT
!> \verbatim
!>          SHIFT is REAL
!>          The shift, used as described under Purpose.
!> \endverbatim
!>
!> \param[out] NUM
!> \verbatim
!>          NUM is INTEGER
!>          The number of eigenvalues of T less than or equal
!>          to SHIFT.
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SSTECT( N, A, B, SHIFT, NUM )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            N, NUM
   REAL               SHIFT
!     ..
!     .. Array Arguments ..
   REAL               A( * ), B( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I
   REAL               M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, &
                      TOM, U, UNFL
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. Executable Statements ..
!
!     Get machine constants
!
   UNFL = SLAMCH( 'Safe minimum' )
   OVFL = SLAMCH( 'Overflow' )
!
!     Find largest entry
!
   MX = ABS( A( 1 ) )
   DO I = 1, N - 1
      MX = MAX( MX, ABS( A( I+1 ) ), ABS( B( I ) ) )
   ENDDO
!
!     Handle easy cases, including zero matrix
!
   IF( SHIFT >= 3.0E+0*MX ) THEN
      NUM = N
      RETURN
   END IF
   IF( SHIFT < -3.0E+0*MX ) THEN
      NUM = 0
      RETURN
   END IF
!
!     Compute scale factors as in Kahan's report
!     At this point, MX  /=  0 so we can divide by it
!
   SUN = SQRT( UNFL )
   SSUN = SQRT( SUN )
   SOV = SQRT( OVFL )
   TOM = SSUN*SOV
   IF( MX <= 1.0E+0 ) THEN
      M1 = 1.0E+0 / MX
      M2 = TOM
   ELSE
      M1 = 1.0E+0
      M2 = TOM / MX
   END IF
!
!     Begin counting
!
   NUM = 0
   SSHIFT = ( SHIFT*M1 )*M2
   U = ( A( 1 )*M1 )*M2 - SSHIFT
   IF( U <= SUN ) THEN
      IF( U <= 0.0E+0 ) THEN
         NUM = NUM + 1
         IF( U > -SUN ) U = -SUN
      ELSE
         U = SUN
      END IF
   END IF
   DO I = 2, N
      TMP = ( B( I-1 )*M1 )*M2
      U = ( ( A( I )*M1 )*M2-TMP*( TMP / U ) ) - SSHIFT
      IF( U <= SUN ) THEN
         IF( U <= 0.0E+0 ) THEN
            NUM = NUM + 1
            IF( U > -SUN ) U = -SUN
         ELSE
            U = SUN
         END IF
      END IF
   ENDDO
   RETURN
!
!     End of SSTECT
!
END


