!> \brief \b DSVDCT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSVDCT( N, S, E, SHIFT, NUM )
!
!       .. Scalar Arguments ..
!       INTEGER            N, NUM
!       DOUBLE PRECISION   SHIFT
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   E( * ), S( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSVDCT counts the number NUM of eigenvalues of a 2*N by 2*N
!> tridiagonal matrix T which are less than or equal to SHIFT.  T is
!> formed by putting zeros on the diagonal and making the off-diagonals
!> equal to S(1), E(1), S(2), E(2), ... , E(N-1), S(N).  If SHIFT is
!> positive, NUM is equal to N plus the number of singular values of a
!> bidiagonal matrix B less than or equal to SHIFT.  Here B has diagonal
!> entries S(1), ..., S(N) and superdiagonal entries E(1), ... E(N-1).
!> If SHIFT is negative, NUM is equal to the number of singular values
!> of B greater than or equal to -SHIFT.
!>
!> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!> Matrix", Report CS41, Computer Science Dept., Stanford University,
!> July 21, 1966
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>          The diagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array of dimension (N-1)
!>          The superdiagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] SHIFT
!> \verbatim
!>          SHIFT is DOUBLE PRECISION
!>          The shift, used as described under Purpose.
!> \endverbatim
!>
!> \param[out] NUM
!> \verbatim
!>          NUM is INTEGER
!>          The number of eigenvalues of T less than or equal to SHIFT.
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
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DSVDCT( N, S, E, SHIFT, NUM )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            N, NUM
   DOUBLE PRECISION   SHIFT
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   E( * ), S( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I
   DOUBLE PRECISION   M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, &
                      TOM, U, UNFL
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. Executable Statements ..
!
!     Get machine constants
!
   UNFL = 2*DLAMCH( 'Safe minimum' )
   OVFL = 1.0D0 / UNFL
!
!     Find largest entry
!
   MX = ABS( S( 1 ) )
   DO I = 1, N - 1
      MX = MAX( MX, ABS( S( I+1 ) ), ABS( E( I ) ) )
   ENDDO
!
   IF( MX == 0.0D+0 ) THEN
      IF( SHIFT < 0.0D+0 ) THEN
         NUM = 0
      ELSE
         NUM = 2*N
      END IF
      RETURN
   END IF
!
!     Compute scale factors as in Kahan's report
!
   SUN = SQRT( UNFL )
   SSUN = SQRT( SUN )
   SOV = SQRT( OVFL )
   TOM = SSUN*SOV
   IF( MX <= 1.0D0 ) THEN
      M1 = 1.0D0 / MX
      M2 = TOM
   ELSE
      M1 = 1.0D0
      M2 = TOM / MX
   END IF
!
!     Begin counting
!
   U = 1.0D0
   NUM = 0
   SSHIFT = ( SHIFT*M1 )*M2
   U = -SSHIFT
   IF( U <= SUN ) THEN
      IF( U <= 0.0D+0 ) THEN
         NUM = NUM + 1
         IF( U > -SUN ) &
            U = -SUN
      ELSE
         U = SUN
      END IF
   END IF
   TMP = ( S( 1 )*M1 )*M2
   U = -TMP*( TMP / U ) - SSHIFT
   IF( U <= SUN ) THEN
      IF( U <= 0.0D+0 ) THEN
         NUM = NUM + 1
         IF( U > -SUN ) &
            U = -SUN
      ELSE
         U = SUN
      END IF
   END IF
   DO I = 1, N - 1
      TMP = ( E( I )*M1 )*M2
      U = -TMP*( TMP / U ) - SSHIFT
      IF( U <= SUN ) THEN
         IF( U <= 0.0D+0 ) THEN
            NUM = NUM + 1
            IF( U > -SUN ) &
               U = -SUN
         ELSE
            U = SUN
         END IF
      END IF
      TMP = ( S( I+1 )*M1 )*M2
      U = -TMP*( TMP / U ) - SSHIFT
      IF( U <= SUN ) THEN
         IF( U <= 0.0D+0 ) THEN
            NUM = NUM + 1
            IF( U > -SUN ) &
               U = -SUN
         ELSE
            U = SUN
         END IF
      END IF
   ENDDO
   RETURN
!
!     End of DSVDCT
!
END

