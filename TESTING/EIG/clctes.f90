!> \brief \b CLCTES
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION CLCTES( Z, D )
!
!       .. Scalar Arguments ..
!       COMPLEX            D, Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLCTES returns .TRUE. if the eigenvalue Z/D is to be selected
!> (specifically, in this subroutine, if the real part of the
!> eigenvalue is negative), and otherwise it returns .FALSE..
!>
!> It is used by the test routine CDRGES to test whether the driver
!> routine CGGES successfully sorts eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX
!>          The numerator part of a complex eigenvalue Z/D.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX
!>          The denominator part of a complex eigenvalue Z/D.
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
!> \ingroup complex_eig
!
!  =====================================================================
   LOGICAL          FUNCTION CLCTES( Z, D )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX            D, Z
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   REAL               ZMAX
!     ..
!     .. Executable Statements ..
!
   IF( D == (0.0E+0,0.0E+0) ) THEN
      CLCTES = ( REAL( Z ) < 0.0E+0 )
   ELSE
      IF( REAL( Z ) == 0.0E+0 .OR. REAL( D ) == 0.0E+0 ) THEN
         CLCTES = ( SIGN( 1.0E+0, AIMAG( Z ) ) /= SIGN( 1.0E+0, AIMAG( D ) ) )
      ELSE IF( AIMAG( Z ) == 0.0E+0 .OR. AIMAG( D ) == 0.0E+0 ) THEN
         CLCTES = ( SIGN( 1.0E+0, REAL( Z ) ) /= SIGN( 1.0E+0, REAL( D ) ) )
      ELSE
         ZMAX = MAX( ABS( REAL( Z ) ), ABS( AIMAG( Z ) ) )
         CLCTES = ( ( REAL( Z ) / ZMAX )*REAL( D )+ &
                  ( AIMAG( Z ) / ZMAX )*AIMAG( D ) < 0.0E+0 )
      END IF
   END IF
!
   RETURN
!
!     End of CLCTES
!
END



