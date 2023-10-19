!> \brief \b ZLCTES
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION ZLCTES( Z, D )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         D, Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLCTES returns .TRUE. if the eigenvalue Z/D is to be selected
!> (specifically, in this subroutine, if the real part of the
!> eigenvalue is negative), and otherwise it returns .FALSE..
!>
!> It is used by the test routine ZDRGES to test whether the driver
!> routine ZGGES successfully sorts eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX*16
!>          The numerator part of a complex eigenvalue Z/D.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX*16
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
!> \ingroup complex16_eig
!
!  =====================================================================
   LOGICAL          FUNCTION ZLCTES( Z, D )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16         D, Z
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION   ZMAX
!     ..
!     .. Executable Statements ..
!
   IF( D == (0.0D+0,0.0D+0) ) THEN
      ZLCTES = ( DBLE( Z ) < 0.0D0 )
   ELSE
      IF( DBLE( Z ) == 0.0D0 .OR. DBLE( D ) == 0.0D0 ) THEN
         ZLCTES = ( SIGN( 1.0D0, DIMAG( Z ) ) /= &
                  SIGN( 1.0D0, DIMAG( D ) ) )
      ELSE IF( DIMAG( Z ) == 0.0D0 .OR. DIMAG( D ) == 0.0D0 ) THEN
         ZLCTES = ( SIGN( 1.0D0, DBLE( Z ) ) /= &
                  SIGN( 1.0D0, DBLE( D ) ) )
      ELSE
         ZMAX = MAX( ABS( DBLE( Z ) ), ABS( DIMAG( Z ) ) )
         ZLCTES = ( ( DBLE( Z ) / ZMAX )*DBLE( D )+ &
                  ( DIMAG( Z ) / ZMAX )*DIMAG( D ) < 0.0D0 )
      END IF
   END IF
!
   RETURN
!
!     End of ZLCTES
!
END



