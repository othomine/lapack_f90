!> \brief \b SLCTES
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION SLCTES( ZR, ZI, D )
!
!       .. Scalar Arguments ..
!       REAL               D, ZI, ZR
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLCTES returns .TRUE. if the eigenvalue (ZR/D) + sqrt(-1)*(ZI/D)
!> is to be selected (specifically, in this subroutine, if the real
!> part of the eigenvalue is negative), and otherwise it returns
!> .FALSE..
!>
!> It is used by the test routine SDRGES to test whether the driver
!> routine SGGES successfully sorts eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ZR
!> \verbatim
!>          ZR is REAL
!>          The numerator of the real part of a complex eigenvalue
!>          (ZR/D) + i*(ZI/D).
!> \endverbatim
!>
!> \param[in] ZI
!> \verbatim
!>          ZI is REAL
!>          The numerator of the imaginary part of a complex eigenvalue
!>          (ZR/D) + i*(ZI).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL
!>          The denominator part of a complex eigenvalue
!>          (ZR/D) + i*(ZI/D).
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
   LOGICAL          FUNCTION SLCTES( ZR, ZI, D )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   REAL               D, ZI, ZR
!     ..
!
!  =====================================================================
!     ..
!     .. Executable Statements ..
!
   IF( D == 0.0E+0 ) THEN
      SLCTES = ( ZR < 0.0E+0 )
   ELSE
      SLCTES = ( SIGN( 1.0E+0, ZR ) /= SIGN( 1.0E+0, D ) )
   END IF
!
   RETURN
!
!     End of SLCTES
!
END


