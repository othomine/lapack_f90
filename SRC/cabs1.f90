!  Authors:
!  ========
!
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup herfs
!
!  =====================================================================
   REAL FUNCTION CABS1( X )
!     .. Parameters ..
   COMPLEX, INTENT(IN) :: X
!     ..
!     .. Statement Function definitions ..
   CABS1 = ABS( REAL( X ) ) + ABS( AIMAG( X ) )
   RETURN
!
!     End of CABS1
!
END
