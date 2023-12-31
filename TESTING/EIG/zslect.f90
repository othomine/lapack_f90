!> \brief \b ZSLECT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION ZSLECT( Z )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSLECT returns .TRUE. if the eigenvalue Z is to be selected,
!> otherwise it returns .FALSE.
!> It is used by ZCHK41 to test if ZGEES successfully sorts eigenvalues,
!> and by ZCHK43 to test if ZGEESX successfully sorts eigenvalues.
!>
!> The common block /SSLCT/ controls how eigenvalues are selected.
!> If SELOPT = 0, then ZSLECT return .TRUE. when real(Z) is less than
!> zero, and .FALSE. otherwise.
!> If SELOPT is at least 1, ZSLECT returns SELVAL(SELOPT) and adds 1
!> to SELOPT, cycling back to 1 at SELMAX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX*16
!>          The eigenvalue Z.
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
   LOGICAL          FUNCTION ZSLECT( Z )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16         Z
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I
   DOUBLE PRECISION   RMIN, X
!     ..
!     .. Scalars in Common ..
   INTEGER            SELDIM, SELOPT
!     ..
!     .. Arrays in Common ..
   LOGICAL            SELVAL( 20 )
   DOUBLE PRECISION   SELWI( 20 ), SELWR( 20 )
!     ..
!     .. Common blocks ..
   COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
!     ..
!     .. Executable Statements ..
!
   IF( SELOPT == 0 ) THEN
      ZSLECT = ( DBLE( Z ) < 0.0D0 )
   ELSE
      RMIN = ABS( Z-DCMPLX( SELWR( 1 ), SELWI( 1 ) ) )
      ZSLECT = SELVAL( 1 )
      DO I = 2, SELDIM
         X = ABS( Z-DCMPLX( SELWR( I ), SELWI( I ) ) )
         IF( X <= RMIN ) THEN
            RMIN = X
            ZSLECT = SELVAL( I )
         END IF
      ENDDO
   END IF
   RETURN
!
!     End of ZSLECT
!
END



