!> \brief \b SSLECT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION SSLECT( ZR, ZI )
!
!       .. Scalar Arguments ..
!       REAL               ZI, ZR
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSLECT returns .TRUE. if the eigenvalue ZR+sqrt(-1)*ZI is to be
!> selected, and otherwise it returns .FALSE.
!> It is used by SCHK41 to test if SGEES successfully sorts eigenvalues,
!> and by SCHK43 to test if SGEESX successfully sorts eigenvalues.
!>
!> The common block /SSLCT/ controls how eigenvalues are selected.
!> If SELOPT = 0, then SSLECT return .TRUE. when ZR is less than zero,
!> and .FALSE. otherwise.
!> If SELOPT is at least 1, SSLECT returns SELVAL(SELOPT) and adds 1
!> to SELOPT, cycling back to 1 at SELMAX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ZR
!> \verbatim
!>          ZR is REAL
!>          The real part of a complex eigenvalue ZR + i*ZI.
!> \endverbatim
!>
!> \param[in] ZI
!> \verbatim
!>          ZI is REAL
!>          The imaginary part of a complex eigenvalue ZR + i*ZI.
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
   LOGICAL          FUNCTION SSLECT( ZR, ZI )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   REAL               ZI, ZR
!     ..
!
!  =====================================================================
!
!     .. Arrays in Common ..
   LOGICAL            SELVAL( 20 )
   REAL               SELWI( 20 ), SELWR( 20 )
!     ..
!     .. Scalars in Common ..
   INTEGER            SELDIM, SELOPT
!     ..
!     .. Common blocks ..
   COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
!     ..
!     .. Local Scalars ..
   INTEGER            I
   REAL               RMIN, X
!     ..
!     .. External Functions ..
   REAL               SLAPY2
   EXTERNAL           SLAPY2
!     ..
!     .. Executable Statements ..
!
   IF( SELOPT == 0 ) THEN
      SSLECT = ( ZR < 0.0E+0 )
   ELSE
      RMIN = SLAPY2( ZR-SELWR( 1 ), ZI-SELWI( 1 ) )
      SSLECT = SELVAL( 1 )
      DO I = 2, SELDIM
         X = SLAPY2( ZR-SELWR( I ), ZI-SELWI( I ) )
         IF( X <= RMIN ) THEN
            RMIN = X
            SSLECT = SELVAL( I )
         END IF
      ENDDO
   END IF
   RETURN
!
!     End of SSLECT
!
END



