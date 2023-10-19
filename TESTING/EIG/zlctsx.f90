!> \brief \b ZLCTSX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION ZLCTSX( ALPHA, BETA )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         ALPHA, BETA
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This function is used to determine what eigenvalues will be
!> selected.  If this is part of the test driver ZDRGSX, do not
!> change the code UNLESS you are testing input examples and not
!> using the built-in examples.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16
!>
!>          parameters to decide whether the pair (ALPHA, BETA) is
!>          selected.
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
   LOGICAL          FUNCTION ZLCTSX( ALPHA, BETA )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   COMPLEX*16         ALPHA, BETA
!     ..
!
!  =====================================================================
!     ..
!     .. Scalars in Common ..
   LOGICAL            FS
   INTEGER            I, M, MPLUSN, N
!     ..
!     .. Common blocks ..
   COMMON             / MN / M, N, MPLUSN, I, FS
!     ..
!     .. Save statement ..
   SAVE
!     ..
!     .. Executable Statements ..
!
   IF( FS ) THEN
      I = I + 1
      ZLCTSX = ( I <= M )
      IF( I == MPLUSN ) THEN
         FS = .FALSE.
         I = 0
      END IF
   ELSE
      I = I + 1
      ZLCTSX = ( I <= N )
      IF( I == MPLUSN ) THEN
         FS = .TRUE.
         I = 0
      END IF
   END IF
!
!      IF( BETA == CZERO ) THEN
!         ZLCTSX = ( DBLE( ALPHA ) > ZERO )
!      ELSE
!         ZLCTSX = ( DBLE( ALPHA/BETA ) > ZERO )
!      END IF
!
   RETURN
!
!     End of ZLCTSX
!
END



