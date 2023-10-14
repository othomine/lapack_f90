!> \brief \b DCHKBK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKBK( NIN, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NIN, NOUT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKBK tests DGEBAK, a routine for backward transformation of
!> the computed right or left eigenvectors if the original matrix
!> was preprocessed by balance subroutine DGEBAL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          The logical unit number for input.  NIN > 0.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The logical unit number for output.  NOUT > 0.
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
   SUBROUTINE DCHKBK( NIN, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NIN, NOUT
!     ..
!
! ======================================================================
!
!     .. Parameters ..
   INTEGER            LDE
   PARAMETER          ( LDE = 20 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IHI, ILO, INFO, J, KNT, N, NINFO
   DOUBLE PRECISION   EPS, RMAX, SAFMIN, VMAX, X
!     ..
!     .. Local Arrays ..
   INTEGER            LMAX( 2 )
   DOUBLE PRECISION   E( LDE, LDE ), EIN( LDE, LDE ), SCALE( LDE )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEBAK
!     ..
!     .. Executable Statements ..
!
   LMAX( 1 ) = 0
   LMAX( 2 ) = 0
   NINFO = 0
   KNT = 0
   RMAX = 0.0D0
   EPS = DLAMCH( 'E' )
   SAFMIN = DLAMCH( 'S' )
!
   DO
!
   READ(NIN,*) N, ILO, IHI
   IF( N == 0 ) EXIT
!
   READ(NIN,*) SCALE(1:N)
   DO I = 1, N
      READ(NIN,*) E(I,1:N)
   ENDDO
!
   DO I = 1, N
      READ(NIN,*) EIN(I,1:N)
   ENDDO
!
   KNT = KNT + 1
   CALL DGEBAK( 'B', 'R', N, ILO, IHI, SCALE, N, E, LDE, INFO )
!
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 1 ) = KNT
   END IF
!
   VMAX = 0.0D0
   DO I = 1, N
      DO J = 1, N
         X = ABS( E( I, J )-EIN( I, J ) ) / EPS
         IF( ABS( E( I, J ) ) > SAFMIN ) X = X / ABS( E( I, J ) )
         VMAX = MAX( VMAX, X )
      ENDDO
   ENDDO
!
   IF( VMAX > RMAX ) THEN
      LMAX( 2 ) = KNT
      RMAX = VMAX
   END IF
!
   ENDDO
!
!
   WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of DGEBAK .. ' )
!
   WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( 1X, 'value of largest test error             = ', D12.3 )
   WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( 1X, 'example number where info is not zero   = ', I4 )
   WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( 1X, 'example number having largest error     = ', I4 )
   WRITE( NOUT, FMT = 9995 )NINFO
 9995 FORMAT( 1X, 'number of examples where info is not 0  = ', I4 )
   WRITE( NOUT, FMT = 9994 )KNT
 9994 FORMAT( 1X, 'total number of examples tested         = ', I4 )
!
   RETURN
!
!     End of DCHKBK
!
END

