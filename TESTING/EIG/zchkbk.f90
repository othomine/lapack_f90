!> \brief \b ZCHKBK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKBK( NIN, NOUT )
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
!> ZCHKBK tests ZGEBAK, a routine for backward transformation of
!> the computed right or left eigenvectors if the original matrix
!> was preprocessed by balance subroutine ZGEBAL.
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup complex16_eig
!
!  =====================================================================
   SUBROUTINE ZCHKBK( NIN, NOUT )
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
   COMPLEX*16         CDUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            LMAX( 2 )
   DOUBLE PRECISION   SCALE( LDE )
   COMPLEX*16         E( LDE, LDE ), EIN( LDE, LDE )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGEBAK
!     ..
!     .. Statement Functions ..
   DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
   LMAX( 1 ) = 0
   LMAX( 2 ) = 0
   NINFO = 0
   KNT = 0
   RMAX = 0.0D+0
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
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZGEBAK( 'B', 'R', N, ILO, IHI, SCALE, N, E, LDE, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZGEBAK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 1 ) = KNT
   END IF
!
   VMAX = 0.0D+0
   DO I = 1, N
      DO J = 1, N
         X = CABS1( E( I, J )-EIN( I, J ) ) / EPS
         IF( CABS1( E( I, J ) ) > SAFMIN ) X = X / CABS1( E( I, J ) )
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
 9999 FORMAT( 1X, '.. test output of ZGEBAK .. ' )
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
!     End of ZCHKBK
!
END




