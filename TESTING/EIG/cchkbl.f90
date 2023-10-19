!> \brief \b CCHKBL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKBL( NIN, NOUT )
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
!> CCHKBL tests CGEBAL, a routine for balancing a general complex
!> matrix and isolating some of its eigenvalues.
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
!> \ingroup complex_eig
!
!  =====================================================================
   SUBROUTINE CCHKBL( NIN, NOUT )
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
   INTEGER            LDA
   PARAMETER          ( LDA = 20 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IHI, IHIIN, ILO, ILOIN, INFO, J, KNT, N, &
                      NINFO
   REAL               ANORM, MEPS, RMAX, SFMIN, TEMP, VMAX
   COMPLEX            CDUM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            LMAX( 3 )
   REAL               DUMMY( 1 ), SCALE( LDA ), SCALIN( LDA )
   COMPLEX            A( LDA, LDA ), AIN( LDA, LDA )
!     ..
!     .. External Functions ..
   REAL               CLANGE, SLAMCH
   EXTERNAL           CLANGE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEBAL
!     ..
!     .. Statement Functions ..
   REAL               CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
   LMAX( 1 ) = 0
   LMAX( 2 ) = 0
   LMAX( 3 ) = 0
   NINFO = 0
   KNT = 0
   RMAX = 0.0E+0
   VMAX = 0.0E+0
   SFMIN = SLAMCH( 'S' )
   MEPS = SLAMCH( 'E' )
!
   DO
!
   READ( NIN, FMT = * )N
   IF( N == 0 ) EXIT
   DO I = 1, N
      READ( NIN, FMT = * ) A( I,1:N)
   ENDDO
!
   READ( NIN, FMT = * )ILOIN, IHIIN
   DO I = 1, N
      READ( NIN, FMT = * ) AIN( I,1:N)
   ENDDO
   READ( NIN, FMT = * ) SCALIN(1:N)
!
   ANORM = CLANGE( 'M', N, N, A, LDA, DUMMY )
   KNT = KNT + 1
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL CGEBAL( 'B', N, A, LDA, ILO, IHI, SCALE, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : CGEBAL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 1 ) = KNT
   END IF
!
   IF( ILO /= ILOIN .OR. IHI /= IHIIN ) THEN
      NINFO = NINFO + 1
      LMAX( 2 ) = KNT
   END IF
!
   DO I = 1, N
      DO J = 1, N
         TEMP = MAX( CABS1( A( I, J ) ), CABS1( AIN( I, J ) ) , SFMIN )
         VMAX = MAX( VMAX, CABS1( A( I, J )-AIN( I, J ) ) / TEMP )
      ENDDO
   ENDDO
!
   DO I = 1, N
      TEMP = MAX( SCALE( I ), SCALIN( I ) , SFMIN )
      VMAX = MAX( VMAX, ABS( SCALE( I )-SCALIN( I ) ) / TEMP )
   ENDDO
!
   IF( VMAX > RMAX ) THEN
      LMAX( 3 ) = KNT
      RMAX = VMAX
   END IF
!
   ENDDO
!
   WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of CGEBAL .. ' )
!
   WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( 1X, 'value of largest test error            = ', E12.3 )
   WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( 1X, 'example number where info is not zero  = ', I4 )
   WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( 1X, 'example number where ILO or IHI wrong  = ', I4 )
   WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( 1X, 'example number having largest error    = ', I4 )
   WRITE( NOUT, FMT = 9994 )NINFO
 9994 FORMAT( 1X, 'number of examples where info is not 0 = ', I4 )
   WRITE( NOUT, FMT = 9993 )KNT
 9993 FORMAT( 1X, 'total number of examples tested        = ', I4 )
!
   RETURN
!
!     End of CCHKBL
!
END




