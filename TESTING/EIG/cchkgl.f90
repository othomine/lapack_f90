!> \brief \b CCHKGL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKGL( NIN, NOUT )
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
!> CCHKGL tests CGGBAL, a routine for balancing a matrix pair (A, B).
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
!> \ingroup complex_eig
!
!  =====================================================================
   SUBROUTINE CCHKGL( NIN, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NIN, NOUT
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            LDA, LDB, LWORK
   PARAMETER          ( LDA = 20, LDB = 20, LWORK = 6*LDA )
   REAL               ZERO
   PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IHI, IHIIN, ILO, ILOIN, INFO, J, KNT, N, &
                      NINFO
   REAL               ANORM, BNORM, EPS, RMAX, VMAX
!     ..
!     .. Local Arrays ..
   INTEGER            LMAX( 3 )
   REAL               LSCALE( LDA ), LSCLIN( LDA ), RSCALE( LDA ), &
                      RSCLIN( LDA ), WORK( LWORK )
   COMPLEX            A( LDA, LDA ), AIN( LDA, LDA ), B( LDB, LDB ), &
                      BIN( LDB, LDB )
!     ..
!     .. External Functions ..
   REAL               CLANGE, SLAMCH
   EXTERNAL           CLANGE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGGBAL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
   LMAX( 1 ) = 0
   LMAX( 2 ) = 0
   LMAX( 3 ) = 0
   NINFO = 0
   KNT = 0
   RMAX = ZERO
!
   EPS = SLAMCH( 'Precision' )
!
10 CONTINUE
!
   READ( NIN, FMT = * )N
   IF( N == 0 ) &
      GO TO 90
   DO I = 1, N
      READ( NIN, FMT = * )( A( I, J ), J = 1, N )
   ENDDO
!
   DO I = 1, N
      READ( NIN, FMT = * )( B( I, J ), J = 1, N )
   ENDDO
!
   READ( NIN, FMT = * )ILOIN, IHIIN
   DO I = 1, N
      READ( NIN, FMT = * )( AIN( I, J ), J = 1, N )
   ENDDO
   DO I = 1, N
      READ( NIN, FMT = * )( BIN( I, J ), J = 1, N )
   ENDDO
!
   READ( NIN, FMT = * )( LSCLIN( I ), I = 1, N )
   READ( NIN, FMT = * )( RSCLIN( I ), I = 1, N )
!
   ANORM = CLANGE( 'M', N, N, A, LDA, WORK )
   BNORM = CLANGE( 'M', N, N, B, LDB, WORK )
!
   KNT = KNT + 1
!
   CALL CGGBAL( 'B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, &
                WORK, INFO )
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
   VMAX = ZERO
   DO I = 1, N
      DO J = 1, N
         VMAX = MAX( VMAX, ABS( A( I, J )-AIN( I, J ) ) )
         VMAX = MAX( VMAX, ABS( B( I, J )-BIN( I, J ) ) )
      ENDDO
   ENDDO
!
   DO I = 1, N
      VMAX = MAX( VMAX, ABS( LSCALE( I )-LSCLIN( I ) ) )
      VMAX = MAX( VMAX, ABS( RSCALE( I )-RSCLIN( I ) ) )
   ENDDO
!
   VMAX = VMAX / ( EPS*MAX( ANORM, BNORM ) )
!
   IF( VMAX > RMAX ) THEN
      LMAX( 3 ) = KNT
      RMAX = VMAX
   END IF
!
   GO TO 10
!
90 CONTINUE
!
   WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( ' .. test output of CGGBAL .. ' )
!
   WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( ' ratio of largest test error              = ', E12.3 )
   WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( ' example number where info is not zero    = ', I4 )
   WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( ' example number where ILO or IHI is wrong = ', I4 )
   WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( ' example number having largest error      = ', I4 )
   WRITE( NOUT, FMT = 9994 )NINFO
 9994 FORMAT( ' number of examples where info is not 0   = ', I4 )
   WRITE( NOUT, FMT = 9993 )KNT
 9993 FORMAT( ' total number of examples tested          = ', I4 )
!
   RETURN
!
!     End of CCHKGL
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
