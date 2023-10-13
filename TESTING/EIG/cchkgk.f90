!> \brief \b CCHKGK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKGK( NIN, NOUT )
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
!> CCHKGK tests CGGBAK, a routine for backward balancing  of
!> a matrix pair (A, B).
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
   SUBROUTINE CCHKGK( NIN, NOUT )
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
   INTEGER            LDA, LDB, LDVL, LDVR
   PARAMETER          ( LDA = 50, LDB = 50, LDVL = 50, LDVR = 50 )
   INTEGER            LDE, LDF, LDWORK, LRWORK
   PARAMETER          ( LDE = 50, LDF = 50, LDWORK = 50, &
                      LRWORK = 6*50 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IHI, ILO, INFO, J, KNT, M, N, NINFO
   REAL               ANORM, BNORM, EPS, RMAX, VMAX
   COMPLEX            CDUM
!     ..
!     .. Local Arrays ..
   INTEGER            LMAX( 4 )
   REAL               LSCALE( LDA ), RSCALE( LDA ), RWORK( LRWORK )
   COMPLEX            A( LDA, LDA ), AF( LDA, LDA ), B( LDB, LDB ), &
                      BF( LDB, LDB ), E( LDE, LDE ), F( LDF, LDF ), &
                      VL( LDVL, LDVL ), VLF( LDVL, LDVL ), &
                      VR( LDVR, LDVR ), VRF( LDVR, LDVR ), &
                      WORK( LDWORK, LDWORK )
!     ..
!     .. External Functions ..
   REAL               CLANGE, SLAMCH
   EXTERNAL           CLANGE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CGGBAK, CGGBAL, CLACPY
!     ..
!     .. Statement Functions ..
   REAL               CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
   LMAX(1:4) = 0
   NINFO = 0
   KNT = 0
   RMAX = 0.0E+0
!
   EPS = SLAMCH( 'Precision' )
!
   DO
   READ( NIN, FMT = * )N, M
   IF( N == 0 ) GO TO 100
   DO I = 1, N
      READ( NIN, FMT = * ) A(I,1:N)
   ENDDO
   DO I = 1, N
      READ( NIN, FMT = * ) B(I,1:N)
   ENDDO
   DO I = 1, N
     READ( NIN, FMT = * ) VL( I,1:M)
   ENDDO
   DO I = 1, N
     READ( NIN, FMT = * ) VR( I,1:M)
   ENDDO
!
   KNT = KNT + 1
!
   ANORM = CLANGE( 'M', N, N, A, LDA, RWORK )
   BNORM = CLANGE( 'M', N, N, B, LDB, RWORK )
!
   CALL CLACPY( 'FULL', N, N, A, LDA, AF, LDA )
   CALL CLACPY( 'FULL', N, N, B, LDB, BF, LDB )
!
   CALL CGGBAL( 'B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, &
                RWORK, INFO )
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 1 ) = KNT
   END IF
!
   CALL CLACPY( 'FULL', N, M, VL, LDVL, VLF, LDVL )
   CALL CLACPY( 'FULL', N, M, VR, LDVR, VRF, LDVR )
!
   CALL CGGBAK( 'B', 'L', N, ILO, IHI, LSCALE, RSCALE, M, VL, LDVL, &
                INFO )
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 2 ) = KNT
   END IF
!
   CALL CGGBAK( 'B', 'R', N, ILO, IHI, LSCALE, RSCALE, M, VR, LDVR, &
                INFO )
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 3 ) = KNT
   END IF
!
!     Test of CGGBAK
!
!     Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
!     where tilde(A) denotes the transformed matrix.
!
   CALL CGEMM( 'N', 'N', N, M, N, (1.0E+0,0.0E+0), AF, LDA, VR, LDVR, (0.0E+0,0.0E+0), &
               WORK, LDWORK )
   CALL CGEMM( 'C', 'N', M, M, N, (1.0E+0,0.0E+0), VL, LDVL, WORK, LDWORK, &
               (0.0E+0,0.0E+0), E, LDE )
!
   CALL CGEMM( 'N', 'N', N, M, N, (1.0E+0,0.0E+0), A, LDA, VRF, LDVR, (0.0E+0,0.0E+0), &
               WORK, LDWORK )
   CALL CGEMM( 'C', 'N', M, M, N, (1.0E+0,0.0E+0), VLF, LDVL, WORK, LDWORK, &
               (0.0E+0,0.0E+0), F, LDF )
!
   VMAX = 0.0E+0
   DO I = 1, M
      DO J = 1, M
         VMAX = MAX(VMAX,CABS1(E(I,J)-F(I,J)))
      ENDDO
   ENDDO
   VMAX = VMAX/(EPS*MAX(ANORM,BNORM))
   IF( VMAX > RMAX ) THEN
      LMAX( 4 ) = KNT
      RMAX = VMAX
   END IF
!
!     Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR
!
   CALL CGEMM( 'N', 'N', N, M, N, (1.0E+0,0.0E+0), BF, LDB, VR, LDVR, (0.0E+0,0.0E+0), &
               WORK, LDWORK )
   CALL CGEMM( 'C', 'N', M, M, N, (1.0E+0,0.0E+0), VL, LDVL, WORK, LDWORK, &
               (0.0E+0,0.0E+0), E, LDE )
!
   CALL CGEMM( 'n', 'n', N, M, N, (1.0E+0,0.0E+0), B, LDB, VRF, LDVR, (0.0E+0,0.0E+0), &
               WORK, LDWORK )
   CALL CGEMM( 'C', 'N', M, M, N, (1.0E+0,0.0E+0), VLF, LDVL, WORK, LDWORK, &
               (0.0E+0,0.0E+0), F, LDF )
!
   VMAX = 0.0E+0
   DO I = 1, M
      DO J = 1, M
         VMAX = MAX(VMAX,CABS1(E(I,J)-F(I,J)))
      ENDDO
   ENDDO
   VMAX = VMAX/(EPS*MAX(ANORM,BNORM))
   IF( VMAX > RMAX ) THEN
      LMAX( 4 ) = KNT
      RMAX = VMAX
   END IF
!
   ENDDO
!
  100 CONTINUE
!
   WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of CGGBAK .. ' )
!
   WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( ' value of largest test error                  =', E12.3 )
   WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( ' example number where CGGBAL info is not 0    =', I4 )
   WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( ' example number where CGGBAK(L) info is not 0 =', I4 )
   WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( ' example number where CGGBAK(R) info is not 0 =', I4 )
   WRITE( NOUT, FMT = 9994 )LMAX( 4 )
 9994 FORMAT( ' example number having largest error          =', I4 )
   WRITE( NOUT, FMT = 9992 )NINFO
 9992 FORMAT( ' number of examples where info is not 0       =', I4 )
   WRITE( NOUT, FMT = 9991 )KNT
 9991 FORMAT( ' total number of examples tested              =', I4 )
!
   RETURN
!
!     End of CCHKGK
!
END

