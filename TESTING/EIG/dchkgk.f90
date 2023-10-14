!> \brief \b DCHKGK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKGK( NIN, NOUT )
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
!> DCHKGK tests DGGBAK, a routine for backward balancing  of
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
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DCHKGK( NIN, NOUT )
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
   INTEGER            LDE, LDF, LDWORK
   PARAMETER          ( LDE = 50, LDF = 50, LDWORK = 50 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IHI, ILO, INFO, J, KNT, M, N, NINFO
   DOUBLE PRECISION   ANORM, BNORM, EPS, RMAX, VMAX
!     ..
!     .. Local Arrays ..
   INTEGER            LMAX( 4 )
   DOUBLE PRECISION   A( LDA, LDA ), AF( LDA, LDA ), B( LDB, LDB ), &
                      BF( LDB, LDB ), E( LDE, LDE ), F( LDF, LDF ), &
                      LSCALE( LDA ), RSCALE( LDA ), VL( LDVL, LDVL ), &
                      VLF( LDVL, LDVL ), VR( LDVR, LDVR ), &
                      VRF( LDVR, LDVR ), WORK( LDWORK, LDWORK )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLANGE
   EXTERNAL           DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMM, DGGBAK, DGGBAL, DLACPY
!     ..
!     .. Executable Statements ..
!
   LMAX(1:4) = 0
   NINFO = 0
   KNT = 0
   RMAX = 0.0D+0
!
   EPS = DLAMCH( 'Precision' )
!
   DO
   READ(NIN,*) N, M
   IF( N == 0 ) GO TO 100
   DO I = 1, N
      READ(NIN,*) A(I,1:N)
   ENDDO
   DO I = 1, N
      READ(NIN,*) B(I,1:N)
   ENDDO
   DO I = 1, N
     READ(NIN,*) VL( I,1:M)
   ENDDO
   DO I = 1, N
     READ(NIN,*) VR( I,1:M)
   ENDDO
!
   KNT = KNT + 1
!
   ANORM = DLANGE( 'M', N, N, A, LDA, WORK )
   BNORM = DLANGE( 'M', N, N, B, LDB, WORK )
!
   CALL DLACPY( 'FULL', N, N, A, LDA, AF, LDA )
   CALL DLACPY( 'FULL', N, N, B, LDB, BF, LDB )
!
   CALL DGGBAL( 'B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, &
                WORK, INFO )
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 1 ) = KNT
   END IF
!
   CALL DLACPY( 'FULL', N, M, VL, LDVL, VLF, LDVL )
   CALL DLACPY( 'FULL', N, M, VR, LDVR, VRF, LDVR )
!
   CALL DGGBAK( 'B', 'L', N, ILO, IHI, LSCALE, RSCALE, M, VL, LDVL, &
                INFO )
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 2 ) = KNT
   END IF
!
   CALL DGGBAK( 'B', 'R', N, ILO, IHI, LSCALE, RSCALE, M, VR, LDVR, &
                INFO )
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 3 ) = KNT
   END IF
!
!     Test of DGGBAK
!
!     Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
!     where tilde(A) denotes the transformed matrix.
!
   CALL DGEMM( 'N', 'N', N, M, N, 1.0D+0, AF, LDA, VR, LDVR, 0.0D+0, WORK, &
               LDWORK )
   CALL DGEMM( 'T', 'N', M, M, N, 1.0D+0, VL, LDVL, WORK, LDWORK, 0.0D+0, &
               E, LDE )
!
   CALL DGEMM( 'N', 'N', N, M, N, 1.0D+0, A, LDA, VRF, LDVR, 0.0D+0, WORK, &
               LDWORK )
   CALL DGEMM( 'T', 'N', M, M, N, 1.0D+0, VLF, LDVL, WORK, LDWORK, 0.0D+0, &
               F, LDF )
!
   VMAX = MAXVAL(ABS(E(1:M,1:M)-F(1:M,1:M))) / ( EPS*MAX( ANORM, BNORM ) )
   IF( VMAX > RMAX ) THEN
      LMAX( 4 ) = KNT
      RMAX = VMAX
   END IF
!
!     Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR
!
   CALL DGEMM( 'N', 'N', N, M, N, 1.0D+0, BF, LDB, VR, LDVR, 0.0D+0, WORK, &
               LDWORK )
   CALL DGEMM( 'T', 'N', M, M, N, 1.0D+0, VL, LDVL, WORK, LDWORK, 0.0D+0, &
               E, LDE )
!
   CALL DGEMM( 'N', 'N', N, M, N, 1.0D+0, B, LDB, VRF, LDVR, 0.0D+0, WORK, &
               LDWORK )
   CALL DGEMM( 'T', 'N', M, M, N, 1.0D+0, VLF, LDVL, WORK, LDWORK, 0.0D+0, &
               F, LDF )
!
   VMAX = MAXVAL(ABS(E(1:M,1:M)-F(1:M,1:M))) / ( EPS*MAX( ANORM, BNORM ) )
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
 9999 FORMAT( 1X, '.. test output of DGGBAK .. ' )
!
   WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( ' value of largest test error                  =', D12.3 )
   WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( ' example number where DGGBAL info is not 0    =', I4 )
   WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( ' example number where DGGBAK(L) info is not 0 =', I4 )
   WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( ' example number where DGGBAK(R) info is not 0 =', I4 )
   WRITE( NOUT, FMT = 9994 )LMAX( 4 )
 9994 FORMAT( ' example number having largest error          =', I4 )
   WRITE( NOUT, FMT = 9993 )NINFO
 9993 FORMAT( ' number of examples where info is not 0       =', I4 )
   WRITE( NOUT, FMT = 9992 )KNT
 9992 FORMAT( ' total number of examples tested              =', I4 )
!
   RETURN
!
!     End of DCHKGK
!
END
