!> \brief \b SCHKGK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKGK( NIN, NOUT )
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
!> SCHKGK tests SGGBAK, a routine for backward balancing  of
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SCHKGK( NIN, NOUT )
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
   REAL               ANORM, BNORM, EPS, RMAX, VMAX
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            LMAX( 4 )
   REAL               A( LDA, LDA ), AF( LDA, LDA ), B( LDB, LDB ), &
                      BF( LDB, LDB ), E( LDE, LDE ), F( LDF, LDF ), &
                      LSCALE( LDA ), RSCALE( LDA ), VL( LDVL, LDVL ), &
                      VLF( LDVL, LDVL ), VR( LDVR, LDVR ), &
                      VRF( LDVR, LDVR ), WORK( LDWORK, LDWORK )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SLANGE
   EXTERNAL           SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMM, SGGBAK, SGGBAL, SLACPY
!     ..
!     .. Executable Statements ..
!
!     Initialization
!
   LMAX( 1:4 ) = 0
   NINFO = 0
   KNT = 0
   RMAX = 0.0E+0
!
   EPS = SLAMCH( 'Precision' )
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
   ANORM = SLANGE( 'M', N, N, A, LDA, WORK )
   BNORM = SLANGE( 'M', N, N, B, LDB, WORK )
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SLACPY( 'FULL', N, N, A, LDA, AF, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SLACPY( 'FULL', N, N, B, LDB, BF, LDB )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGGBAL( 'B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, &
                WORK, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGGBAL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 1 ) = KNT
   END IF
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SLACPY( 'FULL', N, M, VL, LDVL, VLF, LDVL )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SLACPY( 'FULL', N, M, VR, LDVR, VRF, LDVR )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGGBAK( 'B', 'L', N, ILO, IHI, LSCALE, RSCALE, M, VL, LDVL, &
                INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGGBAK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 2 ) = KNT
   END IF
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGGBAK( 'B', 'R', N, ILO, IHI, LSCALE, RSCALE, M, VR, LDVR, &
                INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGGBAK : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( INFO /= 0 ) THEN
      NINFO = NINFO + 1
      LMAX( 3 ) = KNT
   END IF
!
!     Test of SGGBAK
!
!     Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
!     where tilde(A) denotes the transformed matrix.
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'N', 'N', N, M, N, 1.0E+0, AF, LDA, VR, LDVR, 0.0E+0, WORK, &
               LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'T', 'N', M, M, N, 1.0E+0, VL, LDVL, WORK, LDWORK, 0.0E+0, &
               E, LDE )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'N', 'N', N, M, N, 1.0E+0, A, LDA, VRF, LDVR, 0.0E+0, WORK, &
               LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'T', 'N', M, M, N, 1.0E+0, VLF, LDVL, WORK, LDWORK, 0.0E+0, &
               F, LDF )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   VMAX = MAXVAL(ABS(E(1:M,1:M)-F(1:M,1:M))) / ( EPS*MAX( ANORM, BNORM ) )
   IF( VMAX > RMAX ) THEN
      LMAX( 4 ) = KNT
      RMAX = VMAX
   END IF
!
!     Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'N', 'N', N, M, N, 1.0E+0, BF, LDB, VR, LDVR, 0.0E+0, WORK, &
               LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'T', 'N', M, M, N, 1.0E+0, VL, LDVL, WORK, LDWORK, 0.0E+0, &
               E, LDE )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'N', 'N', N, M, N, 1.0E+0, B, LDB, VRF, LDVR, 0.0E+0, WORK, &
               LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEMM( 'T', 'N', M, M, N, 1.0E+0, VLF, LDVL, WORK, LDWORK, 0.0E+0, &
               F, LDF )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEMM : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
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
 9999 FORMAT( 1X, '.. test output of SGGBAK .. ' )
!
   WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( ' value of largest test error                  =', E12.3 )
   WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( ' example number where SGGBAL info is not 0    =', I4 )
   WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( ' example number where SGGBAK(L) info is not 0 =', I4 )
   WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( ' example number where SGGBAK(R) info is not 0 =', I4 )
   WRITE( NOUT, FMT = 9994 )LMAX( 4 )
 9994 FORMAT( ' example number having largest error          =', I4 )
   WRITE( NOUT, FMT = 9992 )NINFO
 9992 FORMAT( ' number of examples where info is not 0       =', I4 )
   WRITE( NOUT, FMT = 9991 )KNT
 9991 FORMAT( ' total number of examples tested              =', I4 )
!
   RETURN
!
!     End of SCHKGK
!
END



