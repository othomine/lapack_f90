!> \brief \b ZGET36
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET36( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NIN, NINFO
!       DOUBLE PRECISION   RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGET36 tests ZTREXC, a routine for reordering diagonal entries of a
!> matrix in complex Schur form. Thus, ZLAEXC computes a unitary matrix
!> Q such that
!>
!>    Q' * T1 * Q  = T2
!>
!> and where one of the diagonal blocks of T1 (the one at row IFST) has
!> been moved to position ILST.
!>
!> The test code verifies that the residual Q'*T1*Q-T2 is small, that T2
!> is in Schur form, and that the final position of the IFST block is
!> ILST.
!>
!> The test matrices are read from a file with logical unit number NIN.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is DOUBLE PRECISION
!>          Value of the largest test ratio.
!> \endverbatim
!>
!> \param[out] LMAX
!> \verbatim
!>          LMAX is INTEGER
!>          Example number where largest test ratio achieved.
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER
!>          Number of examples where INFO is nonzero.
!> \endverbatim
!>
!> \param[out] KNT
!> \verbatim
!>          KNT is INTEGER
!>          Total number of examples tested.
!> \endverbatim
!>
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          Input logical unit number.
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
   SUBROUTINE ZGET36( RMAX, LMAX, NINFO, KNT, NIN )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX, NIN, NINFO
   DOUBLE PRECISION   RMAX
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            LDT, LWORK
   PARAMETER          ( LDT = 10, LWORK = 2*LDT*LDT )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IFST, ILST, INFO1, INFO2, J, N
   DOUBLE PRECISION   EPS, RES
   COMPLEX*16         CTEMP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   RESULT( 2 ), RWORK( LDT )
   COMPLEX*16         DIAG( LDT ), Q( LDT, LDT ), T1( LDT, LDT ), &
                      T2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZCOPY, ZHST01, ZLACPY, ZLASET, ZTREXC
!     ..
!     .. Executable Statements ..
!
   EPS = DLAMCH( 'P' )
   RMAX = 0.0D0
   LMAX = 0
   KNT = 0
   NINFO = 0
!
!     Read input data until N=0
!
   DO
   READ (NIN,*)N, IFST, ILST
   IF( N == 0 ) RETURN
   KNT = KNT + 1
   DO I = 1, N
      READ (NIN,*) TMP( I,1:N)
   ENDDO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'F', N, N, TMP, LDT, T1, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'F', N, N, TMP, LDT, T2, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   RES = 0.0D0
!
!     Test without accumulating Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, (0.0D+0,0.0D+0), (1.0D0,0.0D0), Q, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZTREXC( 'N', N, T1, LDT, Q, LDT, IFST, ILST, INFO1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZTREXC : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DO I = 1, N
      DO J = 1, N
         IF( I == J .AND. Q( I, J ) /= (1.0D0,0.0D0) ) &
            RES = RES + 1.0D0 / EPS
         IF( I /= J .AND. Q( I, J ) /= (0.0D+0,0.0D+0) ) &
            RES = RES + 1.0D0 / EPS
      ENDDO
   ENDDO
!
!     Test with accumulating Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, (0.0D+0,0.0D+0), (1.0D0,0.0D0), Q, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZTREXC( 'V', N, T2, LDT, Q, LDT, IFST, ILST, INFO2 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZTREXC : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compare T1 with T2
!
   DO I = 1, N
      DO J = 1, N
         IF( T1( I, J ) /= T2( I, J ) ) &
            RES = RES + 1.0D0 / EPS
      ENDDO
   ENDDO
   IF( INFO1 /= 0 .OR. INFO2 /= 0 ) &
      NINFO = NINFO + 1
   IF( INFO1 /= INFO2 ) &
      RES = RES + 1.0D0 / EPS
!
!     Test for successful reordering of T2
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZCOPY( N, TMP, LDT+1, DIAG, 1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( IFST < ILST ) THEN
      DO I = IFST + 1, ILST
         CTEMP = DIAG( I )
         DIAG( I ) = DIAG( I-1 )
         DIAG( I-1 ) = CTEMP
      ENDDO
   ELSE IF( IFST > ILST ) THEN
      DO I = IFST - 1, ILST, -1
         CTEMP = DIAG( I+1 )
         DIAG( I+1 ) = DIAG( I )
         DIAG( I ) = CTEMP
      ENDDO
   END IF
   DO I = 1, N
      IF( T2( I, I ) /= DIAG( I ) ) &
         RES = RES + 1.0D0 / EPS
   ENDDO
!
!     Test for small residual, and orthogonality of Q
!
   CALL ZHST01( N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, &
                RWORK, RESULT )
   RES = RES + RESULT( 1 ) + RESULT( 2 )
!
!     Test for T2 being in Schur form
!
   DO J = 1, N - 1
      DO I = J + 1, N
         IF( T2( I, J ) /= (0.0D+0,0.0D+0) ) &
            RES = RES + 1.0D0 / EPS
         ENDDO
      ENDDO
   IF( RES > RMAX ) THEN
      RMAX = RES
      LMAX = KNT
   END IF
   ENDDO
!
!     End of ZGET36
!
END




