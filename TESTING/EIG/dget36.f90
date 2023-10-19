!> \brief \b DGET36
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET36( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NIN
!       DOUBLE PRECISION   RMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            NINFO( 3 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET36 tests DTREXC, a routine for moving blocks (either 1 by 1 or
!> 2 by 2) on the diagonal of a matrix in real Schur form.  Thus, DLAEXC
!> computes an orthogonal matrix Q such that
!>
!>    Q' * T1 * Q  = T2
!>
!> and where one of the diagonal blocks of T1 (the one at row IFST) has
!> been moved to position ILST.
!>
!> The test code verifies that the residual Q'*T1*Q-T2 is small, that T2
!> is in Schur form, and that the final position of the IFST block is
!> ILST (within +-1).
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
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(J) is the number of examples where INFO=J.
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
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DGET36( RMAX, LMAX, NINFO, KNT, NIN )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX, NIN
   DOUBLE PRECISION   RMAX
!     ..
!     .. Array Arguments ..
   INTEGER            NINFO( 3 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            LDT, LWORK
   PARAMETER          ( LDT = 10, LWORK = 2*LDT*LDT )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1, &
                      ILST2, ILSTSV, INFO1, INFO2, J, LOC, N
   DOUBLE PRECISION   EPS, RES
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   Q( LDT, LDT ), RESULT( 2 ), T1( LDT, LDT ), &
                      T2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DHST01, DLACPY, DLASET, DTREXC
!     ..
!     .. Executable Statements ..
!
   EPS = DLAMCH( 'P' )
   RMAX = 0.0D0
   LMAX = 0
   KNT = 0
   NINFO( 1:3 ) = 0
!
!     Read input data until N=0
!
   DO
   READ(NIN,*) N, IFST, ILST
   IF( N == 0 ) RETURN
   KNT = KNT + 1
   DO I = 1, N
      READ(NIN,*) TMP( I,1:N)
   ENDDO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'F', N, N, TMP, LDT, T1, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'F', N, N, TMP, LDT, T2, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IFSTSV = IFST
   ILSTSV = ILST
   IFST1 = IFST
   ILST1 = ILST
   IFST2 = IFST
   ILST2 = ILST
   RES = 0.0D0
!
!     Test without accumulating Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', N, N, 0.0D0, 1.0D0, Q, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DTREXC( 'N', N, T1, LDT, Q, LDT, IFST1, ILST1, WORK, INFO1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DTREXC : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DO I = 1, N
      DO J = 1, N
         IF( I == J .AND. Q( I, J ) /= 1.0D0 ) RES = RES + 1.0D0 / EPS
         IF( I /= J .AND. Q( I, J ) /= 0.0D0 ) RES = RES + 1.0D0 / EPS
      ENDDO
   ENDDO
!
!     Test with accumulating Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', N, N, 0.0D0, 1.0D0, Q, LDT )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DTREXC( 'V', N, T2, LDT, Q, LDT, IFST2, ILST2, WORK, INFO2 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DTREXC : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compare T1 with T2
!
   DO I = 1, N
      DO J = 1, N
         IF( T1( I, J ) /= T2( I, J ) ) RES = RES + 1.0D0 / EPS
      ENDDO
   ENDDO
   IF( IFST1 /= IFST2 ) RES = RES + 1.0D0 / EPS
   IF( ILST1 /= ILST2 ) RES = RES + 1.0D0 / EPS
   IF( INFO1 /= INFO2 ) RES = RES + 1.0D0 / EPS
!
!     Test for successful reordering of T2
!
   IF( INFO2 /= 0 ) THEN
      NINFO( INFO2 ) = NINFO( INFO2 ) + 1
   ELSE
      IF( ABS( IFST2-IFSTSV ) > 1 ) RES = RES + 1.0D0 / EPS
      IF( ABS( ILST2-ILSTSV ) > 1 ) RES = RES + 1.0D0 / EPS
   END IF
!
!     Test for small residual, and orthogonality of Q
!
   CALL DHST01( N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RESULT )
   RES = RES + RESULT( 1 ) + RESULT( 2 )
!
!     Test for T2 being in Schur form
!
   LOC = 1
70 CONTINUE
   IF( T2( LOC+1, LOC ) /= 0.0D0 ) THEN
!
!        2 by 2 block
!
      IF( T2( LOC, LOC+1 ) == 0.0D0 .OR. T2( LOC, LOC ) /= &
          T2( LOC+1, LOC+1 ) .OR. SIGN( 1.0D0, T2( LOC, LOC+1 ) ) == &
          SIGN( 1.0D0, T2( LOC+1, LOC ) ) )RES = RES + 1.0D0 / EPS
      DO I = LOC + 2, N
         IF( T2( I, LOC ) /= 0.0D0 ) RES = RES + 1.0D0 / RES
         IF( T2( I, LOC+1 ) /= 0.0D0 ) RES = RES + 1.0D0 / RES
      ENDDO
      LOC = LOC + 2
   ELSE
!
!        1 by 1 block
!
      DO I = LOC + 1, N
         IF( T2( I, LOC ) /= 0.0D0 ) RES = RES + 1.0D0 / RES
      ENDDO
      LOC = LOC + 1
   END IF
   IF( LOC < N ) GO TO 70
   IF( RES > RMAX ) THEN
      RMAX = RES
      LMAX = KNT
   END IF
   ENDDO
!
!     End of DGET36
!
END




