!> \brief \b DGET40
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET40( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!      INTEGER            KNT, LMAX, NIN
!      DOUBLE PRECISION   RMAX
!      ..
!       .. Array Arguments ..
!      INTEGER            NINFO( 2 )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET40 tests DTGEXC, a routine for swapping adjacent blocks (either
!> 1 by 1 or 2 by 2) on the diagonal of a pencil in real generalized Schur form.
!> Thus, DTGEXC computes an orthogonal matrices Q and Z such that
!>
!>     Q' * ( [ A B ], [ D E ] ) * Z  = ( [ C1 B1 ], [ F1 E1 ] )
!>          ( [ 0 C ]  [   F ] )        ( [ 0  A1 ]  [    D1]  )
!>
!> where (C1,F1) is similar to (C,F) and (A1,D1) is similar to (A,D).
!> Both (A,D) and (C,F) are assumed to be in standard form
!> and (A1,D1) and (C1,F1) are returned with the
!> same properties.
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
!>          NINFO is INTEGER array, dimension (2)
!>          NINFO( 1 ) = DTGEXC without accumulation returned INFO nonzero
!>          NINFO( 2 ) = DTGEXC with accumulation returned INFO nonzero
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
   SUBROUTINE DGET40( RMAX, LMAX, NINFO, KNT, NIN )
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
   INTEGER            NINFO( 2 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            LDT, LWORK
   PARAMETER          ( LDT = 10, LWORK = 100 + 4*LDT + 16 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1, &
                      ILST2, ILSTSV, J, LOC, N
   DOUBLE PRECISION   EPS, RES
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   Q( LDT, LDT ), Z( LDT, LDT ), RESULT( 4 ), &
                      T( LDT, LDT ), T1( LDT, LDT ), T2( LDT, LDT ), &
                      S( LDT, LDT ), S1( LDT, LDT ), S2( LDT, LDT ), &
                      TMP( LDT, LDT ), WORK( LWORK )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DHST01, DLACPY, DLASET, DTGEXC
!     ..
!     .. Executable Statements ..
!
   EPS = DLAMCH( 'P' )
   RMAX = 0.0D0
   LMAX = 0
   KNT = 0
   NINFO( 1 ) = 0
   NINFO( 2 ) = 0
!
!     Read input data until N=0
!
   DO
   READ(NIN,*) N, IFST, ILST
   IF( N == 0 ) &
      RETURN
   KNT = KNT + 1
   DO I = 1, N
      READ(NIN,*) TMP( I,1:N)
   ENDDO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'F', N, N, TMP, LDT, T, LDT )
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
   DO I = 1, N
      READ(NIN,*) TMP( I,1:N)
   ENDDO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLACPY( 'F', N, N, TMP, LDT, S, LDT )
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
   CALL DLACPY( 'F', N, N, TMP, LDT, S1, LDT )
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
   CALL DLACPY( 'F', N, N, TMP, LDT, S2, LDT )
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
!     Test without accumulating Q and Z
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
   CALL DLASET( 'Full', N, N, 0.0D0, 1.0D0, Z, LDT )
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
   CALL DTGEXC( .FALSE., .FALSE., N, T1, LDT, S1, LDT, Q, LDT, &
                Z, LDT, IFST1, ILST1, WORK, LWORK, NINFO ( 1 ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DTGEXC : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   DO I = 1, N
      DO J = 1, N
         IF( I == J .AND. Q( I, J ) /= 1.0D0 ) RES = RES + 1.0D0 / EPS
         IF( I /= J .AND. Q( I, J ) /= 0.0D0 ) RES = RES + 1.0D0 / EPS
         IF( I == J .AND. Z( I, J ) /= 1.0D0 ) RES = RES + 1.0D0 / EPS
         IF( I /= J .AND. Z( I, J ) /= 0.0D0 ) RES = RES + 1.0D0 / EPS
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
   CALL DLASET( 'Full', N, N, 0.0D0, 1.0D0, Z, LDT )
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
   CALL DTGEXC( .TRUE., .TRUE., N, T2, LDT, S2, LDT, Q, LDT, &
                Z, LDT, IFST2, ILST2, WORK, LWORK, NINFO ( 2 ) )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DTGEXC : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Compare T1 with T2 and S1 with S2
!
   DO I = 1, N
      DO J = 1, N
         IF( T1( I, J ) /= T2( I, J ) ) RES = RES + 1.0D0 / EPS
         IF( S1( I, J ) /= S2( I, J ) ) RES = RES + 1.0D0 / EPS
      ENDDO
   ENDDO
   IF( IFST1 /= IFST2 ) RES = RES + 1.0D0 / EPS
   IF( ILST1 /= ILST2 ) RES = RES + 1.0D0 / EPS
   IF( NINFO( 1 ) /= NINFO( 2 ) ) RES = RES + 1.0D0 / EPS
!
!     Test orthogonality of Q and Z and backward error on T2 and S2
!
   CALL DGET51( 1, N, T, LDT, T2, LDT, Q, LDT, Z, LDT, WORK, RESULT( 1 ) )
   CALL DGET51( 1, N, S, LDT, S2, LDT, Q, LDT, Z, LDT, WORK, RESULT( 2 ) )
   CALL DGET51( 3, N, T, LDT, T2, LDT, Q, LDT, Q, LDT, WORK, RESULT( 3 ) )
   CALL DGET51( 3, N, T, LDT, T2, LDT, Z, LDT, Z, LDT, WORK, RESULT( 4 ) )
   RES = RES + RESULT( 1 ) + RESULT( 2 ) + RESULT( 3 ) + RESULT( 4 )
!
!     Read next matrix pair
!
   ENDDO
!
!     End of DGET40
!
END




