*> \brief \b SGET40
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGET40( RMAX, LMAX, NINFO, KNT, NIN )
*
*       .. Scalar Arguments ..
*       INTEGER            KNT, LMAX, NIN
*       REAL               RMAX
*       ..
*       .. Array Arguments ..
*       INTEGER            NINFO( 2 )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGET40 tests STGEXC, a routine for swapping adjacent blocks (either
*> 1 by 1 or 2 by 2) on the diagonal of a pencil in real generalized Schur form.
*> Thus, STGEXC computes an orthogonal matrices Q and Z such that
*>
*>     Q' * ( [ A B ], [ D E ] ) * Z  = ( [ C1 B1 ], [ F1 E1 ] )
*>          ( [ 0 C ]  [   F ] )        ( [ 0  A1 ]  [    D1]  )
*>
*> where (C1,F1) is similar to (C,F) and (A1,D1) is similar to (A,D).
*> Both (A,D) and (C,F) are assumed to be in standard form
*> and (A1,D1) and (C1,F1) are returned with the
*> same properties.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[out] RMAX
*> \verbatim
*>          RMAX is REAL
*>          Value of the largest test ratio.
*> \endverbatim
*>
*> \param[out] LMAX
*> \verbatim
*>          LMAX is INTEGER
*>          Example number where largest test ratio achieved.
*> \endverbatim
*>
*> \param[out] NINFO
*> \verbatim
*>          NINFO is INTEGER array, dimension (2)
*>          NINFO( 1 ) = STGEXC without accumulation returned INFO nonzero
*>          NINFO( 2 ) = STGEXC with accumulation returned INFO nonzero
*> \endverbatim
*>
*> \param[out] KNT
*> \verbatim
*>          KNT is INTEGER
*>          Total number of examples tested.
*> \endverbatim
*>
*> \param[in] NIN
*> \verbatim
*>          NIN is INTEGER
*>          Input logical unit number.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE SGET40( RMAX, LMAX, NINFO, KNT, NIN )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            KNT, LMAX, NIN
      REAL               RMAX
*     ..
*     .. Array Arguments ..
      INTEGER            NINFO( 2 )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0, ONE = 1.0 )
      INTEGER            LDT, LWORK
      PARAMETER          ( LDT = 10, LWORK = 100 + 4*LDT + 16 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1,
     $                   ILST2, ILSTSV, J, LOC, N
      REAL               EPS, RES
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
*     ..
*     .. Local Arrays ..
      REAL               Q( LDT, LDT ), Z( LDT, LDT ), RESULT( 4 ),
     $                   T( LDT, LDT ), T1( LDT, LDT ), T2( LDT, LDT ),
     $                   S( LDT, LDT ), S1( LDT, LDT ), S2( LDT, LDT ),
     $                   TMP( LDT, LDT ), WORK( LWORK )
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGET51, SLACPY, SLASET, STGEXC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'P' )
      RMAX = ZERO
      LMAX = 0
      KNT = 0
      NINFO( 1 ) = 0
      NINFO( 2 ) = 0
*
*     Read input data until N=0
*
   10 CONTINUE
      READ( NIN, FMT = * )N, IFST, ILST
      IF( N.EQ.0 )
     $   RETURN
      KNT = KNT + 1
      DO 20 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   20 CONTINUE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLACPY( 'F', N, N, TMP, LDT, T, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLACPY( 'F', N, N, TMP, LDT, T1, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLACPY( 'F', N, N, TMP, LDT, T2, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      DO 25 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   25 CONTINUE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLACPY( 'F', N, N, TMP, LDT, S, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLACPY( 'F', N, N, TMP, LDT, S1, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLACPY( 'F', N, N, TMP, LDT, S2, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IFSTSV = IFST
      ILSTSV = ILST
      IFST1 = IFST
      ILST1 = ILST
      IFST2 = IFST
      ILST2 = ILST
      RES = ZERO
*
*     Test without accumulating Q and Z
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLASET : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLASET : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL STGEXC( .FALSE., .FALSE., N, T1, LDT, S1, LDT, Q, LDT,
     $             Z, LDT, IFST1, ILST1, WORK, LWORK, NINFO( 1 ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : STGEXC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      DO 40 I = 1, N
         DO 30 J = 1, N
            IF( I.EQ.J .AND. Q( I, J ).NE.ONE )
     $         RES = RES + ONE / EPS
            IF( I.NE.J .AND. Q( I, J ).NE.ZERO )
     $         RES = RES + ONE / EPS
            IF( I.EQ.J .AND. Z( I, J ).NE.ONE )
     $         RES = RES + ONE / EPS
            IF( I.NE.J .AND. Z( I, J ).NE.ZERO )
     $         RES = RES + ONE / EPS
   30    CONTINUE
   40 CONTINUE
*
*     Test with accumulating Q
*
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLASET : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : SLASET : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S1_time)
#endif
      CALL STGEXC( .TRUE., .TRUE., N, T2, LDT, S2, LDT, Q, LDT,
     $             Z, LDT, IFST2, ILST2, WORK, LWORK, NINFO( 2 ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,
     $count=S2_time)
      open(file='results.out',unit=10,access='append')
      write(10,'(A,F16.10,A)') 'Total time : STGEXC : ',
     $real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
*
*     Compare T1 with T2 and S1 with S2
*
      DO 60 I = 1, N
         DO 50 J = 1, N
            IF( T1( I, J ).NE.T2( I, J ) )
     $         RES = RES + ONE / EPS
            IF( S1( I, J ).NE.S2( I, J ) )
     $         RES = RES + ONE / EPS
   50    CONTINUE
   60 CONTINUE
      IF( IFST1.NE.IFST2 )
     $   RES = RES + ONE / EPS
      IF( ILST1.NE.ILST2 )
     $   RES = RES + ONE / EPS
      IF( NINFO( 1 ).NE.NINFO( 2 ) )
     $   RES = RES + ONE / EPS
*
*     Test orthogonality of Q and Z and backward error on T2 and S2
*
      CALL SGET51( 1, N, T, LDT, T2, LDT, Q, LDT, Z, LDT, WORK,
     $             RESULT( 1 ) )
      CALL SGET51( 1, N, S, LDT, S2, LDT, Q, LDT, Z, LDT, WORK,
     $             RESULT( 2 ) )
      CALL SGET51( 3, N, T, LDT, T2, LDT, Q, LDT, Q, LDT, WORK,
     $             RESULT( 3 ) )
      CALL SGET51( 3, N, T, LDT, T2, LDT, Z, LDT, Z, LDT, WORK,
     $             RESULT( 4 ) )
      RES = RES + RESULT( 1 ) + RESULT( 2 ) + RESULT( 3 ) + RESULT( 4 )
*
*     Read next matrix pair
*
      GO TO 10
*
*     End of SGET40
*
      END

