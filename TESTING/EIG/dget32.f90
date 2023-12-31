!> \brief \b DGET32
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET32( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NINFO
!       DOUBLE PRECISION   RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET32 tests DLASY2, a routine for solving
!>
!>         op(TL)*X + ISGN*X*op(TR) = SCALE*B
!>
!> where TL is N1 by N1, TR is N2 by N2, and N1,N2 =1 or 2 only.
!> X and B are N1 by N2, op() is an optional transpose, an
!> ISGN = 1 or -1. SCALE is chosen less than or equal to 1 to
!> avoid overflow in X.
!>
!> The test condition is that the scaled residual
!>
!> norm( op(TL)*X + ISGN*X*op(TR) = SCALE*B )
!>      / ( max( ulp*norm(TL), ulp*norm(TR)) * norm(X), SMLNUM )
!>
!> should be on the order of 1. Here, ulp is the machine precision.
!> Also, it is verified that SCALE is less than or equal to 1, and
!> that XNORM = infinity-norm(X).
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
!>          Number of examples returned with INFO /= 0.
!> \endverbatim
!>
!> \param[out] KNT
!> \verbatim
!>          KNT is INTEGER
!>          Total number of examples tested.
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
   SUBROUTINE DGET32( RMAX, LMAX, NINFO, KNT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX, NINFO
   DOUBLE PRECISION   RMAX
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            LTRANL, LTRANR
   INTEGER            IB, IB1, IB2, IB3, INFO, ISGN, ITL, ITLSCL, &
                      ITR, ITRANL, ITRANR, ITRSCL, N1, N2
   DOUBLE PRECISION   BIGNUM, DEN, EPS, RES, SCALE, SGN, SMLNUM, TMP, &
                      TNRM, XNORM, XNRM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            ITVAL( 2, 2, 8 )
   DOUBLE PRECISION   B( 2, 2 ), TL( 2, 2 ), TR( 2, 2 ), VAL( 3 ), &
                      X( 2, 2 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLASY2
!     ..
!     .. Data statements ..
   DATA               ITVAL / 8, 4, 2, 1, 4, 8, 1, 2, 2, 1, 8, 4, 1, &
                      2, 4, 8, 9, 4, 2, 1, 4, 9, 1, 2, 2, 1, 9, 4, 1, &
                      2, 4, 9 /
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = DLAMCH( 'P' )
   SMLNUM = DLAMCH( 'S' ) / EPS
   BIGNUM = 1.0D0 / SMLNUM
!
!     Set up test case parameters
!
   VAL( 1 ) = SQRT( SMLNUM )
   VAL( 2 ) = 1.0D0
   VAL( 3 ) = SQRT( BIGNUM )
!
   KNT = 0
   NINFO = 0
   LMAX = 0
   RMAX = 0.0D0
!
!     Begin test loop
!
   DO ITRANL = 0, 1
      DO ITRANR = 0, 1
         DO ISGN = -1, 1, 2
            SGN = ISGN
            LTRANL = ITRANL == 1
            LTRANR = ITRANR == 1
!
            N1 = 1
            N2 = 1
            DO ITL = 1, 3
               DO ITR = 1, 3
                  DO IB = 1, 3
                     TL( 1, 1 ) = VAL( ITL )
                     TR( 1, 1 ) = VAL( ITR )
                     B( 1, 1 ) = VAL( IB )
                     KNT = KNT + 1
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                     CALL DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, &
                                  2, TR, 2, B, 2, SCALE, X, 2, XNORM, &
                                  INFO )
#ifdef _TIMER
                     call system_clock(count_rate=nb_periods_sec,count=S2_time)
                     open(file='results.out', unit=10, position = 'append')
                     write(10,'(A,F16.10,A)') 'Total time : DLASY2 : ',&
                           real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                     close(10)
#endif
                     IF( INFO /= 0 ) &
                        NINFO = NINFO + 1
                     RES = ABS( ( TL( 1, 1 )+SGN*TR( 1, 1 ) )* &
                           X( 1, 1 )-SCALE*B( 1, 1 ) )
                     IF( INFO == 0 ) THEN
                        DEN = MAX( EPS*( ( ABS( TR( 1, &
                              1 ) )+ABS( TL( 1, 1 ) ) )*ABS( X( 1, &
                              1 ) ) ), SMLNUM )
                     ELSE
                        DEN = SMLNUM*MAX( ABS( X( 1, 1 ) ), 1.0D0 )
                     END IF
                     RES = RES / DEN
                     IF( SCALE > 1.0D0 ) &
                        RES = RES + 1.0D0 / EPS
                     RES = RES + ABS( XNORM-ABS( X( 1, 1 ) ) ) / &
                           MAX( SMLNUM, XNORM ) / EPS
                     IF( INFO /= 0 .AND. INFO /= 1 ) &
                        RES = RES + 1.0D0 / EPS
                     IF( RES > RMAX ) THEN
                        LMAX = KNT
                        RMAX = RES
                     END IF
                  ENDDO
               ENDDO
            ENDDO
!
            N1 = 2
            N2 = 1
            DO ITL = 1, 8
               DO ITLSCL = 1, 3
                  DO ITR = 1, 3
                     DO IB1 = 1, 3
                        DO IB2 = 1, 3
                           B( 1, 1 ) = VAL( IB1 )
                           B( 2, 1 ) = -4.0D0*VAL( IB2 )
                           TL( 1, 1 ) = ITVAL( 1, 1, ITL )* &
                                        VAL( ITLSCL )
                           TL( 2, 1 ) = ITVAL( 2, 1, ITL )* &
                                        VAL( ITLSCL )
                           TL( 1, 2 ) = ITVAL( 1, 2, ITL )* &
                                        VAL( ITLSCL )
                           TL( 2, 2 ) = ITVAL( 2, 2, ITL )* &
                                        VAL( ITLSCL )
                           TR( 1, 1 ) = VAL( ITR )
                           KNT = KNT + 1
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL DLASY2( LTRANL, LTRANR, ISGN, N1, N2, &
                                        TL, 2, TR, 2, B, 2, SCALE, X, &
                                        2, XNORM, INFO )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : DLASY2 : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                           IF( INFO /= 0 ) &
                              NINFO = NINFO + 1
                           IF( LTRANL ) THEN
                              TMP = TL( 1, 2 )
                              TL( 1, 2 ) = TL( 2, 1 )
                              TL( 2, 1 ) = TMP
                           END IF
                           RES = ABS( ( TL( 1, 1 )+SGN*TR( 1, 1 ) )* &
                                 X( 1, 1 )+TL( 1, 2 )*X( 2, 1 )- &
                                 SCALE*B( 1, 1 ) )
                           RES = RES + ABS( ( TL( 2, 2 )+SGN*TR( 1, &
                                 1 ) )*X( 2, 1 )+TL( 2, 1 )* &
                                 X( 1, 1 )-SCALE*B( 2, 1 ) )
                           TNRM = ABS( TR( 1, 1 ) ) + &
                                  ABS( TL( 1, 1 ) ) + &
                                  ABS( TL( 1, 2 ) ) + &
                                  ABS( TL( 2, 1 ) ) + &
                                  ABS( TL( 2, 2 ) )
                           XNRM = MAX( ABS( X( 1, 1 ) ), &
                                  ABS( X( 2, 1 ) ) )
                           DEN = MAX( SMLNUM, SMLNUM*XNRM, &
                                 ( TNRM*EPS )*XNRM )
                           RES = RES / DEN
                           IF( SCALE > 1.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           RES = RES + ABS( XNORM-XNRM ) / &
                                 MAX( SMLNUM, XNORM ) / EPS
                           IF( RES > RMAX ) THEN
                              LMAX = KNT
                              RMAX = RES
                           END IF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
!
            N1 = 1
            N2 = 2
            DO ITR = 1, 8
               DO ITRSCL = 1, 3
                  DO ITL = 1, 3
                     DO IB1 = 1, 3
                        DO IB2 = 1, 3
                           B( 1, 1 ) = VAL( IB1 )
                           B( 1, 2 ) = -2.0D0*VAL( IB2 )
                           TR( 1, 1 ) = ITVAL( 1, 1, ITR )* &
                                        VAL( ITRSCL )
                           TR( 2, 1 ) = ITVAL( 2, 1, ITR )* &
                                        VAL( ITRSCL )
                           TR( 1, 2 ) = ITVAL( 1, 2, ITR )* &
                                        VAL( ITRSCL )
                           TR( 2, 2 ) = ITVAL( 2, 2, ITR )* &
                                        VAL( ITRSCL )
                           TL( 1, 1 ) = VAL( ITL )
                           KNT = KNT + 1
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                           CALL DLASY2( LTRANL, LTRANR, ISGN, N1, N2, &
                                        TL, 2, TR, 2, B, 2, SCALE, X, &
                                        2, XNORM, INFO )
#ifdef _TIMER
                           call system_clock(count_rate=nb_periods_sec,count=S2_time)
                           open(file='results.out', unit=10, position = 'append')
                           write(10,'(A,F16.10,A)') 'Total time : DLASY2 : ',&
                                 real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                           close(10)
#endif
                           IF( INFO /= 0 ) &
                              NINFO = NINFO + 1
                           IF( LTRANR ) THEN
                              TMP = TR( 1, 2 )
                              TR( 1, 2 ) = TR( 2, 1 )
                              TR( 2, 1 ) = TMP
                           END IF
                           TNRM = ABS( TL( 1, 1 ) ) + &
                                  ABS( TR( 1, 1 ) ) + &
                                  ABS( TR( 1, 2 ) ) + &
                                  ABS( TR( 2, 2 ) ) + &
                                  ABS( TR( 2, 1 ) )
                           XNRM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
                           RES = ABS( ( ( TL( 1, 1 )+SGN*TR( 1, &
                                 1 ) ) )*( X( 1, 1 ) )+ &
                                 ( SGN*TR( 2, 1 ) )*( X( 1, 2 ) )- &
                                 ( SCALE*B( 1, 1 ) ) )
                           RES = RES + ABS( ( ( TL( 1, 1 )+SGN*TR( 2, &
                                 2 ) ) )*( X( 1, 2 ) )+ &
                                 ( SGN*TR( 1, 2 ) )*( X( 1, 1 ) )- &
                                 ( SCALE*B( 1, 2 ) ) )
                           DEN = MAX( SMLNUM, SMLNUM*XNRM, &
                                 ( TNRM*EPS )*XNRM )
                           RES = RES / DEN
                           IF( SCALE > 1.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           RES = RES + ABS( XNORM-XNRM ) / &
                                 MAX( SMLNUM, XNORM ) / EPS
                           IF( RES > RMAX ) THEN
                              LMAX = KNT
                              RMAX = RES
                           END IF
                        ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
!
            N1 = 2
            N2 = 2
            DO ITR = 1, 8
               DO ITRSCL = 1, 3
                  DO ITL = 1, 8
                     DO ITLSCL = 1, 3
                        DO IB1 = 1, 3
                           DO IB2 = 1, 3
                              DO IB3 = 1, 3
                                 B( 1, 1 ) = VAL( IB1 )
                                 B( 2, 1 ) = -4.0D0*VAL( IB2 )
                                 B( 1, 2 ) = -2.0D0*VAL( IB3 )
                                 B( 2, 2 ) = 8.0D0* &
                                             MIN( VAL( IB1 ), VAL &
                                             ( IB2 ), VAL( IB3 ) )
                                 TR( 1, 1 ) = ITVAL( 1, 1, ITR )* &
                                              VAL( ITRSCL )
                                 TR( 2, 1 ) = ITVAL( 2, 1, ITR )* &
                                              VAL( ITRSCL )
                                 TR( 1, 2 ) = ITVAL( 1, 2, ITR )* &
                                              VAL( ITRSCL )
                                 TR( 2, 2 ) = ITVAL( 2, 2, ITR )* &
                                              VAL( ITRSCL )
                                 TL( 1, 1 ) = ITVAL( 1, 1, ITL )* &
                                              VAL( ITLSCL )
                                 TL( 2, 1 ) = ITVAL( 2, 1, ITL )* &
                                              VAL( ITLSCL )
                                 TL( 1, 2 ) = ITVAL( 1, 2, ITL )* &
                                              VAL( ITLSCL )
                                 TL( 2, 2 ) = ITVAL( 2, 2, ITL )* &
                                              VAL( ITLSCL )
                                 KNT = KNT + 1
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                                 CALL DLASY2( LTRANL, LTRANR, ISGN, &
                                              N1, N2, TL, 2, TR, 2, &
                                              B, 2, SCALE, X, 2, &
                                              XNORM, INFO )
#ifdef _TIMER
                                 call system_clock(count_rate=nb_periods_sec,count=S2_time)
                                 open(file='results.out', unit=10, position = 'append')
                                 write(10,'(A,F16.10,A)') 'Total time : DLASY2 : ',&
                                       real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                                 close(10)
#endif
                                 IF( INFO /= 0 ) &
                                    NINFO = NINFO + 1
                                 IF( LTRANR ) THEN
                                    TMP = TR( 1, 2 )
                                    TR( 1, 2 ) = TR( 2, 1 )
                                    TR( 2, 1 ) = TMP
                                 END IF
                                 IF( LTRANL ) THEN
                                    TMP = TL( 1, 2 )
                                    TL( 1, 2 ) = TL( 2, 1 )
                                    TL( 2, 1 ) = TMP
                                 END IF
                                 TNRM = ABS( TR( 1, 1 ) ) + &
                                        ABS( TR( 2, 1 ) ) + &
                                        ABS( TR( 1, 2 ) ) + &
                                        ABS( TR( 2, 2 ) ) + &
                                        ABS( TL( 1, 1 ) ) + &
                                        ABS( TL( 2, 1 ) ) + &
                                        ABS( TL( 1, 2 ) ) + &
                                        ABS( TL( 2, 2 ) )
                                 XNRM = MAX( ABS( X( 1, 1 ) )+ &
                                        ABS( X( 1, 2 ) ), &
                                        ABS( X( 2, 1 ) )+ &
                                        ABS( X( 2, 2 ) ) )
                                 RES = ABS( ( ( TL( 1, 1 )+SGN*TR( 1, &
                                       1 ) ) )*( X( 1, 1 ) )+ &
                                       ( SGN*TR( 2, 1 ) )* &
                                       ( X( 1, 2 ) )+( TL( 1, 2 ) )* &
                                       ( X( 2, 1 ) )- &
                                       ( SCALE*B( 1, 1 ) ) )
                                 RES = RES + ABS( ( TL( 1, 1 ) )* &
                                       ( X( 1, 2 ) )+ &
                                       ( SGN*TR( 1, 2 ) )* &
                                       ( X( 1, 1 ) )+ &
                                       ( SGN*TR( 2, 2 ) )* &
                                       ( X( 1, 2 ) )+( TL( 1, 2 ) )* &
                                       ( X( 2, 2 ) )- &
                                       ( SCALE*B( 1, 2 ) ) )
                                 RES = RES + ABS( ( TL( 2, 1 ) )* &
                                       ( X( 1, 1 ) )+ &
                                       ( SGN*TR( 1, 1 ) )* &
                                       ( X( 2, 1 ) )+ &
                                       ( SGN*TR( 2, 1 ) )* &
                                       ( X( 2, 2 ) )+( TL( 2, 2 ) )* &
                                       ( X( 2, 1 ) )- &
                                       ( SCALE*B( 2, 1 ) ) )
                                 RES = RES + ABS( ( ( TL( 2, &
                                       2 )+SGN*TR( 2, 2 ) ) )* &
                                       ( X( 2, 2 ) )+ &
                                       ( SGN*TR( 1, 2 ) )* &
                                       ( X( 2, 1 ) )+( TL( 2, 1 ) )* &
                                       ( X( 1, 2 ) )- &
                                       ( SCALE*B( 2, 2 ) ) )
                                 DEN = MAX( SMLNUM, SMLNUM*XNRM, &
                                       ( TNRM*EPS )*XNRM )
                                 RES = RES / DEN
                                 IF( SCALE > 1.0D0 ) &
                                    RES = RES + 1.0D0 / EPS
                                 RES = RES + ABS( XNORM-XNRM ) / &
                                       MAX( SMLNUM, XNORM ) / EPS
                                 IF( RES > RMAX ) THEN
                                    LMAX = KNT
                                    RMAX = RES
                                 END IF
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
   RETURN
!
!     End of DGET32
!
END




