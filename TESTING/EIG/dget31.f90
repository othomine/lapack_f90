!> \brief \b DGET31
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET31( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX
!       DOUBLE PRECISION   RMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            NINFO( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET31 tests DLALN2, a routine for solving
!>
!>    (ca A - w D)X = sB
!>
!> where A is an NA by NA matrix (NA=1 or 2 only), w is a real (NW=1) or
!> complex (NW=2) constant, ca is a real constant, D is an NA by NA real
!> diagonal matrix, and B is an NA by NW matrix (when NW=2 the second
!> column of B contains the imaginary part of the solution).  The code
!> returns X and s, where s is a scale factor, less than or equal to 1,
!> which is chosen to avoid overflow in X.
!>
!> If any singular values of ca A-w D are less than another input
!> parameter SMIN, they are perturbed up to SMIN.
!>
!> The test condition is that the scaled residual
!>
!>     norm( (ca A-w D)*X - s*B ) /
!>           ( max( ulp*norm(ca A-w D), SMIN )*norm(X) )
!>
!> should be on the order of 1.  Here, ulp is the machine precision.
!> Also, it is verified that SCALE is less than or equal to 1, and that
!> XNORM = infinity-norm(X).
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
!>          NINFO(1) = number of examples with INFO less than 0
!>          NINFO(2) = number of examples with INFO greater than 0
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
!
!> \ingroup double_eig
!
!  =====================================================================
   SUBROUTINE DGET31( RMAX, LMAX, NINFO, KNT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX
   DOUBLE PRECISION   RMAX
!     ..
!     .. Array Arguments ..
   INTEGER            NINFO( 2 )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            IA, IB, ICA, ID1, ID2, INFO, ISMIN, ITRANS, &
                      IWI, IWR, NA, NW
   DOUBLE PRECISION   BIGNUM, CA, D1, D2, DEN, EPS, RES, SCALE, SMIN, &
                      SMLNUM, TMP, UNFL, WI, WR, XNORM
!     ..
!     .. Local Arrays ..
   LOGICAL            LTRANS( 0: 1 )
   DOUBLE PRECISION   A( 2, 2 ), B( 2, 2 ), VAB( 3 ), VCA( 5 ), &
                      VDD( 4 ), VSMIN( 4 ), VWI( 4 ), VWR( 4 ), &
                      X( 2, 2 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLALN2
!     ..
!     .. Data statements ..
   DATA               LTRANS / .FALSE., .TRUE. /
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = DLAMCH( 'P' )
   UNFL = DLAMCH( 'U' )
   SMLNUM = DLAMCH( 'S' ) / EPS
   BIGNUM = 1.0D0 / SMLNUM
!
!     Set up test case parameters
!
   VSMIN( 1 ) = SMLNUM
   VSMIN( 2 ) = EPS
   VSMIN( 3 ) = 1.0D0 / ( 10.0D0*10.0D0 )
   VSMIN( 4 ) = 1.0D0 / EPS
   VAB( 1 ) = SQRT( SMLNUM )
   VAB( 2 ) = 1.0D0
   VAB( 3 ) = SQRT( BIGNUM )
   VWR( 1 ) = 0.0D0
   VWR( 2 ) = 0.5D0
   VWR( 3 ) = 2.0D0
   VWR( 4 ) = 1.0D0
   VWI( 1 ) = SMLNUM
   VWI( 2 ) = EPS
   VWI( 3 ) = 1.0D0
   VWI( 4 ) = 2.0D0
   VDD( 1 ) = SQRT( SMLNUM )
   VDD( 2 ) = 1.0D0
   VDD( 3 ) = 2.0D0
   VDD( 4 ) = SQRT( BIGNUM )
   VCA( 1 ) = 0.0D0
   VCA( 2 ) = SQRT( SMLNUM )
   VCA( 3 ) = EPS
   VCA( 4 ) = 0.5D0
   VCA( 5 ) = 1.0D0
!
   KNT = 0
   NINFO( 1:2 ) = 0
   LMAX = 0
   RMAX = 0.0D0
!
!     Begin test loop
!
   DO ID1 = 1, 4
      D1 = VDD( ID1 )
      DO ID2 = 1, 4
         D2 = VDD( ID2 )
         DO ICA = 1, 5
            CA = VCA( ICA )
            DO ITRANS = 0, 1
               DO ISMIN = 1, 4
                  SMIN = VSMIN( ISMIN )
                  NA = 1
                  NW = 1
                  DO IA = 1, 3
                     A( 1, 1 ) = VAB( IA )
                     DO IB = 1, 3
                        B( 1, 1 ) = VAB( IB )
                        DO IWR = 1, 4
                           IF( D1 == 1.0D0 .AND. D2 == 1.0D0 .AND. CA == &
                               1.0D0 ) THEN
                              WR = VWR( IWR )*A( 1, 1 )
                           ELSE
                              WR = VWR( IWR )
                           END IF
                           WI = 0.0D0
                           CALL DLALN2( LTRANS( ITRANS ), NA, NW, &
                                        SMIN, CA, A, 2, D1, D2, B, 2, &
                                        WR, WI, X, 2, SCALE, XNORM, &
                                        INFO )
                           IF( INFO < 0 ) &
                              NINFO( 1 ) = NINFO( 1 ) + 1
                           IF( INFO > 0 ) &
                              NINFO( 2 ) = NINFO( 2 ) + 1
                           RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* &
                                 X( 1, 1 )-SCALE*B( 1, 1 ) )
                           IF( INFO == 0 ) THEN
                              DEN = MAX( EPS*( ABS( ( CA*A( 1, &
                                    1 )-WR*D1 )*X( 1, 1 ) ) ), &
                                    SMLNUM )
                           ELSE
                              DEN = MAX( SMIN*ABS( X( 1, 1 ) ), &
                                    SMLNUM )
                           END IF
                           RES = RES / DEN
                           IF( ABS( X( 1, 1 ) ) < UNFL .AND. &
                               ABS( B( 1, 1 ) ) <= SMLNUM* &
                               ABS( CA*A( 1, 1 )-WR*D1 ) )RES = 0.0D0
                           IF( SCALE > 1.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           RES = RES + ABS( XNORM-ABS( X( 1, 1 ) ) ) &
                                  / MAX( SMLNUM, XNORM ) / EPS
                           IF( INFO /= 0 .AND. INFO /= 1 ) &
                              RES = RES + 1.0D0 / EPS
                           KNT = KNT + 1
                           IF( RES > RMAX ) THEN
                              LMAX = KNT
                              RMAX = RES
                           END IF
                        ENDDO
                     ENDDO
                  ENDDO
!
                  NA = 1
                  NW = 2
                  DO IA = 1, 3
                     A( 1, 1 ) = VAB( IA )
                     DO IB = 1, 3
                        B( 1, 1 ) = VAB( IB )
                        B( 1, 2 ) = -0.5D0*VAB( IB )
                        DO IWR = 1, 4
                           IF( D1 == 1.0D0 .AND. D2 == 1.0D0 .AND. CA == &
                               1.0D0 ) THEN
                              WR = VWR( IWR )*A( 1, 1 )
                           ELSE
                              WR = VWR( IWR )
                           END IF
                           DO IWI = 1, 4
                              IF( D1 == 1.0D0 .AND. D2 == 1.0D0 .AND. &
                                  CA == 1.0D0 ) THEN
                                 WI = VWI( IWI )*A( 1, 1 )
                              ELSE
                                 WI = VWI( IWI )
                              END IF
                              CALL DLALN2( LTRANS( ITRANS ), NA, NW, &
                                           SMIN, CA, A, 2, D1, D2, B, &
                                           2, WR, WI, X, 2, SCALE, &
                                           XNORM, INFO )
                              IF( INFO < 0 ) NINFO( 1 ) = NINFO( 1 ) + 1
                              IF( INFO > 0 ) NINFO( 2 ) = NINFO( 2 ) + 1
                              RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* &
                                    X( 1, 1 )+( WI*D1 )*X( 1, 2 )- &
                                    SCALE*B( 1, 1 ) )
                              RES = RES + ABS( ( -WI*D1 )*X( 1, 1 )+ &
                                    ( CA*A( 1, 1 )-WR*D1 )*X( 1, 2 )- &
                                    SCALE*B( 1, 2 ) )
                              IF( INFO == 0 ) THEN
                                 DEN = MAX( EPS*( MAX( ABS( CA*A( 1, &
                                       1 )-WR*D1 ), ABS( D1*WI ) )* &
                                       ( ABS( X( 1, 1 ) )+ABS( X( 1, &
                                       2 ) ) ) ), SMLNUM )
                              ELSE
                                 DEN = MAX( SMIN*( ABS( X( 1, &
                                       1 ) )+ABS( X( 1, 2 ) ) ), &
                                       SMLNUM )
                              END IF
                              RES = RES / DEN
                              IF( ABS( X( 1, 1 ) ) < UNFL .AND. &
                                  ABS( X( 1, 2 ) ) < UNFL .AND. &
                                  ABS( B( 1, 1 ) ) <= SMLNUM* &
                                  ABS( CA*A( 1, 1 )-WR*D1 ) ) &
                                  RES = 0.0D0
                              IF( SCALE > 1.0D0 ) &
                                 RES = RES + 1.0D0 / EPS
                              RES = RES + ABS( XNORM- &
                                    ABS( X( 1, 1 ) )- &
                                    ABS( X( 1, 2 ) ) ) / &
                                    MAX( SMLNUM, XNORM ) / EPS
                              IF( INFO /= 0 .AND. INFO /= 1 ) &
                                 RES = RES + 1.0D0 / EPS
                              KNT = KNT + 1
                              IF( RES > RMAX ) THEN
                                 LMAX = KNT
                                 RMAX = RES
                              END IF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
!
                  NA = 2
                  NW = 1
                  DO IA = 1, 3
                     A( 1, 1 ) = VAB( IA )
                     A( 1, 2 ) = -3.0D0*VAB( IA )
                     A( 2, 1 ) = -7.0D0*VAB( IA )
                     A( 2, 2 ) = 21.0D0*VAB( IA )
                     DO IB = 1, 3
                        B( 1, 1 ) = VAB( IB )
                        B( 2, 1 ) = -2.0D0*VAB( IB )
                        DO IWR = 1, 4
                           IF( D1 == 1.0D0 .AND. D2 == 1.0D0 .AND. CA == &
                               1.0D0 ) THEN
                              WR = VWR( IWR )*A( 1, 1 )
                           ELSE
                              WR = VWR( IWR )
                           END IF
                           WI = 0.0D0
                           CALL DLALN2( LTRANS( ITRANS ), NA, NW, &
                                        SMIN, CA, A, 2, D1, D2, B, 2, &
                                        WR, WI, X, 2, SCALE, XNORM, &
                                        INFO )
                           IF( INFO < 0 ) &
                              NINFO( 1 ) = NINFO( 1 ) + 1
                           IF( INFO > 0 ) &
                              NINFO( 2 ) = NINFO( 2 ) + 1
                           IF( ITRANS == 1 ) THEN
                              TMP = A( 1, 2 )
                              A( 1, 2 ) = A( 2, 1 )
                              A( 2, 1 ) = TMP
                           END IF
                           RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* &
                                 X( 1, 1 )+( CA*A( 1, 2 ) )* &
                                 X( 2, 1 )-SCALE*B( 1, 1 ) )
                           RES = RES + ABS( ( CA*A( 2, 1 ) )* &
                                 X( 1, 1 )+( CA*A( 2, 2 )-WR*D2 )* &
                                 X( 2, 1 )-SCALE*B( 2, 1 ) )
                           IF( INFO == 0 ) THEN
                              DEN = MAX( EPS*( MAX( ABS( CA*A( 1, &
                                    1 )-WR*D1 )+ABS( CA*A( 1, 2 ) ), &
                                    ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, &
                                    2 )-WR*D2 ) )*MAX( ABS( X( 1, &
                                    1 ) ), ABS( X( 2, 1 ) ) ) ), &
                                    SMLNUM )
                           ELSE
                              DEN = MAX( EPS*( MAX( SMIN / EPS, &
                                    MAX( ABS( CA*A( 1, &
                                    1 )-WR*D1 )+ABS( CA*A( 1, 2 ) ), &
                                    ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, &
                                    2 )-WR*D2 ) ) )*MAX( ABS( X( 1, &
                                    1 ) ), ABS( X( 2, 1 ) ) ) ), &
                                    SMLNUM )
                           END IF
                           RES = RES / DEN
                           IF( ABS( X( 1, 1 ) ) < UNFL .AND. &
                               ABS( X( 2, 1 ) ) < UNFL .AND. &
                               ABS( B( 1, 1 ) )+ABS( B( 2, 1 ) ) <= &
                               SMLNUM*( ABS( CA*A( 1, &
                               1 )-WR*D1 )+ABS( CA*A( 1, &
                               2 ) )+ABS( CA*A( 2, &
                               1 ) )+ABS( CA*A( 2, 2 )-WR*D2 ) ) ) &
                               RES = 0.0D0
                           IF( SCALE > 1.0D0 ) &
                              RES = RES + 1.0D0 / EPS
                           RES = RES + ABS( XNORM- &
                                 MAX( ABS( X( 1, 1 ) ), ABS( X( 2, &
                                 1 ) ) ) ) / MAX( SMLNUM, XNORM ) / &
                                 EPS
                           IF( INFO /= 0 .AND. INFO /= 1 ) &
                              RES = RES + 1.0D0 / EPS
                           KNT = KNT + 1
                           IF( RES > RMAX ) THEN
                              LMAX = KNT
                              RMAX = RES
                           END IF
                        ENDDO
                     ENDDO
                     ENDDO
!
                  NA = 2
                  NW = 2
                  DO IA = 1, 3
                     A( 1, 1 ) = VAB( IA )*2.0D0
                     A( 1, 2 ) = -3.0D0*VAB( IA )
                     A( 2, 1 ) = -7.0D0*VAB( IA )
                     A( 2, 2 ) = 21.0D0*VAB( IA )
                     DO IB = 1, 3
                        B( 1, 1 ) = VAB( IB )
                        B( 2, 1 ) = -2.0D0*VAB( IB )
                        B( 1, 2 ) = 4.0D0*VAB( IB )
                        B( 2, 2 ) = -7.0D0*VAB( IB )
                        DO IWR = 1, 4
                           IF( D1 == 1.0D0 .AND. D2 == 1.0D0 .AND. CA == &
                               1.0D0 ) THEN
                              WR = VWR( IWR )*A( 1, 1 )
                           ELSE
                              WR = VWR( IWR )
                           END IF
                           DO IWI = 1, 4
                              IF( D1 == 1.0D0 .AND. D2 == 1.0D0 .AND. &
                                  CA == 1.0D0 ) THEN
                                 WI = VWI( IWI )*A( 1, 1 )
                              ELSE
                                 WI = VWI( IWI )
                              END IF
                              CALL DLALN2( LTRANS( ITRANS ), NA, NW, &
                                           SMIN, CA, A, 2, D1, D2, B, &
                                           2, WR, WI, X, 2, SCALE, &
                                           XNORM, INFO )
                              IF( INFO < 0 ) &
                                 NINFO( 1 ) = NINFO( 1 ) + 1
                              IF( INFO > 0 ) &
                                 NINFO( 2 ) = NINFO( 2 ) + 1
                              IF( ITRANS == 1 ) THEN
                                 TMP = A( 1, 2 )
                                 A( 1, 2 ) = A( 2, 1 )
                                 A( 2, 1 ) = TMP
                              END IF
                              RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* &
                                    X( 1, 1 )+( CA*A( 1, 2 ) )* &
                                    X( 2, 1 )+( WI*D1 )*X( 1, 2 )- &
                                    SCALE*B( 1, 1 ) )
                              RES = RES + ABS( ( CA*A( 1, &
                                    1 )-WR*D1 )*X( 1, 2 )+ &
                                    ( CA*A( 1, 2 ) )*X( 2, 2 )- &
                                    ( WI*D1 )*X( 1, 1 )-SCALE* &
                                    B( 1, 2 ) )
                              RES = RES + ABS( ( CA*A( 2, 1 ) )* &
                                    X( 1, 1 )+( CA*A( 2, 2 )-WR*D2 )* &
                                    X( 2, 1 )+( WI*D2 )*X( 2, 2 )- &
                                    SCALE*B( 2, 1 ) )
                              RES = RES + ABS( ( CA*A( 2, 1 ) )* &
                                    X( 1, 2 )+( CA*A( 2, 2 )-WR*D2 )* &
                                    X( 2, 2 )-( WI*D2 )*X( 2, 1 )- &
                                    SCALE*B( 2, 2 ) )
                              IF( INFO == 0 ) THEN
                                 DEN = MAX( EPS*( MAX( ABS( CA*A( 1, &
                                       1 )-WR*D1 )+ABS( CA*A( 1, &
                                       2 ) )+ABS( WI*D1 ), &
                                       ABS( CA*A( 2, &
                                       1 ) )+ABS( CA*A( 2, &
                                       2 )-WR*D2 )+ABS( WI*D2 ) )* &
                                       MAX( ABS( X( 1, &
                                       1 ) )+ABS( X( 2, 1 ) ), &
                                       ABS( X( 1, 2 ) )+ABS( X( 2, &
                                       2 ) ) ) ), SMLNUM )
                              ELSE
                                 DEN = MAX( EPS*( MAX( SMIN / EPS, &
                                       MAX( ABS( CA*A( 1, &
                                       1 )-WR*D1 )+ABS( CA*A( 1, &
                                       2 ) )+ABS( WI*D1 ), &
                                       ABS( CA*A( 2, &
                                       1 ) )+ABS( CA*A( 2, &
                                       2 )-WR*D2 )+ABS( WI*D2 ) ) )* &
                                       MAX( ABS( X( 1, &
                                       1 ) )+ABS( X( 2, 1 ) ), &
                                       ABS( X( 1, 2 ) )+ABS( X( 2, &
                                       2 ) ) ) ), SMLNUM )
                              END IF
                              RES = RES / DEN
                              IF( ABS( X( 1, 1 ) ) < UNFL .AND. &
                                  ABS( X( 2, 1 ) ) < UNFL .AND. &
                                  ABS( X( 1, 2 ) ) < UNFL .AND. &
                                  ABS( X( 2, 2 ) ) < UNFL .AND. &
                                  ABS( B( 1, 1 ) )+ &
                                  ABS( B( 2, 1 ) ) <= SMLNUM* &
                                  ( ABS( CA*A( 1, 1 )-WR*D1 )+ &
                                  ABS( CA*A( 1, 2 ) )+ABS( CA*A( 2, &
                                  1 ) )+ABS( CA*A( 2, &
                                  2 )-WR*D2 )+ABS( WI*D2 )+ABS( WI* &
                                  D1 ) ) )RES = 0.0D0
                              IF( SCALE > 1.0D0 ) &
                                 RES = RES + 1.0D0 / EPS
                              RES = RES + ABS( XNORM- &
                                    MAX( ABS( X( 1, 1 ) )+ABS( X( 1, &
                                    2 ) ), ABS( X( 2, &
                                    1 ) )+ABS( X( 2, 2 ) ) ) ) / &
                                    MAX( SMLNUM, XNORM ) / EPS
                              IF( INFO /= 0 .AND. INFO /= 1 ) &
                                 RES = RES + 1.0D0 / EPS
                              KNT = KNT + 1
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
!
   RETURN
!
!     End of DGET31
!
END

