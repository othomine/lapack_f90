!> \brief \b SGET35
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET35( RMAX, LMAX, NINFO, KNT )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NINFO
!       REAL               RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET35 tests STRSYL, a routine for solving the Sylvester matrix
!> equation
!>
!>    op(A)*X + ISGN*X*op(B) = scale*C,
!>
!> A and B are assumed to be in Schur canonical form, op() represents an
!> optional transpose, and ISGN can be -1 or +1.  Scale is an output
!> less than or equal to 1, chosen to avoid overflow in X.
!>
!> The test code verifies that the following residual is order 1:
!>
!>    norm(op(A)*X + ISGN*X*op(B) - scale*C) /
!>        (EPS*max(norm(A),norm(B))*norm(X))
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is REAL
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
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SGET35( RMAX, LMAX, NINFO, KNT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX, NINFO
   REAL               RMAX
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   CHARACTER          TRANA, TRANB
   INTEGER            I, IMA, IMB, IMLDA1, IMLDA2, IMLDB1, IMLOFF, &
                      INFO, ISGN, ITRANA, ITRANB, J, M, N
   REAL               BIGNUM, CNRM, EPS, RES, RES1, RMUL, SCALE, &
                      SMLNUM, TNRM, XNRM
!     ..
!     .. Local Arrays ..
   INTEGER            IDIM( 8 ), IVAL( 6, 6, 8 )
   REAL               A( 6, 6 ), B( 6, 6 ), C( 6, 6 ), CC( 6, 6 ), &
                      DUM( 1 ), VM1( 3 ), VM2( 3 )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SLANGE
   EXTERNAL           SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGEMM, STRSYL
!     ..
!     .. Data statements ..
   DATA               IDIM / 1, 2, 3, 4, 3, 3, 6, 4 /
   DATA               IVAL / 1, 35*0, 1, 2, 4*0, -2, 0, 28*0, 1, 5*0, &
                      5, 1, 2, 3*0, -8, -2, 1, 21*0, 3, 4, 4*0, -5, &
                      3, 4*0, 1, 2, 1, 4, 2*0, -3, -9, -1, 1, 14*0, &
                      1, 5*0, 2, 3, 4*0, 5, 6, 7, 21*0, 1, 5*0, 1, 3, &
                      -4, 3*0, 2, 5, 2, 21*0, 1, 2, 4*0, -2, 0, 4*0, &
                      5, 6, 3, 4, 2*0, -1, -9, -5, 2, 2*0, 4*8, 5, 6, &
                      4*9, -7, 5, 1, 5*0, 1, 5, 2, 3*0, 2, -21, 5, &
                      3*0, 1, 2, 3, 4, 14*0 /
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = SLAMCH( 'P' )
   SMLNUM = SLAMCH( 'S' )*4.0E+0 / EPS
   BIGNUM = 1.0E+0 / SMLNUM
!
!     Set up test case parameters
!
   VM1( 1 ) = SQRT( SMLNUM )
   VM1( 2 ) = 1.0E+0
   VM1( 3 ) = SQRT( BIGNUM )
   VM2( 1 ) = 1.0E+0
   VM2( 2 ) = 1.0E+0 + 2.0E+0*EPS
   VM2( 3 ) = 2.0E+0
!
   KNT = 0
   NINFO = 0
   LMAX = 0
   RMAX = 0.0E+0
!
!     Begin test loop
!
   DO ITRANA = 1, 2
      DO ITRANB = 1, 2
         DO ISGN = -1, 1, 2
            DO IMA = 1, 8
               DO IMLDA1 = 1, 3
                  DO IMLDA2 = 1, 3
                     DO IMLOFF = 1, 2
                        DO IMB = 1, 8
                           DO IMLDB1 = 1, 3
                              IF( ITRANA == 1 ) &
                                 TRANA = 'N'
                              IF( ITRANA == 2 ) &
                                 TRANA = 'T'
                              IF( ITRANB == 1 ) &
                                 TRANB = 'N'
                              IF( ITRANB == 2 ) &
                                 TRANB = 'T'
                              M = IDIM( IMA )
                              N = IDIM( IMB )
                              TNRM = 0.0E+0
                              DO I = 1, M
                                 DO J = 1, M
                                    A( I, J ) = IVAL( I, J, IMA )
                                    IF( ABS( I-J ) <= 1 ) THEN
                                       A( I, J ) = A( I, J )* &
                                                   VM1( IMLDA1 )
                                       A( I, J ) = A( I, J )* &
                                                   VM2( IMLDA2 )
                                    ELSE
                                       A( I, J ) = A( I, J )* &
                                                   VM1( IMLOFF )
                                    END IF
                                    TNRM = MAX( TNRM, &
                                           ABS( A( I, J ) ) )
                                 ENDDO
                              ENDDO
                              DO I = 1, N
                                 DO J = 1, N
                                    B( I, J ) = IVAL( I, J, IMB )
                                    IF( ABS( I-J ) <= 1 ) THEN
                                       B( I, J ) = B( I, J )* &
                                                   VM1( IMLDB1 )
                                    ELSE
                                       B( I, J ) = B( I, J )* &
                                                   VM1( IMLOFF )
                                    END IF
                                    TNRM = MAX( TNRM, &
                                           ABS( B( I, J ) ) )
                                 ENDDO
                              ENDDO
                              CNRM = 0.0E+0
                              DO I = 1, M
                                 DO J = 1, N
                                    C( I, J ) = SIN( REAL( I*J ) )
                                    CNRM = MAX( CNRM, C( I, J ) )
                                    CC( I, J ) = C( I, J )
                                 ENDDO
                              ENDDO
                              KNT = KNT + 1
                              CALL STRSYL( TRANA, TRANB, ISGN, M, N, &
                                           A, 6, B, 6, C, 6, SCALE, &
                                           INFO )
                              IF( INFO /= 0 ) &
                                 NINFO = NINFO + 1
                              XNRM = SLANGE( 'M', M, N, C, 6, DUM )
                              RMUL = 1.0E+0
                              IF( XNRM > 1.0E+0 .AND. TNRM > 1.0E+0 ) &
                                   THEN
                                 IF( XNRM > BIGNUM / TNRM ) THEN
                                    RMUL = 1.0E+0 / MAX( XNRM, TNRM )
                                 END IF
                              END IF
                              CALL SGEMM( TRANA, 'N', M, N, M, RMUL, &
                                          A, 6, C, 6, -SCALE*RMUL, &
                                          CC, 6 )
                              CALL SGEMM( 'N', TRANB, M, N, N, &
                                          REAL( ISGN )*RMUL, C, 6, B, &
                                          6, 1.0E+0, CC, 6 )
                              RES1 = SLANGE( 'M', M, N, CC, 6, DUM )
                              RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, &
                                    ( ( RMUL*TNRM )*EPS )*XNRM )
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
!     End of SGET35
!
END

