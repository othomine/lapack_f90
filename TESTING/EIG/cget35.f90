!> \brief \b CGET35
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGET35( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NIN, NINFO
!       REAL               RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGET35 tests CTRSYL, a routine for solving the Sylvester matrix
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
!
!> \ingroup complex_eig
!
!  =====================================================================
   SUBROUTINE CGET35( RMAX, LMAX, NINFO, KNT, NIN )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, LMAX, NIN, NINFO
   REAL               RMAX
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            LDT
   PARAMETER          ( LDT = 10 )
!     ..
!     .. Local Scalars ..
   CHARACTER          TRANA, TRANB
   INTEGER            I, IMLA, IMLAD, IMLB, IMLC, INFO, ISGN, ITRANA, &
                      ITRANB, J, M, N
   REAL               BIGNUM, EPS, RES, RES1, SCALE, SMLNUM, TNRM, &
                      XNRM
   COMPLEX            RMUL
!     ..
!     .. Local Arrays ..
   REAL               DUM( 1 ), VM1( 3 ), VM2( 3 )
   COMPLEX            A( LDT, LDT ), ATMP( LDT, LDT ), B( LDT, LDT ), &
                      BTMP( LDT, LDT ), C( LDT, LDT ), &
                      CSAV( LDT, LDT ), CTMP( LDT, LDT )
!     ..
!     .. External Functions ..
   REAL               CLANGE, SLAMCH
   EXTERNAL           CLANGE, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           CGEMM, CTRSYL
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
   EPS = SLAMCH( 'P' )
   SMLNUM = SLAMCH( 'S' ) / EPS
   BIGNUM = 1.0E0 / SMLNUM
!
!     Set up test case parameters
!
   VM1( 1 ) = SQRT( SMLNUM )
   VM1( 2 ) = 1.0E0
   VM1( 3 ) = 1.0E6
   VM2( 1 ) = 1.0E0
   VM2( 2 ) = 1.0E0 + 2.0E0*EPS
   VM2( 3 ) = 2.0E0
!
   KNT = 0
   NINFO = 0
   LMAX = 0
   RMAX = 0.0E0
!
!     Begin test loop
!
   DO
   READ(NIN,*)M, N
   IF( N == 0 ) RETURN
   DO I = 1, M
      READ(NIN,*) ATMP(I,1:M)
   ENDDO
   DO I = 1, N
      READ(NIN,*) BTMP(I,1:N)
   ENDDO
   DO I = 1, M
      READ(NIN,*) CTMP(I,1:N)
   ENDDO
   DO IMLA = 1, 3
      DO IMLAD = 1, 3
         DO IMLB = 1, 3
            DO IMLC = 1, 3
               DO ITRANA = 1, 2
                  DO ITRANB = 1, 2
                     DO ISGN = -1, 1, 2
                        IF( ITRANA == 1 ) TRANA = 'N'
                        IF( ITRANA == 2 ) TRANA = 'C'
                        IF( ITRANB == 1 ) TRANB = 'N'
                        IF( ITRANB == 2 ) TRANB = 'C'
                        TNRM = 0.0E0
                        DO I = 1, M
                           A(I,1:M) = ATMP(I,1:M)*VM1(IMLA)
                           A(I,I) = A(I,I)*VM2(IMLAD)
                        ENDDO
                        TNRM = MAXVAL(ABS(A(1:M,1:M)))
                        B(1:N,1:N) = BTMP(1:N,1:N)*VM1( IMLB )
                        TNRM = MAXVAL(ABS(B(1:N,1:N)))
                        IF( TNRM == 0.0E0 ) TNRM = 1.0E0
                        C(1:M,1:N) = CTMP(1:M,1:N)*VM1( IMLC )
                        CSAV(1:M,1:N) = C(1:M,1:N)
                        KNT = KNT + 1
                        CALL CTRSYL( TRANA, TRANB, ISGN, M, N, A, &
                                     LDT, B, LDT, C, LDT, SCALE, INFO )
                        IF( INFO /= 0 ) NINFO = NINFO + 1
                        XNRM = CLANGE( 'M', M, N, C, LDT, DUM )
                        RMUL = 1.0E0
                        IF( XNRM > 1.0E0 .AND. TNRM > 1.0E0 ) THEN
                           IF( XNRM > BIGNUM / TNRM ) THEN
                              RMUL = MAX( XNRM, TNRM )
                              RMUL = 1.0E0 / RMUL
                           END IF
                        END IF
                        CALL CGEMM( TRANA, 'N', M, N, M, RMUL, A, &
                                    LDT, C, LDT, -SCALE*RMUL, CSAV, LDT )
                        CALL CGEMM( 'N', TRANB, M, N, N, &
                                    REAL( ISGN )*RMUL, C, LDT, B, &
                                    LDT, (1.0E0,0.0E0), CSAV, LDT )
                        RES1 = CLANGE( 'M', M, N, CSAV, LDT, DUM )
                        RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, &
                              ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )
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
!
!     End of CGET35
!
END

