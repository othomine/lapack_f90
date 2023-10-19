!> \brief \b ZGET35
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET35( RMAX, LMAX, NINFO, KNT, NIN )
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
!> ZGET35 tests ZTRSYL, a routine for solving the Sylvester matrix
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
   SUBROUTINE ZGET35( RMAX, LMAX, NINFO, KNT, NIN )
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
   INTEGER            LDT
   PARAMETER          ( LDT = 10 )
   DOUBLE PRECISION   LARGE
   PARAMETER          ( LARGE = 1.0D6 )
!     ..
!     .. Local Scalars ..
   CHARACTER          TRANA, TRANB
   INTEGER            I, IMLA, IMLAD, IMLB, IMLC, INFO, ISGN, ITRANA, &
                      ITRANB, J, M, N
   DOUBLE PRECISION   BIGNUM, EPS, RES, RES1, SCALE, SMLNUM, TNRM, &
                      XNRM
   COMPLEX*16         RMUL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   DUM( 1 ), VM1( 3 ), VM2( 3 )
   COMPLEX*16         A( LDT, LDT ), ATMP( LDT, LDT ), B( LDT, LDT ), &
                      BTMP( LDT, LDT ), C( LDT, LDT ), &
                      CSAV( LDT, LDT ), CTMP( LDT, LDT )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, ZLANGE
   EXTERNAL           DLAMCH, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZGEMM, ZTRSYL
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
   VM1( 1 ) = SQRT( SMLNUM )
   VM1( 2 ) = 1.0D0
   VM1( 3 ) = LARGE
   VM2( 1 ) = 1.0D0
   VM2( 2 ) = 1.0D0 + 2.0D+0*EPS
   VM2( 3 ) = 2.0D+0
!
   KNT = 0
   NINFO = 0
   LMAX = 0
   RMAX = 0.0D0
!
!     Begin test loop
!
   DO
   READ(NIN,*)M, N
   IF( N == 0 ) RETURN
   DO I = 1, M
      READ(NIN,*) ATMP( I,1:M)
   ENDDO
   DO I = 1, N
      READ(NIN,*) BTMP( I,1:N)
   ENDDO
   DO I = 1, M
      READ(NIN,*) CTMP( I,1:N)
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
                        TNRM = 0.0D0
                        DO I = 1, M
                           DO J = 1, M
                              A( I, J ) = ATMP( I, J )*VM1( IMLA )
                              TNRM = MAX( TNRM, ABS( A( I, J ) ) )
                           ENDDO
                           A( I, I ) = A( I, I )*VM2( IMLAD )
                           TNRM = MAX( TNRM, ABS( A( I, I ) ) )
                        ENDDO
                        DO I = 1, N
                           DO J = 1, N
                              B( I, J ) = BTMP( I, J )*VM1( IMLB )
                              TNRM = MAX( TNRM, ABS( B( I, J ) ) )
                           ENDDO
                        ENDDO
                        IF( TNRM == 0.0D0 ) &
                           TNRM = 1.0D0
                        DO I = 1, M
                           DO J = 1, N
                              C( I, J ) = CTMP( I, J )*VM1( IMLC )
                              CSAV( I, J ) = C( I, J )
                           ENDDO
                           ENDDO
                        KNT = KNT + 1
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL ZTRSYL( TRANA, TRANB, ISGN, M, N, A, &
                                     LDT, B, LDT, C, LDT, SCALE, &
                                     INFO )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : ZTRSYL : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        IF( INFO /= 0 ) &
                           NINFO = NINFO + 1
                        XNRM = ZLANGE( 'M', M, N, C, LDT, DUM )
                        RMUL = (1.0D0,0.0D0)
                        IF( XNRM > 1.0D0 .AND. TNRM > 1.0D0 ) THEN
                           IF( XNRM > BIGNUM / TNRM ) THEN
                              RMUL = MAX( XNRM, TNRM )
                              RMUL = (1.0D0,0.0D0) / RMUL
                           END IF
                        END IF
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL ZGEMM( TRANA, 'N', M, N, M, RMUL, A, &
                                    LDT, C, LDT, -SCALE*RMUL, CSAV, &
                                    LDT )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                        CALL ZGEMM( 'N', TRANB, M, N, N, &
                                    DBLE( ISGN )*RMUL, C, LDT, B, &
                                    LDT, (1.0D0,0.0D0), CSAV, LDT )
#ifdef _TIMER
                        call system_clock(count_rate=nb_periods_sec,count=S2_time)
                        open(file='results.out', unit=10, position = 'append')
                        write(10,'(A,F16.10,A)') 'Total time : ZGEMM : ',&
                              real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                        close(10)
#endif
                        RES1 = ZLANGE( 'M', M, N, CSAV, LDT, DUM )
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
!     End of ZGET35
!
END




