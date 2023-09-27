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
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
   COMPLEX*16         CZERO, CONE
   PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), &
                      CONE = ( 1.0D+0, 0.0D+0 ) )
   INTEGER            LDT, LWORK
   PARAMETER          ( LDT = 10, LWORK = 2*LDT*LDT )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IFST, ILST, INFO1, INFO2, J, N
   DOUBLE PRECISION   EPS, RES
   COMPLEX*16         CTEMP
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
   RMAX = ZERO
   LMAX = 0
   KNT = 0
   NINFO = 0
!
!     Read input data until N=0
!
10 CONTINUE
   READ( NIN, FMT = * )N, IFST, ILST
   IF( N == 0 ) &
      RETURN
   KNT = KNT + 1
   DO I = 1, N
      READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   ENDDO
   CALL ZLACPY( 'F', N, N, TMP, LDT, T1, LDT )
   CALL ZLACPY( 'F', N, N, TMP, LDT, T2, LDT )
   RES = ZERO
!
!     Test without accumulating Q
!
   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDT )
   CALL ZTREXC( 'N', N, T1, LDT, Q, LDT, IFST, ILST, INFO1 )
   DO I = 1, N
      DO J = 1, N
         IF( I == J .AND. Q( I, J ) /= CONE ) &
            RES = RES + ONE / EPS
         IF( I /= J .AND. Q( I, J ) /= CZERO ) &
            RES = RES + ONE / EPS
      ENDDO
   ENDDO
!
!     Test with accumulating Q
!
   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDT )
   CALL ZTREXC( 'V', N, T2, LDT, Q, LDT, IFST, ILST, INFO2 )
!
!     Compare T1 with T2
!
   DO I = 1, N
      DO J = 1, N
         IF( T1( I, J ) /= T2( I, J ) ) &
            RES = RES + ONE / EPS
      ENDDO
   ENDDO
   IF( INFO1 /= 0 .OR. INFO2 /= 0 ) &
      NINFO = NINFO + 1
   IF( INFO1 /= INFO2 ) &
      RES = RES + ONE / EPS
!
!     Test for successful reordering of T2
!
   CALL ZCOPY( N, TMP, LDT+1, DIAG, 1 )
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
         RES = RES + ONE / EPS
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
         IF( T2( I, J ) /= CZERO ) &
            RES = RES + ONE / EPS
         ENDDO
      ENDDO
   IF( RES > RMAX ) THEN
      RMAX = RES
      LMAX = KNT
   END IF
   GO TO 10
!
!     End of ZGET36
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
