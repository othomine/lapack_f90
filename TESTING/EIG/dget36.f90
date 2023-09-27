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
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
   INTEGER            LDT, LWORK
   PARAMETER          ( LDT = 10, LWORK = 2*LDT*LDT )
!     ..
!     .. Local Scalars ..
   INTEGER            I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1, &
                      ILST2, ILSTSV, INFO1, INFO2, J, LOC, N
   DOUBLE PRECISION   EPS, RES
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
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, SIGN
!     ..
!     .. Executable Statements ..
!
   EPS = DLAMCH( 'P' )
   RMAX = ZERO
   LMAX = 0
   KNT = 0
   NINFO( 1 ) = 0
   NINFO( 2 ) = 0
   NINFO( 3 ) = 0
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
   CALL DLACPY( 'F', N, N, TMP, LDT, T1, LDT )
   CALL DLACPY( 'F', N, N, TMP, LDT, T2, LDT )
   IFSTSV = IFST
   ILSTSV = ILST
   IFST1 = IFST
   ILST1 = ILST
   IFST2 = IFST
   ILST2 = ILST
   RES = ZERO
!
!     Test without accumulating Q
!
   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
   CALL DTREXC( 'N', N, T1, LDT, Q, LDT, IFST1, ILST1, WORK, INFO1 )
   DO I = 1, N
      DO J = 1, N
         IF( I == J .AND. Q( I, J ) /= ONE ) &
            RES = RES + ONE / EPS
         IF( I /= J .AND. Q( I, J ) /= ZERO ) &
            RES = RES + ONE / EPS
      ENDDO
   ENDDO
!
!     Test with accumulating Q
!
   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
   CALL DTREXC( 'V', N, T2, LDT, Q, LDT, IFST2, ILST2, WORK, INFO2 )
!
!     Compare T1 with T2
!
   DO I = 1, N
      DO J = 1, N
         IF( T1( I, J ) /= T2( I, J ) ) &
            RES = RES + ONE / EPS
      ENDDO
   ENDDO
   IF( IFST1 /= IFST2 ) &
      RES = RES + ONE / EPS
   IF( ILST1 /= ILST2 ) &
      RES = RES + ONE / EPS
   IF( INFO1 /= INFO2 ) &
      RES = RES + ONE / EPS
!
!     Test for successful reordering of T2
!
   IF( INFO2 /= 0 ) THEN
      NINFO( INFO2 ) = NINFO( INFO2 ) + 1
   ELSE
      IF( ABS( IFST2-IFSTSV ) > 1 ) &
         RES = RES + ONE / EPS
      IF( ABS( ILST2-ILSTSV ) > 1 ) &
         RES = RES + ONE / EPS
   END IF
!
!     Test for small residual, and orthogonality of Q
!
   CALL DHST01( N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, &
                RESULT )
   RES = RES + RESULT( 1 ) + RESULT( 2 )
!
!     Test for T2 being in Schur form
!
   LOC = 1
70 CONTINUE
   IF( T2( LOC+1, LOC ) /= ZERO ) THEN
!
!        2 by 2 block
!
      IF( T2( LOC, LOC+1 ) == ZERO .OR. T2( LOC, LOC ) /= &
          T2( LOC+1, LOC+1 ) .OR. SIGN( ONE, T2( LOC, LOC+1 ) ) == &
          SIGN( ONE, T2( LOC+1, LOC ) ) )RES = RES + ONE / EPS
      DO I = LOC + 2, N
         IF( T2( I, LOC ) /= ZERO ) &
            RES = RES + ONE / RES
         IF( T2( I, LOC+1 ) /= ZERO ) &
            RES = RES + ONE / RES
      ENDDO
      LOC = LOC + 2
   ELSE
!
!        1 by 1 block
!
      DO I = LOC + 1, N
         IF( T2( I, LOC ) /= ZERO ) &
            RES = RES + ONE / RES
      ENDDO
      LOC = LOC + 1
   END IF
   IF( LOC < N ) &
      GO TO 70
   IF( RES > RMAX ) THEN
      RMAX = RES
      LMAX = KNT
   END IF
   GO TO 10
!
!     End of DGET36
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
