!> \brief \b SGET38
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET38( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, NIN
!       ..
!       .. Array Arguments ..
!       INTEGER            LMAX( 3 ), NINFO( 3 )
!       REAL               RMAX( 3 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET38 tests STRSEN, a routine for estimating condition numbers of a
!> cluster of eigenvalues and/or its associated right invariant subspace
!>
!> The test matrices are read from a file with logical unit number NIN.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is REAL array, dimension (3)
!>          Values of the largest test ratios.
!>          RMAX(1) = largest residuals from SHST01 or comparing
!>                    different calls to STRSEN
!>          RMAX(2) = largest error in reciprocal condition
!>                    numbers taking their conditioning into account
!>          RMAX(3) = largest error in reciprocal condition
!>                    numbers not taking their conditioning into
!>                    account (may be larger than RMAX(2))
!> \endverbatim
!>
!> \param[out] LMAX
!> \verbatim
!>          LMAX is INTEGER array, dimension (3)
!>          LMAX(i) is example number where largest test ratio
!>          RMAX(i) is achieved. Also:
!>          If SGEHRD returns INFO nonzero on example i, LMAX(1)=i
!>          If SHSEQR returns INFO nonzero on example i, LMAX(2)=i
!>          If STRSEN returns INFO nonzero on example i, LMAX(3)=i
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = No. of times SGEHRD returned INFO nonzero
!>          NINFO(2) = No. of times SHSEQR returned INFO nonzero
!>          NINFO(3) = No. of times STRSEN returned INFO nonzero
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SGET38( RMAX, LMAX, NINFO, KNT, NIN )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KNT, NIN
!     ..
!     .. Array Arguments ..
   INTEGER            LMAX( 3 ), NINFO( 3 )
   REAL               RMAX( 3 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE, TWO
   PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
   REAL               EPSIN
   PARAMETER          ( EPSIN = 5.9605E-8 )
   INTEGER            LDT, LWORK
   PARAMETER          ( LDT = 20, LWORK = 2*LDT*( 10+LDT ) )
   INTEGER            LIWORK
   PARAMETER          ( LIWORK = LDT*LDT )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, ISCL, ITMP, J, KMIN, M, N, NDIM
   REAL               BIGNUM, EPS, S, SEP, SEPIN, SEPTMP, SIN, &
                      SMLNUM, STMP, TNRM, TOL, TOLIN, V, VIMIN, VMAX, &
                      VMUL, VRMIN
!     ..
!     .. Local Arrays ..
   LOGICAL            SELECT( LDT )
   INTEGER            IPNT( LDT ), ISELEC( LDT ), IWORK( LIWORK )
   REAL               Q( LDT, LDT ), QSAV( LDT, LDT ), &
                      QTMP( LDT, LDT ), RESULT( 2 ), T( LDT, LDT ), &
                      TMP( LDT, LDT ), TSAV( LDT, LDT ), &
                      TSAV1( LDT, LDT ), TTMP( LDT, LDT ), VAL( 3 ), &
                      WI( LDT ), WITMP( LDT ), WORK( LWORK ), &
                      WR( LDT ), WRTMP( LDT )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SLANGE
   EXTERNAL           SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SCOPY, SGEHRD, SHSEQR, SHST01, SLACPY, SORGHR, &
                      SSCAL, STRSEN
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, REAL, SQRT
!     ..
!     .. Executable Statements ..
!
   EPS = SLAMCH( 'P' )
   SMLNUM = SLAMCH( 'S' ) / EPS
   BIGNUM = ONE / SMLNUM
!
!     EPSIN = 2**(-24) = precision to which input data computed
!
   EPS = MAX( EPS, EPSIN )
   RMAX( 1 ) = ZERO
   RMAX( 2 ) = ZERO
   RMAX( 3 ) = ZERO
   LMAX( 1 ) = 0
   LMAX( 2 ) = 0
   LMAX( 3 ) = 0
   KNT = 0
   NINFO( 1 ) = 0
   NINFO( 2 ) = 0
   NINFO( 3 ) = 0
!
   VAL( 1 ) = SQRT( SMLNUM )
   VAL( 2 ) = ONE
   VAL( 3 ) = SQRT( SQRT( BIGNUM ) )
!
!     Read input data until N=0.  Assume input eigenvalues are sorted
!     lexicographically (increasing by real part, then decreasing by
!     imaginary part)
!
10 CONTINUE
   READ( NIN, FMT = * )N, NDIM
   IF( N == 0 ) &
      RETURN
   READ( NIN, FMT = * )( ISELEC( I ), I = 1, NDIM )
   DO I = 1, N
      READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   ENDDO
   READ( NIN, FMT = * )SIN, SEPIN
!
   TNRM = SLANGE( 'M', N, N, TMP, LDT, WORK )
   DO ISCL = 1, 3
!
!        Scale input matrix
!
      KNT = KNT + 1
      CALL SLACPY( 'F', N, N, TMP, LDT, T, LDT )
      VMUL = VAL( ISCL )
      DO I = 1, N
         CALL SSCAL( N, VMUL, T( 1, I ), 1 )
      ENDDO
      IF( TNRM == ZERO ) &
         VMUL = ONE
      CALL SLACPY( 'F', N, N, T, LDT, TSAV, LDT )
!
!        Compute Schur form
!
      CALL SGEHRD( N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, &
                   INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 1 ) = KNT
         NINFO( 1 ) = NINFO( 1 ) + 1
         GO TO 160
      END IF
!
!        Generate orthogonal matrix
!
      CALL SLACPY( 'L', N, N, T, LDT, Q, LDT )
      CALL SORGHR( N, 1, N, Q, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, &
                   INFO )
!
!        Compute Schur form
!
      CALL SHSEQR( 'S', 'V', N, 1, N, T, LDT, WR, WI, Q, LDT, WORK, &
                   LWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 2 ) = KNT
         NINFO( 2 ) = NINFO( 2 ) + 1
         GO TO 160
      END IF
!
!        Sort, select eigenvalues
!
      DO I = 1, N
         IPNT( I ) = I
         SELECT( I ) = .FALSE.
      ENDDO
      CALL SCOPY( N, WR, 1, WRTMP, 1 )
      CALL SCOPY( N, WI, 1, WITMP, 1 )
      DO I = 1, N - 1
         KMIN = I
         VRMIN = WRTMP( I )
         VIMIN = WITMP( I )
         DO J = I + 1, N
            IF( WRTMP( J ) < VRMIN ) THEN
               KMIN = J
               VRMIN = WRTMP( J )
               VIMIN = WITMP( J )
            END IF
         ENDDO
         WRTMP( KMIN ) = WRTMP( I )
         WITMP( KMIN ) = WITMP( I )
         WRTMP( I ) = VRMIN
         WITMP( I ) = VIMIN
         ITMP = IPNT( I )
         IPNT( I ) = IPNT( KMIN )
         IPNT( KMIN ) = ITMP
      ENDDO
      DO I = 1, NDIM
         SELECT( IPNT( ISELEC( I ) ) ) = .TRUE.
      ENDDO
!
!        Compute condition numbers
!
      CALL SLACPY( 'F', N, N, Q, LDT, QSAV, LDT )
      CALL SLACPY( 'F', N, N, T, LDT, TSAV1, LDT )
      CALL STRSEN( 'B', 'V', SELECT, N, T, LDT, Q, LDT, WRTMP, WITMP, &
                   M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      SEPTMP = SEP / VMUL
      STMP = S
!
!        Compute residuals
!
      CALL SHST01( N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, &
                   RESULT )
      VMAX = MAX( RESULT( 1 ), RESULT( 2 ) )
      IF( VMAX > RMAX( 1 ) ) THEN
         RMAX( 1 ) = VMAX
         IF( NINFO( 1 ) == 0 ) &
            LMAX( 1 ) = KNT
      END IF
!
!        Compare condition number for eigenvalue cluster
!        taking its condition number into account
!
      V = MAX( TWO*REAL( N )*EPS*TNRM, SMLNUM )
      IF( TNRM == ZERO ) &
         V = ONE
      IF( V > SEPTMP ) THEN
         TOL = ONE
      ELSE
         TOL = V / SEPTMP
      END IF
      IF( V > SEPIN ) THEN
         TOLIN = ONE
      ELSE
         TOLIN = V / SEPIN
      END IF
      TOL = MAX( TOL, SMLNUM / EPS )
      TOLIN = MAX( TOLIN, SMLNUM / EPS )
      IF( EPS*( SIN-TOLIN ) > STMP+TOL ) THEN
         VMAX = ONE / EPS
      ELSE IF( SIN-TOLIN > STMP+TOL ) THEN
         VMAX = ( SIN-TOLIN ) / ( STMP+TOL )
      ELSE IF( SIN+TOLIN < EPS*( STMP-TOL ) ) THEN
         VMAX = ONE / EPS
      ELSE IF( SIN+TOLIN < STMP-TOL ) THEN
         VMAX = ( STMP-TOL ) / ( SIN+TOLIN )
      ELSE
         VMAX = ONE
      END IF
      IF( VMAX > RMAX( 2 ) ) THEN
         RMAX( 2 ) = VMAX
         IF( NINFO( 2 ) == 0 ) &
            LMAX( 2 ) = KNT
      END IF
!
!        Compare condition numbers for invariant subspace
!        taking its condition number into account
!
      IF( V > SEPTMP*STMP ) THEN
         TOL = SEPTMP
      ELSE
         TOL = V / STMP
      END IF
      IF( V > SEPIN*SIN ) THEN
         TOLIN = SEPIN
      ELSE
         TOLIN = V / SIN
      END IF
      TOL = MAX( TOL, SMLNUM / EPS )
      TOLIN = MAX( TOLIN, SMLNUM / EPS )
      IF( EPS*( SEPIN-TOLIN ) > SEPTMP+TOL ) THEN
         VMAX = ONE / EPS
      ELSE IF( SEPIN-TOLIN > SEPTMP+TOL ) THEN
         VMAX = ( SEPIN-TOLIN ) / ( SEPTMP+TOL )
      ELSE IF( SEPIN+TOLIN < EPS*( SEPTMP-TOL ) ) THEN
         VMAX = ONE / EPS
      ELSE IF( SEPIN+TOLIN < SEPTMP-TOL ) THEN
         VMAX = ( SEPTMP-TOL ) / ( SEPIN+TOLIN )
      ELSE
         VMAX = ONE
      END IF
      IF( VMAX > RMAX( 2 ) ) THEN
         RMAX( 2 ) = VMAX
         IF( NINFO( 2 ) == 0 ) &
            LMAX( 2 ) = KNT
      END IF
!
!        Compare condition number for eigenvalue cluster
!        without taking its condition number into account
!
      IF( SIN <= REAL( 2*N )*EPS .AND. STMP <= REAL( 2*N )*EPS ) THEN
         VMAX = ONE
      ELSE IF( EPS*SIN > STMP ) THEN
         VMAX = ONE / EPS
      ELSE IF( SIN > STMP ) THEN
         VMAX = SIN / STMP
      ELSE IF( SIN < EPS*STMP ) THEN
         VMAX = ONE / EPS
      ELSE IF( SIN < STMP ) THEN
         VMAX = STMP / SIN
      ELSE
         VMAX = ONE
      END IF
      IF( VMAX > RMAX( 3 ) ) THEN
         RMAX( 3 ) = VMAX
         IF( NINFO( 3 ) == 0 ) &
            LMAX( 3 ) = KNT
      END IF
!
!        Compare condition numbers for invariant subspace
!        without taking its condition number into account
!
      IF( SEPIN <= V .AND. SEPTMP <= V ) THEN
         VMAX = ONE
      ELSE IF( EPS*SEPIN > SEPTMP ) THEN
         VMAX = ONE / EPS
      ELSE IF( SEPIN > SEPTMP ) THEN
         VMAX = SEPIN / SEPTMP
      ELSE IF( SEPIN < EPS*SEPTMP ) THEN
         VMAX = ONE / EPS
      ELSE IF( SEPIN < SEPTMP ) THEN
         VMAX = SEPTMP / SEPIN
      ELSE
         VMAX = ONE
      END IF
      IF( VMAX > RMAX( 3 ) ) THEN
         RMAX( 3 ) = VMAX
         IF( NINFO( 3 ) == 0 ) &
            LMAX( 3 ) = KNT
      END IF
!
!        Compute eigenvalue condition number only and compare
!        Update Q
!
      VMAX = ZERO
      CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -ONE
      STMP = -ONE
      CALL STRSEN( 'E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( S /= STMP ) &
         VMAX = ONE / EPS
      IF( -ONE /= SEPTMP ) &
         VMAX = ONE / EPS
      DO I = 1, N
         DO J = 1, N
            IF( TTMP( I, J ) /= T( I, J ) ) &
               VMAX = ONE / EPS
            IF( QTMP( I, J ) /= Q( I, J ) ) &
               VMAX = ONE / EPS
         ENDDO
      ENDDO
!
!        Compute invariant subspace condition number only and compare
!        Update Q
!
      CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -ONE
      STMP = -ONE
      CALL STRSEN( 'V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( -ONE /= STMP ) &
         VMAX = ONE / EPS
      IF( SEP /= SEPTMP ) &
         VMAX = ONE / EPS
      DO I = 1, N
         DO J = 1, N
            IF( TTMP( I, J ) /= T( I, J ) ) &
               VMAX = ONE / EPS
            IF( QTMP( I, J ) /= Q( I, J ) ) &
               VMAX = ONE / EPS
            ENDDO
         ENDDO
!
!        Compute eigenvalue condition number only and compare
!        Do not update Q
!
      CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -ONE
      STMP = -ONE
      CALL STRSEN( 'E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( S /= STMP ) &
         VMAX = ONE / EPS
      IF( -ONE /= SEPTMP ) &
         VMAX = ONE / EPS
      DO I = 1, N
         DO J = 1, N
            IF( TTMP( I, J ) /= T( I, J ) ) &
               VMAX = ONE / EPS
            IF( QTMP( I, J ) /= QSAV( I, J ) ) &
               VMAX = ONE / EPS
            ENDDO
         ENDDO
!
!        Compute invariant subspace condition number only and compare
!        Do not update Q
!
      CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -ONE
      STMP = -ONE
      CALL STRSEN( 'V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( -ONE /= STMP ) &
         VMAX = ONE / EPS
      IF( SEP /= SEPTMP ) &
         VMAX = ONE / EPS
      DO I = 1, N
         DO J = 1, N
            IF( TTMP( I, J ) /= T( I, J ) ) &
               VMAX = ONE / EPS
            IF( QTMP( I, J ) /= QSAV( I, J ) ) &
               VMAX = ONE / EPS
            ENDDO
         ENDDO
      IF( VMAX > RMAX( 1 ) ) THEN
         RMAX( 1 ) = VMAX
         IF( NINFO( 1 ) == 0 ) &
            LMAX( 1 ) = KNT
      END IF
  160 CONTINUE
      ENDDO
   GO TO 10
!
!     End of SGET38
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
