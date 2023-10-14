!> \brief \b DGET38
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET38( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, NIN
!       ..
!       .. Array Arguments ..
!       INTEGER            LMAX( 3 ), NINFO( 3 )
!       DOUBLE PRECISION   RMAX( 3 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET38 tests DTRSEN, a routine for estimating condition numbers of a
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
!>          RMAX is DOUBLE PRECISION array, dimension (3)
!>          Values of the largest test ratios.
!>          RMAX(1) = largest residuals from DHST01 or comparing
!>                    different calls to DTRSEN
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
!>          If DGEHRD returns INFO nonzero on example i, LMAX(1)=i
!>          If DHSEQR returns INFO nonzero on example i, LMAX(2)=i
!>          If DTRSEN returns INFO nonzero on example i, LMAX(3)=i
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = No. of times DGEHRD returned INFO nonzero
!>          NINFO(2) = No. of times DHSEQR returned INFO nonzero
!>          NINFO(3) = No. of times DTRSEN returned INFO nonzero
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
   SUBROUTINE DGET38( RMAX, LMAX, NINFO, KNT, NIN )
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
   DOUBLE PRECISION   RMAX( 3 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   EPSIN
   PARAMETER          ( EPSIN = 5.9605D-8 )
   INTEGER            LDT, LWORK
   PARAMETER          ( LDT = 20, LWORK = 2*LDT*( 10+LDT ) )
   INTEGER            LIWORK
   PARAMETER          ( LIWORK = LDT*LDT )
!     ..
!     .. Local Scalars ..
   INTEGER            I, INFO, ISCL, ITMP, J, KMIN, M, N, NDIM
   DOUBLE PRECISION   BIGNUM, EPS, S, SEP, SEPIN, SEPTMP, SIN, &
                      SMLNUM, STMP, TNRM, TOL, TOLIN, V, VIMIN, VMAX, &
                      VMUL, VRMIN
!     ..
!     .. Local Arrays ..
   LOGICAL            SELECT( LDT )
   INTEGER            IPNT( LDT ), ISELEC( LDT ), IWORK( LIWORK )
   DOUBLE PRECISION   Q( LDT, LDT ), QSAV( LDT, LDT ), &
                      QTMP( LDT, LDT ), RESULT( 2 ), T( LDT, LDT ), &
                      TMP( LDT, LDT ), TSAV( LDT, LDT ), &
                      TSAV1( LDT, LDT ), TTMP( LDT, LDT ), VAL( 3 ), &
                      WI( LDT ), WITMP( LDT ), WORK( LWORK ), &
                      WR( LDT ), WRTMP( LDT )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLANGE
   EXTERNAL           DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DGEHRD, DHSEQR, DHST01, DLACPY, DORGHR, &
                      DSCAL, DTRSEN
!     ..
!     .. Executable Statements ..
!
   EPS = DLAMCH( 'P' )
   SMLNUM = DLAMCH( 'S' ) / EPS
   BIGNUM = 1.0D0 / SMLNUM
!
!     EPSIN = 2**(-24) = precision to which input data computed
!
   EPS = MAX( EPS, EPSIN )
   RMAX( 1 ) = 0.0D0
   RMAX( 2 ) = 0.0D0
   RMAX( 3 ) = 0.0D0
   LMAX( 1 ) = 0
   LMAX( 2 ) = 0
   LMAX( 3 ) = 0
   KNT = 0
   NINFO( 1 ) = 0
   NINFO( 2 ) = 0
   NINFO( 3 ) = 0
!
   VAL( 1 ) = SQRT( SMLNUM )
   VAL( 2 ) = 1.0D0
   VAL( 3 ) = SQRT( SQRT( BIGNUM ) )
!
!     Read input data until N=0.  Assume input eigenvalues are sorted
!     lexicographically (increasing by real part, then decreasing by
!     imaginary part)
!
   DO
   READ(NIN,*) N, NDIM
   IF( N == 0 ) &
      RETURN
   READ(NIN,*) ISELEC(1:NDIM)
   DO I = 1, N
      READ(NIN,*) TMP(I,1:N)
   ENDDO
   READ(NIN,*) SIN, SEPIN
!
   TNRM = DLANGE( 'M', N, N, TMP, LDT, WORK )
   DO ISCL = 1, 3
!
!        Scale input matrix
!
      KNT = KNT + 1
      CALL DLACPY( 'F', N, N, TMP, LDT, T, LDT )
      VMUL = VAL( ISCL )
      DO I = 1, N
         CALL DSCAL( N, VMUL, T( 1, I ), 1 )
      ENDDO
      IF( TNRM == 0.0D0 ) &
         VMUL = 1.0D0
      CALL DLACPY( 'F', N, N, T, LDT, TSAV, LDT )
!
!        Compute Schur form
!
      CALL DGEHRD( N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, &
                   INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 1 ) = KNT
         NINFO( 1 ) = NINFO( 1 ) + 1
         GO TO 160
      END IF
!
!        Generate orthogonal matrix
!
      CALL DLACPY( 'L', N, N, T, LDT, Q, LDT )
      CALL DORGHR( N, 1, N, Q, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, &
                   INFO )
!
!        Compute Schur form
!
      CALL DHSEQR( 'S', 'V', N, 1, N, T, LDT, WR, WI, Q, LDT, WORK, &
                   LWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 2 ) = KNT
         NINFO( 2 ) = NINFO( 2 ) + 1
         GO TO 160
      END IF
!
!        Sort, select eigenvalues
!
      DO I = 1,N
         IPNT(I) = I
      ENDDO
      SELECT(1:N) = .FALSE.
      CALL DCOPY( N, WR, 1, WRTMP, 1 )
      CALL DCOPY( N, WI, 1, WITMP, 1 )
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
      SELECT( IPNT( ISELEC(1:NDIM) ) ) = .TRUE.
!
!        Compute condition numbers
!
      CALL DLACPY( 'F', N, N, Q, LDT, QSAV, LDT )
      CALL DLACPY( 'F', N, N, T, LDT, TSAV1, LDT )
      CALL DTRSEN( 'B', 'V', SELECT, N, T, LDT, Q, LDT, WRTMP, WITMP, &
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
      CALL DHST01( N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, &
                   RESULT )
      VMAX = MAX( RESULT( 1 ), RESULT( 2 ) )
      IF( VMAX > RMAX( 1 ) ) THEN
         RMAX( 1 ) = VMAX
         IF( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT
      END IF
!
!        Compare condition number for eigenvalue cluster
!        taking its condition number into account
!
      V = MAX( 2.0D0*DBLE( N )*EPS*TNRM, SMLNUM )
      IF( TNRM == 0.0D0 ) V = 1.0D0
      IF( V > SEPTMP ) THEN
         TOL = 1.0D0
      ELSE
         TOL = V / SEPTMP
      END IF
      IF( V > SEPIN ) THEN
         TOLIN = 1.0D0
      ELSE
         TOLIN = V / SEPIN
      END IF
      TOL = MAX( TOL, SMLNUM / EPS )
      TOLIN = MAX( TOLIN, SMLNUM / EPS )
      IF( EPS*( SIN-TOLIN ) > STMP+TOL ) THEN
         VMAX = 1.0D0 / EPS
      ELSE IF( SIN-TOLIN > STMP+TOL ) THEN
         VMAX = ( SIN-TOLIN ) / ( STMP+TOL )
      ELSE IF( SIN+TOLIN < EPS*( STMP-TOL ) ) THEN
         VMAX = 1.0D0 / EPS
      ELSE IF( SIN+TOLIN < STMP-TOL ) THEN
         VMAX = ( STMP-TOL ) / ( SIN+TOLIN )
      ELSE
         VMAX = 1.0D0
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
         VMAX = 1.0D0 / EPS
      ELSE IF( SEPIN-TOLIN > SEPTMP+TOL ) THEN
         VMAX = ( SEPIN-TOLIN ) / ( SEPTMP+TOL )
      ELSE IF( SEPIN+TOLIN < EPS*( SEPTMP-TOL ) ) THEN
         VMAX = 1.0D0 / EPS
      ELSE IF( SEPIN+TOLIN < SEPTMP-TOL ) THEN
         VMAX = ( SEPTMP-TOL ) / ( SEPIN+TOLIN )
      ELSE
         VMAX = 1.0D0
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
      IF( SIN <= DBLE( 2*N )*EPS .AND. STMP <= DBLE( 2*N )*EPS ) THEN
         VMAX = 1.0D0
      ELSE IF( EPS*SIN > STMP ) THEN
         VMAX = 1.0D0 / EPS
      ELSE IF( SIN > STMP ) THEN
         VMAX = SIN / STMP
      ELSE IF( SIN < EPS*STMP ) THEN
         VMAX = 1.0D0 / EPS
      ELSE IF( SIN < STMP ) THEN
         VMAX = STMP / SIN
      ELSE
         VMAX = 1.0D0
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
         VMAX = 1.0D0
      ELSE IF( EPS*SEPIN > SEPTMP ) THEN
         VMAX = 1.0D0 / EPS
      ELSE IF( SEPIN > SEPTMP ) THEN
         VMAX = SEPIN / SEPTMP
      ELSE IF( SEPIN < EPS*SEPTMP ) THEN
         VMAX = 1.0D0 / EPS
      ELSE IF( SEPIN < SEPTMP ) THEN
         VMAX = SEPTMP / SEPIN
      ELSE
         VMAX = 1.0D0
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
      VMAX = 0.0D0
      CALL DLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL DLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -1.0D0
      STMP = -1.0D0
      CALL DTRSEN( 'E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( S /= STMP ) VMAX = 1.0D0 / EPS
      IF( -1.0D0 /= SEPTMP ) VMAX = 1.0D0 / EPS
      IF (ANY(TTMP(1:N,1:N) /= T(1:N,1:N))) VMAX = 1.0D0 / EPS
      IF (ANY(QTMP(1:N,1:N) /= Q(1:N,1:N))) VMAX = 1.0D0 / EPS
!
!        Compute invariant subspace condition number only and compare
!        Update Q
!
      CALL DLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL DLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -1.0D0
      STMP = -1.0D0
      CALL DTRSEN( 'V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( -1.0D0 /= STMP ) VMAX = 1.0D0 / EPS
      IF( SEP /= SEPTMP ) VMAX = 1.0D0 / EPS
      IF (ANY(TTMP(1:N,1:N) /= T(1:N,1:N))) VMAX = 1.0D0 / EPS
      IF (ANY(QTMP(1:N,1:N) /= Q(1:N,1:N))) VMAX = 1.0D0 / EPS
!
!        Compute eigenvalue condition number only and compare
!        Do not update Q
!
      CALL DLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL DLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -1.0D0
      STMP = -1.0D0
      CALL DTRSEN( 'E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( S /= STMP ) &
         VMAX = 1.0D0 / EPS
      IF( -1.0D0 /= SEPTMP ) VMAX = 1.0D0 / EPS
      IF (ANY(TTMP(1:N,1:N) /= T(1:N,1:N))) VMAX = 1.0D0 / EPS
      IF (ANY(QTMP(1:N,1:N) /= QSAV(1:N,1:N))) VMAX = 1.0D0 / EPS
!
!        Compute invariant subspace condition number only and compare
!        Do not update Q
!
      CALL DLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
      CALL DLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
      SEPTMP = -1.0D0
      STMP = -1.0D0
      CALL DTRSEN( 'V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, &
                   WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, &
                   LIWORK, INFO )
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 160
      END IF
      IF( -1.0D0 /= STMP ) VMAX = 1.0D0 / EPS
      IF( SEP /= SEPTMP ) VMAX = 1.0D0 / EPS
      IF (ANY(TTMP(1:N,1:N) /= T(1:N,1:N))) VMAX = 1.0D0 / EPS
      IF (ANY(QTMP(1:N,1:N) /= QSAV(1:N,1:N))) VMAX = 1.0D0 / EPS
      IF( VMAX > RMAX( 1 ) ) THEN
         RMAX( 1 ) = VMAX
         IF( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT
      END IF
  160 CONTINUE
      ENDDO
   ENDDO
!
!     End of DGET38
!
END

