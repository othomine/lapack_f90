!> \brief \b DGET37
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET37( RMAX, LMAX, NINFO, KNT, NIN )
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
!> DGET37 tests DTRSNA, a routine for estimating condition numbers of
!> eigenvalues and/or right eigenvectors of a matrix.
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
!>          Value of the largest test ratio.
!>          RMAX(1) = largest ratio comparing different calls to DTRSNA
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
!>          If DTRSNA returns INFO nonzero on example i, LMAX(3)=i
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = No. of times DGEHRD returned INFO nonzero
!>          NINFO(2) = No. of times DHSEQR returned INFO nonzero
!>          NINFO(3) = No. of times DTRSNA returned INFO nonzero
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
!>          Input logical unit number
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
   SUBROUTINE DGET37( RMAX, LMAX, NINFO, KNT, NIN )
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
!     ..
!     .. Local Scalars ..
   INTEGER            I, ICMP, IFND, INFO, ISCL, J, KMIN, M, N
   DOUBLE PRECISION   BIGNUM, EPS, SMLNUM, TNRM, TOL, TOLIN, V, &
                      VIMIN, VMAX, VMUL, VRMIN
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   LOGICAL            SELECT( LDT )
   INTEGER            IWORK( 2*LDT ), LCMP( 3 )
   DOUBLE PRECISION   DUM( 1 ), LE( LDT, LDT ), RE( LDT, LDT ), &
                      S( LDT ), SEP( LDT ), SEPIN( LDT ), &
                      SEPTMP( LDT ), SIN( LDT ), STMP( LDT ), &
                      T( LDT, LDT ), TMP( LDT, LDT ), VAL( 3 ), &
                      WI( LDT ), WIIN( LDT ), WITMP( LDT ), &
                      WORK( LWORK ), WR( LDT ), WRIN( LDT ), &
                      WRTMP( LDT )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLANGE
   EXTERNAL           DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DGEHRD, DHSEQR, DLACPY, DSCAL, DTREVC, &
                      DTRSNA
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
   VAL( 3 ) = SQRT( BIGNUM )
!
!     Read input data until N=0.  Assume input eigenvalues are sorted
!     lexicographically (increasing by real part, then decreasing by
!     imaginary part)
!
   DO
   READ(NIN,*) N
   IF( N == 0 ) &
      RETURN
   DO I = 1, N
      READ(NIN,*) TMP(I,1:N)
   ENDDO
   DO I = 1, N
      READ(NIN,*) WRIN( I ), WIIN( I ), SIN( I ), SEPIN( I )
   ENDDO
   TNRM = DLANGE( 'M', N, N, TMP, LDT, WORK )
!
!     Begin test
!
   DO ISCL = 1, 3
!
!        Scale input matrix
!
      KNT = KNT + 1
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLACPY( 'F', N, N, TMP, LDT, T, LDT )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      VMUL = VAL( ISCL )
      DO I = 1, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSCAL( N, VMUL, T( 1, I ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
      IF( TNRM == 0.0D0 ) VMUL = 1.0D0
!
!        Compute eigenvalues and eigenvectors
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGEHRD( N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGEHRD : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 1 ) = KNT
         NINFO( 1 ) = NINFO( 1 ) + 1
         GO TO 240
      END IF
      DO J = 1, N - 2
         T(J+2:N,J) = 0.0D0
      ENDDO
!
!        Compute Schur form
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DHSEQR( 'S', 'N', N, 1, N, T, LDT, WR, WI, DUM, 1, WORK, &
                   LWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DHSEQR : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 2 ) = KNT
         NINFO( 2 ) = NINFO( 2 ) + 1
         GO TO 240
      END IF
!
!        Compute eigenvectors
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTREVC( 'Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, N, M, WORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTREVC : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        Compute condition numbers
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, S, SEP, N, M, WORK, N, IWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
!
!        Sort eigenvalues and condition numbers lexicographically
!        to compare with inputs
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, WR, 1, WRTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, WI, 1, WITMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, S, 1, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, SEP, 1, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSCAL( N, 1.0D0 / VMUL, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
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
         VRMIN = STMP( KMIN )
         STMP( KMIN ) = STMP( I )
         STMP( I ) = VRMIN
         VRMIN = SEPTMP( KMIN )
         SEPTMP( KMIN ) = SEPTMP( I )
         SEPTMP( I ) = VRMIN
      ENDDO
!
!        Compare condition numbers for eigenvalues
!        taking their condition numbers into account
!
      V = MAX( 2.0D0*DBLE( N )*EPS*TNRM, SMLNUM )
      IF( TNRM == 0.0D0 ) V = 1.0D0
      DO I = 1, N
         IF( V > SEPTMP( I ) ) THEN
            TOL = 1.0D0
         ELSE
            TOL = V / SEPTMP( I )
         END IF
         IF( V > SEPIN( I ) ) THEN
            TOLIN = 1.0D0
         ELSE
            TOLIN = V / SEPIN( I )
         END IF
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         IF( EPS*( SIN( I )-TOLIN ) > STMP( I )+TOL ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SIN( I )-TOLIN > STMP( I )+TOL ) THEN
            VMAX = ( SIN( I )-TOLIN ) / ( STMP( I )+TOL )
         ELSE IF( SIN( I )+TOLIN < EPS*( STMP( I )-TOL ) ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SIN( I )+TOLIN < STMP( I )-TOL ) THEN
            VMAX = ( STMP( I )-TOL ) / ( SIN( I )+TOLIN )
         ELSE
            VMAX = 1.0D0
         END IF
         IF( VMAX > RMAX( 2 ) ) THEN
            RMAX( 2 ) = VMAX
            IF( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT
         END IF
      ENDDO
!
!        Compare condition numbers for eigenvectors
!        taking their condition numbers into account
!
      DO I = 1, N
         IF( V > SEPTMP( I )*STMP( I ) ) THEN
            TOL = SEPTMP( I )
         ELSE
            TOL = V / STMP( I )
         END IF
         IF( V > SEPIN( I )*SIN( I ) ) THEN
            TOLIN = SEPIN( I )
         ELSE
            TOLIN = V / SIN( I )
         END IF
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         IF( EPS*( SEPIN( I )-TOLIN ) > SEPTMP( I )+TOL ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SEPIN( I )-TOLIN > SEPTMP( I )+TOL ) THEN
            VMAX = ( SEPIN( I )-TOLIN ) / ( SEPTMP( I )+TOL )
         ELSE IF( SEPIN( I )+TOLIN < EPS*( SEPTMP( I )-TOL ) ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SEPIN( I )+TOLIN < SEPTMP( I )-TOL ) THEN
            VMAX = ( SEPTMP( I )-TOL ) / ( SEPIN( I )+TOLIN )
         ELSE
            VMAX = 1.0D0
         END IF
         IF( VMAX > RMAX( 2 ) ) THEN
            RMAX( 2 ) = VMAX
            IF( NINFO( 2 ) == 0 ) &
               LMAX( 2 ) = KNT
         END IF
         ENDDO
!
!        Compare condition numbers for eigenvalues
!        without taking their condition numbers into account
!
      DO I = 1, N
         IF( SIN( I ) <= DBLE( 2*N )*EPS .AND. STMP( I ) <= &
             DBLE( 2*N )*EPS ) THEN
            VMAX = 1.0D0
         ELSE IF( EPS*SIN( I ) > STMP( I ) ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SIN( I ) > STMP( I ) ) THEN
            VMAX = SIN( I ) / STMP( I )
         ELSE IF( SIN( I ) < EPS*STMP( I ) ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SIN( I ) < STMP( I ) ) THEN
            VMAX = STMP( I ) / SIN( I )
         ELSE
            VMAX = 1.0D0
         END IF
         IF( VMAX > RMAX( 3 ) ) THEN
            RMAX( 3 ) = VMAX
            IF( NINFO( 3 ) == 0 ) &
               LMAX( 3 ) = KNT
         END IF
         ENDDO
!
!        Compare condition numbers for eigenvectors
!        without taking their condition numbers into account
!
      DO I = 1, N
         IF( SEPIN( I ) <= V .AND. SEPTMP( I ) <= V ) THEN
            VMAX = 1.0D0
         ELSE IF( EPS*SEPIN( I ) > SEPTMP( I ) ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SEPIN( I ) > SEPTMP( I ) ) THEN
            VMAX = SEPIN( I ) / SEPTMP( I )
         ELSE IF( SEPIN( I ) < EPS*SEPTMP( I ) ) THEN
            VMAX = 1.0D0 / EPS
         ELSE IF( SEPIN( I ) < SEPTMP( I ) ) THEN
            VMAX = SEPTMP( I ) / SEPIN( I )
         ELSE
            VMAX = 1.0D0
         END IF
         IF( VMAX > RMAX( 3 ) ) THEN
            RMAX( 3 ) = VMAX
            IF( NINFO( 3 ) == 0 ) &
               LMAX( 3 ) = KNT
         END IF
         ENDDO
!
!        Compute eigenvalue condition numbers only and compare
!
      VMAX = 0.0D0
      DUM( 1 ) = -1.0D0
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Eigcond', 'All', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, N
         IF( STMP( I ) /= S( I ) ) &
            VMAX = 1.0D0 / EPS
         IF( SEPTMP( I ) /= DUM( 1 ) ) &
            VMAX = 1.0D0 / EPS
         ENDDO
!
!        Compute eigenvector condition numbers only and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Veccond', 'All', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, N
         IF( STMP( I ) /= DUM( 1 ) ) &
            VMAX = 1.0D0 / EPS
         IF( SEPTMP( I ) /= SEP( I ) ) &
            VMAX = 1.0D0 / EPS
         ENDDO
!
!        Compute all condition numbers using SELECT and compare
!
      DO I = 1, N
         SELECT( I ) = .TRUE.
         ENDDO
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, &
                   RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, N
         IF( SEPTMP( I ) /= SEP( I ) ) &
            VMAX = 1.0D0 / EPS
         IF( STMP( I ) /= S( I ) ) &
            VMAX = 1.0D0 / EPS
         ENDDO
!
!        Compute eigenvalue condition numbers using SELECT and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, N
         IF( STMP( I ) /= S( I ) ) &
            VMAX = 1.0D0 / EPS
         IF( SEPTMP( I ) /= DUM( 1 ) ) &
            VMAX = 1.0D0 / EPS
         ENDDO
!
!        Compute eigenvector condition numbers using SELECT and compare
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( N, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, N
         IF( STMP( I ) /= DUM( 1 ) ) &
            VMAX = 1.0D0 / EPS
         IF( SEPTMP( I ) /= SEP( I ) ) &
            VMAX = 1.0D0 / EPS
         ENDDO
      IF( VMAX > RMAX( 1 ) ) THEN
         RMAX( 1 ) = VMAX
         IF( NINFO( 1 ) == 0 ) &
            LMAX( 1 ) = KNT
      END IF
!
!        Select first real and first complex eigenvalue
!
      IF( WI( 1 ) == 0.0D0 ) THEN
         LCMP( 1 ) = 1
         IFND = 0
         DO I = 2, N
            IF( IFND == 1 .OR. WI( I ) == 0.0D0 ) THEN
               SELECT( I ) = .FALSE.
            ELSE
               IFND = 1
               LCMP( 2 ) = I
               LCMP( 3 ) = I + 1
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( N, RE( 1, I ), 1, RE( 1, 2 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( N, RE( 1, I+1 ), 1, RE( 1, 3 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( N, LE( 1, I ), 1, LE( 1, 2 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( N, LE( 1, I+1 ), 1, LE( 1, 3 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            END IF
            ENDDO
         IF( IFND == 0 ) THEN
            ICMP = 1
         ELSE
            ICMP = 3
         END IF
      ELSE
         LCMP( 1 ) = 1
         LCMP( 2 ) = 2
         IFND = 0
         DO I = 3, N
            IF( IFND == 1 .OR. WI( I ) /= 0.0D0 ) THEN
               SELECT( I ) = .FALSE.
            ELSE
               LCMP( 3 ) = I
               IFND = 1
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( N, RE( 1, I ), 1, RE( 1, 3 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DCOPY( N, LE( 1, I ), 1, LE( 1, 3 ), 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            END IF
            ENDDO
         IF( IFND == 0 ) THEN
            ICMP = 2
         ELSE
            ICMP = 3
         END IF
      END IF
!
!        Compute all selected condition numbers
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( ICMP, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( ICMP, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, &
                   RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, ICMP
         J = LCMP( I )
         IF( SEPTMP( I ) /= SEP( J ) ) &
            VMAX = 1.0D0 / EPS
         IF( STMP( I ) /= S( J ) ) &
            VMAX = 1.0D0 / EPS
         ENDDO
!
!        Compute selected eigenvalue condition numbers
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( ICMP, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( ICMP, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, ICMP
         J = LCMP( I )
         IF( STMP( I ) /= S( J ) ) &
            VMAX = 1.0D0 / EPS
         IF( SEPTMP( I ) /= DUM( 1 ) ) &
            VMAX = 1.0D0 / EPS
         ENDDO
!
!        Compute selected eigenvector condition numbers
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( ICMP, DUM, 0, STMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DCOPY( ICMP, DUM, 0, SEPTMP, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DTRSNA( 'Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, &
                   LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DTRSNA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      IF( INFO /= 0 ) THEN
         LMAX( 3 ) = KNT
         NINFO( 3 ) = NINFO( 3 ) + 1
         GO TO 240
      END IF
      DO I = 1, ICMP
         J = LCMP( I )
         IF( STMP( I ) /= DUM( 1 ) ) VMAX = 1.0D0 / EPS
         IF( SEPTMP( I ) /= SEP( J ) ) VMAX = 1.0D0 / EPS
         ENDDO
      IF( VMAX > RMAX( 1 ) ) THEN
         RMAX( 1 ) = VMAX
         IF( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT
      END IF
  240 CONTINUE
      ENDDO
   ENDDO
!
!     End of DGET37
!
END



