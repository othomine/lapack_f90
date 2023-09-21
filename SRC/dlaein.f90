!> \brief \b DLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B,
!                          LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            NOINIT, RIGHTV
!       INTEGER            INFO, LDB, LDH, N
!       DOUBLE PRECISION   BIGNUM, EPS3, SMLNUM, WI, WR
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   B( LDB, * ), H( LDH, * ), VI( * ), VR( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEIN uses inverse iteration to find a right or left eigenvector
!> corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
!> matrix H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RIGHTV
!> \verbatim
!>          RIGHTV is LOGICAL
!>          = .TRUE. : compute right eigenvector;
!>          = .FALSE.: compute left eigenvector.
!> \endverbatim
!>
!> \param[in] NOINIT
!> \verbatim
!>          NOINIT is LOGICAL
!>          = .TRUE. : no initial vector supplied in (VR,VI).
!>          = .FALSE.: initial vector supplied in (VR,VI).
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          The upper Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max(1,N).
!> \endverbatim
!>
!> \param[in] WR
!> \verbatim
!>          WR is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] WI
!> \verbatim
!>          WI is DOUBLE PRECISION
!>          The real and imaginary parts of the eigenvalue of H whose
!>          corresponding right or left eigenvector is to be computed.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[in,out] VI
!> \verbatim
!>          VI is DOUBLE PRECISION array, dimension (N)
!>          On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain
!>          a real starting vector for inverse iteration using the real
!>          eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI
!>          must contain the real and imaginary parts of a complex
!>          starting vector for inverse iteration using the complex
!>          eigenvalue (WR,WI); otherwise VR and VI need not be set.
!>          On exit, if WI = 0.0 (real eigenvalue), VR contains the
!>          computed real eigenvector; if WI.ne.0.0 (complex eigenvalue),
!>          VR and VI contain the real and imaginary parts of the
!>          computed complex eigenvector. The eigenvector is normalized
!>          so that the component of largest magnitude has magnitude 1;
!>          here the magnitude of a complex number (x,y) is taken to be
!>          |x| + |y|.
!>          VI is not referenced if WI = 0.0.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= N+1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[in] EPS3
!> \verbatim
!>          EPS3 is DOUBLE PRECISION
!>          A small machine-dependent value which is used to perturb
!>          close eigenvalues, and to replace zero pivots.
!> \endverbatim
!>
!> \param[in] SMLNUM
!> \verbatim
!>          SMLNUM is DOUBLE PRECISION
!>          A machine-dependent value close to the underflow threshold.
!> \endverbatim
!>
!> \param[in] BIGNUM
!> \verbatim
!>          BIGNUM is DOUBLE PRECISION
!>          A machine-dependent value close to the overflow threshold.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          = 1:  inverse iteration did not converge; VR is set to the
!>                last iterate, and so is VI if WI.ne.0.0.
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
!> \ingroup laein
!
!  =====================================================================
   SUBROUTINE DLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B, &
                      LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            NOINIT, RIGHTV
   INTEGER            INFO, LDB, LDH, N
   DOUBLE PRECISION   BIGNUM, EPS3, SMLNUM, WI, WR
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   B( LDB, * ), H( LDH, * ), VI( * ), VR( * ), &
                      WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE, TENTH
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TENTH = 1.0D-1 )
!     ..
!     .. Local Scalars ..
   CHARACTER          NORMIN, TRANS
   INTEGER            I, I1, I2, I3, IERR, ITS, J
   DOUBLE PRECISION   ABSBII, ABSBJJ, EI, EJ, GROWTO, NORM, NRMSML, &
                      REC, ROOTN, SCALE, TEMP, VCRIT, VMAX, VNORM, W, &
                      W1, X, XI, XR, Y
!     ..
!     .. External Functions ..
   INTEGER            IDAMAX
   DOUBLE PRECISION   DASUM, DLAPY2, DNRM2
   EXTERNAL           IDAMAX, DASUM, DLAPY2, DNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLADIV, DLATRS, DSCAL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
   INFO = 0
!
!     GROWTO is the threshold used in the acceptance test for an
!     eigenvector.
!
   ROOTN = SQRT( DBLE( N ) )
   GROWTO = TENTH / ROOTN
   NRMSML = MAX( ONE, EPS3*ROOTN )*SMLNUM
!
!     Form B = H - (WR,WI)*I (except that the subdiagonal elements and
!     the imaginary parts of the diagonal elements are not stored).
!
   DO J = 1, N
      DO I = 1, J - 1
         B( I, J ) = H( I, J )
      ENDDO
      B( J, J ) = H( J, J ) - WR
   ENDDO
!
   IF( WI == ZERO ) THEN
!
!        Real eigenvalue.
!
      IF( NOINIT ) THEN
!
!           Set initial vector.
!
         DO I = 1, N
            VR( I ) = EPS3
         ENDDO
      ELSE
!
!           Scale supplied initial vector.
!
         VNORM = DNRM2( N, VR, 1 )
         CALL DSCAL( N, ( EPS3*ROOTN ) / MAX( VNORM, NRMSML ), VR, &
                     1 )
      END IF
!
      IF( RIGHTV ) THEN
!
!           LU decomposition with partial pivoting of B, replacing zero
!           pivots by EPS3.
!
         DO I = 1, N - 1
            EI = H( I+1, I )
            IF( ABS( B( I, I ) ) < ABS( EI ) ) THEN
!
!                 Interchange rows and eliminate.
!
               X = B( I, I ) / EI
               B( I, I ) = EI
               DO J = I + 1, N
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - X*TEMP
                  B( I, J ) = TEMP
               ENDDO
            ELSE
!
!                 Eliminate without interchange.
!
               IF( B( I, I ) == ZERO ) &
                  B( I, I ) = EPS3
               X = EI / B( I, I )
               IF( X /= ZERO ) THEN
                  DO J = I + 1, N
                     B( I+1, J ) = B( I+1, J ) - X*B( I, J )
                  ENDDO
               END IF
            END IF
         ENDDO
         IF( B( N, N ) == ZERO ) &
            B( N, N ) = EPS3
!
         TRANS = 'N'
!
      ELSE
!
!           UL decomposition with partial pivoting of B, replacing zero
!           pivots by EPS3.
!
         DO J = N, 2, -1
            EJ = H( J, J-1 )
            IF( ABS( B( J, J ) ) < ABS( EJ ) ) THEN
!
!                 Interchange columns and eliminate.
!
               X = B( J, J ) / EJ
               B( J, J ) = EJ
               DO I = 1, J - 1
                  TEMP = B( I, J-1 )
                  B( I, J-1 ) = B( I, J ) - X*TEMP
                  B( I, J ) = TEMP
               ENDDO
            ELSE
!
!                 Eliminate without interchange.
!
               IF( B( J, J ) == ZERO ) &
                  B( J, J ) = EPS3
               X = EJ / B( J, J )
               IF( X /= ZERO ) THEN
                  DO I = 1, J - 1
                     B( I, J-1 ) = B( I, J-1 ) - X*B( I, J )
                  ENDDO
               END IF
            END IF
         ENDDO
         IF( B( 1, 1 ) == ZERO ) &
            B( 1, 1 ) = EPS3
!
         TRANS = 'T'
!
      END IF
!
      NORMIN = 'N'
      DO ITS = 1, N
!
!           Solve U*x = scale*v for a right eigenvector
!             or U**T*x = scale*v for a left eigenvector,
!           overwriting x on v.
!
         CALL DLATRS( 'Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, &
                      VR, SCALE, WORK, IERR )
         NORMIN = 'Y'
!
!           Test for sufficient growth in the norm of v.
!
         VNORM = DASUM( N, VR, 1 )
         IF( VNORM >= GROWTO*SCALE ) &
            GO TO 120
!
!           Choose new orthogonal starting vector and try again.
!
         TEMP = EPS3 / ( ROOTN+ONE )
         VR( 1 ) = EPS3
         DO I = 2, N
            VR( I ) = TEMP
            ENDDO
         VR( N-ITS+1 ) = VR( N-ITS+1 ) - EPS3*ROOTN
         ENDDO
!
!        Failure to find eigenvector in N iterations.
!
      INFO = 1
!
  120    CONTINUE
!
!        Normalize eigenvector.
!
      I = IDAMAX( N, VR, 1 )
      CALL DSCAL( N, ONE / ABS( VR( I ) ), VR, 1 )
   ELSE
!
!        Complex eigenvalue.
!
      IF( NOINIT ) THEN
!
!           Set initial vector.
!
         DO I = 1, N
            VR( I ) = EPS3
            VI( I ) = ZERO
            ENDDO
      ELSE
!
!           Scale supplied initial vector.
!
         NORM = DLAPY2( DNRM2( N, VR, 1 ), DNRM2( N, VI, 1 ) )
         REC = ( EPS3*ROOTN ) / MAX( NORM, NRMSML )
         CALL DSCAL( N, REC, VR, 1 )
         CALL DSCAL( N, REC, VI, 1 )
      END IF
!
      IF( RIGHTV ) THEN
!
!           LU decomposition with partial pivoting of B, replacing zero
!           pivots by EPS3.
!
!           The imaginary part of the (i,j)-th element of U is stored in
!           B(j+1,i).
!
         B( 2, 1 ) = -WI
         DO I = 2, N
            B( I+1, 1 ) = ZERO
            ENDDO
!
         DO I = 1, N - 1
            ABSBII = DLAPY2( B( I, I ), B( I+1, I ) )
            EI = H( I+1, I )
            IF( ABSBII < ABS( EI ) ) THEN
!
!                 Interchange rows and eliminate.
!
               XR = B( I, I ) / EI
               XI = B( I+1, I ) / EI
               B( I, I ) = EI
               B( I+1, I ) = ZERO
               DO J = I + 1, N
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - XR*TEMP
                  B( J+1, I+1 ) = B( J+1, I ) - XI*TEMP
                  B( I, J ) = TEMP
                  B( J+1, I ) = ZERO
                  ENDDO
               B( I+2, I ) = -WI
               B( I+1, I+1 ) = B( I+1, I+1 ) - XI*WI
               B( I+2, I+1 ) = B( I+2, I+1 ) + XR*WI
            ELSE
!
!                 Eliminate without interchanging rows.
!
               IF( ABSBII == ZERO ) THEN
                  B( I, I ) = EPS3
                  B( I+1, I ) = ZERO
                  ABSBII = EPS3
               END IF
               EI = ( EI / ABSBII ) / ABSBII
               XR = B( I, I )*EI
               XI = -B( I+1, I )*EI
               DO J = I + 1, N
                  B( I+1, J ) = B( I+1, J ) - XR*B( I, J ) + &
                                XI*B( J+1, I )
                  B( J+1, I+1 ) = -XR*B( J+1, I ) - XI*B( I, J )
                  ENDDO
               B( I+2, I+1 ) = B( I+2, I+1 ) - WI
            END IF
!
!              Compute 1-norm of offdiagonal elements of i-th row.
!
            WORK( I ) = DASUM( N-I, B( I, I+1 ), LDB ) + &
                        DASUM( N-I, B( I+2, I ), 1 )
            ENDDO
         IF( B( N, N ) == ZERO .AND. B( N+1, N ) == ZERO ) &
            B( N, N ) = EPS3
         WORK( N ) = ZERO
!
         I1 = N
         I2 = 1
         I3 = -1
      ELSE
!
!           UL decomposition with partial pivoting of conjg(B),
!           replacing zero pivots by EPS3.
!
!           The imaginary part of the (i,j)-th element of U is stored in
!           B(j+1,i).
!
         B( N+1, N ) = WI
         DO J = 1, N - 1
            B( N+1, J ) = ZERO
            ENDDO
!
         DO J = N, 2, -1
            EJ = H( J, J-1 )
            ABSBJJ = DLAPY2( B( J, J ), B( J+1, J ) )
            IF( ABSBJJ < ABS( EJ ) ) THEN
!
!                 Interchange columns and eliminate
!
               XR = B( J, J ) / EJ
               XI = B( J+1, J ) / EJ
               B( J, J ) = EJ
               B( J+1, J ) = ZERO
               DO I = 1, J - 1
                  TEMP = B( I, J-1 )
                  B( I, J-1 ) = B( I, J ) - XR*TEMP
                  B( J, I ) = B( J+1, I ) - XI*TEMP
                  B( I, J ) = TEMP
                  B( J+1, I ) = ZERO
                  ENDDO
               B( J+1, J-1 ) = WI
               B( J-1, J-1 ) = B( J-1, J-1 ) + XI*WI
               B( J, J-1 ) = B( J, J-1 ) - XR*WI
            ELSE
!
!                 Eliminate without interchange.
!
               IF( ABSBJJ == ZERO ) THEN
                  B( J, J ) = EPS3
                  B( J+1, J ) = ZERO
                  ABSBJJ = EPS3
               END IF
               EJ = ( EJ / ABSBJJ ) / ABSBJJ
               XR = B( J, J )*EJ
               XI = -B( J+1, J )*EJ
               DO I = 1, J - 1
                  B( I, J-1 ) = B( I, J-1 ) - XR*B( I, J ) + &
                                XI*B( J+1, I )
                  B( J, I ) = -XR*B( J+1, I ) - XI*B( I, J )
                  ENDDO
               B( J, J-1 ) = B( J, J-1 ) + WI
            END IF
!
!              Compute 1-norm of offdiagonal elements of j-th column.
!
            WORK( J ) = DASUM( J-1, B( 1, J ), 1 ) + &
                        DASUM( J-1, B( J+1, 1 ), LDB )
            ENDDO
         IF( B( 1, 1 ) == ZERO .AND. B( 2, 1 ) == ZERO ) &
            B( 1, 1 ) = EPS3
         WORK( 1 ) = ZERO
!
         I1 = 1
         I2 = N
         I3 = 1
      END IF
!
      DO ITS = 1, N
         SCALE = ONE
         VMAX = ONE
         VCRIT = BIGNUM
!
!           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
!             or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector,
!           overwriting (xr,xi) on (vr,vi).
!
         DO I = I1, I2, I3
!
            IF( WORK( I ) > VCRIT ) THEN
               REC = ONE / VMAX
               CALL DSCAL( N, REC, VR, 1 )
               CALL DSCAL( N, REC, VI, 1 )
               SCALE = SCALE*REC
               VMAX = ONE
               VCRIT = BIGNUM
            END IF
!
            XR = VR( I )
            XI = VI( I )
            IF( RIGHTV ) THEN
               DO J = I + 1, N
                  XR = XR - B( I, J )*VR( J ) + B( J+1, I )*VI( J )
                  XI = XI - B( I, J )*VI( J ) - B( J+1, I )*VR( J )
                  ENDDO
            ELSE
               DO J = 1, I - 1
                  XR = XR - B( J, I )*VR( J ) + B( I+1, J )*VI( J )
                  XI = XI - B( J, I )*VI( J ) - B( I+1, J )*VR( J )
                  ENDDO
            END IF
!
            W = ABS( B( I, I ) ) + ABS( B( I+1, I ) )
            IF( W > SMLNUM ) THEN
               IF( W < ONE ) THEN
                  W1 = ABS( XR ) + ABS( XI )
                  IF( W1 > W*BIGNUM ) THEN
                     REC = ONE / W1
                     CALL DSCAL( N, REC, VR, 1 )
                     CALL DSCAL( N, REC, VI, 1 )
                     XR = VR( I )
                     XI = VI( I )
                     SCALE = SCALE*REC
                     VMAX = VMAX*REC
                  END IF
               END IF
!
!                 Divide by diagonal element of B.
!
               CALL DLADIV( XR, XI, B( I, I ), B( I+1, I ), VR( I ), &
                            VI( I ) )
               VMAX = MAX( ABS( VR( I ) )+ABS( VI( I ) ), VMAX )
               VCRIT = BIGNUM / VMAX
            ELSE
               DO J = 1, N
                  VR( J ) = ZERO
                  VI( J ) = ZERO
                  ENDDO
               VR( I ) = ONE
               VI( I ) = ONE
               SCALE = ZERO
               VMAX = ONE
               VCRIT = BIGNUM
            END IF
            ENDDO
!
!           Test for sufficient growth in the norm of (VR,VI).
!
         VNORM = DASUM( N, VR, 1 ) + DASUM( N, VI, 1 )
         IF( VNORM >= GROWTO*SCALE ) &
            GO TO 280
!
!           Choose a new orthogonal starting vector and try again.
!
         Y = EPS3 / ( ROOTN+ONE )
         VR( 1 ) = EPS3
         VI( 1 ) = ZERO
!
         DO I = 2, N
            VR( I ) = Y
            VI( I ) = ZERO
            ENDDO
         VR( N-ITS+1 ) = VR( N-ITS+1 ) - EPS3*ROOTN
         ENDDO
!
!        Failure to find eigenvector in N iterations
!
      INFO = 1
!
  280    CONTINUE
!
!        Normalize eigenvector.
!
      VNORM = ZERO
      DO I = 1, N
         VNORM = MAX( VNORM, ABS( VR( I ) )+ABS( VI( I ) ) )
         ENDDO
      CALL DSCAL( N, ONE / VNORM, VR, 1 )
      CALL DSCAL( N, ONE / VNORM, VI, 1 )
!
   END IF
!
   RETURN
!
!     End of DLAEIN
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
