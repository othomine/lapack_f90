!> \brief \b CLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK,
!                          EPS3, SMLNUM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            NOINIT, RIGHTV
!       INTEGER            INFO, LDB, LDH, N
!       REAL               EPS3, SMLNUM
!       COMPLEX            W
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            B( LDB, * ), H( LDH, * ), V( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAEIN uses inverse iteration to find a right or left eigenvector
!> corresponding to the eigenvalue W of a complex upper Hessenberg
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
!>          = .TRUE. : no initial vector supplied in V
!>          = .FALSE.: initial vector supplied in V.
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
!>          H is COMPLEX array, dimension (LDH,N)
!>          The upper Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max(1,N).
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is COMPLEX
!>          The eigenvalue of H whose corresponding right or left
!>          eigenvector is to be computed.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX array, dimension (N)
!>          On entry, if NOINIT = .FALSE., V must contain a starting
!>          vector for inverse iteration; otherwise V need not be set.
!>          On exit, V contains the computed eigenvector, normalized so
!>          that the component of largest magnitude has magnitude 1; here
!>          the magnitude of a complex number (x,y) is taken to be
!>          |x| + |y|.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[in] EPS3
!> \verbatim
!>          EPS3 is REAL
!>          A small machine-dependent value which is used to perturb
!>          close eigenvalues, and to replace zero pivots.
!> \endverbatim
!>
!> \param[in] SMLNUM
!> \verbatim
!>          SMLNUM is REAL
!>          A machine-dependent value close to the underflow threshold.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          = 1:  inverse iteration did not converge; V is set to the
!>                last iterate.
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
!> \ingroup laein
!
!  =====================================================================
   SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK, &
                      EPS3, SMLNUM, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            NOINIT, RIGHTV
   INTEGER            INFO, LDB, LDH, N
   REAL               EPS3, SMLNUM
   COMPLEX            W
!     ..
!     .. Array Arguments ..
   REAL               RWORK( * )
   COMPLEX            B( LDB, * ), H( LDH, * ), V( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   CHARACTER          NORMIN, TRANS
   INTEGER            I, IERR, ITS, J
   REAL               GROWTO, NRMSML, ROOTN, RTEMP, SCALE, VNORM
   COMPLEX            CDUM, EI, EJ, TEMP, X, B_TMP( LDB )
!     ..
!     .. External Functions ..
   INTEGER            ICAMAX
   REAL               SCNRM2
   COMPLEX            CLADIV
   EXTERNAL           ICAMAX, SCNRM2, CLADIV
!     ..
!     .. External Subroutines ..
   EXTERNAL           CLATRS
!     ..
!     .. Statement Functions ..
   REAL               CABS1
!     ..
!     .. Statement Function definitions ..
   CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
   INFO = 0
!
!     GROWTO is the threshold used in the acceptance test for an
!     eigenvector.
!
   ROOTN = SQRT( REAL( N ) )
   GROWTO = 0.1E+0 / ROOTN
   NRMSML = MAX( 1.0E+0, EPS3*ROOTN )*SMLNUM
!
!     Form B = H - W*I (except that the subdiagonal elements are not
!     stored).
!
   DO J = 1, N
      B(1:J-1,J) = H(1:J-1,J)
      B( J, J ) = H( J, J ) - W
   ENDDO
!
   IF( NOINIT ) THEN
!
!        Initialize V.
!
      V(1:N) = EPS3
   ELSE
!
!        Scale supplied initial vector.
!
      VNORM = SCNRM2( N, V, 1 )
      V(1:N) = V(1:N)*( EPS3*ROOTN ) / MAX( VNORM, NRMSML )
   END IF
!
   IF( RIGHTV ) THEN
!
!        LU decomposition with partial pivoting of B, replacing zero
!        pivots by EPS3.
!
      DO I = 1, N - 1
         EI = H( I+1, I )
         IF( CABS1( B( I, I ) ) < CABS1( EI ) ) THEN
!
!              Interchange rows and eliminate.
!
            X = CLADIV( B( I, I ), EI )
            B( I, I ) = EI
            B_TMP(I+1:N) = B(I+1,I+1:N)
            B(I+1,I+1:N) = B(I,I+1:N) - X*B_TMP(I+1:N)
            B(I,I+1:N) = B_TMP(I+1:N)
         ELSE
!
!              Eliminate without interchange.
!
            IF( B( I, I ) == (0.0E+0,0.0E+0) ) B( I, I ) = EPS3
            X = CLADIV( EI, B( I, I ) )
            IF( X /= (0.0E+0,0.0E+0) ) B(I+1,I+1:N) = B( I+1,I+1:N) - X*B( I,I+1:N)
         END IF
      ENDDO
      IF( B( N, N ) == (0.0E+0,0.0E+0) ) B( N, N ) = EPS3
!
      TRANS = 'N'
!
   ELSE
!
!        UL decomposition with partial pivoting of B, replacing zero
!        pivots by EPS3.
!
      DO J = N, 2, -1
         EJ = H( J, J-1 )
         IF( CABS1( B( J, J ) ) < CABS1( EJ ) ) THEN
!
!              Interchange columns and eliminate.
!
            X = CLADIV( B( J, J ), EJ )
            B( J, J ) = EJ
            B_TMP(1:J-1) = B(1:J-1, J-1 )
            B(1:J-1, J-1 ) = B(1:J-1, J ) - X*B_TMP(1:J-1)
            B(1:J-1, J ) = B_TMP(1:J-1)
         ELSE
!
!              Eliminate without interchange.
!
            IF( B( J, J ) == (0.0E+0,0.0E+0) ) B( J, J ) = EPS3
            X = CLADIV( EJ, B( J, J ) )
            IF( X /= (0.0E+0,0.0E+0) ) B(1:J-1,J-1) = B(1:J-1,J-1) - X*B(1:J-1,J)
         END IF
      ENDDO
      IF( B( 1, 1 ) == (0.0E+0,0.0E+0) ) B( 1, 1 ) = EPS3
!
      TRANS = 'C'
!
   END IF
!
   NORMIN = 'N'
   DO ITS = 1, N
!
!        Solve U*x = scale*v for a right eigenvector
!          or U**H *x = scale*v for a left eigenvector,
!        overwriting x on v.
!
      CALL CLATRS( 'Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, V, &
                   SCALE, RWORK, IERR )
      NORMIN = 'Y'
!
!        Test for sufficient growth in the norm of v.
!
      VNORM = sum(ABS(REAL(V(1:N))) + ABS(AIMAG(V(1:N))))
      IF( VNORM >= GROWTO*SCALE ) GO TO 120
!
!        Choose new orthogonal starting vector and try again.
!
      RTEMP = EPS3 / ( ROOTN+1.0E+0 )
      V( 1 ) = EPS3
      V(2:N) = RTEMP
      V( N-ITS+1 ) = V( N-ITS+1 ) - EPS3*ROOTN
   ENDDO
!
!     Failure to find eigenvector in N iterations.
!
   INFO = 1
!
  120 CONTINUE
!
!     Normalize eigenvector.
!
   I = ICAMAX( N, V, 1 )
   V(1:N) = V(1:N) / CABS1(V(I))
!
   RETURN
!
!     End of CLAEIN
!
END
