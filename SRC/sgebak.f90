!> \brief \b SGEBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEBAK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgebak.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgebak.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgebak.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       REAL               V( LDV, * ), SCALE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEBAK forms the right or left eigenvectors of a real general matrix
!> by backward transformation on the computed eigenvectors of the
!> balanced matrix output by SGEBAL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the type of backward transformation required:
!>          = 'N': do nothing, return immediately;
!>          = 'P': do backward transformation for permutation only;
!>          = 'S': do backward transformation for scaling only;
!>          = 'B': do backward transformations for both permutation and
!>                 scaling.
!>          JOB must be the same as the argument JOB supplied to SGEBAL.
!> \endverbatim
!>
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  V contains right eigenvectors;
!>          = 'L':  V contains left eigenvectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrix V.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          The integers ILO and IHI determined by SGEBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is REAL array, dimension (N)
!>          Details of the permutation and scaling factors, as returned
!>          by SGEBAL.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns of the matrix V.  M >= 0.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is REAL array, dimension (LDV,M)
!>          On entry, the matrix of right or left eigenvectors to be
!>          transformed, as returned by SHSEIN or STREVC.
!>          On exit, V is overwritten by the transformed eigenvectors.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup gebak
!
!  =====================================================================
   SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
                      INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          JOB, SIDE
   INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
   REAL               V( LDV, * ), SCALE( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ONE
   PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LEFTV, RIGHTV
   INTEGER            I, II, K
   REAL               S
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           SSCAL, SSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
   RIGHTV = LSAME( SIDE, 'R' )
   LEFTV = LSAME( SIDE, 'L' )
!
   INFO = 0
   IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
       .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
      INFO = -2
   ELSE IF( N < 0 ) THEN
      INFO = -3
   ELSE IF( ILO < 1 .OR. ILO > MAX( 1, N ) ) THEN
      INFO = -4
   ELSE IF( IHI < MIN( ILO, N ) .OR. IHI > N ) THEN
      INFO = -5
   ELSE IF( M < 0 ) THEN
      INFO = -7
   ELSE IF( LDV < MAX( 1, N ) ) THEN
      INFO = -9
   END IF
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'SGEBAK', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) &
      RETURN
   IF( M == 0 ) &
      RETURN
   IF( LSAME( JOB, 'N' ) ) &
      RETURN
!
   IF( ILO == IHI ) &
      GO TO 30
!
!     Backward balance
!
   IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
!
      IF( RIGHTV ) THEN
         DO I = ILO, IHI
            S = SCALE( I )
            CALL SSCAL( M, S, V( I, 1 ), LDV )
         ENDDO
      END IF
!
      IF( LEFTV ) THEN
         DO I = ILO, IHI
            S = ONE / SCALE( I )
            CALL SSCAL( M, S, V( I, 1 ), LDV )
         ENDDO
      END IF
!
   END IF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
30 CONTINUE
   IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
      IF( RIGHTV ) THEN
         DO II = 1, N
            I = II
            IF( I >= ILO .AND. I <= IHI ) &
               GO TO 40
            IF( I < ILO ) &
               I = ILO - II
            K = INT( SCALE( I ) )
            IF( K == I ) &
               GO TO 40
            CALL SSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
40       CONTINUE
         ENDDO
      END IF
!
      IF( LEFTV ) THEN
         DO II = 1, N
            I = II
            IF( I >= ILO .AND. I <= IHI ) &
               GO TO 50
            IF( I < ILO ) &
               I = ILO - II
            K = INT( SCALE( I ) )
            IF( K == I ) &
               GO TO 50
            CALL SSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
50       CONTINUE
         ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of SGEBAK
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

