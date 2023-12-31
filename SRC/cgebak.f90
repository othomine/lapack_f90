!> \brief \b CGEBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEBAK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebak.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebak.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebak.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       REAL               SCALE( * )
!       COMPLEX            V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEBAK forms the right or left eigenvectors of a complex general
!> matrix by backward transformation on the computed eigenvectors of the
!> balanced matrix output by CGEBAL.
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
!>          JOB must be the same as the argument JOB supplied to CGEBAL.
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
!>          The integers ILO and IHI determined by CGEBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is REAL array, dimension (N)
!>          Details of the permutation and scaling factors, as returned
!>          by CGEBAL.
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
!>          V is COMPLEX array, dimension (LDV,M)
!>          On entry, the matrix of right or left eigenvectors to be
!>          transformed, as returned by CHSEIN or CTREVC.
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
   SUBROUTINE CGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
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
   REAL               SCALE( * )
   COMPLEX            V( LDV, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   COMPLEX            V_TMP( M )
   LOGICAL            LEFTV, RIGHTV
   INTEGER            I, II, K
   REAL               S
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA
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
      CALL XERBLA( 'CGEBAK', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
   IF( M == 0 ) RETURN
   IF( LSAME( JOB, 'N' ) ) RETURN
!
   IF( ILO /= IHI ) THEN
!
!     Backward balance
!
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
!
         IF( RIGHTV ) THEN
            DO I = ILO, IHI
               V(I,1:M) = V(I,1:M) * SCALE( I )
            ENDDO
         END IF
!
         IF( LEFTV ) THEN
            DO I = ILO, IHI
               V(I,1:M) = V(I,1:M) / SCALE( I )
            ENDDO
         END IF
!
      END IF
   END IF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
   IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
      IF( RIGHTV ) THEN
         DO II = 1, N
            I = II
            IF( I < ILO .OR. I > IHI ) THEN
               IF( I < ILO ) I = ILO - II
               K = INT( SCALE( I ) )
               IF( K /= I ) THEN
                  V_TMP(1:M) = V(I,1:M)
                  V(I,1:M) = V(K,1:M)
                  V(K,1:M) = V_TMP(1:M)
               ENDIF
            ENDIF
         ENDDO
      END IF
!
      IF( LEFTV ) THEN
         DO II = 1, N
            I = II
            IF( I < ILO .OR. I > IHI ) THEN
               IF( I < ILO ) I = ILO - II
               K = INT( SCALE( I ) )
               IF( K /= I ) THEN
                  V_TMP(1:M) = V(I,1:M)
                  V(I,1:M) = V(K,1:M)
                  V(K,1:M) = V_TMP(1:M)
               ENDIF
            ENDIF
         ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of CGEBAK
!
END

