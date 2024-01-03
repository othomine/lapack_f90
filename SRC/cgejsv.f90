!> \brief \b CGEJSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEJSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgejsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgejsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgejsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!     SUBROUTINE CGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP,
!                         M, N, A, LDA, SVA, U, LDU, V, LDV,
!                         CWORK, LWORK, RWORK, LRWORK, IWORK, INFO )
!
!     .. Scalar Arguments ..
!     IMPLICIT    N1.0E+0
!     INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N
!     ..
!     .. Array Arguments ..
!     COMPLEX     A( LDA, * ),  U( LDU, * ), V( LDV, * ), CWORK( LWORK )
!     REAL        SVA( N ), RWORK( LRWORK )
!     INTEGER     IWORK( * )
!     CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEJSV computes the singular value decomposition (SVD) of a complex M-by-N
!> matrix [A], where M >= N. The SVD of [A] is written as
!>
!>              [A] = [U] * [SIGMA] * [V]^*,
!>
!> where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N
!> diagonal elements, [U] is an M-by-N (or M-by-M) unitary matrix, and
!> [V] is an N-by-N unitary matrix. The diagonal elements of [SIGMA] are
!> the singular values of [A]. The columns of [U] and [V] are the left and
!> the right singular vectors of [A], respectively. The matrices [U] and [V]
!> are computed and stored in the arrays U and V, respectively. The diagonal
!> of [SIGMA] is computed and stored in the array SVA.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBA
!> \verbatim
!>          JOBA is CHARACTER*1
!>         Specifies the level of accuracy:
!>       = 'C': This option works well (high relative accuracy) if A = B * D,
!>              with well-conditioned B and arbitrary diagonal matrix D.
!>              The accuracy cannot be spoiled by COLUMN scaling. The
!>              accuracy of the computed output depends on the condition of
!>              B, and the procedure aims at the best theoretical accuracy.
!>              The relative error max_{i=1:N}|d sigma_i| / sigma_i is
!>              bounded by f(M,N)*epsilon* cond(B), independent of D.
!>              The input matrix is preprocessed with the QRF with column
!>              pivoting. This initial preprocessing and preconditioning by
!>              a rank revealing QR factorization is common for all values of
!>              JOBA. Additional actions are specified as follows:
!>       = 'E': Computation as with 'C' with an additional estimate of the
!>              condition number of B. It provides a realistic error bound.
!>       = 'F': If A = D1 * C * D2 with ill-conditioned diagonal scalings
!>              D1, D2, and well-conditioned matrix C, this option gives
!>              higher accuracy than the 'C' option. If the structure of the
!>              input matrix is not known, and relative accuracy is
!>              desirable, then this option is advisable. The input matrix A
!>              is preprocessed with QR factorization with FULL (row and
!>              column) pivoting.
!>       = 'G': Computation as with 'F' with an additional estimate of the
!>              condition number of B, where A=B*D. If A has heavily weighted
!>              rows, then using this condition number gives too pessimistic
!>              error bound.
!>       = 'A': Small singular values are not well determined by the data
!>              and are considered as noisy; the matrix is treated as
!>              numerically rank deficient. The error in the computed
!>              singular values is bounded by f(m,n)*epsilon*||A||.
!>              The computed SVD A = U * S * V^* restores A up to
!>              f(m,n)*epsilon*||A||.
!>              This gives the procedure the licence to discard (set to zero)
!>              all singular values below N*epsilon*||A||.
!>       = 'R': Similar as in 'A'. Rank revealing property of the initial
!>              QR factorization is used do reveal (using triangular factor)
!>              a gap sigma_{r+1} < epsilon * sigma_r in which case the
!>              numerical RANK is declared to be r. The SVD is computed with
!>              absolute error bounds, but more accurately than with 'A'.
!> \endverbatim
!>
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>         Specifies whether to compute the columns of U:
!>       = 'U': N columns of U are returned in the array U.
!>       = 'F': full set of M left sing. vectors is returned in the array U.
!>       = 'W': U may be used as workspace of length M*N. See the description
!>              of U.
!>       = 'N': U is not computed.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>         Specifies whether to compute the matrix V:
!>       = 'V': N columns of V are returned in the array V; Jacobi rotations
!>              are not explicitly accumulated.
!>       = 'J': N columns of V are returned in the array V, but they are
!>              computed as the product of Jacobi rotations, if JOBT = 'N'.
!>       = 'W': V may be used as workspace of length N*N. See the description
!>              of V.
!>       = 'N': V is not computed.
!> \endverbatim
!>
!> \param[in] JOBR
!> \verbatim
!>          JOBR is CHARACTER*1
!>         Specifies the RANGE for the singular values. Issues the licence to
!>         set to zero small positive singular values if they are outside
!>         specified range. If A  /=  0 is scaled so that the largest singular
!>         value of c*A is around SQRT(BIG), BIG=SLAMCH('O'), then JOBR issues
!>         the licence to kill columns of A whose norm in c*A is less than
!>         SQRT(SFMIN) (for JOBR = 'R'), or less than SMALL=SFMIN/EPSLN,
!>         where SFMIN=SLAMCH('S'), EPSLN=SLAMCH('E').
!>       = 'N': Do not kill small columns of c*A. This option assumes that
!>              BLAS and QR factorizations and triangular solvers are
!>              implemented to work in that range. If the condition of A
!>              is greater than BIG, use CGESVJ.
!>       = 'R': RESTRICTED range for sigma(c*A) is [SQRT(SFMIN), SQRT(BIG)]
!>              (roughly, as described above). This option is recommended.
!>                                             ===========================
!>         For computing the singular values in the FULL range [SFMIN,BIG]
!>         use CGESVJ.
!> \endverbatim
!>
!> \param[in] JOBT
!> \verbatim
!>          JOBT is CHARACTER*1
!>         If the matrix is square then the procedure may determine to use
!>         transposed A if A^* seems to be better with respect to convergence.
!>         If the matrix is not square, JOBT is ignored.
!>         The decision is based on two values of entropy over the adjoint
!>         orbit of A^* * A. See the descriptions of RWORK(6) and RWORK(7).
!>       = 'T': transpose if entropy test indicates possibly faster
!>         convergence of Jacobi process if A^* is taken as input. If A is
!>         replaced with A^*, then the row pivoting is included automatically.
!>       = 'N': do not speculate.
!>         The option 'T' can be used to compute only the singular values, or
!>         the full SVD (U, SIGMA and V). For only one set of singular vectors
!>         (U or V), the caller should provide both U and V, as one of the
!>         matrices is used as workspace if the matrix A is transposed.
!>         The implementer can easily remove this constraint and make the
!>         code more complicated. See the descriptions of U and V.
!>         In general, this option is considered experimental, and 'N'; should
!>         be preferred. This is subject to changes in the future.
!> \endverbatim
!>
!> \param[in] JOBP
!> \verbatim
!>          JOBP is CHARACTER*1
!>         Issues the licence to introduce structured perturbations to drown
!>         denormalized numbers. This licence should be active if the
!>         denormals are poorly implemented, causing slow computation,
!>         especially in cases of fast convergence (!). For details see [1,2].
!>         For the sake of simplicity, this perturbations are included only
!>         when the full SVD or only the singular values are requested. The
!>         implementer/user can easily add the perturbation for the cases of
!>         computing one set of singular vectors.
!>       = 'P': introduce perturbation
!>       = 'N': do not perturb
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>         The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The number of columns of the input matrix A. M >= N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] SVA
!> \verbatim
!>          SVA is REAL array, dimension (N)
!>          On exit,
!>          - For RWORK(1)/RWORK(2) = 1.0E+0: The singular values of A. During
!>            the computation SVA contains Euclidean column norms of the
!>            iterated matrices in the array A.
!>          - For RWORK(1)  /=  RWORK(2): The singular values of A are
!>            (RWORK(1)/RWORK(2)) * SVA(1:N). This factored form is used if
!>            sigma_max(A) overflows or if small singular values have been
!>            saved from underflow by scaling the input matrix A.
!>          - If JOBR='R' then some of the singular values may be returned
!>            as exact zeros obtained by "set to zero" because they are
!>            below the numerical rank threshold or are denormalized numbers.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array, dimension ( LDU, N ) or ( LDU, M )
!>          If JOBU = 'U', then U contains on exit the M-by-N matrix of
!>                         the left singular vectors.
!>          If JOBU = 'F', then U contains on exit the M-by-M matrix of
!>                         the left singular vectors, including an ONB
!>                         of the orthogonal complement of the Range(A).
!>          If JOBU = 'W'  .AND. (JOBV = 'V' .AND. JOBT = 'T' .AND. M = N),
!>                         then U is used as workspace if the procedure
!>                         replaces A with A^*. In that case, [V] is computed
!>                         in U as left singular vectors of A^* and then
!>                         copied back to the V array. This 'W' option is just
!>                         a reminder to the caller that in this case U is
!>                         reserved as workspace of length N*N.
!>          If JOBU = 'N'  U is not referenced, unless JOBT='T'.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U,  LDU >= 1.
!>          IF  JOBU = 'U' or 'F' or 'W',  then LDU >= M.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension ( LDV, N )
!>          If JOBV = 'V', 'J' then V contains on exit the N-by-N matrix of
!>                         the right singular vectors;
!>          If JOBV = 'W', AND (JOBU = 'U' AND JOBT = 'T' AND M = N),
!>                         then V is used as workspace if the procedure
!>                         replaces A with A^*. In that case, [U] is computed
!>                         in V as right singular vectors of A^* and then
!>                         copied back to the U array. This 'W' option is just
!>                         a reminder to the caller that in this case V is
!>                         reserved as workspace of length N*N.
!>          If JOBV = 'N'  V is not referenced, unless JOBT='T'.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V,  LDV >= 1.
!>          If JOBV = 'V' or 'J' or 'W', then LDV >= N.
!> \endverbatim
!>
!> \param[out] CWORK
!> \verbatim
!>          CWORK is COMPLEX array, dimension (MAX(2,LWORK))
!>          If the call to CGEJSV is a workspace query (indicated by LWORK=-1 or
!>          LRWORK=-1), then on exit CWORK(1) contains the required length of
!>          CWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          Length of CWORK to confirm proper allocation of workspace.
!>          LWORK depends on the job:
!>
!>          1. If only SIGMA is needed ( JOBU = 'N', JOBV = 'N' ) and
!>            1.1 .. no scaled condition estimate required (JOBA /= 'E'.AND.JOBA /= 'G'):
!>               LWORK >= 2*N+1. This is the minimal requirement.
!>               ->> For optimal performance (blocked code) the optimal value
!>               is LWORK >= N + (N+1)*NB. Here NB is the optimal
!>               block size for CGEQP3 and CGEQRF.
!>               In general, optimal LWORK is computed as
!>               LWORK >= max(N+LWORK(CGEQP3),N+LWORK(CGEQRF), LWORK(CGESVJ)).
!>            1.2. .. an estimate of the scaled condition number of A is
!>               required (JOBA='E', or 'G'). In this case, LWORK the minimal
!>               requirement is LWORK >= N*N + 2*N.
!>               ->> For optimal performance (blocked code) the optimal value
!>               is LWORK >= max(N+(N+1)*NB, N*N+2*N)=N**2+2*N.
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(CGEQP3),N+LWORK(CGEQRF), LWORK(CGESVJ),
!>                            N*N+LWORK(CPOCON)).
!>          2. If SIGMA and the right singular vectors are needed (JOBV = 'V'),
!>             (JOBU = 'N')
!>            2.1   .. no scaled condition estimate requested (JOBE = 'N'):
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance,
!>               LWORK >= max(N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for CGEQP3, CGEQRF, CGELQF,
!>               CUNMLQ. In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(CGEQP3), N+LWORK(CGESVJ),
!>                       N+LWORK(CGELQF), 2*N+LWORK(CGEQRF), N+LWORK(CUNMLQ)).
!>            2.2 .. an estimate of the scaled condition number of A is
!>               required (JOBA='E', or 'G').
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance,
!>               LWORK >= max(N+(N+1)*NB, 2*N,2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for CGEQP3, CGEQRF, CGELQF,
!>               CUNMLQ. In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(CGEQP3), LWORK(CPOCON), N+LWORK(CGESVJ),
!>                       N+LWORK(CGELQF), 2*N+LWORK(CGEQRF), N+LWORK(CUNMLQ)).
!>          3. If SIGMA and the left singular vectors are needed
!>            3.1  .. no scaled condition estimate requested (JOBE = 'N'):
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance:
!>               if JOBU = 'U' :: LWORK >= max(3*N, N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for CGEQP3, CGEQRF, CUNMQR.
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(CGEQP3), 2*N+LWORK(CGEQRF), N+LWORK(CUNMQR)).
!>            3.2  .. an estimate of the scaled condition number of A is
!>               required (JOBA='E', or 'G').
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance:
!>               if JOBU = 'U' :: LWORK >= max(3*N, N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for CGEQP3, CGEQRF, CUNMQR.
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(CGEQP3),N+LWORK(CPOCON),
!>                        2*N+LWORK(CGEQRF), N+LWORK(CUNMQR)).
!>
!>          4. If the full SVD is needed: (JOBU = 'U' or JOBU = 'F') and
!>            4.1. if JOBV = 'V'
!>               the minimal requirement is LWORK >= 5*N+2*N*N.
!>            4.2. if JOBV = 'J' the minimal requirement is
!>               LWORK >= 4*N+N*N.
!>            In both cases, the allocated CWORK can accommodate blocked runs
!>            of CGEQP3, CGEQRF, CGELQF, CUNMQR, CUNMLQ.
!>
!>          If the call to CGEJSV is a workspace query (indicated by LWORK=-1 or
!>          LRWORK=-1), then on exit CWORK(1) contains the optimal and CWORK(2) contains the
!>          minimal length of CWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(7,LRWORK))
!>          On exit,
!>          RWORK(1) = Determines the scaling factor SCALE = RWORK(2) / RWORK(1)
!>                    such that SCALE*SVA(1:N) are the computed singular values
!>                    of A. (See the description of SVA().)
!>          RWORK(2) = See the description of RWORK(1).
!>          RWORK(3) = SCONDA is an estimate for the condition number of
!>                    column equilibrated A. (If JOBA = 'E' or 'G')
!>                    SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1).
!>                    It is computed using CPOCON. It holds
!>                    N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
!>                    where R is the triangular factor from the QRF of A.
!>                    However, if R is truncated and the numerical rank is
!>                    determined to be strictly smaller than N, SCONDA is
!>                    returned as -1, thus indicating that the smallest
!>                    singular values might be lost.
!>
!>          If full SVD is needed, the following two condition numbers are
!>          useful for the analysis of the algorithm. They are provided for
!>          a developer/implementer who is familiar with the details of
!>          the method.
!>
!>          RWORK(4) = an estimate of the scaled condition number of the
!>                    triangular factor in the first QR factorization.
!>          RWORK(5) = an estimate of the scaled condition number of the
!>                    triangular factor in the second QR factorization.
!>          The following two parameters are computed if JOBT = 'T'.
!>          They are provided for a developer/implementer who is familiar
!>          with the details of the method.
!>          RWORK(6) = the entropy of A^* * A :: this is the Shannon entropy
!>                    of diag(A^* * A) / Trace(A^* * A) taken as point in the
!>                    probability simplex.
!>          RWORK(7) = the entropy of A * A^*. (See the description of RWORK(6).)
!>          If the call to CGEJSV is a workspace query (indicated by LWORK=-1 or
!>          LRWORK=-1), then on exit RWORK(1) contains the required length of
!>          RWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          Length of RWORK to confirm proper allocation of workspace.
!>          LRWORK depends on the job:
!>
!>       1. If only the singular values are requested i.e. if
!>          LSAME(JOBU,'N') .AND. LSAME(JOBV,'N')
!>          then:
!>          1.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>               then: LRWORK = max( 7, 2 * M ).
!>          1.2. Otherwise, LRWORK  = max( 7,  N ).
!>       2. If singular values with the right singular vectors are requested
!>          i.e. if
!>          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J')) .AND.
!>          .NOT.(LSAME(JOBU,'U').OR.LSAME(JOBU,'F'))
!>          then:
!>          2.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>          then LRWORK = max( 7, 2 * M ).
!>          2.2. Otherwise, LRWORK  = max( 7,  N ).
!>       3. If singular values with the left singular vectors are requested, i.e. if
!>          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND.
!>          .NOT.(LSAME(JOBV,'V').OR.LSAME(JOBV,'J'))
!>          then:
!>          3.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>          then LRWORK = max( 7, 2 * M ).
!>          3.2. Otherwise, LRWORK  = max( 7,  N ).
!>       4. If singular values with both the left and the right singular vectors
!>          are requested, i.e. if
!>          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND.
!>          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J'))
!>          then:
!>          4.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>          then LRWORK = max( 7, 2 * M ).
!>          4.2. Otherwise, LRWORK  = max( 7, N ).
!>
!>          If, on entry, LRWORK = -1 or LWORK=-1, a workspace query is assumed and
!>          the length of RWORK is returned in RWORK(1).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, of dimension at least 4, that further depends
!>          on the job:
!>
!>          1. If only the singular values are requested then:
!>             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>             then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>          2. If the singular values and the right singular vectors are requested then:
!>             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>             then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>          3. If the singular values and the left singular vectors are requested then:
!>             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>             then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>          4. If the singular values with both the left and the right singular vectors
!>             are requested, then:
!>             4.1. If LSAME(JOBV,'J') the length of IWORK is determined as follows:
!>                  If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>                  then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>             4.2. If LSAME(JOBV,'V') the length of IWORK is determined as follows:
!>                  If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>                  then the length of IWORK is 2*N+M; otherwise the length of IWORK is 2*N.
!>
!>          On exit,
!>          IWORK(1) = the numerical rank determined after the initial
!>                     QR factorization with pivoting. See the descriptions
!>                     of JOBA and JOBR.
!>          IWORK(2) = the number of the computed nonzero singular values
!>          IWORK(3) = if nonzero, a warning message:
!>                     If IWORK(3) = 1 then some of the column norms of A
!>                     were denormalized floats. The requested high accuracy
!>                     is not warranted by the data.
!>          IWORK(4) = 1 or -1. If IWORK(4) = 1, then the procedure used A^* to
!>                     do the job as specified by the JOB parameters.
!>          If the call to CGEJSV is a workspace query (indicated by LWORK = -1 and
!>          LRWORK = -1), then on exit IWORK(1) contains the required length of
!>          IWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           < 0:  if INFO = -i, then the i-th argument had an illegal value.
!>           = 0:  successful exit;
!>           > 0:  CGEJSV  did not converge in the maximal allowed number
!>                 of sweeps. The computed values may be inaccurate.
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
!> \ingroup gejsv
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>  CGEJSV implements a preconditioned Jacobi SVD algorithm. It uses CGEQP3,
!>  CGEQRF, and CGELQF as preprocessors and preconditioners. Optionally, an
!>  additional row pivoting can be used as a preprocessor, which in some
!>  cases results in much higher accuracy. An example is matrix A with the
!>  structure A = D1 * C * D2, where D1, D2 are arbitrarily ill-conditioned
!>  diagonal matrices and C is well-conditioned matrix. In that case, complete
!>  pivoting in the first QR factorizations provides accuracy dependent on the
!>  condition number of C, and independent of D1, D2. Such higher accuracy is
!>  not completely understood theoretically, but it works well in practice.
!>  Further, if A can be written as A = B*D, with well-conditioned B and some
!>  diagonal D, then the high accuracy is guaranteed, both theoretically and
!>  in software, independent of D. For more details see [1], [2].
!>     The computational range for the singular values can be the full range
!>  ( UNDERFLOW,OVERFLOW ), provided that the machine arithmetic and the BLAS
!>  & LAPACK routines called by CGEJSV are implemented to work in that range.
!>  If that is not the case, then the restriction for safe computation with
!>  the singular values in the range of normalized IEEE numbers is that the
!>  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not
!>  overflow. This code (CGEJSV) is best used in this restricted range,
!>  meaning that singular values of magnitude below ||A||_2 / SLAMCH('O') are
!>  returned as zeros. See JOBR for details on this.
!>     Further, this implementation is somewhat slower than the one described
!>  in [1,2] due to replacement of some non-LAPACK components, and because
!>  the choice of some tuning parameters in the iterative part (CGESVJ) is
!>  left to the implementer on a particular machine.
!>     The rank revealing QR factorization (in this code: CGEQP3) should be
!>  implemented as in [3]. We have a new version of CGEQP3 under development
!>  that is more robust than the current one in LAPACK, with a cleaner cut in
!>  rank deficient cases. It will be available in the SIGMA library [4].
!>  If M is much larger than N, it is obvious that the initial QRF with
!>  column pivoting can be preprocessed by the QRF without pivoting. That
!>  well known trick is not used in CGEJSV because in some cases heavy row
!>  weighting can be treated with complete pivoting. The overhead in cases
!>  M much larger than N is then only due to pivoting, but the benefits in
!>  terms of accuracy have prevailed. The implementer/user can incorporate
!>  this extra QRF step easily. The implementer can also improve data movement
!>  (matrix transpose, matrix copy, matrix transposed copy) - this
!>  implementation of CGEJSV uses only the simplest, naive data movement.
!> \endverbatim
!
!> \par Contributor:
!  ==================
!>
!>  Zlatko Drmac (Zagreb, Croatia)
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!> [1] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I.
!>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342.
!>     LAPACK Working note 169.
!> [2] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II.
!>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362.
!>     LAPACK Working note 170.
!> [3] Z. Drmac and Z. Bujanovic: On the failure of rank-revealing QR
!>     factorization software - a case study.
!>     ACM Trans. Math. Softw. Vol. 35, No 2 (2008), pp. 1-28.
!>     LAPACK Working note 176.
!> [4] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,
!>     QSVD, (H,K)-SVD computations.
!>     Department of Mathematics, University of Zagreb, 2008, 2016.
!> \endverbatim
!
!>  \par Bugs, examples and comments:
!   =================================
!>
!>  Please report all bugs and send interesting examples and/or comments to
!>  drmac@math.hr. Thank you.
!>
!  =====================================================================
   SUBROUTINE CGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, &
                      M, N, A, LDA, SVA, U, LDU, V, LDV, &
                      CWORK, LWORK, RWORK, LRWORK, IWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   IMPLICIT    NONE
   INTEGER     INFO, LDA, LDU, LDV, LWORK, LRWORK, M, N
!     ..
!     .. Array Arguments ..
   COMPLEX     A( LDA, * ), U( LDU, * ), V( LDV, * ), CWORK( LWORK )
   REAL        SVA( N ), RWORK( LRWORK )
   INTEGER     IWORK( * )
   CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV
!     ..
!
!  ===========================================================================
!     ..
!     .. Local Scalars ..
   COMPLEX CTEMP, U_TMP(N,N)
   REAL    AAPP,   AAQQ,   AATMAX, AATMIN, BIG,    BIG1,   COND_OK, &
           CONDR1, CONDR2, ENTRA,  ENTRAT, EPSLN,  MAXPRJ, SCALEM, &
           SCONDA, SFMIN,  SMALL,  TEMP1,  USCAL1, USCAL2, XSC
   INTEGER IERR,   N1,     NR,     NUMRANK,        p, q,   WARNING, I
   LOGICAL ALMORT, DEFR,   ERREST, GOSCAL,  JRACC,  KILL,   LQUERY, &
           LSVEC,  L2ABER, L2KILL, L2PERT,  L2RANK, L2TRAN, NOSCAL, &
           ROWPIV, RSVEC,  TRANSP
!
   INTEGER OPTWRK, MINWRK, MINRWRK, MINIWRK
   INTEGER LWCON,  LWLQF, LWQP3, LWQRF, LWUNMLQ, LWUNMQR, LWUNMQRM, &
           LWSVDJ, LWSVDJV, LRWQP3, LRWCON, LRWSVDJ, IWOFF
   INTEGER LWRK_CGELQF, LWRK_CGEQP3,  LWRK_CGEQP3N, LWRK_CGEQRF, &
           LWRK_CGESVJ, LWRK_CGESVJV, LWRK_CGESVJU, LWRK_CUNMLQ, &
           LWRK_CUNMQR, LWRK_CUNMQRM
!     ..
!     .. Local Arrays
   COMPLEX CDUMMY(1)
   REAL    RDUMMY(1)
!     ..
!     .. External Functions ..
   REAL      SLAMCH, SCNRM2
   INTEGER   ISAMAX, ICAMAX
   LOGICAL   LSAME
   EXTERNAL  ISAMAX, ICAMAX, LSAME, SLAMCH, SCNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL  SLASSQ, CCOPY,  CGELQF, CGEQP3, CGEQRF, CLACPY, CLAPMR, &
             CLASCL, SLASCL, CLASET, CLASSQ, CLASWP, CUNGQR, CUNMLQ, &
             CUNMQR, CPOCON, SSCAL,  CSSCAL, CTRSM,  XERBLA
!
   EXTERNAL  CGESVJ
!     ..
!
!     Test the input arguments
!
   LSVEC  = LSAME( JOBU, 'U' ) .OR. LSAME( JOBU, 'F' )
   JRACC  = LSAME( JOBV, 'J' )
   RSVEC  = LSAME( JOBV, 'V' ) .OR. JRACC
   ROWPIV = LSAME( JOBA, 'F' ) .OR. LSAME( JOBA, 'G' )
   L2RANK = LSAME( JOBA, 'R' )
   L2ABER = LSAME( JOBA, 'A' )
   ERREST = LSAME( JOBA, 'E' ) .OR. LSAME( JOBA, 'G' )
   L2TRAN = LSAME( JOBT, 'T' ) .AND. ( M  ==  N )
   L2KILL = LSAME( JOBR, 'R' )
   DEFR   = LSAME( JOBR, 'N' )
   L2PERT = LSAME( JOBP, 'P' )
!
   LQUERY = ( LWORK  ==  -1 ) .OR. ( LRWORK  ==  -1 )
!
   IF ( .NOT.(ROWPIV .OR. L2RANK .OR. L2ABER .OR. &
        ERREST .OR. LSAME( JOBA, 'C' ) )) THEN
      INFO = - 1
   ELSE IF ( .NOT.( LSVEC .OR. LSAME( JOBU, 'N' ) .OR. &
      ( LSAME( JOBU, 'W' ) .AND. RSVEC .AND. L2TRAN ) ) ) THEN
      INFO = - 2
   ELSE IF ( .NOT.( RSVEC .OR. LSAME( JOBV, 'N' ) .OR. &
      ( LSAME( JOBV, 'W' ) .AND. LSVEC .AND. L2TRAN ) ) ) THEN
      INFO = - 3
   ELSE IF ( .NOT. ( L2KILL .OR. DEFR ) )    THEN
      INFO = - 4
   ELSE IF ( .NOT. ( LSAME(JOBT,'T') .OR. LSAME(JOBT,'N') ) ) THEN
      INFO = - 5
   ELSE IF ( .NOT. ( L2PERT .OR. LSAME( JOBP, 'N' ) ) ) THEN
      INFO = - 6
   ELSE IF ( M  <  0 ) THEN
      INFO = - 7
   ELSE IF ( ( N  <  0 ) .OR. ( N  >  M ) ) THEN
      INFO = - 8
   ELSE IF ( LDA  <  M ) THEN
      INFO = - 10
   ELSE IF ( LSVEC .AND. ( LDU  <  M ) ) THEN
      INFO = - 13
   ELSE IF ( RSVEC .AND. ( LDV  <  N ) ) THEN
      INFO = - 15
   ELSE
!        #:)
      INFO = 0
   END IF
!
   IF ( INFO  ==  0 ) THEN
!         .. compute the minimal and the optimal workspace lengths
!         [[The expressions for computing the minimal and the optimal
!         values of LCWORK, LRWORK are written with a lot of redundancy and
!         can be simplified. However, this verbose form is useful for
!         maintenance and modifications of the code.]]
!
!        .. minimal workspace length for CGEQP3 of an M x N matrix,
!         CGEQRF of an N x N matrix, CGELQF of an N x N matrix,
!         CUNMLQ for computing N x N matrix, CUNMQR for computing N x N
!         matrix, CUNMQR for computing M x N matrix, respectively.
       LWQP3 = N+1
       LWQRF = MAX( 1, N )
       LWLQF = MAX( 1, N )
       LWUNMLQ  = MAX( 1, N )
       LWUNMQR  = MAX( 1, N )
       LWUNMQRM = MAX( 1, M )
!        .. minimal workspace length for CPOCON of an N x N matrix
       LWCON = 2 * N
!        .. minimal workspace length for CGESVJ of an N x N matrix,
!         without and with explicit accumulation of Jacobi rotations
       LWSVDJ  = MAX( 2 * N, 1 )
       LWSVDJV = MAX( 2 * N, 1 )
!         .. minimal REAL workspace length for CGEQP3, CPOCON, CGESVJ
       LRWQP3  = 2 * N
       LRWCON  = N
       LRWSVDJ = N
       IF ( LQUERY ) THEN
           CALL CGEQP3( M, N, A, LDA, IWORK, CDUMMY, CDUMMY, -1, &
                RDUMMY, IERR )
           LWRK_CGEQP3 = INT( CDUMMY(1) )
           CALL CGEQRF( N, N, A, LDA, CDUMMY, CDUMMY,-1, IERR )
           LWRK_CGEQRF = INT( CDUMMY(1) )
           CALL CGELQF( N, N, A, LDA, CDUMMY, CDUMMY,-1, IERR )
           LWRK_CGELQF = INT( CDUMMY(1) )
       END IF
       MINWRK  = 2
       OPTWRK  = 2
       MINIWRK = N
       IF ( .NOT. (LSVEC .OR. RSVEC ) ) THEN
!             .. minimal and optimal sizes of the complex workspace if
!             only the singular values are requested
           IF ( ERREST ) THEN
               MINWRK = MAX( N+LWQP3, N**2+LWCON, N+LWQRF, LWSVDJ )
           ELSE
               MINWRK = MAX( N+LWQP3, N+LWQRF, LWSVDJ )
           END IF
           IF ( LQUERY ) THEN
               CALL CGESVJ( 'L', 'N', 'N', N, N, A, LDA, SVA, N, V, &
                    LDV, CDUMMY, -1, RDUMMY, -1, IERR )
               LWRK_CGESVJ = INT( CDUMMY(1) )
               IF ( ERREST ) THEN
                   OPTWRK = MAX( N+LWRK_CGEQP3, N**2+LWCON, &
                                 N+LWRK_CGEQRF, LWRK_CGESVJ )
               ELSE
                   OPTWRK = MAX( N+LWRK_CGEQP3, N+LWRK_CGEQRF, &
                                 LWRK_CGESVJ )
               END IF
           END IF
           IF ( L2TRAN .OR. ROWPIV ) THEN
               IF ( ERREST ) THEN
                  MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWCON, LRWSVDJ )
               ELSE
                  MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ )
               END IF
           ELSE
               IF ( ERREST ) THEN
                  MINRWRK = MAX( 7, LRWQP3, LRWCON, LRWSVDJ )
               ELSE
                  MINRWRK = MAX( 7, LRWQP3, LRWSVDJ )
               END IF
           END IF
           IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
       ELSE IF ( RSVEC .AND. (.NOT.LSVEC) ) THEN
!            .. minimal and optimal sizes of the complex workspace if the
!            singular values and the right singular vectors are requested
          IF ( ERREST ) THEN
              MINWRK = MAX( N+LWQP3, LWCON, LWSVDJ, N+LWLQF, &
                            2*N+LWQRF, N+LWSVDJ, N+LWUNMLQ )
          ELSE
              MINWRK = MAX( N+LWQP3, LWSVDJ, N+LWLQF, 2*N+LWQRF, &
                            N+LWSVDJ, N+LWUNMLQ )
          END IF
          IF ( LQUERY ) THEN
              CALL CGESVJ( 'L', 'U', 'N', N,N, U, LDU, SVA, N, A, &
                   LDA, CDUMMY, -1, RDUMMY, -1, IERR )
              LWRK_CGESVJ = INT( CDUMMY(1) )
              CALL CUNMLQ( 'L', 'C', N, N, N, A, LDA, CDUMMY, &
                   V, LDV, CDUMMY, -1, IERR )
              LWRK_CUNMLQ = INT( CDUMMY(1) )
              IF ( ERREST ) THEN
              OPTWRK = MAX( N+LWRK_CGEQP3, LWCON, LWRK_CGESVJ, &
                            N+LWRK_CGELQF, 2*N+LWRK_CGEQRF, &
                            N+LWRK_CGESVJ,  N+LWRK_CUNMLQ )
              ELSE
              OPTWRK = MAX( N+LWRK_CGEQP3, LWRK_CGESVJ,N+LWRK_CGELQF, &
                            2*N+LWRK_CGEQRF, N+LWRK_CGESVJ, &
                            N+LWRK_CUNMLQ )
              END IF
          END IF
          IF ( L2TRAN .OR. ROWPIV ) THEN
               IF ( ERREST ) THEN
                  MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON )
               ELSE
                  MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ )
               END IF
          ELSE
               IF ( ERREST ) THEN
                  MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON )
               ELSE
                  MINRWRK = MAX( 7, LRWQP3, LRWSVDJ )
               END IF
          END IF
          IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
       ELSE IF ( LSVEC .AND. (.NOT.RSVEC) ) THEN
!            .. minimal and optimal sizes of the complex workspace if the
!            singular values and the left singular vectors are requested
          IF ( ERREST ) THEN
              MINWRK = N + MAX( LWQP3,LWCON,N+LWQRF,LWSVDJ,LWUNMQRM )
          ELSE
              MINWRK = N + MAX( LWQP3, N+LWQRF, LWSVDJ, LWUNMQRM )
          END IF
          IF ( LQUERY ) THEN
              CALL CGESVJ( 'L', 'U', 'N', N,N, U, LDU, SVA, N, A, &
                   LDA, CDUMMY, -1, RDUMMY, -1, IERR )
              LWRK_CGESVJ = INT( CDUMMY(1) )
              CALL CUNMQR( 'L', 'N', M, N, N, A, LDA, CDUMMY, U, &
                  LDU, CDUMMY, -1, IERR )
              LWRK_CUNMQRM = INT( CDUMMY(1) )
              IF ( ERREST ) THEN
              OPTWRK = N + MAX( LWRK_CGEQP3, LWCON, N+LWRK_CGEQRF, &
                                LWRK_CGESVJ, LWRK_CUNMQRM )
              ELSE
              OPTWRK = N + MAX( LWRK_CGEQP3, N+LWRK_CGEQRF, &
                                LWRK_CGESVJ, LWRK_CUNMQRM )
              END IF
          END IF
          IF ( L2TRAN .OR. ROWPIV ) THEN
              IF ( ERREST ) THEN
                 MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON )
              ELSE
                 MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ )
              END IF
          ELSE
              IF ( ERREST ) THEN
                 MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON )
              ELSE
                 MINRWRK = MAX( 7, LRWQP3, LRWSVDJ )
              END IF
          END IF
          IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
       ELSE
!            .. minimal and optimal sizes of the complex workspace if the
!            full SVD is requested
          IF ( .NOT. JRACC ) THEN
              IF ( ERREST ) THEN
                 MINWRK = MAX( N+LWQP3, N+LWCON,  2*N+N**2+LWCON, &
                            2*N+LWQRF,         2*N+LWQP3, &
                            2*N+N**2+N+LWLQF,  2*N+N**2+N+N**2+LWCON, &
                            2*N+N**2+N+LWSVDJ, 2*N+N**2+N+LWSVDJV, &
                            2*N+N**2+N+LWUNMQR,2*N+N**2+N+LWUNMLQ, &
                            N+N**2+LWSVDJ,   N+LWUNMQRM )
              ELSE
                 MINWRK = MAX( N+LWQP3, 3*N+N**2+MAX(&
                            LWLQF, N**2+LWCON, LWSVDJ, LWSVDJV, LWUNMQR,&
                            LWUNMLQ), 2*N+MAX(LWQRF, LWQP3, N**2+LWCON), &
                            N+N**2+LWSVDJ,      N+LWUNMQRM )
              END IF
              MINIWRK = MINIWRK + N
              IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
          ELSE
              IF ( ERREST ) THEN
                 MINWRK = MAX( N+LWQP3, N+LWCON, 2*N+LWQRF, &
                            2*N+N**2+LWSVDJV, 2*N+N**2+N+LWUNMQR, &
                            N+LWUNMQRM )
              ELSE
                 MINWRK = MAX( N+LWQP3, 2*N+LWQRF, &
                            2*N+N**2+LWSVDJV, 2*N+N**2+N+LWUNMQR, &
                            N+LWUNMQRM )
              END IF
              IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
          END IF
          IF ( LQUERY ) THEN
              CALL CUNMQR( 'L', 'N', M, N, N, A, LDA, CDUMMY, U, &
                   LDU, CDUMMY, -1, IERR )
              LWRK_CUNMQRM = INT( CDUMMY(1) )
              CALL CUNMQR( 'L', 'N', N, N, N, A, LDA, CDUMMY, U, &
                   LDU, CDUMMY, -1, IERR )
              LWRK_CUNMQR = INT( CDUMMY(1) )
              IF ( .NOT. JRACC ) THEN
                  CALL CGEQP3( N,N, A, LDA, IWORK, CDUMMY,CDUMMY, -1, &
                       RDUMMY, IERR )
                  LWRK_CGEQP3N = INT( CDUMMY(1) )
                  CALL CGESVJ( 'L', 'U', 'N', N, N, U, LDU, SVA, &
                       N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                  LWRK_CGESVJ = INT( CDUMMY(1) )
                  CALL CGESVJ( 'U', 'U', 'N', N, N, U, LDU, SVA, &
                       N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                  LWRK_CGESVJU = INT( CDUMMY(1) )
                  CALL CGESVJ( 'L', 'U', 'V', N, N, U, LDU, SVA, &
                       N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                  LWRK_CGESVJV = INT( CDUMMY(1) )
                  CALL CUNMLQ( 'L', 'C', N, N, N, A, LDA, CDUMMY, &
                       V, LDV, CDUMMY, -1, IERR )
                  LWRK_CUNMLQ = INT( CDUMMY(1) )
                  IF ( ERREST ) THEN
                    OPTWRK = MAX( N+LWRK_CGEQP3, N+LWCON, &
                             2*N+N**2+LWCON, 2*N+LWRK_CGEQRF, &
                             2*N+LWRK_CGEQP3N, &
                             2*N+N**2+N+LWRK_CGELQF, &
                             2*N+N**2+N+N**2+LWCON, &
                             2*N+N**2+N+LWRK_CGESVJ, &
                             2*N+N**2+N+LWRK_CGESVJV, &
                             2*N+N**2+N+LWRK_CUNMQR, &
                             2*N+N**2+N+LWRK_CUNMLQ, &
                             N+N**2+LWRK_CGESVJU, &
                             N+LWRK_CUNMQRM )
                  ELSE
                    OPTWRK = MAX( N+LWRK_CGEQP3, &
                             2*N+N**2+LWCON, 2*N+LWRK_CGEQRF, &
                             2*N+LWRK_CGEQP3N, &
                             2*N+N**2+N+LWRK_CGELQF, &
                             2*N+N**2+N+N**2+LWCON, &
                             2*N+N**2+N+LWRK_CGESVJ, &
                             2*N+N**2+N+LWRK_CGESVJV, &
                             2*N+N**2+N+LWRK_CUNMQR, &
                             2*N+N**2+N+LWRK_CUNMLQ, &
                             N+N**2+LWRK_CGESVJU, &
                             N+LWRK_CUNMQRM )
                  END IF
              ELSE
                  CALL CGESVJ( 'L', 'U', 'V', N, N, U, LDU, SVA, &
                       N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                  LWRK_CGESVJV = INT( CDUMMY(1) )
                  CALL CUNMQR( 'L', 'N', N, N, N, CDUMMY, N, CDUMMY, &
                       V, LDV, CDUMMY, -1, IERR )
                  LWRK_CUNMQR = INT( CDUMMY(1) )
                  CALL CUNMQR( 'L', 'N', M, N, N, A, LDA, CDUMMY, U, &
                       LDU, CDUMMY, -1, IERR )
                  LWRK_CUNMQRM = INT( CDUMMY(1) )
                  IF ( ERREST ) THEN
                     OPTWRK = MAX( N+LWRK_CGEQP3, N+LWCON, &
                              2*N+LWRK_CGEQRF, 2*N+N**2, &
                              2*N+N**2+LWRK_CGESVJV, &
                              2*N+N**2+N+LWRK_CUNMQR,N+LWRK_CUNMQRM )
                  ELSE
                     OPTWRK = MAX( N+LWRK_CGEQP3, 2*N+LWRK_CGEQRF, &
                              2*N+N**2, 2*N+N**2+LWRK_CGESVJV, &
                              2*N+N**2+N+LWRK_CUNMQR, &
                              N+LWRK_CUNMQRM )
                  END IF
              END IF
          END IF
          IF ( L2TRAN .OR. ROWPIV ) THEN
              MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON )
          ELSE
              MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON )
          END IF
       END IF
       MINWRK = MAX( 2, MINWRK )
       OPTWRK = MAX( OPTWRK, MINWRK )
       IF ( LWORK   <  MINWRK  .AND. (.NOT.LQUERY) ) INFO = - 17
       IF ( LRWORK  <  MINRWRK .AND. (.NOT.LQUERY) ) INFO = - 19
   END IF
!
   IF ( INFO  /=  0 ) THEN
!       #:(
      CALL XERBLA( 'CGEJSV', - INFO )
      RETURN
   ELSE IF ( LQUERY ) THEN
       CWORK(1) = OPTWRK
       CWORK(2) = MINWRK
       RWORK(1) = MINRWRK
       IWORK(1) = MAX( 4, MINIWRK )
       RETURN
   END IF
!
!     Quick return for void matrix (Y3K safe)
! #:)
   IF ( ( M  ==  0 ) .OR. ( N  ==  0 ) ) THEN
      IWORK(1:4) = 0
      RWORK(1:7) = 0
      RETURN
   ENDIF
!
!     Determine whether the matrix U should be M x N or M x M
!
   IF ( LSVEC ) THEN
      N1 = N
      IF ( LSAME( JOBU, 'F' ) ) N1 = M
   END IF
!
!     Set numerical parameters
!
!!    NOTE: Make sure SLAMCH() does not fail on the target architecture.
!
   EPSLN = SLAMCH('Epsilon')
   SFMIN = SLAMCH('SafeMinimum')
   SMALL = SFMIN / EPSLN
   BIG   = SLAMCH('O')
!     BIG   = 1.0E+0 / SFMIN
!
!     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
!
!(!)  If necessary, scale SVA() to protect the largest norm from
!     overflow. It is possible that this scaling pushes the smallest
!     column norm left from the underflow threshold (extreme case).
!
   SCALEM  = 1.0E+0 / SQRT(REAL(M)*REAL(N))
   NOSCAL  = .TRUE.
   GOSCAL  = .TRUE.
   DO p = 1, N
      AAPP = 0.0E+0
      AAQQ = 1.0E+0
      CALL CLASSQ( M, A(1,p), 1, AAPP, AAQQ )
      IF ( AAPP  >  BIG ) THEN
         INFO = - 9
         CALL XERBLA( 'CGEJSV', -INFO )
         RETURN
      END IF
      AAQQ = SQRT(AAQQ)
      IF ( ( AAPP  <  (BIG / AAQQ) ) .AND. NOSCAL  ) THEN
         SVA(p)  = AAPP * AAQQ
      ELSE
         NOSCAL  = .FALSE.
         SVA(p)  = AAPP * ( AAQQ * SCALEM )
         IF ( GOSCAL ) THEN
            GOSCAL = .FALSE.
            SVA(1:P-1) = SCALEM * SVA(1:P-1)
         END IF
      END IF
      ENDDO
!
   IF ( NOSCAL ) SCALEM = 1.0E+0
!
   AAPP = 0.0E+0
   AAQQ = BIG
   DO p = 1, N
      AAPP = MAX( AAPP, SVA(p) )
      IF ( SVA(p)  /=  0.0E+0 ) AAQQ = MIN( AAQQ, SVA(p) )
   ENDDO
!
!     Quick return for zero M x N matrix
! #:)
   IF ( AAPP  ==  0.0E+0 ) THEN
      IF ( LSVEC ) CALL CLASET( 'G', M, N1, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0), U, LDU )
      IF ( RSVEC ) THEN
         V(1:N,1:N) = (0.0E+0,0.0E+0)
         DO I = 1, N
           V(I,I) = (1.0E+0,0.0E+0)
         ENDDO
      ENDIF
      RWORK(1) = 1.0E+0
      RWORK(2) = 1.0E+0
      IF ( ERREST ) RWORK(3) = 1.0E+0
      IF ( LSVEC .AND. RSVEC ) THEN
         RWORK(4) = 1.0E+0
         RWORK(5) = 1.0E+0
      END IF
      IF ( L2TRAN ) THEN
         RWORK(6) = 0.0E+0
         RWORK(7) = 0.0E+0
      END IF
      IWORK(1) = 0
      IWORK(2) = 0
      IWORK(3) = 0
      IWORK(4) = -1
      RETURN
   END IF
!
!     Issue warning if denormalized column norms detected. Override the
!     high relative accuracy request. Issue licence to kill nonzero columns
!     (set them to zero) whose norm is less than sigma_max / BIG (roughly).
! #:(
   WARNING = 0
   IF ( AAQQ  <=  SFMIN ) THEN
      L2RANK = .TRUE.
      L2KILL = .TRUE.
      WARNING = 1
   END IF
!
!     Quick return for one-column matrix
! #:)
   IF ( N  ==  1 ) THEN
!
      IF ( LSVEC ) THEN
         CALL CLASCL( 'G',0,0,SVA(1),SCALEM, M,1,A(1,1),LDA,IERR )
         U(1:M,1) = A(1:M,1)
!           computing all M left singular vectors of the M x 1 matrix
         IF ( N1  /=  N  ) THEN
           CALL CGEQRF( M, N, U,LDU, CWORK, CWORK(N+1),LWORK-N,IERR )
           CALL CUNGQR( M,N1,1, U,LDU,CWORK,CWORK(N+1),LWORK-N,IERR )
           A(1:M,1) = U(1:M,1)
         END IF
      END IF
      IF ( RSVEC ) THEN
          V(1,1) = (1.0E+0,0.0E+0)
      END IF
      IF ( SVA(1)  <  (BIG*SCALEM) ) THEN
         SVA(1)  = SVA(1) / SCALEM
         SCALEM  = 1.0E+0
      END IF
      RWORK(1) = 1.0E+0 / SCALEM
      RWORK(2) = 1.0E+0
      IF ( SVA(1)  /=  0.0E+0 ) THEN
         IWORK(1) = 1
         IF ( ( SVA(1) / SCALEM)  >=  SFMIN ) THEN
            IWORK(2) = 1
         ELSE
            IWORK(2) = 0
         END IF
      ELSE
         IWORK(1) = 0
         IWORK(2) = 0
      END IF
      IWORK(3) = 0
      IWORK(4) = -1
      IF ( ERREST ) RWORK(3) = 1.0E+0
      IF ( LSVEC .AND. RSVEC ) THEN
         RWORK(4) = 1.0E+0
         RWORK(5) = 1.0E+0
      END IF
      IF ( L2TRAN ) THEN
         RWORK(6) = 0.0E+0
         RWORK(7) = 0.0E+0
      END IF
      RETURN
!
   END IF
!
   TRANSP = .FALSE.
!
   AATMAX = -1.0E+0
   AATMIN =  BIG
   IF ( ROWPIV .OR. L2TRAN ) THEN
!
!     Compute the row norms, needed to determine row pivoting sequence
!     (in the case of heavily row weighted A, row pivoting is strongly
!     advised) and to collect information needed to compare the
!     structures of A * A^* and A^* * A (in the case L2TRAN == .TRUE.).
!
      IF ( L2TRAN ) THEN
         DO p = 1, M
            XSC   = 0.0E+0
            TEMP1 = 1.0E+0
            CALL CLASSQ( N, A(p,1), LDA, XSC, TEMP1 )
!              CLASSQ gets both the ell_2 and the ell_infinity norm
!              in one pass through the vector
            RWORK(M+p)  = XSC * SCALEM
            RWORK(p)    = XSC * (SCALEM*SQRT(TEMP1))
            AATMAX = MAX( AATMAX, RWORK(p) )
            IF (RWORK(p)  /=  0.0E+0) &
               AATMIN = MIN(AATMIN,RWORK(p))
            ENDDO
      ELSE
         DO p = 1, M
            RWORK(M+p) = SCALEM*ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
            AATMAX = MAX( AATMAX, RWORK(M+p) )
            AATMIN = MIN( AATMIN, RWORK(M+p) )
            ENDDO
      END IF
!
   END IF
!
!     For square matrix A try to determine whether A^*  would be better
!     input for the preconditioned Jacobi SVD, with faster convergence.
!     The decision is based on an O(N) function of the vector of column
!     and row norms of A, based on the Shannon entropy. This should give
!     the right choice in most cases when the difference actually matters.
!     It may fail and pick the slower converging side.
!
   ENTRA  = 0.0E+0
   ENTRAT = 0.0E+0
   IF ( L2TRAN ) THEN
!
      XSC   = 0.0E+0
      TEMP1 = 1.0E+0
      CALL SLASSQ( N, SVA, 1, XSC, TEMP1 )
      TEMP1 = 1.0E+0 / TEMP1
!
      ENTRA = 0.0E+0
      DO p = 1, N
         BIG1  = ( ( SVA(p) / XSC )**2 ) * TEMP1
         IF ( BIG1  /=  0.0E+0 ) ENTRA = ENTRA + BIG1 * ALOG(BIG1)
      ENDDO
      ENTRA = - ENTRA / ALOG(REAL(N))
!
!        Now, SVA().^2/Trace(A^* * A) is a point in the probability simplex.
!        It is derived from the diagonal of  A^* * A.  Do the same with the
!        diagonal of A * A^*, compute the entropy of the corresponding
!        probability distribution. Note that A * A^* and A^* * A have the
!        same trace.
!
      ENTRAT = 0.0E+0
      DO p = 1, M
         BIG1 = ( ( RWORK(p) / XSC )**2 ) * TEMP1
         IF ( BIG1  /=  0.0E+0 ) ENTRAT = ENTRAT + BIG1 * ALOG(BIG1)
      ENDDO
      ENTRAT = - ENTRAT / ALOG(REAL(M))
!
!        Analyze the entropies and decide A or A^*. Smaller entropy
!        usually means better input for the algorithm.
!
      TRANSP = ( ENTRAT  <  ENTRA )
!
!        If A^* is better than A, take the adjoint of A. This is allowed
!        only for square matrices, M=N.
      IF ( TRANSP ) THEN
!           In an optimal implementation, this trivial transpose
!           should be replaced with faster transpose.
         DO p = 1, N - 1
            A(p,p) = CONJG(A(p,p))
            DO q = p + 1, N
                CTEMP = CONJG(A(q,p))
               A(q,p) = CONJG(A(p,q))
               A(p,q) = CTEMP
            ENDDO
         ENDDO
         A(N,N) = CONJG(A(N,N))
         DO p = 1, N
            RWORK(M+p) = SVA(p)
            SVA(p) = RWORK(p)
!              previously computed row 2-norms are now column 2-norms
!              of the transposed matrix
         ENDDO
         TEMP1  = AAPP
         AAPP   = AATMAX
         AATMAX = TEMP1
         TEMP1  = AAQQ
         AAQQ   = AATMIN
         AATMIN = TEMP1
         KILL   = LSVEC
         LSVEC  = RSVEC
         RSVEC  = KILL
         IF ( LSVEC ) N1 = N
!
         ROWPIV = .TRUE.
      END IF
!
   END IF
!     END IF L2TRAN
!
!     Scale the matrix so that its maximal singular value remains less
!     than SQRT(BIG) -- the matrix is scaled so that its maximal column
!     has Euclidean norm equal to SQRT(BIG/N). The only reason to keep
!     SQRT(BIG) instead of BIG is the fact that CGEJSV uses LAPACK and
!     BLAS routines that, in some implementations, are not capable of
!     working in the full interval [SFMIN,BIG] and that they may provoke
!     overflows in the intermediate results. If the singular values spread
!     from SFMIN to BIG, then CGESVJ will compute them. So, in that case,
!     one should use CGESVJ instead of CGEJSV.
   BIG1   = SQRT( BIG )
   TEMP1  = SQRT( BIG / REAL(N) )
!     >> for future updates: allow bigger range, i.e. the largest column
!     will be allowed up to BIG/N and CGESVJ will do the rest. However, for
!     this all other (LAPACK) components must allow such a range.
!     TEMP1  = BIG/REAL(N)
!     TEMP1  = BIG * EPSLN  this should 'almost' work with current LAPACK components
   CALL SLASCL( 'G', 0, 0, AAPP, TEMP1, N, 1, SVA, N, IERR )
   IF ( AAQQ  >  (AAPP * SFMIN) ) THEN
       AAQQ = ( AAQQ / AAPP ) * TEMP1
   ELSE
       AAQQ = ( AAQQ * TEMP1 ) / AAPP
   END IF
   TEMP1 = TEMP1 * SCALEM
   CALL CLASCL( 'G', 0, 0, AAPP, TEMP1, M, N, A, LDA, IERR )
!
!     To undo scaling at the end of this procedure, multiply the
!     computed singular values with USCAL2 / USCAL1.
!
   USCAL1 = TEMP1
   USCAL2 = AAPP
!
   IF ( L2KILL ) THEN
!        L2KILL enforces computation of nonzero singular values in
!        the restricted range of condition number of the initial A,
!        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN).
      XSC = SQRT( SFMIN )
   ELSE
      XSC = SMALL
!
!        Now, if the condition number of A is too big,
!        sigma_max(A) / sigma_min(A)  >  SQRT(BIG/N) * EPSLN / SFMIN,
!        as a precaution measure, the full SVD is computed using CGESVJ
!        with accumulated Jacobi rotations. This provides numerically
!        more robust computation, at the cost of slightly increased run
!        time. Depending on the concrete implementation of BLAS and LAPACK
!        (i.e. how they behave in presence of extreme ill-conditioning) the
!        implementor may decide to remove this switch.
      IF ( ( AAQQ < SQRT(SFMIN) ) .AND. LSVEC .AND. RSVEC ) THEN
         JRACC = .TRUE.
      END IF
!
   END IF
   IF ( AAQQ  <  XSC ) THEN
      DO p = 1, N
         IF ( SVA(p)  <  XSC ) THEN
            A(1:M,p) = (0.0E+0,0.0E+0)
            SVA(p) = 0.0E+0
         END IF
      ENDDO
   END IF
!
!     Preconditioning using QR factorization with pivoting
!
   IF ( ROWPIV ) THEN
!        Optional row permutation (Bjoerck row pivoting):
!        A result by Cox and Higham shows that the Bjoerck's
!        row pivoting combined with standard column pivoting
!        has similar effect as Powell-Reid complete pivoting.
!        The ell-infinity norms of A are made nonincreasing.
      IF ( ( LSVEC .AND. RSVEC ) .AND. .NOT.( JRACC ) ) THEN
           IWOFF = 2*N
      ELSE
           IWOFF = N
      END IF
      DO p = 1, M - 1
         q = ISAMAX( M-p+1, RWORK(M+p), 1 ) + p - 1
         IWORK(IWOFF+p) = q
         IF ( p  /=  q ) THEN
            TEMP1      = RWORK(M+p)
            RWORK(M+p) = RWORK(M+q)
            RWORK(M+q) = TEMP1
         END IF
      ENDDO
      CALL CLASWP( N, A, LDA, 1, M-1, IWORK(IWOFF+1), 1 )
   END IF
!
!     End of the preparation phase (scaling, optional sorting and
!     transposing, optional flushing of small columns).
!
!     Preconditioning
!
!     If the full SVD is needed, the right singular vectors are computed
!     from a matrix equation, and for that we need theoretical analysis
!     of the Businger-Golub pivoting. So we use CGEQP3 as the first RR QRF.
!     In all other cases the first RR QRF can be chosen by other criteria
!     (eg speed by replacing global with restricted window pivoting, such
!     as in xGEQPX from TOMS # 782). Good results will be obtained using
!     xGEQPX with properly (!) chosen numerical parameters.
!     Any improvement of CGEQP3 improves overall performance of CGEJSV.
!
!     A * P1 = Q1 * [ R1^* 0]^*:
!        .. all columns are free columns
   IWORK(1:N) = 0
   CALL CGEQP3( M, N, A, LDA, IWORK, CWORK, CWORK(N+1), LWORK-N, &
                RWORK, IERR )
!
!     The upper triangular matrix R1 from the first QRF is inspected for
!     rank deficiency and possibilities for deflation, or possible
!     ill-conditioning. Depending on the user specified flag L2RANK,
!     the procedure explores possibilities to reduce the numerical
!     rank by inspecting the computed upper triangular factor. If
!     L2RANK or L2ABER are up, then CGEJSV will compute the SVD of
!     A + dA, where ||dA|| <= f(M,N)*EPSLN.
!
   NR = 1
   IF ( L2ABER ) THEN
!        Standard absolute error bound suffices. All sigma_i with
!        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
!        aggressive enforcement of lower numerical rank by introducing a
!        backward error of the order of N*EPSLN*||A||.
      TEMP1 = SQRT(REAL(N))*EPSLN
      DO p = 2, N
         IF ( ABS(A(p,p))  >=  (TEMP1*ABS(A(1,1))) ) THEN
            NR = NR + 1
         ELSE
            EXIT
         END IF
      ENDDO
   ELSE IF ( L2RANK ) THEN
!        .. similarly as above, only slightly more gentle (less aggressive).
!        Sudden drop on the diagonal of R1 is used as the criterion for
!        close-to-rank-deficient.
      TEMP1 = SQRT(SFMIN)
      DO p = 2, N
         IF ( ( ABS(A(p,p))  <  (EPSLN*ABS(A(p-1,p-1))) ) .OR. &
              ( ABS(A(p,p))  <  SMALL ) .OR. &
              ( L2KILL .AND. (ABS(A(p,p))  <  TEMP1) ) ) EXIT
         NR = NR + 1
      ENDDO
!
   ELSE
!        The goal is high relative accuracy. However, if the matrix
!        has high scaled condition number the relative accuracy is in
!        general not feasible. Later on, a condition number estimator
!        will be deployed to estimate the scaled condition number.
!        Here we just remove the underflowed part of the triangular
!        factor. This prevents the situation in which the code is
!        working hard to get the accuracy not warranted by the data.
      TEMP1  = SQRT(SFMIN)
      DO p = 2, N
         IF ( ( ABS(A(p,p))  <  SMALL ) .OR. &
              ( L2KILL .AND. (ABS(A(p,p))  <  TEMP1) ) ) EXIT
         NR = NR + 1
      ENDDO
!
   END IF
!
   ALMORT = .FALSE.
   IF ( NR  ==  N ) THEN
      MAXPRJ = 1.0E+0
      DO p = 2, N
         TEMP1  = ABS(A(p,p)) / SVA(IWORK(p))
         MAXPRJ = MIN( MAXPRJ, TEMP1 )
      ENDDO
      IF ( MAXPRJ**2  >=  1.0E+0 - REAL(N)*EPSLN ) ALMORT = .TRUE.
   END IF
!
!
   SCONDA = - 1.0E+0
   CONDR1 = - 1.0E+0
   CONDR2 = - 1.0E+0
!
   IF ( ERREST ) THEN
      IF ( N  ==  NR ) THEN
         IF ( RSVEC ) THEN
!              .. V is available as workspace
            CALL CLACPY( 'U', N, N, A, LDA, V, LDV )
            DO p = 1, N
               V(1:p,p) = V(1:p,p)/SVA(IWORK(p))
            ENDDO
            IF ( LSVEC )THEN
                CALL CPOCON( 'U', N, V, LDV, 1.0E+0, TEMP1, &
                     CWORK(N+1), RWORK, IERR )
            ELSE
                CALL CPOCON( 'U', N, V, LDV, 1.0E+0, TEMP1, &
                     CWORK, RWORK, IERR )
            END IF
!
         ELSE IF ( LSVEC ) THEN
!              .. U is available as workspace
            CALL CLACPY( 'U', N, N, A, LDA, U, LDU )
            DO p = 1, N
               U(1:p,p) = U(1:p,p)/SVA(IWORK(p))
            ENDDO
            CALL CPOCON( 'U', N, U, LDU, 1.0E+0, TEMP1, &
                 CWORK(N+1), RWORK, IERR )
         ELSE
            CALL CLACPY( 'U', N, N, A, LDA, CWORK, N )
![]            CALL CLACPY( 'U', N, N, A, LDA, CWORK(N+1), N )
!              Change: here index shifted by N to the left, CWORK(1:N)
!              not needed for SIGMA only computation
            DO p = 1, N
!                TEMP1 = SVA(IWORK(p))
![]               CALL CSSCAL( p, 1.0E+0/TEMP1, CWORK(N+(p-1)*N+1), 1 )
               CWORK((p-1)*N+1:(p-1)*N+p) = CWORK((p-1)*N+1:(p-1)*N+p)/SVA(IWORK(p))
!                CALL CSSCAL( p, 1.0E+0/TEMP1, CWORK((p-1)*N+1), 1 )
            ENDDO
!           .. the columns of R are scaled to have unit Euclidean lengths.
![]               CALL CPOCON( 'U', N, CWORK(N+1), N, 1.0E+0, TEMP1,
![]     $              CWORK(N+N*N+1), RWORK, IERR )
            CALL CPOCON( 'U', N, CWORK, N, 1.0E+0, TEMP1, CWORK(N*N+1), RWORK, IERR )
!
         END IF
         IF ( TEMP1  /=  0.0E+0 ) THEN
            SCONDA = 1.0E+0 / SQRT(TEMP1)
         ELSE
            SCONDA = - 1.0E+0
         END IF
!           SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1).
!           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
      ELSE
         SCONDA = - 1.0E+0
      END IF
   END IF
!
   L2PERT = L2PERT .AND. ( ABS( A(1,1)/A(NR,NR) )  >  SQRT(BIG1) )
!     If there is no violent scaling, artificial perturbation is not needed.
!
!     Phase 3:
!
   IF ( .NOT. ( RSVEC .OR. LSVEC ) ) THEN
!
!         Singular Values only
!
!         .. transpose A(1:NR,1:N)
      DO p = 1, MIN( N-1, NR )
         A(p+1:N,p) = A(p,p+1:N)
         A(p:N,p) = CONJG(A(p:N,p))
      ENDDO
      IF ( NR  ==  N ) A(N,N) = CONJG(A(N,N))
!
!        The following two DO-loops introduce small relative perturbation
!        into the strict upper triangle of the lower triangular matrix.
!        Small entries below the main diagonal are also changed.
!        This modification is useful if the computing environment does not
!        provide/allow FLUSH TO 0.0E+0 underflow, for it prevents many
!        annoying denormalized numbers in case of strongly scaled matrices.
!        The perturbation is structured so that it does not introduce any
!        new perturbation of the singular values, and it does not destroy
!        the job done by the preconditioner.
!        The licence for this perturbation is in the variable L2PERT, which
!        should be .FALSE. if FLUSH TO 0.0E+0 underflow is active.
!
      IF ( .NOT. ALMORT ) THEN
!
         IF ( L2PERT ) THEN
!              XSC = SQRT(SMALL)
            XSC = EPSLN / REAL(N)
            DO q = 1, NR
               CTEMP = CMPLX(XSC*ABS(A(q,q)),0.0E+0)
               DO p = 1, N
                  IF ( ( (p > q) .AND. (ABS(A(p,q)) <= TEMP1) ) &
                       .OR. ( p  <  q ) ) A(p,q) = CTEMP
!     $                     A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) &
               ENDDO
            ENDDO
         ELSE
            CALL CLASET( 'U', NR-1,NR-1, (0.0E+0,0.0E+0),(0.0E+0,0.0E+0), A(1,2),LDA )
         END IF
!
!            .. second preconditioning using the QR factorization
!
         CALL CGEQRF( N,NR, A,LDA, CWORK, CWORK(N+1),LWORK-N, IERR )
!
!           .. and transpose upper to lower triangular
         DO p = 1, NR - 1
            A(p+1:NR,p) = A(p,p+1:NR)
            A(p:NR,p) = CONJG(A(p:NR,p))
         ENDDO
!
      END IF
!
!           Row-cyclic Jacobi SVD algorithm with column pivoting
!
!           .. again some perturbation (a "background noise") is added
!           to drown denormals
         IF ( L2PERT ) THEN
!              XSC = SQRT(SMALL)
            XSC = EPSLN / REAL(N)
            DO q = 1, NR
               CTEMP = CMPLX(XSC*ABS(A(q,q)),0.0E+0)
               DO p = 1, NR
                  IF ( ( (p > q) .AND. (ABS(A(p,q)) <= TEMP1) ) &
                          .OR. ( p  <  q ) ) A(p,q) = CTEMP
!     $                   A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) &
               ENDDO
            ENDDO
         ELSE
            CALL CLASET( 'U', NR-1, NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), A(1,2), LDA )
         END IF
!
!           .. and one-sided Jacobi rotations are started on a lower
!           triangular matrix (plus perturbation which is ignored in
!           the part which destroys triangular form (confusing?!))
!
         CALL CGESVJ( 'L', 'N', 'N', NR, NR, A, LDA, SVA, &
                   N, V, LDV, CWORK, LWORK, RWORK, LRWORK, INFO )
!
         SCALEM  = RWORK(1)
         NUMRANK = NINT(RWORK(2))
!
!
   ELSE IF ( ( RSVEC .AND. ( .NOT. LSVEC ) .AND. ( .NOT. JRACC ) ) .OR. &
      ( JRACC .AND. ( .NOT. LSVEC ) .AND. ( NR  /=  N ) ) ) THEN
!
!        -> Singular Values and Right Singular Vectors <-
!
      IF ( ALMORT ) THEN
!
!           .. in this case NR equals N
         DO p = 1, NR
            V(p,p:N) = A(p,p:N)
            V(p:N,p) = CONJG(V(p:N,p))
         ENDDO
         CALL CLASET( 'U', NR-1,NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), V(1,2), LDV )
!
         CALL CGESVJ( 'L','U','N', N, NR, V, LDV, SVA, NR, A, LDA, &
                     CWORK, LWORK, RWORK, LRWORK, INFO )
         SCALEM  = RWORK(1)
         NUMRANK = NINT(RWORK(2))

      ELSE
!
!        .. two more QR factorizations ( one QRF is not enough, two require
!        accumulated product of Jacobi rotations, three are perfect )
!
         CALL CLASET( 'L', NR-1,NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), A(2,1), LDA )
         CALL CGELQF( NR,N, A, LDA, CWORK, CWORK(N+1), LWORK-N, IERR)
         CALL CLACPY( 'L', NR, NR, A, LDA, V, LDV )
         CALL CLASET( 'U', NR-1,NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), V(1,2), LDV )
         CALL CGEQRF( NR, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), &
                      LWORK-2*N, IERR )
         DO p = 1, NR
            V(p:NR,p) = CONJG(V(p:NR,p))
         ENDDO
         CALL CLASET('U', NR-1, NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), V(1,2), LDV)
!
         CALL CGESVJ( 'L', 'U','N', NR, NR, V,LDV, SVA, NR, U, &
                     LDU, CWORK(N+1), LWORK-N, RWORK, LRWORK, INFO )
         SCALEM  = RWORK(1)
         NUMRANK = NINT(RWORK(2))
         IF ( NR  <  N ) THEN
            V(NR+1:N,1:NR) = (0.0E+0,0.0E+0)
            V(1:NR,NR+1:N) = (0.0E+0,0.0E+0)
            V(NR+1:N,NR+1:N) = (0.0E+0,0.0E+0)
            DO I = NR+1, N
               V( I, I ) = (1.0E+0,0.0E+0)
            ENDDO
         END IF
!
      CALL CUNMLQ( 'L', 'C', N, N, NR, A, LDA, CWORK, V, LDV, CWORK(N+1), LWORK-N, IERR )
!
      END IF
!         .. permute the rows of V
!         DO 8991 p = 1, N
!            CALL CCOPY( N, V(p,1), LDV, A(IWORK(p),1), LDA )
! 8991    CONTINUE
!         CALL CLACPY( 'All', N, N, A, LDA, V, LDV )
      CALL CLAPMR( .FALSE., N, N, V, LDV, IWORK )
!
      IF ( TRANSP ) U(1:N,1:N) = V(1:N,1:N)
!
   ELSE IF ( JRACC .AND. (.NOT. LSVEC) .AND. ( NR ==  N ) ) THEN
!
      CALL CLASET( 'L', N-1,N-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), A(2,1), LDA )
!
      CALL CGESVJ( 'U','N','V', N, N, A, LDA, SVA, N, V, LDV, &
                  CWORK, LWORK, RWORK, LRWORK, INFO )
       SCALEM  = RWORK(1)
       NUMRANK = NINT(RWORK(2))
       CALL CLAPMR( .FALSE., N, N, V, LDV, IWORK )
!
   ELSE IF ( LSVEC .AND. ( .NOT. RSVEC ) ) THEN
!
!        .. Singular Values and Left Singular Vectors                 ..
!
!        .. second preconditioning step to avoid need to accumulate
!        Jacobi rotations in the Jacobi iterations.
      DO p = 1, NR
         U(p,p:N) = A(p,p:N)
         U(p:N,p) = CONJG(U(p:N,p))
      ENDDO
      CALL CLASET( 'U', NR-1, NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), U(1,2), LDU )
!
      CALL CGEQRF( N, NR, U, LDU, CWORK(N+1), CWORK(2*N+1), &
                 LWORK-2*N, IERR )
!
      DO p = 1, NR - 1
         U(p+1:NR,p) = U(p,p+1:NR)
         U(p:N,p) = CONJG(U(p:N,p))
      ENDDO
      CALL CLASET( 'U', NR-1, NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), U(1,2), LDU )
!
      CALL CGESVJ( 'L', 'U', 'N', NR,NR, U, LDU, SVA, NR, A, &
           LDA, CWORK(N+1), LWORK-N, RWORK, LRWORK, INFO )
      SCALEM  = RWORK(1)
      NUMRANK = NINT(RWORK(2))
!
      IF ( NR  <  M ) THEN
         U(NR+1:M,1:NR) = (0.0E+0,0.0E+0)
         IF ( NR  <  N1 ) THEN
            U(1:NR,NR+1:N1) = (0.0E+0,0.0E+0)
            U(NR+1:M,NR+1:N1) = (0.0E+0,0.0E+0)
         END IF
      END IF
!
      CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR )
!
      IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
!
      DO p = 1, N1
         U(1:M,p) = U(1:M,p) / SCNRM2( M, U(1,p), 1 )
      ENDDO
!
      IF ( TRANSP ) V(1:N,1:N) = U(1:N,1:N)
!
   ELSE
!
!        .. Full SVD ..
!
      IF ( .NOT. JRACC ) THEN
!
      IF ( .NOT. ALMORT ) THEN
!
!           Second Preconditioning Step (QRF [with pivoting])
!           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
!           equivalent to an LQF CALL. Since in many libraries the QRF
!           seems to be better optimized than the LQF, we do explicit
!           transpose and use the QRF. This is subject to changes in an
!           optimized implementation of CGEJSV.
!
         DO p = 1, NR
            V(p:N,p) = A(p,p:N)
            V(p:N,p) = CONJG(V(p:N,p))
         ENDDO
!
!           .. the following two loops perturb small entries to avoid
!           denormals in the second QR factorization, where they are
!           as good as zeros. This is done to avoid painfully slow
!           computation with denormals. The relative size of the perturbation
!           is a parameter that can be changed by the implementer.
!           This perturbation device will be obsolete on machines with
!           properly implemented arithmetic.
!           To switch it off, set L2PERT=.FALSE. To remove it from  the
!           code, remove the action under L2PERT=.TRUE., leave the ELSE part.
!           The following two loops should be blocked and fused with the
!           transposed copy above.
!
         IF ( L2PERT ) THEN
            XSC = SQRT(SMALL)
            DO q = 1, NR
               CTEMP = CMPLX(XSC*ABS( V(q,q) ),0.0E+0)
               DO p = 1, N
                  IF ( ( p  >  q ) .AND. ( ABS(V(p,q))  <=  TEMP1 ) &
                      .OR. ( p  <  q ) ) V(p,q) = CTEMP
!     $                   V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) &
                  IF ( p  <  q ) V(p,q) = - V(p,q)
               ENDDO
            ENDDO
         ELSE
            CALL CLASET( 'U', NR-1, NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), V(1,2), LDV )
         END IF
!
!           Estimate the row scaled condition number of R1
!           (If R1 is rectangular, N > NR, then the condition number
!           of the leading NR x NR submatrix is estimated.)
!
         CALL CLACPY( 'L', NR, NR, V, LDV, CWORK(2*N+1), NR )
         DO p = 1, NR
            CWORK(2*N+(p-1)*NR+p:2*N+p*NR) = CWORK(2*N+(p-1)*NR+p:2*N+p*NR)/&
              SCNRM2(NR-p+1,CWORK(2*N+(p-1)*NR+p),1)
         ENDDO
         CALL CPOCON('L',NR,CWORK(2*N+1),NR,1.0E+0,TEMP1, &
                      CWORK(2*N+NR*NR+1),RWORK,IERR)
         CONDR1 = 1.0E+0 / SQRT(TEMP1)
!           .. here need a second opinion on the condition number
!           .. then assume worst case scenario
!           R1 is OK for inverse <=> CONDR1  <  REAL(N)
!           more conservative    <=> CONDR1  <  SQRT(REAL(N))
!
         COND_OK = SQRT(SQRT(REAL(NR)))
![TP]       COND_OK is a tuning parameter.
!
         IF ( CONDR1  <  COND_OK ) THEN
!              .. the second QRF without pivoting. Note: in an optimized
!              implementation, this QRF should be implemented as the QRF
!              of a lower triangular matrix.
!              R1^* = Q2 * R2
            CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), &
                 LWORK-2*N, IERR )
!
            IF ( L2PERT ) THEN
               XSC = SQRT(SMALL)/EPSLN
               DO p = 2, NR
                  DO q = 1, p - 1
                     CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), &
                                 0.0E+0)
                     IF ( ABS(V(q,p))  <=  TEMP1 ) V(q,p) = CTEMP
!     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) &
                  ENDDO
               ENDDO
            END IF
!
            IF ( NR  /=  N ) CALL CLACPY( 'A', N, NR, V, LDV, CWORK(2*N+1), N )
!              .. save ...
!
!           .. this transposed copy should be better than naive
            DO p = 1, NR - 1
               V(p+1:NR,p) = V(p,p+1:NR)
               V(p:NR,p) = CONJG(V(p:NR,p))
            ENDDO
            V(NR,NR)=CONJG(V(NR,NR))
!
            CONDR2 = CONDR1
!
         ELSE
!
!              .. ill-conditioned case: second QRF with pivoting
!              Note that windowed pivoting would be equally good
!              numerically, and more run-time efficient. So, in
!              an optimal implementation, the next call to CGEQP3
!              should be replaced with eg. CALL CGEQPX (ACM TOMS #782)
!              with properly (carefully) chosen parameters.
!
!              R1^* * P2 = Q2 * R2
            IWORK(N+1:N+NR) = 0
            CALL CGEQP3( N, NR, V, LDV, IWORK(N+1), CWORK(N+1), &
                     CWORK(2*N+1), LWORK-2*N, RWORK, IERR )
!*               CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1),
!*     $              LWORK-2*N, IERR )
            IF ( L2PERT ) THEN
               XSC = SQRT(SMALL)
               DO p = 2, NR
                  DO q = 1, p - 1
                     CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), 0.0E+0)
                     IF ( ABS(V(q,p))  <=  TEMP1 ) V(q,p) = CTEMP
!     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) &
                     ENDDO
                  ENDDO
            END IF
!
            CALL CLACPY( 'A', N, NR, V, LDV, CWORK(2*N+1), N )
!
            IF ( L2PERT ) THEN
               XSC = SQRT(SMALL)
               DO p = 2, NR
                  DO q = 1, p - 1
                     CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), 0.0E+0)
!                        V(p,q) = - TEMP1*( V(q,p) / ABS(V(q,p)) )
                     V(p,q) = - CTEMP
                  ENDDO
               ENDDO
            ELSE
               CALL CLASET( 'L',NR-1,NR-1,(0.0E+0,0.0E+0),(0.0E+0,0.0E+0),V(2,1),LDV )
            END IF
!              Now, compute R2 = L3 * Q3, the LQ factorization.
            CALL CGELQF( NR, NR, V, LDV, CWORK(2*N+N*NR+1), &
                  CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR )
!              .. and estimate the condition number
            CALL CLACPY( 'L',NR,NR,V,LDV,CWORK(2*N+N*NR+NR+1),NR )
            DO p = 1, NR
               TEMP1 = 1.0E+0/SCNRM2( p, CWORK(2*N+N*NR+NR+p), NR )
               CALL CSSCAL( p, TEMP1, CWORK(2*N+N*NR+NR+p), NR )
            ENDDO
            CALL CPOCON( 'L',NR,CWORK(2*N+N*NR+NR+1),NR,1.0E+0,TEMP1, &
                 CWORK(2*N+N*NR+NR+NR*NR+1),RWORK,IERR )
            CONDR2 = 1.0E+0 / SQRT(TEMP1)
!
!
            IF ( CONDR2  >=  COND_OK ) THEN
!                 .. save the Householder vectors used for Q3
!                 (this overwrites the copy of R2, as it will not be
!                 needed in this branch, but it does not overwrite the
!                 Huseholder vectors of Q2.).
               CALL CLACPY( 'U', NR, NR, V, LDV, CWORK(2*N+1), N )
!                 .. and the rest of the information on Q3 is in
!                 WORK(2*N+N*NR+1:2*N+N*NR+N)
            END IF
!
         END IF
!
         IF ( L2PERT ) THEN
            XSC = SQRT(SMALL)
            DO q = 2, NR
               CTEMP = XSC * V(q,q)
!                     V(p,q) = - TEMP1*( V(p,q) / ABS(V(p,q)) )
               V(1:q-1,q) = - CTEMP
            ENDDO
         ELSE
            CALL CLASET( 'U', NR-1,NR-1, (0.0E+0,0.0E+0),(0.0E+0,0.0E+0), V(1,2), LDV )
         END IF
!
!        Second preconditioning finished; continue with Jacobi SVD
!        The input matrix is lower triangular.
!
!        Recover the right singular vectors as solution of a well
!        conditioned triangular matrix equation.
!
         IF ( CONDR1  <  COND_OK ) THEN
!
            CALL CGESVJ( 'L','U','N',NR,NR,V,LDV,SVA,NR,U, LDU, &
                 CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,RWORK, &
                 LRWORK, INFO )
            SCALEM  = RWORK(1)
            NUMRANK = NINT(RWORK(2))
            DO p = 1, NR
               U(1:NR,p) = V(1:NR,p)
               V(1:NR,p) = SVA(p)*V(1:NR,p)
            ENDDO

!        .. pick the right matrix equation and solve it
!
            IF ( NR  ==  N ) THEN
! :))             .. best case, R1 is inverted. The solution of this matrix
!                 equation is Q2*V2 = the product of the Jacobi rotations
!                 used in CGESVJ, premultiplied with the orthogonal matrix
!                 from the second QR factorization.
               CALL CTRSM('L','U','N','N', NR,NR,(1.0E+0,0.0E+0), A,LDA, V,LDV)
            ELSE
!                 .. R1 is well conditioned, but non-square. Adjoint of R2
!                 is inverted to get the product of the Jacobi rotations
!                 used in CGESVJ. The Q-factor from the second QR
!                 factorization is then built in explicitly.
               CALL CTRSM('L','U','C','N',NR,NR,(1.0E+0,0.0E+0),CWORK(2*N+1), &
                    N,V,LDV)
               IF ( NR  <  N ) THEN
                  V(NR+1:N,1:NR) = (0.0E+0,0.0E+0)
                  V(1:NR,1:N) = (0.0E+0,0.0E+0)
                  V(NR+1:N,NR+1:N) = (0.0E+0,0.0E+0)
                  DO I = NR+1, N
                     V( I, I ) = (1.0E+0,0.0E+0)
                  ENDDO
               END IF
               CALL CUNMQR('L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), &
                   V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR)
            END IF
!
         ELSE IF ( CONDR2  <  COND_OK ) THEN
!
!              The matrix R2 is inverted. The solution of the matrix equation
!              is Q3^* * V3 = the product of the Jacobi rotations (applied to
!              the lower triangular L3 from the LQ factorization of
!              R2=L3*Q3), pre-multiplied with the transposed Q3.
            CALL CGESVJ( 'L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, &
                 LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, &
                 RWORK, LRWORK, INFO )
            SCALEM  = RWORK(1)
            NUMRANK = NINT(RWORK(2))
            DO p = 1, NR
               U(1:NR,p) = SVA(p)*V(1:NR,p)
            ENDDO
            CALL CTRSM('L','U','N','N',NR,NR,(1.0E+0,0.0E+0),CWORK(2*N+1),N, &
                       U,LDU)
!              .. apply the permutation from the second QR factorization
            DO q = 1, NR
               DO p = 1, NR
                  CWORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q)
               ENDDO
               DO p = 1, NR
                  U(p,q) = CWORK(2*N+N*NR+NR+p)
               ENDDO
            ENDDO
            IF ( NR  <  N ) THEN
               V(NR+1:N,1:NR) = (0.0E+0,0.0E+0)
               V(1:NR,NR+1:N) = (0.0E+0,0.0E+0)
               V(NR+1:N,NR+1:N) = (0.0E+0,0.0E+0)
               DO I = NR+1, N
                  V( I, I ) = (1.0E+0,0.0E+0)
               ENDDO
            END IF
            CALL CUNMQR( 'L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), &
                 V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
         ELSE
!              Last line of defense.
! #:(          This is a rather pathological case: no scaled condition
!              improvement after two pivoted QR factorizations. Other
!              possibility is that the rank revealing QR factorization
!              or the condition estimator has failed, or the COND_OK
!              is set very close to 1.0E+0 (which is unnecessary). Normally,
!              this branch should never be executed, but in rare cases of
!              failure of the RRQR or condition estimator, the last line of
!              defense ensures that CGEJSV completes the task.
!              Compute the full SVD of L3 using CGESVJ with explicit
!              accumulation of Jacobi rotations.
            CALL CGESVJ( 'L', 'U', 'V', NR, NR, V, LDV, SVA, NR, U, &
                 LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, &
                            RWORK, LRWORK, INFO )
            SCALEM  = RWORK(1)
            NUMRANK = NINT(RWORK(2))
            IF ( NR  <  N ) THEN
               V(NR+1:N,1:NR) = (0.0E+0,0.0E+0)
               V(1:NR,NR+1:N) = (0.0E+0,0.0E+0)
               V(NR+1:N,NR+1:N) = (0.0E+0,0.0E+0)
               DO I = NR+1, N
                  V( I, I ) = (1.0E+0,0.0E+0)
               ENDDO
            END IF
            CALL CUNMQR( 'L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), &
                 V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
!
            CALL CUNMLQ( 'L', 'C', NR, NR, NR, CWORK(2*N+1), N, &
                 CWORK(2*N+N*NR+1), U, LDU, CWORK(2*N+N*NR+NR+1), &
                 LWORK-2*N-N*NR-NR, IERR )
            DO q = 1, NR
               DO p = 1, NR
                  CWORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q)
               ENDDO
               U(1:NR,q) = CWORK((2+NR)*N+NR+1:(2+NR)*N+2*NR)
            ENDDO
!
         END IF
!
!           Permute the rows of V using the (column) permutation from the
!           first QRF. Also, scale the columns to make them unit in
!           Euclidean norm. This applies to all cases.
!
         TEMP1 = SQRT(REAL(N)) * EPSLN
         DO q = 1, N
            DO p = 1, N
               CWORK(2*N+N*NR+NR+IWORK(p)) = V(p,q)
            ENDDO
            V(1:N,q) = CWORK((2+NR)*N+NR+1:(2+NR)*N+NR+N)
            XSC = 1.0E+0 / SCNRM2( N, V(1,q), 1 )
            IF ( (XSC  <  (1.0E+0-TEMP1)) .OR. (XSC  >  (1.0E+0+TEMP1)) ) &
              V(1:N,q) = XSC*V(1:N,q)
            ENDDO
!           At this moment, V contains the right singular vectors of A.
!           Next, assemble the left singular vector matrix U (M x N).
         IF ( NR  <  M ) THEN
            U(NR+1:M,1:NR) = (0.0E+0,0.0E+0)
            IF ( NR  <  N1 ) THEN
               U(1:NR,NR+1:N1) = (0.0E+0,0.0E+0)
               CALL CLASET('A',M-NR,N1-NR,(0.0E+0,0.0E+0),(1.0E+0,0.0E+0), &
                           U(NR+1,NR+1),LDU)
            END IF
         END IF
!
!           The Q matrix from the first QRF is built into the left singular
!           matrix U. This applies to all cases.
!
         CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, &
              LDU, CWORK(N+1), LWORK-N, IERR )

!           The columns of U are normalized. The cost is O(M*N) flops.
         TEMP1 = SQRT(REAL(M)) * EPSLN
         DO p = 1, NR
            XSC = 1.0E+0 / SCNRM2( M, U(1,p), 1 )
            IF ( (XSC  <  (1.0E+0-TEMP1)) .OR. (XSC  >  (1.0E+0+TEMP1)) ) U(1:M,p) = U(1:M,p)*XSC
         ENDDO
!
!           If the initial QRF is computed with row pivoting, the left
!           singular vectors must be adjusted.
!
         IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
!
      ELSE
!
!        .. the initial matrix A has almost orthogonal columns and
!        the second QRF is not needed
!
         CALL CLACPY( 'U', N, N, A, LDA, CWORK(N+1), N )
         IF ( L2PERT ) THEN
            XSC = SQRT(SMALL)
            DO p = 2, N
               CTEMP = XSC * CWORK( N + (p-1)*N + p )
               DO q = 1, p - 1
!                     CWORK(N+(q-1)*N+p)=-TEMP1 * ( CWORK(N+(p-1)*N+q) /
!     $                                        ABS(CWORK(N+(p-1)*N+q)) )
                  CWORK(N+(q-1)*N+p)=-CTEMP
               ENDDO
            ENDDO
         ELSE
            CALL CLASET( 'L',N-1,N-1,(0.0E+0,0.0E+0),(0.0E+0,0.0E+0),CWORK(N+2),N )
         END IF
!
         CALL CGESVJ( 'U', 'U', 'N', N, N, CWORK(N+1), N, SVA, &
              N, U, LDU, CWORK(N+N*N+1), LWORK-N-N*N, RWORK, LRWORK, &
          INFO )
!
         SCALEM  = RWORK(1)
         NUMRANK = NINT(RWORK(2))
         DO p = 1, N
            U(1:N,p) = CWORK(p*N+1:N*(p+1))
            CWORK(p*N+1:(p+1)*N) = SVA(p)*CWORK(p*N+1:(p+1)*N)
         ENDDO
!
         CALL CTRSM( 'L', 'U', 'N', 'N', N, N, &
              (1.0E+0,0.0E+0), A, LDA, CWORK(N+1), N )
         DO p = 1, N
            CALL CCOPY( N, CWORK(N+p), N, V(IWORK(p),1), LDV )
         ENDDO
         TEMP1 = SQRT(REAL(N))*EPSLN
         DO p = 1, N
            XSC = 1.0E+0 / SCNRM2( N, V(1,p), 1 )
            IF ( (XSC  <  (1.0E+0-TEMP1)) .OR. (XSC  >  (1.0E+0+TEMP1)) ) V(1:N,p) = XSC*V(1:N,p)
         ENDDO
!
!           Assemble the left singular vector matrix U (M x N).
!
         IF ( N  <  M ) THEN
            U(N+1:M,1:N) = (0.0E+0,0.0E+0)
            IF ( N  <  N1 ) THEN
               U(1:N,N+1:N1) = (0.0E+0,0.0E+0)
               CALL CLASET( 'A',M-N,N1-N, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0),U(N+1,N+1),LDU)
            END IF
         END IF
         CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, &
              LDU, CWORK(N+1), LWORK-N, IERR )
         TEMP1 = SQRT(REAL(M))*EPSLN
         DO p = 1, N1
            XSC = 1.0E+0 / SCNRM2( M, U(1,p), 1 )
            IF ( (XSC  <  (1.0E+0-TEMP1)) .OR. (XSC  >  (1.0E+0+TEMP1)) ) U(1:M,p) = XSC * U(1:M,p)
         ENDDO
!
         IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
!
      END IF
!
!        end of the  >> almost orthogonal case <<  in the full SVD
!
      ELSE
!
!        This branch deploys a preconditioned Jacobi SVD with explicitly
!        accumulated rotations. It is included as optional, mainly for
!        experimental purposes. It does perform well, and can also be used.
!        In this implementation, this branch will be automatically activated
!        if the  condition number sigma_max(A) / sigma_min(A) is predicted
!        to be greater than the overflow threshold. This is because the
!        a posteriori computation of the singular vectors assumes robust
!        implementation of BLAS and some LAPACK procedures, capable of working
!        in presence of extreme values, e.g. when the singular values spread from
!        the underflow to the overflow threshold.
!
      DO p = 1, NR
         V(p:N,p) = CONJG(A(p,p:N))
      ENDDO
!
      IF ( L2PERT ) THEN
         XSC = SQRT(SMALL/EPSLN)
         DO q = 1, NR
            CTEMP = CMPLX(XSC*ABS( V(q,q) ),0.0E+0)
            DO p = 1, N
               IF ( ( p  >  q ) .AND. ( ABS(V(p,q))  <=  TEMP1 ) &
                   .OR. ( p  <  q ) ) V(p,q) = CTEMP
!     $                V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) &
               IF ( p  <  q ) V(p,q) = - V(p,q)
            ENDDO
         ENDDO
      ELSE
         CALL CLASET( 'U', NR-1, NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), V(1,2), LDV )
      END IF

      CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR )
      CALL CLACPY( 'L', N, NR, V, LDV, CWORK(2*N+1), N )
!
      DO p = 1, NR
         U(p:NR,p) = CONJG(V(p,p:NR))
      ENDDO

      IF ( L2PERT ) THEN
         XSC = SQRT(SMALL/EPSLN)
         DO q = 2, NR
            DO p = 1, q - 1
               CTEMP = CMPLX(XSC * MIN(ABS(U(p,p)),ABS(U(q,q))), 0.0E+0)
!                  U(p,q) = - TEMP1 * ( U(q,p) / ABS(U(q,p)) )
               U(p,q) = - CTEMP
            ENDDO
         ENDDO
      ELSE
         CALL CLASET('U', NR-1, NR-1, (0.0E+0,0.0E+0), (0.0E+0,0.0E+0), U(1,2), LDU )
      END IF

      CALL CGESVJ( 'L', 'U', 'V', NR, NR, U, LDU, SVA, &
           N, V, LDV, CWORK(2*N+N*NR+1), LWORK-2*N-N*NR, &
            RWORK, LRWORK, INFO )
      SCALEM  = RWORK(1)
      NUMRANK = NINT(RWORK(2))

      IF ( NR  <  N ) THEN
         V(NR+1:N,1:NR) = (0.0E+0,0.0E+0)
         V(1:NR,NR+1:N) = (0.0E+0,0.0E+0)
         V(NR+1:N,NR+1:N) = (0.0E+0,0.0E+0)
         DO I = NR+1, N
            V( I, I ) = (1.0E+0,0.0E+0)
         ENDDO
      END IF

      CALL CUNMQR( 'L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), &
           V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
!
!           Permute the rows of V using the (column) permutation from the
!           first QRF. Also, scale the columns to make them unit in
!           Euclidean norm. This applies to all cases.
!
         TEMP1 = SQRT(REAL(N)) * EPSLN
         DO q = 1, N
            DO p = 1, N
               CWORK(2*N+N*NR+NR+IWORK(p)) = V(p,q)
            ENDDO
            V(1:N,q) = CWORK((2+NR)*N+NR+1:(3+NR)*N+NR)
            XSC = 1.0E+0 / SCNRM2( N, V(1,q), 1 )
            IF ( (XSC  <  (1.0E+0-TEMP1)) .OR. (XSC  >  (1.0E+0+TEMP1)) ) V(1:N,q) = XSC*V(1:N,q)
         ENDDO
!
!           At this moment, V contains the right singular vectors of A.
!           Next, assemble the left singular vector matrix U (M x N).
!
      IF ( NR  <  M ) THEN
         U(NR+1:M,1:NR) = (0.0E+0,0.0E+0)
         IF ( NR  <  N1 ) THEN
            U(1:NR,NR+1:N1) = (0.0E+0,0.0E+0)
            CALL CLASET('A',M-NR,N1-NR, (0.0E+0,0.0E+0), (1.0E+0,0.0E+0),U(NR+1,NR+1),LDU)
         END IF
      END IF
!
      CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR )
!
      IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
!
!
      END IF
      IF ( TRANSP ) THEN
!           .. swap U and V because the procedure worked on A^*
         U_TMP(1:N,1:N) = U(1:N,1:N)
         U(1:N,1:N) = V(1:N,1:N)
         V(1:N,1:N) = U_TMP(1:N,1:N)
      END IF
!
   END IF
!     end of the full SVD
!
!     Undo scaling, if necessary (and possible)
!
   IF ( USCAL2  <=  (BIG/SVA(1))*USCAL1 ) THEN
      CALL SLASCL( 'G', 0, 0, USCAL1, USCAL2, NR, 1, SVA, N, IERR )
      USCAL1 = 1.0E+0
      USCAL2 = 1.0E+0
   END IF
!
   IF ( NR  <  N ) SVA(NR+1:N) = 0.0E+0
!
   RWORK(1) = USCAL2 * SCALEM
   RWORK(2) = USCAL1
   IF ( ERREST ) RWORK(3) = SCONDA
   IF ( LSVEC .AND. RSVEC ) THEN
      RWORK(4) = CONDR1
      RWORK(5) = CONDR2
   END IF
   IF ( L2TRAN ) THEN
      RWORK(6) = ENTRA
      RWORK(7) = ENTRAT
   END IF
!
   IWORK(1) = NR
   IWORK(2) = NUMRANK
   IWORK(3) = WARNING
   IF ( TRANSP ) THEN
       IWORK(4) =  1
   ELSE
       IWORK(4) = -1
   END IF

!
   RETURN
!     ..
!     .. END OF CGEJSV
!     ..
   END

