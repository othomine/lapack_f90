!> \brief \b DLATME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATME( N, DIST, ISEED, D, MODE, COND, DMAX, EI,
!         RSIGN,
!                          UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM,
!         A,
!                          LDA, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIST, RSIGN, SIM, UPPER
!       INTEGER            INFO, KL, KU, LDA, MODE, MODES, N
!       DOUBLE PRECISION   ANORM, COND, CONDS, DMAX
!       ..
!       .. Array Arguments ..
!       CHARACTER          EI( * )
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * ), D( * ), DS( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLATME generates random non-symmetric square matrices with
!>    specified eigenvalues for testing LAPACK programs.
!>
!>    DLATME operates by applying the following sequence of
!>    operations:
!>
!>    1. Set the diagonal to D, where D may be input or
!>         computed according to MODE, COND, DMAX, and RSIGN
!>         as described below.
!>
!>    2. If complex conjugate pairs are desired (MODE=0 and EI(1)='R',
!>         or MODE=5), certain pairs of adjacent elements of D are
!>         interpreted as the real and complex parts of a complex
!>         conjugate pair; A thus becomes block diagonal, with 1x1
!>         and 2x2 blocks.
!>
!>    3. If UPPER='T', the upper triangle of A is set to random values
!>         out of distribution DIST.
!>
!>    4. If SIM='T', A is multiplied on the left by a random matrix
!>         X, whose singular values are specified by DS, MODES, and
!>         CONDS, and on the right by X inverse.
!>
!>    5. If KL < N-1, the lower bandwidth is reduced to KL using
!>         Householder transformations.  If KU < N-1, the upper
!>         bandwidth is reduced to KU.
!>
!>    6. If ANORM is not negative, the matrix is scaled to have
!>         maximum-element-norm ANORM.
!>
!>    (Note: since the matrix cannot be reduced beyond Hessenberg form,
!>     no packing options are available.)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The number of columns (or rows) of A. Not modified.
!> \endverbatim
!>
!> \param[in] DIST
!> \verbatim
!>          DIST is CHARACTER*1
!>           On entry, DIST specifies the type of distribution to be used
!>           to generate the random eigen-/singular values, and for the
!>           upper triangle (see UPPER).
!>           'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )
!>           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )
!>           'N' => NORMAL( 0, 1 )   ( 'N' for normal )
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. They should lie between 0 and 4095 inclusive,
!>           and ISEED(4) should be odd. The random number generator
!>           uses a linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to DLATME
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension ( N )
!>           This array is used to specify the eigenvalues of A.  If
!>           MODE=0, then D is assumed to contain the eigenvalues (but
!>           see the description of EI), otherwise they will be
!>           computed according to MODE, COND, DMAX, and RSIGN and
!>           placed in D.
!>           Modified if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry this describes how the eigenvalues are to
!>           be specified:
!>           MODE = 0 means use D (with EI) as input
!>           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
!>           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
!>           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
!>           MODE = 5 sets D to random numbers in the range
!>                    ( 1/COND , 1 ) such that their logarithms
!>                    are uniformly distributed.  Each odd-even pair
!>                    of elements will be either used as two real
!>                    eigenvalues or as the real and imaginary part
!>                    of a complex conjugate pair of eigenvalues;
!>                    the choice of which is done is random, with
!>                    50-50 probability, for each pair.
!>           MODE = 6 set D to random numbers from same distribution
!>                    as the rest of the matrix.
!>           MODE < 0 has the same meaning as ABS(MODE), except that
!>              the order of the elements of D is reversed.
!>           Thus if MODE is between 1 and 4, D has entries ranging
!>              from 1 to 1/COND, if between -1 and -4, D has entries
!>              ranging from 1/COND to 1,
!>           Not modified.
!> \endverbatim
!>
!> \param[in] COND
!> \verbatim
!>          COND is DOUBLE PRECISION
!>           On entry, this is used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] DMAX
!> \verbatim
!>          DMAX is DOUBLE PRECISION
!>           If MODE is neither -6, 0 nor 6, the contents of D, as
!>           computed according to MODE and COND, will be scaled by
!>           DMAX / max(abs(D(i))).  Note that DMAX need not be
!>           positive: if DMAX is negative (or zero), D will be
!>           scaled by a negative number (or zero).
!>           Not modified.
!> \endverbatim
!>
!> \param[in] EI
!> \verbatim
!>          EI is CHARACTER*1 array, dimension ( N )
!>           If MODE is 0, and EI(1) is not ' ' (space character),
!>           this array specifies which elements of D (on input) are
!>           real eigenvalues and which are the real and imaginary parts
!>           of a complex conjugate pair of eigenvalues.  The elements
!>           of EI may then only have the values 'R' and 'I'.  If
!>           EI(j)='R' and EI(j+1)='I', then the j-th eigenvalue is
!>           CMPLX( D(j) , D(j+1) ), and the (j+1)-th is the complex
!>           conjugate thereof.  If EI(j)=EI(j+1)='R', then the j-th
!>           eigenvalue is D(j) (i.e., real).  EI(1) may not be 'I',
!>           nor may two adjacent elements of EI both have the value 'I'.
!>           If MODE is not 0, then EI is ignored.  If MODE is 0 and
!>           EI(1)=' ', then the eigenvalues will all be real.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] RSIGN
!> \verbatim
!>          RSIGN is CHARACTER*1
!>           If MODE is not 0, 6, or -6, and RSIGN='T', then the
!>           elements of D, as computed according to MODE and COND, will
!>           be multiplied by a random sign (+1 or -1).  If RSIGN='F',
!>           they will not be.  RSIGN may only have the values 'T' or
!>           'F'.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] UPPER
!> \verbatim
!>          UPPER is CHARACTER*1
!>           If UPPER='T', then the elements of A above the diagonal
!>           (and above the 2x2 diagonal blocks, if A has complex
!>           eigenvalues) will be set to random numbers out of DIST.
!>           If UPPER='F', they will not.  UPPER may only have the
!>           values 'T' or 'F'.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] SIM
!> \verbatim
!>          SIM is CHARACTER*1
!>           If SIM='T', then A will be operated on by a "similarity
!>           transform", i.e., multiplied on the left by a matrix X and
!>           on the right by X inverse.  X = U S V, where U and V are
!>           random unitary matrices and S is a (diagonal) matrix of
!>           singular values specified by DS, MODES, and CONDS.  If
!>           SIM='F', then A will not be transformed.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] DS
!> \verbatim
!>          DS is DOUBLE PRECISION array, dimension ( N )
!>           This array is used to specify the singular values of X,
!>           in the same way that D specifies the eigenvalues of A.
!>           If MODE=0, the DS contains the singular values, which
!>           may not be zero.
!>           Modified if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODES
!> \verbatim
!>          MODES is INTEGER
!> \endverbatim
!>
!> \param[in] CONDS
!> \verbatim
!>          CONDS is DOUBLE PRECISION
!>           Same as MODE and COND, but for specifying the diagonal
!>           of S.  MODES=-6 and +6 are not allowed (since they would
!>           result in randomly ill-conditioned eigenvalues.)
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           This specifies the lower bandwidth of the  matrix.  KL=1
!>           specifies upper Hessenberg form.  If KL is at least N-1,
!>           then A will have full lower bandwidth.  KL must be at
!>           least 1.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           This specifies the upper bandwidth of the  matrix.  KU=1
!>           specifies lower Hessenberg form.  If KU is at least N-1,
!>           then A will have full upper bandwidth; if KU and KL
!>           are both at least N-1, then A will be dense.  Only one of
!>           KU and KL may be less than N-1.  KU must be at least 1.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is DOUBLE PRECISION
!>           If ANORM is not negative, then A will be scaled by a non-
!>           negative real number to make the maximum-element-norm of A
!>           to be ANORM.
!>           Not modified.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, N )
!>           On exit A is the desired test matrix.
!>           Modified.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           LDA specifies the first dimension of A as declared in the
!>           calling program.  LDA must be at least N.
!>           Not modified.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension ( 3*N )
!>           Workspace.
!>           Modified.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           Error code.  On exit, INFO will be set to one of the
!>           following values:
!>             0 => normal return
!>            -1 => N negative
!>            -2 => DIST illegal string
!>            -5 => MODE not in range -6 to 6
!>            -6 => COND less than 1.0, and MODE neither -6, 0 nor 6
!>            -8 => EI(1) is not ' ' or 'R', EI(j) is not 'R' or 'I', or
!>                  two adjacent elements of EI are 'I'.
!>            -9 => RSIGN is not 'T' or 'F'
!>           -10 => UPPER is not 'T' or 'F'
!>           -11 => SIM   is not 'T' or 'F'
!>           -12 => MODES=0 and DS has a zero singular value.
!>           -13 => MODES is not in the range -5 to 5.
!>           -14 => MODES is nonzero and CONDS is less than 1.
!>           -15 => KL is less than 1.
!>           -16 => KU is less than 1, or KL and KU are both less than
!>                  N-1.
!>           -19 => LDA is less than N.
!>            1  => Error return from DLATM1 (computing D)
!>            2  => Cannot scale to DMAX (max. eigenvalue is 0)
!>            3  => Error return from DLATM1 (computing DS)
!>            4  => Error return from DLARGE
!>            5  => Zero singular value from DLATM1.
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
!> \ingroup double_matgen
!
!  =====================================================================
   SUBROUTINE DLATME( N, DIST, ISEED, D, MODE, COND, DMAX, EI, &
     RSIGN, &
                      UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM, &
     A, &
                      LDA, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          DIST, RSIGN, SIM, UPPER
   INTEGER            INFO, KL, KU, LDA, MODE, MODES, N
   DOUBLE PRECISION   ANORM, COND, CONDS, DMAX
!     ..
!     .. Array Arguments ..
   CHARACTER          EI( * )
   INTEGER            ISEED( 4 )
   DOUBLE PRECISION   A( LDA, * ), D( * ), DS( * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            BADEI, BADS, USEEI
   INTEGER            I, IC, ICOLS, IDIST, IINFO, IR, IROWS, IRSIGN, &
                      ISIM, IUPPER, J, JC, JCR, JR
   DOUBLE PRECISION   ALPHA, TAU, TEMP, XNORMS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   TEMPA( 1 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLANGE, DLARAN
   EXTERNAL           LSAME, DLANGE, DLARAN
!     ..
!     .. External Subroutines ..
   EXTERNAL           DCOPY, DGEMV, DGER, DLARFG, DLARGE, DLARNV, &
                      DLASET, DLATM1, DSCAL, XERBLA
!     ..
!     .. Executable Statements ..
!
!     1)      Decode and Test the input parameters.
!             Initialize flags & seed.
!
   INFO = 0
!
!     Quick return if possible
!
   IF( N == 0 ) RETURN
!
!     Decode DIST
!
   IF( LSAME( DIST, 'U' ) ) THEN
      IDIST = 1
   ELSE IF( LSAME( DIST, 'S' ) ) THEN
      IDIST = 2
   ELSE IF( LSAME( DIST, 'N' ) ) THEN
      IDIST = 3
   ELSE
      IDIST = -1
   END IF
!
!     Check EI
!
   USEEI = .TRUE.
   BADEI = .FALSE.
   IF( LSAME( EI( 1 ), ' ' ) .OR. MODE /= 0 ) THEN
      USEEI = .FALSE.
   ELSE
      IF( LSAME( EI( 1 ), 'R' ) ) THEN
         DO J = 2, N
            IF( LSAME( EI( J ), 'I' ) ) THEN
               IF( LSAME( EI( J-1 ), 'I' ) ) &
                  BADEI = .TRUE.
            ELSE
               IF( .NOT.LSAME( EI( J ), 'R' ) ) &
                  BADEI = .TRUE.
            END IF
         ENDDO
      ELSE
         BADEI = .TRUE.
      END IF
   END IF
!
!     Decode RSIGN
!
   IF( LSAME( RSIGN, 'T' ) ) THEN
      IRSIGN = 1
   ELSE IF( LSAME( RSIGN, 'F' ) ) THEN
      IRSIGN = 0
   ELSE
      IRSIGN = -1
   END IF
!
!     Decode UPPER
!
   IF( LSAME( UPPER, 'T' ) ) THEN
      IUPPER = 1
   ELSE IF( LSAME( UPPER, 'F' ) ) THEN
      IUPPER = 0
   ELSE
      IUPPER = -1
   END IF
!
!     Decode SIM
!
   IF( LSAME( SIM, 'T' ) ) THEN
      ISIM = 1
   ELSE IF( LSAME( SIM, 'F' ) ) THEN
      ISIM = 0
   ELSE
      ISIM = -1
   END IF
!
!     Check DS, if MODES=0 and ISIM=1
!
   BADS = .FALSE.
   IF( MODES == 0 .AND. ISIM == 1 ) BADS = ANY(DS(1:N) == 0.0D+0 )
!
!     Set INFO if an error
!
   IF( N < 0 ) THEN
      INFO = -1
   ELSE IF( IDIST == -1 ) THEN
      INFO = -2
   ELSE IF( ABS( MODE ) > 6 ) THEN
      INFO = -5
   ELSE IF( ( MODE /= 0 .AND. ABS( MODE ) /= 6 ) .AND. COND < 1.0D+0 ) &
             THEN
      INFO = -6
   ELSE IF( BADEI ) THEN
      INFO = -8
   ELSE IF( IRSIGN == -1 ) THEN
      INFO = -9
   ELSE IF( IUPPER == -1 ) THEN
      INFO = -10
   ELSE IF( ISIM == -1 ) THEN
      INFO = -11
   ELSE IF( BADS ) THEN
      INFO = -12
   ELSE IF( ISIM == 1 .AND. ABS( MODES ) > 5 ) THEN
      INFO = -13
   ELSE IF( ISIM == 1 .AND. MODES /= 0 .AND. CONDS < 1.0D+0 ) THEN
      INFO = -14
   ELSE IF( KL < 1 ) THEN
      INFO = -15
   ELSE IF( KU < 1 .OR. ( KU < N-1 .AND. KL < N-1 ) ) THEN
      INFO = -16
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -19
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DLATME', -INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      RETURN
   END IF
!
!     Initialize random number generator
!
   ISEED(1:4) = MOD( ABS( ISEED(1:4) ), 4096 )
!
   IF( MOD( ISEED( 4 ), 2 ) /= 1 ) ISEED( 4 ) = ISEED( 4 ) + 1
!
!     2)      Set up diagonal of A
!
!             Compute D according to COND and MODE
!
   CALL DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO )
   IF( IINFO /= 0 ) THEN
      INFO = 1
      RETURN
   END IF
   IF( MODE /= 0 .AND. ABS( MODE ) /= 6 ) THEN
!
!        Scale by DMAX
!
      TEMP = MAXVAL(ABS( D(1:N) ) )
!
      IF( TEMP > 0.0D+0 ) THEN
         ALPHA = DMAX / TEMP
      ELSE IF( DMAX /= 0.0D+0 ) THEN
         INFO = 2
         RETURN
      ELSE
         ALPHA = 0.0D+0
      END IF
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSCAL( N, ALPHA, D, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   END IF
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', N, N, 0.0D+0, 0.0D+0, A, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DCOPY( N, D, 1, A, LDA+1 )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Set up complex conjugate pairs
!
   IF( MODE == 0 ) THEN
      IF( USEEI ) THEN
         DO J = 2, N
            IF( LSAME( EI( J ), 'I' ) ) THEN
               A( J-1, J ) = A( J, J )
               A( J, J-1 ) = -A( J, J )
               A( J, J ) = A( J-1, J-1 )
            END IF
         ENDDO
      END IF
!
   ELSE IF( ABS( MODE ) == 5 ) THEN
!
      DO J = 2, N, 2
         IF( DLARAN( ISEED ) > 0.5D+0 ) THEN
            A( J-1, J ) = A( J, J )
            A( J, J-1 ) = -A( J, J )
            A( J, J ) = A( J-1, J-1 )
         END IF
      ENDDO
   END IF
!
!     3)      If UPPER='T', set upper triangle of A to random numbers.
!             (but don't modify the corners of 2x2 blocks.)
!
   IF( IUPPER /= 0 ) THEN
      DO JC = 2, N
         IF( A( JC-1, JC ) /= 0.0D+0 ) THEN
            JR = JC - 2
         ELSE
            JR = JC - 1
         END IF
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLARNV( IDIST, ISEED, JR, A( 1, JC ) )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
!
!     4)      If SIM='T', apply similarity transformation.
!
!                                -1
!             Transform is  X A X  , where X = U S V, thus
!
!             it is  U S V A V' (1/S) U'
!
   IF( ISIM /= 0 ) THEN
!
!        Compute S (singular values of the eigenvector matrix)
!        according to CONDS and MODES
!
      CALL DLATM1( MODES, CONDS, 0, 0, ISEED, DS, N, IINFO )
      IF( IINFO /= 0 ) THEN
         INFO = 3
         RETURN
      END IF
!
!        Multiply by V and V'
!
      CALL DLARGE( N, A, LDA, ISEED, WORK, IINFO )
      IF( IINFO /= 0 ) THEN
         INFO = 4
         RETURN
      END IF
!
!        Multiply by S and (1/S)
!
      DO J = 1, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSCAL( N, DS( J ), A( J, 1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IF( DS( J ) /= 0.0D+0 ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSCAL( N, 1.0D+0 / DS( J ), A( 1, J ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ELSE
            INFO = 5
            RETURN
         END IF
      ENDDO
!
!        Multiply by U and U'
!
      CALL DLARGE( N, A, LDA, ISEED, WORK, IINFO )
      IF( IINFO /= 0 ) THEN
         INFO = 4
         RETURN
      END IF
   END IF
!
!     5)      Reduce the bandwidth.
!
   IF( KL < N-1 ) THEN
!
!        Reduce bandwidth -- kill column
!
      DO JCR = KL + 1, N - 1
         IC = JCR - KL
         IROWS = N + 1 - JCR
         ICOLS = N + KL - JCR
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DCOPY( IROWS, A( JCR, IC ), 1, WORK, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         XNORMS = WORK( 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLARFG( IROWS, XNORMS, WORK( 2 ), 1, TAU )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLARFG : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WORK( 1 ) = 1.0D+0
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEMV( 'T', IROWS, ICOLS, 1.0D+0, A( JCR, IC+1 ), LDA, &
                     WORK, 1, 0.0D+0, WORK( IROWS+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGEMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGER( IROWS, ICOLS, -TAU, WORK, 1, WORK( IROWS+1 ), 1, &
                    A( JCR, IC+1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGER : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEMV( 'N', N, IROWS, 1.0D+0, A( 1, JCR ), LDA, WORK, 1, &
                     0.0D+0, WORK( IROWS+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGEMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGER( N, IROWS, -TAU, WORK( IROWS+1 ), 1, WORK, 1, &
                    A( 1, JCR ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGER : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         A( JCR, IC ) = XNORMS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLASET( 'Full', IROWS-1, 1, 0.0D+0, 0.0D+0, A( JCR+1, IC ), &
                      LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   ELSE IF( KU < N-1 ) THEN
!
!        Reduce upper bandwidth -- kill a row at a time.
!
      DO JCR = KU + 1, N - 1
         IR = JCR - KU
         IROWS = N + KU - JCR
         ICOLS = N + 1 - JCR
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DCOPY( ICOLS, A( IR, JCR ), LDA, WORK, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         XNORMS = WORK( 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLARFG( ICOLS, XNORMS, WORK( 2 ), 1, TAU )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLARFG : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         WORK( 1 ) = 1.0D+0
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEMV( 'N', IROWS, ICOLS, 1.0D+0, A( IR+1, JCR ), LDA, &
                     WORK, 1, 0.0D+0, WORK( ICOLS+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGEMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGER( IROWS, ICOLS, -TAU, WORK( ICOLS+1 ), 1, WORK, 1, &
                    A( IR+1, JCR ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGER : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEMV( 'C', ICOLS, N, 1.0D+0, A( JCR, 1 ), LDA, WORK, 1, &
                     0.0D+0, WORK( ICOLS+1 ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGEMV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGER( ICOLS, N, -TAU, WORK, 1, WORK( ICOLS+1 ), 1, &
                    A( JCR, 1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGER : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         A( IR, JCR ) = XNORMS
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLASET( 'Full', 1, ICOLS-1, 0.0D+0, 0.0D+0, A( IR, JCR+1 ), &
                      LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
!
!     Scale the matrix to have norm ANORM
!
   IF( ANORM >= 0.0D+0 ) THEN
      TEMP = DLANGE( 'M', N, N, A, LDA, TEMPA )
      IF( TEMP > 0.0D+0 ) THEN
         ALPHA = ANORM / TEMP
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DSCAL( N, ALPHA, A( 1, J ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDDO
      END IF
   END IF
!
   RETURN
!
!     End of DLATME
!
END

