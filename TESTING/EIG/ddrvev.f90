!> \brief \b DDRVEV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVEV( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL,
!                          VR, LDVR, LRE, LDLRE, RESULT, WORK, NWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDLRE, LDVL, LDVR, NOUNIT, NSIZES,
!      $                   NTYPES, NWORK
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       DOUBLE PRECISION   A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ),
!      $                   RESULT( 7 ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WI( * ), WI1( * ), WORK( * ), WR( * ), WR1( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DDRVEV  checks the nonsymmetric eigenvalue problem driver DGEEV.
!>
!>    When DDRVEV is called, a number of matrix "sizes" ("n's") and a
!>    number of matrix "types" are specified.  For each size ("n")
!>    and each type of matrix, one matrix will be generated and used
!>    to test the nonsymmetric eigenroutines.  For each matrix, 7
!>    tests will be performed:
!>
!>    (1)     | A * VR - VR * W | / ( n |A| ulp )
!>
!>      Here VR is the matrix of unit right eigenvectors.
!>      W is a block diagonal matrix, with a 1x1 block for each
!>      real eigenvalue and a 2x2 block for each complex conjugate
!>      pair.  If eigenvalues j and j+1 are a complex conjugate pair,
!>      so WR(j) = WR(j+1) = wr and WI(j) = - WI(j+1) = wi, then the
!>      2 x 2 block corresponding to the pair will be:
!>
!>              (  wr  wi  )
!>              ( -wi  wr  )
!>
!>      Such a block multiplying an n x 2 matrix  ( ur ui ) on the
!>      right will be the same as multiplying  ur + i*ui  by  wr + i*wi.
!>
!>    (2)     | A**H * VL - VL * W**H | / ( n |A| ulp )
!>
!>      Here VL is the matrix of unit left eigenvectors, A**H is the
!>      conjugate transpose of A, and W is as above.
!>
!>    (3)     | |VR(i)| - 1 | / ulp and whether largest component real
!>
!>      VR(i) denotes the i-th column of VR.
!>
!>    (4)     | |VL(i)| - 1 | / ulp and whether largest component real
!>
!>      VL(i) denotes the i-th column of VL.
!>
!>    (5)     W(full) = W(partial)
!>
!>      W(full) denotes the eigenvalues computed when both VR and VL
!>      are also computed, and W(partial) denotes the eigenvalues
!>      computed when only W, only W and VR, or only W and VL are
!>      computed.
!>
!>    (6)     VR(full) = VR(partial)
!>
!>      VR(full) denotes the right eigenvectors computed when both VR
!>      and VL are computed, and VR(partial) denotes the result
!>      when only VR is computed.
!>
!>     (7)     VL(full) = VL(partial)
!>
!>      VL(full) denotes the left eigenvectors computed when both VR
!>      and VL are also computed, and VL(partial) denotes the result
!>      when only VL is computed.
!>
!>    The "sizes" are specified by an array NN(1:NSIZES); the value of
!>    each element NN(j) specifies one size.
!>    The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!>    if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!>    Currently, the list of possible types is:
!>
!>    (1)  The zero matrix.
!>    (2)  The identity matrix.
!>    (3)  A (transposed) Jordan block, with 1's on the diagonal.
!>
!>    (4)  A diagonal matrix with evenly spaced entries
!>         1, ..., ULP  and random signs.
!>         (ULP = (first number larger than 1) - 1 )
!>    (5)  A diagonal matrix with geometrically spaced entries
!>         1, ..., ULP  and random signs.
!>    (6)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>         and random signs.
!>
!>    (7)  Same as (4), but multiplied by a constant near
!>         the overflow threshold
!>    (8)  Same as (4), but multiplied by a constant near
!>         the underflow threshold
!>
!>    (9)  A matrix of the form  U' T U, where U is orthogonal and
!>         T has evenly spaced entries 1, ..., ULP with random signs
!>         on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (10) A matrix of the form  U' T U, where U is orthogonal and
!>         T has geometrically spaced entries 1, ..., ULP with random
!>         signs on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (11) A matrix of the form  U' T U, where U is orthogonal and
!>         T has "clustered" entries 1, ULP,..., ULP with random
!>         signs on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (12) A matrix of the form  U' T U, where U is orthogonal and
!>         T has real or complex conjugate paired eigenvalues randomly
!>         chosen from ( ULP, 1 ) and random O(1) entries in the upper
!>         triangle.
!>
!>    (13) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has evenly spaced entries 1, ..., ULP
!>         with random signs on the diagonal and random O(1) entries
!>         in the upper triangle.
!>
!>    (14) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has geometrically spaced entries
!>         1, ..., ULP with random signs on the diagonal and random
!>         O(1) entries in the upper triangle.
!>
!>    (15) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has "clustered" entries 1, ULP,..., ULP
!>         with random signs on the diagonal and random O(1) entries
!>         in the upper triangle.
!>
!>    (16) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has real or complex conjugate paired
!>         eigenvalues randomly chosen from ( ULP, 1 ) and random
!>         O(1) entries in the upper triangle.
!>
!>    (17) Same as (16), but multiplied by a constant
!>         near the overflow threshold
!>    (18) Same as (16), but multiplied by a constant
!>         near the underflow threshold
!>
!>    (19) Nonsymmetric matrix with random entries chosen from (-1,1).
!>         If N is at least 4, all entries in first two rows and last
!>         row, and first column and last two columns are zero.
!>    (20) Same as (19), but multiplied by a constant
!>         near the overflow threshold
!>    (21) Same as (19), but multiplied by a constant
!>         near the underflow threshold
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          DDRVEV does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, DDRVEV
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A.  This
!>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size in NN a
!>          matrix of that size and of type j will be generated.
!>          If NTYPES is smaller than the maximum number of types
!>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>          MAXTYP will not be generated.  If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to DDRVEV to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns INFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, max(NN))
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.  On exit, A contains the last matrix actually used.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, and H. LDA must be at
!>          least 1 and at least max(NN).
!> \endverbatim
!>
!> \param[out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDA, max(NN))
!>          Another copy of the test matrix A, modified by DGEEV.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (max(NN))
!>
!>          The real and imaginary parts of the eigenvalues of A.
!>          On exit, WR + WI*i are the eigenvalues of the matrix in A.
!> \endverbatim
!>
!> \param[out] WR1
!> \verbatim
!>          WR1 is DOUBLE PRECISION array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] WI1
!> \verbatim
!>          WI1 is DOUBLE PRECISION array, dimension (max(NN))
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but those computed when DGEEV only computes a partial
!>          eigendecomposition, i.e. not the eigenvalues and left
!>          and right eigenvectors.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL, max(NN))
!>          VL holds the computed left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          Leading dimension of VL. Must be at least max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (LDVR, max(NN))
!>          VR holds the computed right eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          Leading dimension of VR. Must be at least max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] LRE
!> \verbatim
!>          LRE is DOUBLE PRECISION array, dimension (LDLRE,max(NN))
!>          LRE holds the computed right or left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDLRE
!> \verbatim
!>          LDLRE is INTEGER
!>          Leading dimension of LRE. Must be at least max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (7)
!>          The values computed by the seven tests described above.
!>          The values are currently limited to 1/ulp, to avoid overflow.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (NWORK)
!> \endverbatim
!>
!> \param[in] NWORK
!> \verbatim
!>          NWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          5*NN(j)+2*NN(j)**2 for all j.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -6: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -16: LDVL < 1 or LDVL < NMAX, where NMAX is max( NN(j) ).
!>          -18: LDVR < 1 or LDVR < NMAX, where NMAX is max( NN(j) ).
!>          -20: LDLRE < 1 or LDLRE < NMAX, where NMAX is max( NN(j) ).
!>          -23: NWORK too small.
!>          If  DLATMR, SLATMS, SLATME or DGEEV returns an error code,
!>              the absolute value of it is returned.
!>
!>-----------------------------------------------------------------------
!>
!>     Some Local Variables and Parameters:
!>     ---- ----- --------- --- ----------
!>
!>     0.0D0, 1.0D0       Real 0 and 1.
!>     MAXTYP          The number of types defined.
!>     NMAX            Largest value in NN.
!>     NERRS           The number of tests which have exceeded THRESH
!>     COND, CONDS,
!>     IMODE           Values to be passed to the matrix generators.
!>     ANORM           Norm of A; passed to matrix generators.
!>
!>     OVFL, UNFL      Overflow and underflow thresholds.
!>     ULP, ULPINV     Finest relative precision and its inverse.
!>     RTULP, RTULPI   Square roots of the previous 4 values.
!>
!>             The following four arrays decode JTYPE:
!>     KTYPE(j)        The general type (1-10) for type "j".
!>     KMODE(j)        The MODE value to be passed to the matrix
!>                     generator for type "j".
!>     KMAGN(j)        The order of magnitude ( O(1),
!>                     O(overflow^(1/2) ), O(underflow^(1/2) )
!>     KCONDS(j)       Selectw whether CONDS is to be 1 or
!>                     1/sqrt(ulp).  (0 means irrelevant.)
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
   SUBROUTINE DDRVEV( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, &
                      NOUNIT, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, &
                      VR, LDVR, LRE, LDLRE, RESULT, WORK, NWORK, &
                      IWORK, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDLRE, LDVL, LDVR, NOUNIT, NSIZES, &
                      NTYPES, NWORK
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
   DOUBLE PRECISION   A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), &
                      RESULT( 7 ), VL( LDVL, * ), VR( LDVR, * ), &
                      WI( * ), WI1( * ), WORK( * ), WR( * ), WR1( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 21 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADNN
   CHARACTER*3        PATH
   INTEGER            IINFO, IMODE, ITYPE, IWK, J, JCOL, JJ, JSIZE, &
                      JTYPE, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, &
                      NTEST, NTESTF, NTESTT
   DOUBLE PRECISION   ANORM, COND, CONDS, OVFL, RTULP, RTULPI, TNRM, &
                      ULP, ULPINV, UNFL, VMX, VRMX, VTST
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   CHARACTER          ADUMMA( 1 )
   INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), &
                      KMAGN( MAXTYP ), KMODE( MAXTYP ), &
                      KTYPE( MAXTYP )
   DOUBLE PRECISION   DUM( 1 ), RES( 2 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
   EXTERNAL           DLAMCH, DLAPY2, DNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEEV, DGET22, DLACPY, DLASET, DLASUM, &
                      DLATME, DLATMR, DLATMS, XERBLA
!     ..
!     .. Data statements ..
   DATA               KTYPE / 1, 2, 3, 5*4, 4*6, 6*6, 3*9 /
   DATA               KMAGN / 3*1, 1, 1, 1, 2, 3, 4*1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /
   DATA               KMODE / 3*0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /
   DATA               KCONDS / 3*0, 5*0, 4*1, 6*2, 3*0 /
!     ..
!     .. Executable Statements ..
!
   PATH( 1: 1 ) = 'Double precision'
   PATH( 2: 3 ) = 'EV'
!
!     Check for errors
!
   NTESTT = 0
   NTESTF = 0
   INFO = 0
!
!     Important constants
!
   BADNN = ANY(NN(1:NSIZES) < 0)
   NMAX = MAXVAL(NN(1:NSIZES))
!
!     Check for errors
!
   IF( NSIZES < 0 ) THEN
      INFO = -1
   ELSE IF( BADNN ) THEN
      INFO = -2
   ELSE IF( NTYPES < 0 ) THEN
      INFO = -3
   ELSE IF( THRESH < 0.0D0 ) THEN
      INFO = -6
   ELSE IF( NOUNIT <= 0 ) THEN
      INFO = -7
   ELSE IF( LDA < 1 .OR. LDA < NMAX ) THEN
      INFO = -9
   ELSE IF( LDVL < 1 .OR. LDVL < NMAX ) THEN
      INFO = -16
   ELSE IF( LDVR < 1 .OR. LDVR < NMAX ) THEN
      INFO = -18
   ELSE IF( LDLRE < 1 .OR. LDLRE < NMAX ) THEN
      INFO = -20
   ELSE IF( 5*NMAX+2*NMAX**2 > NWORK ) THEN
      INFO = -23
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DDRVEV', -INFO )
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
!     Quick return if nothing to do
!
   IF( NSIZES == 0 .OR. NTYPES == 0 ) RETURN
!
!     More Important constants
!
   UNFL = DLAMCH( 'Safe minimum' )
   OVFL = 1.0D0 / UNFL
   ULP = DLAMCH( 'Precision' )
   ULPINV = 1.0D0 / ULP
   RTULP = SQRT( ULP )
   RTULPI = 1.0D0 / RTULP
!
!     Loop over sizes, types
!
   NERRS = 0
!
   DO JSIZE = 1, NSIZES
      N = NN( JSIZE )
      IF( NSIZES /= 1 ) THEN
         MTYPES = MIN( MAXTYP, NTYPES )
      ELSE
         MTYPES = MIN( MAXTYP+1, NTYPES )
      END IF
!
      DO JTYPE = 1, MTYPES
         IF( .NOT.DOTYPE( JTYPE ) ) GO TO 260
!
!           Save ISEED in case of an error.
!
         IOLDSD(1:4) = ISEED(1:4)
!
!           Compute "A"
!
!           Control parameters:
!
!           KMAGN  KCONDS  KMODE        KTYPE
!       =1  O(1)   1       clustered 1  zero
!       =2  large  large   clustered 2  identity
!       =3  small          exponential  Jordan
!       =4                 arithmetic   diagonal, (w/ eigenvalues)
!       =5                 random log   symmetric, w/ eigenvalues
!       =6                 random       general, w/ eigenvalues
!       =7                              random diagonal
!       =8                              random symmetric
!       =9                              random general
!       =10                             random triangular
!
         IF( MTYPES > MAXTYP ) GO TO 90
!
         ITYPE = KTYPE( JTYPE )
         IMODE = KMODE( JTYPE )
!
!           Compute norm
!
         SELECT CASE (KMAGN(JTYPE))
          CASE (1)
           ANORM = 1.0D+0
          CASE (2)
           ANORM = OVFL*ULP
          CASE (3)
           ANORM = UNFL*ULPINV
         END SELECT
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLASET( 'Full', LDA, N, 0.0D0, 0.0D0, A, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         IINFO = 0
         COND = ULPINV
!
!           Special Matrices -- Identity & Jordan block
!
!              Zero
!
         IF( ITYPE == 1 ) THEN
            IINFO = 0
!
         ELSE IF( ITYPE == 2 ) THEN
!
!              Identity
!
            FORALL(JCOL = 1:N) A( JCOL, JCOL ) = ANORM
!
         ELSE IF( ITYPE == 3 ) THEN
!
!              Jordan Block
!
            DO JCOL = 1, N
               A( JCOL, JCOL ) = ANORM
               IF( JCOL > 1 ) A( JCOL, JCOL-1 ) = 1.0D0
            ENDDO
!
         ELSE IF( ITYPE == 4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
            CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, &
                         ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 5 ) THEN
!
!              Symmetric, eigenvalues specified
!
            CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, &
                         ANORM, N, N, 'N', A, LDA, WORK( N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 6 ) THEN
!
!              General, eigenvalues specified
!
            IF( KCONDS( JTYPE ) == 1 ) THEN
               CONDS = 1.0D0
            ELSE IF( KCONDS( JTYPE ) == 2 ) THEN
               CONDS = RTULPI
            ELSE
               CONDS = 0.0D0
            END IF
!
            ADUMMA( 1 ) = ' '
            CALL DLATME( N, 'S', ISEED, WORK, IMODE, COND, 1.0D0, &
                         ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4, &
                         CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), &
                         IINFO )
!
         ELSE IF( ITYPE == 7 ) THEN
!
!              Diagonal, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, 0, 0, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 8 ) THEN
!
!              Symmetric, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, N, N, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE IF( ITYPE == 9 ) THEN
!
!              General, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, N, N, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
            IF( N >= 4 ) THEN
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL DLASET( 'Full', 2, N, 0.0D0, 0.0D0, A, LDA )
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
               CALL DLASET( 'Full', N-3, 1, 0.0D0, 0.0D0, A( 3, 1 ), &
                            LDA )
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
               CALL DLASET( 'Full', N-3, 2, 0.0D0, 0.0D0, A( 3, N-1 ), &
                            LDA )
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
               CALL DLASET( 'Full', 1, N, 0.0D0, 0.0D0, A( N, 1 ), &
                            LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
            END IF
!
         ELSE IF( ITYPE == 10 ) THEN
!
!              Triangular, random eigenvalues
!
            CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, 1.0D0, 1.0D0, &
                         'T', 'N', WORK( N+1 ), 1, 1.0D0, &
                         WORK( 2*N+1 ), 1, 1.0D0, 'N', IDUMMA, N, 0, &
                         0.0D0, ANORM, 'NO', A, LDA, IWORK, IINFO )
!
         ELSE
!
            IINFO = 1
         END IF
!
         IF( IINFO /= 0 ) THEN
            WRITE( NOUNIT, FMT = 9993 )'Generator', IINFO, N, JTYPE, &
               IOLDSD
            INFO = ABS( IINFO )
            RETURN
         END IF
!
90       CONTINUE
!
!           Test for minimal and generous workspace
!
         DO IWK = 1, 2
            IF( IWK == 1 ) THEN
               NNWORK = 4*N
            ELSE
               NNWORK = 5*N + 2*N**2
            END IF
            NNWORK = MAX( NNWORK, 1 )
!
!              Initialize RESULT
!
            RESULT(1:7) = -1.0D+0
!
!              Compute eigenvalues and eigenvectors, and test them
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLACPY( 'F', N, N, A, LDA, H, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGEEV( 'V', 'V', N, H, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, NNWORK, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGEEV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9993 )'DGEEV1', IINFO, N, JTYPE, &
                  IOLDSD
               INFO = ABS( IINFO )
               GO TO 220
            END IF
!
!              Do Test (1)
!
            CALL DGET22( 'N', 'N', 'N', N, A, LDA, VR, LDVR, WR, WI, &
                         WORK, RES )
            RESULT( 1 ) = RES( 1 )
!
!              Do Test (2)
!
            CALL DGET22( 'T', 'N', 'T', N, A, LDA, VL, LDVL, WR, WI, &
                         WORK, RES )
            RESULT( 2 ) = RES( 1 )
!
!              Do Test (3)
!
            DO J = 1, N
               TNRM = 1.0D0
               IF( WI( J ) == 0.0D0 ) THEN
                  TNRM = DNRM2( N, VR( 1, J ), 1 )
               ELSE IF( WI( J ) > 0.0D0 ) THEN
                  TNRM = DLAPY2( DNRM2( N, VR( 1, J ), 1 ), &
                         DNRM2( N, VR( 1, J+1 ), 1 ) )
               END IF
               RESULT( 3 ) = MAX( RESULT( 3 ), &
                             MIN( ULPINV, ABS( TNRM-1.0D0 ) / ULP ) )
               IF( WI( J ) > 0.0D0 ) THEN
                  VMX = 0.0D0
                  VRMX = 0.0D0
                  DO JJ = 1, N
                     VTST = DLAPY2( VR( JJ, J ), VR( JJ, J+1 ) )
                     IF( VTST > VMX ) &
                        VMX = VTST
                     IF( VR( JJ, J+1 ) == 0.0D0 .AND. &
                         ABS( VR( JJ, J ) ) > VRMX ) &
                         VRMX = ABS( VR( JJ, J ) )
                     ENDDO
                  IF( VRMX / VMX < 1.0D0-2.0D0*ULP ) &
                     RESULT( 3 ) = ULPINV
               END IF
               ENDDO
!
!              Do Test (4)
!
            DO J = 1, N
               TNRM = 1.0D0
               IF( WI( J ) == 0.0D0 ) THEN
                  TNRM = DNRM2( N, VL( 1, J ), 1 )
               ELSE IF( WI( J ) > 0.0D0 ) THEN
                  TNRM = DLAPY2( DNRM2( N, VL( 1, J ), 1 ), &
                         DNRM2( N, VL( 1, J+1 ), 1 ) )
               END IF
               RESULT( 4 ) = MAX( RESULT( 4 ), &
                             MIN( ULPINV, ABS( TNRM-1.0D0 ) / ULP ) )
               IF( WI( J ) > 0.0D0 ) THEN
                  VMX = 0.0D0
                  VRMX = 0.0D0
                  DO JJ = 1, N
                     VTST = DLAPY2( VL( JJ, J ), VL( JJ, J+1 ) )
                     IF( VTST > VMX ) &
                        VMX = VTST
                     IF( VL( JJ, J+1 ) == 0.0D0 .AND. &
                         ABS( VL( JJ, J ) ) > VRMX ) &
                         VRMX = ABS( VL( JJ, J ) )
                     ENDDO
                  IF( VRMX / VMX < 1.0D0-2.0D0*ULP ) &
                     RESULT( 4 ) = ULPINV
               END IF
               ENDDO
!
!              Compute eigenvalues only, and test them
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLACPY( 'F', N, N, A, LDA, H, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGEEV( 'N', 'N', N, H, LDA, WR1, WI1, DUM, 1, DUM, &
                        1, WORK, NNWORK, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGEEV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9993 )'DGEEV2', IINFO, N, JTYPE, &
                  IOLDSD
               INFO = ABS( IINFO )
               GO TO 220
            END IF
!
!              Do Test (5)
!
            IF (ANY(WR(1:N) /= WR1(1:N) .or. WI(1:N) /= WI1(1:N))) RESULT( 5 ) = ULPINV
!
!              Compute eigenvalues and right eigenvectors, and test them
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLACPY( 'F', N, N, A, LDA, H, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGEEV( 'N', 'V', N, H, LDA, WR1, WI1, DUM, 1, LRE, &
                        LDLRE, WORK, NNWORK, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGEEV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9993 )'DGEEV3', IINFO, N, JTYPE, &
                  IOLDSD
               INFO = ABS( IINFO )
               GO TO 220
            END IF
!
!              Do Test (5) again
!
            IF (ANY(WR(1:N) /= WR1(1:N) .or. WI(1:N) /= WI1(1:N))) RESULT( 5 ) = ULPINV
!
!              Do Test (6)
!
            IF (ANY(VR(1:N,1:N) /= LRE(1:N,1:N))) RESULT( 6 ) = ULPINV
!
!              Compute eigenvalues and left eigenvectors, and test them
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLACPY( 'F', N, N, A, LDA, H, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGEEV( 'V', 'N', N, H, LDA, WR1, WI1, LRE, LDLRE, &
                        DUM, 1, WORK, NNWORK, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGEEV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9993 )'DGEEV4', IINFO, N, JTYPE, &
                  IOLDSD
               INFO = ABS( IINFO )
               GO TO 220
            END IF
!
!              Do Test (5) again
!
            IF (ANY(WR(1:N) /= WR1(1:N) .or. WI(1:N) /= WI1(1:N))) RESULT( 5 ) = ULPINV
!
!              Do Test (7)
!
            IF (ANY(VL(1:N,1:N) /= LRE(1:N,1:N))) RESULT( 7 ) = ULPINV
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
  220          CONTINUE
!
            NTEST = COUNT(RESULT(1:7) >= 0.0E+0)
            NFAIL = COUNT(RESULT(1:7) >= THRESH)
!
            IF( NFAIL > 0 ) NTESTF = NTESTF + 1
            IF( NTESTF == 1 ) THEN
               WRITE( NOUNIT, FMT = 9999 )PATH
               WRITE( NOUNIT, FMT = 9998 )
               WRITE( NOUNIT, FMT = 9997 )
               WRITE( NOUNIT, FMT = 9996 )
               WRITE( NOUNIT, FMT = 9995 )THRESH
               NTESTF = 2
            END IF
!
            DO J = 1, 7
               IF( RESULT( J ) >= THRESH ) THEN
                  WRITE( NOUNIT, FMT = 9994 )N, IWK, IOLDSD, JTYPE, J, RESULT( J )
               END IF
               ENDDO
!
            NERRS = NERRS + NFAIL
            NTESTT = NTESTT + NTEST
!
            ENDDO
  260    CONTINUE
         ENDDO
      ENDDO
!
!     Summary
!
   CALL DLASUM( PATH, NOUNIT, NERRS, NTESTT )
!
 9999 FORMAT( / 1X, A3, ' -- Real Eigenvalue-Eigenvector Decomposition', &
         ' Driver', / ' Matrix types (see DDRVEV for details): ' )
!
 9998 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ', &
         '           ', '  5=Diagonal: geometr. spaced entries.', &
         / '  2=Identity matrix.                    ', '  6=Diagona', &
         'l: clustered entries.', / '  3=Transposed Jordan block.  ', &
         '          ', '  7=Diagonal: large, evenly spaced.', / '  ', &
         '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s', &
         'mall, evenly spaced.' )
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev', &
         'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e', &
         'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ', &
         ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond', &
         'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp', &
         'lex ', / ' 12=Well-cond., random complex ', 6X, '   ', &
         ' 17=Ill-cond., large rand. complx ', / ' 13=Ill-condi', &
         'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.', &
         ' complx ' )
 9996 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ', &
         'with small random entries.', / ' 20=Matrix with large ran', &
         'dom entries.   ', / )
 9995 FORMAT( ' Tests performed with test threshold =', F8.2, &
         / / ' 1 = | A VR - VR W | / ( n |A| ulp ) ', &
         / ' 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) ', &
         / ' 3 = | |VR(i)| - 1 | / ulp ', &
         / ' 4 = | |VL(i)| - 1 | / ulp ', &
         / ' 5 = 0 if W same no matter if VR or VL computed,', &
         ' 1/ulp otherwise', / &
         ' 6 = 0 if VR same no matter if VL computed,', &
         '  1/ulp otherwise', / &
         ' 7 = 0 if VL same no matter if VR computed,', &
         '  1/ulp otherwise', / )
 9994 FORMAT( ' N=', I5, ', IWK=', I2, ', seed=', 4( I4, ',' ), &
         ' type ', I2, ', test(', I2, ')=', G10.3 )
 9993 FORMAT( ' DDRVEV: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', &
         I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
!
   RETURN
!
!     End of DDRVEV
!
END




