!> \brief \b ZCHKHB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKHB( NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED,
!                          THRESH, NOUNIT, A, LDA, SD, SE, U, LDU, WORK,
!                          LWORK, RWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES,
!      $                   NWDTHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), KK( * ), NN( * )
!       DOUBLE PRECISION   RESULT( * ), RWORK( * ), SD( * ), SE( * )
!       COMPLEX*16         A( LDA, * ), U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKHB tests the reduction of a Hermitian band matrix to tridiagonal
!> from, used with the Hermitian eigenvalue problem.
!>
!> ZHBTRD factors a Hermitian band matrix A as  U S U* , where * means
!> conjugate transpose, S is symmetric tridiagonal, and U is unitary.
!> ZHBTRD can use either just the lower or just the upper triangle
!> of A; ZCHKHB checks both cases.
!>
!> When ZCHKHB is called, a number of matrix "sizes" ("n's"), a number
!> of bandwidths ("k's"), and a number of matrix "types" are
!> specified.  For each size ("n"), each bandwidth ("k") less than or
!> equal to "n", and each type of matrix, one matrix will be generated
!> and used to test the hermitian banded reduction routine.  For each
!> matrix, a number of tests will be performed:
!>
!> (1)     | A - V S V* | / ( |A| n ulp )  computed by ZHBTRD with
!>                                         UPLO='U'
!>
!> (2)     | I - UU* | / ( n ulp )
!>
!> (3)     | A - V S V* | / ( |A| n ulp )  computed by ZHBTRD with
!>                                         UPLO='L'
!>
!> (4)     | I - UU* | / ( n ulp )
!>
!> The "sizes" are specified by an array NN(1:NSIZES); the value of
!> each element NN(j) specifies one size.
!> The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!> if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!> Currently, the list of possible types is:
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!>
!> (3)  A diagonal matrix with evenly spaced entries
!>      1, ..., ULP  and random signs.
!>      (ULP = (first number larger than 1) - 1 )
!> (4)  A diagonal matrix with geometrically spaced entries
!>      1, ..., ULP  and random signs.
!> (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>      and random signs.
!>
!> (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!> (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!> (8)  A matrix of the form  U* D U, where U is unitary and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!>
!> (9)  A matrix of the form  U* D U, where U is unitary and
!>      D has geometrically spaced entries 1, ..., ULP with random
!>      signs on the diagonal.
!>
!> (10) A matrix of the form  U* D U, where U is unitary and
!>      D has "clustered" entries 1, ULP,..., ULP with random
!>      signs on the diagonal.
!>
!> (11) Same as (8), but multiplied by SQRT( overflow threshold )
!> (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!> (13) Hermitian matrix with random entries chosen from (-1,1).
!> (14) Same as (13), but multiplied by SQRT( overflow threshold )
!> (15) Same as (13), but multiplied by SQRT( underflow threshold )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          ZCHKHB does nothing.  It must be at least zero.
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
!> \param[in] NWDTHS
!> \verbatim
!>          NWDTHS is INTEGER
!>          The number of bandwidths to use.  If it is zero,
!>          ZCHKHB does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] KK
!> \verbatim
!>          KK is INTEGER array, dimension (NWDTHS)
!>          An array containing the bandwidths to be used for the band
!>          matrices.  The values must be at least zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, ZCHKHB
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
!>          next call to ZCHKHB to continue the same random number
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
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                            (LDA, max(NN))
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 2 (not 1!)
!>          and at least max( KK )+1.
!> \endverbatim
!>
!> \param[out] SD
!> \verbatim
!>          SD is DOUBLE PRECISION array, dimension (max(NN))
!>          Used to hold the diagonal of the tridiagonal matrix computed
!>          by ZHBTRD.
!> \endverbatim
!>
!> \param[out] SE
!> \verbatim
!>          SE is DOUBLE PRECISION array, dimension (max(NN))
!>          Used to hold the off-diagonal of the tridiagonal matrix
!>          computed by ZHBTRD.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU, max(NN))
!>          Used to hold the unitary matrix computed by ZHBTRD.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          max( LDA+1, max(NN)+1 )*max(NN).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (4)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>
!>-----------------------------------------------------------------------
!>
!>       Some Local Variables and Parameters:
!>       ---- ----- --------- --- ----------
!>       0.0D+0, 1.0D+0       Real 0 and 1.
!>       MAXTYP          The number of types defined.
!>       NTEST           The number of tests performed, or which can
!>                       be performed so far, for the current matrix.
!>       NTESTT          The total number of tests performed so far.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far.
!>       COND, IMODE     Values to be passed to the matrix generators.
!>       ANORM           Norm of A; passed to matrix generators.
!>
!>       OVFL, UNFL      Overflow and underflow thresholds.
!>       ULP, ULPINV     Finest relative precision and its inverse.
!>       RTOVFL, RTUNFL  Square roots of the previous 2 values.
!>               The following four arrays decode JTYPE:
!>       KTYPE(j)        The general type (1-10) for type "j".
!>       KMODE(j)        The MODE value to be passed to the matrix
!>                       generator for type "j".
!>       KMAGN(j)        The order of magnitude ( O(1),
!>                       O(overflow^(1/2) ), O(underflow^(1/2) )
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
!> \ingroup complex16_eig
!
!  =====================================================================
   SUBROUTINE ZCHKHB( NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED, &
                      THRESH, NOUNIT, A, LDA, SD, SE, U, LDU, WORK, &
                      LWORK, RWORK, RESULT, INFO )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES, &
                      NWDTHS
   DOUBLE PRECISION   THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), KK( * ), NN( * )
   DOUBLE PRECISION   RESULT( * ), RWORK( * ), SD( * ), SE( * )
   COMPLEX*16         A( LDA, * ), U( LDU, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..

   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 15 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADNN, BADNNB
   INTEGER            I, IINFO, IMODE, ITYPE, J, JC, JCOL, JR, JSIZE, &
                      JTYPE, JWIDTH, K, KMAX, MTYPES, N, NERRS, &
                      NMATS, NMAX, NTEST, NTESTT
   DOUBLE PRECISION   ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, &
                      TEMP1, ULP, ULPINV, UNFL
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), &
                      KMODE( MAXTYP ), KTYPE( MAXTYP )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLAMCH
   EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           DLASUM, XERBLA, ZHBT21, ZHBTRD, ZLACPY, ZLASET, &
                      ZLATMR, ZLATMS
!     ..
!     .. Data statements ..
   DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8 /
   DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 /
   DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 /
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
   NTESTT = 0
   INFO = 0
!
!     Important constants
!
   NMAX = MAXVAL(NN(1:NSIZES))
   BADNN = ANY(NN(1:NSIZES) < 0 )
!
   NMAX = MAXVAL(NN(1:NSIZES))
   BADNN = ANY(NN(1:NSIZES) < 0 )

   KMAX = MAXVAL(KK(1:NSIZES))
   BADNNB = ANY(KK(1:NSIZES) < 0 )
   KMAX = MIN( NMAX-1, KMAX )
!
!     Check for errors
!
   IF( NSIZES < 0 ) THEN
      INFO = -1
   ELSE IF( BADNN ) THEN
      INFO = -2
   ELSE IF( NWDTHS < 0 ) THEN
      INFO = -3
   ELSE IF( BADNNB ) THEN
      INFO = -4
   ELSE IF( NTYPES < 0 ) THEN
      INFO = -5
   ELSE IF( LDA < KMAX+1 ) THEN
      INFO = -11
   ELSE IF( LDU < NMAX ) THEN
      INFO = -15
   ELSE IF( ( MAX( LDA, NMAX )+1 )*NMAX > LWORK ) THEN
      INFO = -17
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'ZCHKHB', -INFO )
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
!     Quick return if possible
!
   IF( NSIZES == 0 .OR. NTYPES == 0 .OR. NWDTHS == 0 ) RETURN
!
!     More Important constants
!
   UNFL = DLAMCH( 'Safe minimum' )
   OVFL = 1.0D+0 / UNFL
   ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
   ULPINV = 1.0D+0 / ULP
   RTUNFL = SQRT( UNFL )
   RTOVFL = SQRT( OVFL )
!
!     Loop over sizes, types
!
   NERRS = 0
   NMATS = 0
!
   DO JSIZE = 1, NSIZES
      N = NN( JSIZE )
      ANINV = 1.0D+0 / DBLE( MAX( 1, N ) )
!
      DO JWIDTH = 1, NWDTHS
         K = KK( JWIDTH )
         IF( K > N ) GO TO 180
         K = MAX( 0, MIN( N-1, K ) )
!
         IF( NSIZES /= 1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
!
         DO JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 170
            NMATS = NMATS + 1
            NTEST = 0
!
            IOLDSD(1:4) = ISEED(1:4)
!
!              Compute "A".
!              Store as "Upper"; later, we will copy to other format.
!
!              Control parameters:
!
!                  KMAGN  KMODE        KTYPE
!              =1  O(1)   clustered 1  zero
!              =2  large  clustered 2  identity
!              =3  small  exponential  (none)
!              =4         arithmetic   diagonal, (w/ eigenvalues)
!              =5         random log   hermitian, w/ eigenvalues
!              =6         random       (none)
!              =7                      random diagonal
!              =8                      random hermitian
!              =9                      positive definite
!              =10                     diagonally dominant tridiagonal
!
            IF( MTYPES > MAXTYP ) GO TO 100
!
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
!
!              Compute norm
!
            SELECT CASE (KMAGN(JTYPE))
             CASE (1)
              ANORM = (1.0D+0,0.0D+0)
             CASE (2)
              ANORM = ( RTOVFL*ULP )*ANINV
             CASE (3)
              ANORM = RTUNFL*N*ULPINV
            END SELECT
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLASET( 'Full', LDA, N, (0.0D+0,0.0D+0), (0.0D+0,0.0D+0), A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IINFO = 0
            IF( JTYPE <= 15 ) THEN
               COND = ULPINV
            ELSE
               COND = ULPINV*ANINV / 10.0D+0
            END IF
!
!              Special Matrices -- Identity & Jordan block
!
!                 Zero
!
            IF( ITYPE == 1 ) THEN
               IINFO = 0
!
            ELSE IF( ITYPE == 2 ) THEN
!
!                 Identity
!
               A(K+1,1:N) = ANORM
!
            ELSE IF( ITYPE == 4 ) THEN
!
!                 Diagonal Matrix, [Eigen]values Specified
!
               CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, &
                            COND, ANORM, 0, 0, 'Q', A( K+1, 1 ), LDA, &
                            WORK, IINFO )
!
            ELSE IF( ITYPE == 5 ) THEN
!
!                 Hermitian, eigenvalues specified
!
               CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, &
                            COND, ANORM, K, K, 'Q', A, LDA, WORK, &
                            IINFO )
!
            ELSE IF( ITYPE == 7 ) THEN
!
!                 Diagonal, random eigenvalues
!
               CALL ZLATMR( N, N, 'S', ISEED, 'H', WORK, 6, 1.0D+0, &
                            (1.0D+0,0.0D+0), 'T', 'N', WORK( N+1 ), 1, 1.0D+0, &
                            WORK( 2*N+1 ), 1, 1.0D+0, 'N', IDUMMA, 0, 0, &
                            0.0D+0, ANORM, 'Q', A( K+1, 1 ), LDA, &
                            IDUMMA, IINFO )
!
            ELSE IF( ITYPE == 8 ) THEN
!
!                 Hermitian, random eigenvalues
!
               CALL ZLATMR( N, N, 'S', ISEED, 'H', WORK, 6, 1.0D+0, &
                            (1.0D+0,0.0D+0), 'T', 'N', WORK( N+1 ), 1, 1.0D+0, &
                            WORK( 2*N+1 ), 1, 1.0D+0, 'N', IDUMMA, K, K, &
                            0.0D+0, ANORM, 'Q', A, LDA, IDUMMA, IINFO )
!
            ELSE IF( ITYPE == 9 ) THEN
!
!                 Positive definite, eigenvalues specified.
!
               CALL ZLATMS( N, N, 'S', ISEED, 'P', RWORK, IMODE, &
                            COND, ANORM, K, K, 'Q', A, LDA, &
                            WORK( N+1 ), IINFO )
!
            ELSE IF( ITYPE == 10 ) THEN
!
!                 Positive definite tridiagonal, eigenvalues specified.
!
               IF( N > 1 ) K = MAX( 1, K )
               CALL ZLATMS( N, N, 'S', ISEED, 'P', RWORK, IMODE, &
                            COND, ANORM, 1, 1, 'Q', A( K, 1 ), LDA, &
                            WORK, IINFO )
               DO I = 2, N
                  TEMP1 = ABS(A(K,I))/SQRT(ABS(A(K+1,I-1)*A(K+1,I)))
                  IF(TEMP1>0.5D+0) A(K,I) = 0.5D+0*SQRT(ABS(A(K+1,I-1)*A(K+1,I)))
               ENDDO
!
            ELSE
!
               IINFO = 1
            END IF
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
  100          CONTINUE
!
!              Call ZHBTRD to compute S and U from upper triangle.
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLACPY( ' ', K+1, N, A, LDA, WORK, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            NTEST = 1
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZHBTRD( 'V', 'U', N, K, WORK, LDA, SD, SE, U, LDU, &
                         WORK( LDA*N+1 ), IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZHBTRD(U)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 1 ) = ULPINV
                  GO TO 150
               END IF
            END IF
!
!              Do tests 1 and 2
!
            CALL ZHBT21( 'Upper', N, K, 1, A, LDA, SD, SE, U, LDU, &
                         WORK, RWORK, RESULT( 1 ) )
!
!              Convert A from Upper-Triangle-Only storage to
!              Lower-Triangle-Only storage.
!
            DO JC = 1, N
               DO JR = 0, MIN( K, N-JC )
                  A( JR+1, JC ) = DCONJG( A( K+1-JR, JC+JR ) )
               ENDDO
            ENDDO
            DO JC = N + 1 - K, N
               DO JR = MIN( K, N-JC ) + 1, K
                  A( JR+1, JC ) = 0.0D+0
               ENDDO
            ENDDO
!
!              Call ZHBTRD to compute S and U from lower triangle
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZLACPY( ' ', K+1, N, A, LDA, WORK, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            NTEST = 3
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL ZHBTRD( 'V', 'L', N, K, WORK, LDA, SD, SE, U, LDU, &
                         WORK( LDA*N+1 ), IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZHBTRD : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
            IF( IINFO /= 0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZHBTRD(L)', IINFO, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO < 0 ) THEN
                  RETURN
               ELSE
                  RESULT( 3 ) = ULPINV
                  GO TO 150
               END IF
            END IF
            NTEST = 4
!
!              Do tests 3 and 4
!
            CALL ZHBT21( 'Lower', N, K, 1, A, LDA, SD, SE, U, LDU, &
                         WORK, RWORK, RESULT( 3 ) )
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
  150          CONTINUE
            NTESTT = NTESTT + NTEST
!
!              Print out tests which fail.
!
            DO JR = 1, NTEST
               IF( RESULT( JR ) >= THRESH ) THEN
!
!                    If this is the first test to fail,
!                    print a header to the data file.
!
                  IF( NERRS == 0 ) THEN
                     WRITE( NOUNIT, FMT = 9998 )'ZHB'
                     WRITE( NOUNIT, FMT = 9997 )
                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )'Hermitian'
                     WRITE( NOUNIT, FMT = 9994 )'unitary', '*', &
                        'conjugate transpose', ( '*', J = 1, 4 )
                  END IF
                  NERRS = NERRS + 1
                  WRITE( NOUNIT, FMT = 9993 )N, K, IOLDSD, JTYPE, &
                     JR, RESULT( JR )
               END IF
               ENDDO
!
  170       CONTINUE
            ENDDO
  180    CONTINUE
         ENDDO
      ENDDO
!
!     Summary
!
   CALL DLASUM( 'ZHB', NOUNIT, NERRS, NTESTT )
   RETURN
!
 9999 FORMAT( ' ZCHKHB: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', &
         I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9998 FORMAT( / 1X, A3, &
        ' -- Complex Hermitian Banded Tridiagonal Reduction Routines' &
          )
 9997 FORMAT( ' Matrix types (see DCHK23 for details): ' )
!
 9996 FORMAT( / ' Special Matrices:', &
         / '  1=Zero matrix.                        ', &
         '  5=Diagonal: clustered entries.', &
         / '  2=Identity matrix.                    ', &
         '  6=Diagonal: large, evenly spaced.', &
         / '  3=Diagonal: evenly spaced entries.    ', &
         '  7=Diagonal: small, evenly spaced.', &
         / '  4=Diagonal: geometr. spaced entries.' )
 9995 FORMAT( ' Dense ', A, ' Banded Matrices:', &
         / '  8=Evenly spaced eigenvals.            ', &
         ' 12=Small, evenly spaced eigenvals.', &
         / '  9=Geometrically spaced eigenvals.     ', &
         ' 13=Matrix with random O(1) entries.', &
         / ' 10=Clustered eigenvalues.              ', &
         ' 14=Matrix with large random entries.', &
         / ' 11=Large, evenly spaced eigenvals.     ', &
         ' 15=Matrix with small random entries.' )
!
 9994 FORMAT( / ' Tests performed:   (S is Tridiag,  U is ', A, ',', &
         / 20X, A, ' means ', A, '.', / ' UPLO=''U'':', &
         / '  1= | A - U S U', A1, ' | / ( |A| n ulp )     ', &
         '  2= | I - U U', A1, ' | / ( n ulp )', / ' UPLO=''L'':', &
         / '  3= | A - U S U', A1, ' | / ( |A| n ulp )     ', &
         '  4= | I - U U', A1, ' | / ( n ulp )' )
 9993 FORMAT( ' N=', I5, ', K=', I4, ', seed=', 4( I4, ',' ), ' type ', &
         I2, ', test(', I2, ')=', G10.3 )
!
!     End of ZCHKHB
!
END


