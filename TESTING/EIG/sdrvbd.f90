!> \brief \b SDRVBD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S,
!                          SSAV, E, WORK, LWORK, IWORK, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUT, NSIZES,
!      $                   NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
!       REAL               A( LDA, * ), ASAV( LDA, * ), E( * ), S( * ),
!      $                   SSAV( * ), U( LDU, * ), USAV( LDU, * ),
!      $                   VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SDRVBD checks the singular value decomposition (SVD) drivers
!> SGESVD, SGESDD, SGESVDQ, SGESVJ, SGEJSV, and DGESVDX.
!>
!> Both SGESVD and SGESDD factor A = U diag(S) VT, where U and VT are
!> orthogonal and diag(S) is diagonal with the entries of the array S
!> on its diagonal. The entries of S are the singular values,
!> nonnegative and stored in decreasing order.  U and VT can be
!> optionally not computed, overwritten on A, or computed partially.
!>
!> A is M by N. Let MNMIN = min( M, N ). S has dimension MNMIN.
!> U can be M by M or M by MNMIN. VT can be N by N or MNMIN by N.
!>
!> When SDRVBD is called, a number of matrix "sizes" (M's and N's)
!> and a number of matrix "types" are specified.  For each size (M,N)
!> and each type of matrix, and for the minimal workspace as well as
!> workspace adequate to permit blocking, an  M x N  matrix "A" will be
!> generated and used to test the SVD routines.  For each matrix, A will
!> be factored as A = U diag(S) VT and the following 12 tests computed:
!>
!> Test for SGESVD:
!>
!> (1)    | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (2)    | I - U'U | / ( M ulp )
!>
!> (3)    | I - VT VT' | / ( N ulp )
!>
!> (4)    S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> (5)    | U - Upartial | / ( M ulp ) where Upartial is a partially
!>        computed U.
!>
!> (6)    | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>        computed VT.
!>
!> (7)    | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>        vector of singular values from the partial SVD
!>
!> Test for SGESDD:
!>
!> (8)    | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (9)    | I - U'U | / ( M ulp )
!>
!> (10)   | I - VT VT' | / ( N ulp )
!>
!> (11)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> (12)   | U - Upartial | / ( M ulp ) where Upartial is a partially
!>        computed U.
!>
!> (13)   | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>        computed VT.
!>
!> (14)   | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>        vector of singular values from the partial SVD
!>
!> Test for SGESVDQ:
!>
!> (36)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (37)   | I - U'U | / ( M ulp )
!>
!> (38)   | I - VT VT' | / ( N ulp )
!>
!> (39)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for SGESVJ:
!>
!> (15)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (16)   | I - U'U | / ( M ulp )
!>
!> (17)   | I - VT VT' | / ( N ulp )
!>
!> (18)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for SGEJSV:
!>
!> (19)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (20)   | I - U'U | / ( M ulp )
!>
!> (21)   | I - VT VT' | / ( N ulp )
!>
!> (22)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for SGESVDX( 'V', 'V', 'A' )/SGESVDX( 'N', 'N', 'A' )
!>
!> (23)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (24)   | I - U'U | / ( M ulp )
!>
!> (25)   | I - VT VT' | / ( N ulp )
!>
!> (26)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> (27)   | U - Upartial | / ( M ulp ) where Upartial is a partially
!>        computed U.
!>
!> (28)   | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>        computed VT.
!>
!> (29)   | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>        vector of singular values from the partial SVD
!>
!> Test for SGESVDX( 'V', 'V', 'I' )
!>
!> (30)   | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
!>
!> (31)   | I - U'U | / ( M ulp )
!>
!> (32)   | I - VT VT' | / ( N ulp )
!>
!> Test for SGESVDX( 'V', 'V', 'V' )
!>
!> (33)   | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
!>
!> (34)   | I - U'U | / ( M ulp )
!>
!> (35)   | I - VT VT' | / ( N ulp )
!>
!> The "sizes" are specified by the arrays MM(1:NSIZES) and
!> NN(1:NSIZES); the value of each element pair (MM(j),NN(j))
!> specifies one size.  The "types" are specified by a logical array
!> DOTYPE( 1:NTYPES ); if DOTYPE(j) is .TRUE., then matrix type "j"
!> will be generated.
!> Currently, the list of possible types is:
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!> (3)  A matrix of the form  U D V, where U and V are orthogonal and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!> (4)  Same as (3), but multiplied by the underflow-threshold / ULP.
!> (5)  Same as (3), but multiplied by the overflow-threshold * ULP.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of matrix sizes (M,N) contained in the vectors
!>          MM and NN.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER array, dimension (NSIZES)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, SDRVBD
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrices are in A and B.
!>          This is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size (m,n), a matrix
!>          of type j will be generated.  If NTYPES is smaller than the
!>          maximum number of types defined (PARAMETER MAXTYP), then
!>          types NTYPES+1 through MAXTYP will not be generated.  If
!>          NTYPES is larger than MAXTYP, DOTYPE(MAXTYP+1) through
!>          DOTYPE(NTYPES) will be ignored.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed of the random number generator.  The array
!>          elements should be between 0 and 4095; if not they will be
!>          reduced mod 4096.  Also, ISEED(4) must be odd.
!>          On exit, ISEED is changed and can be used in the next call to
!>          SDRVBD to continue the same random number sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  The test
!>          ratios are scaled to be O(1), so THRESH should be a small
!>          multiple of 1, e.g., 10 or 100.  To have every test ratio
!>          printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,NMAX)
!>          where NMAX is the maximum value of N in NN.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,MMAX),
!>          where MMAX is the maximum value of M in MM.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array, dimension (LDU,MMAX)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,MMAX).
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT,NMAX)
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is REAL array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[out] USAV
!> \verbatim
!>          USAV is REAL array, dimension (LDU,MMAX)
!> \endverbatim
!>
!> \param[out] VTSAV
!> \verbatim
!>          VTSAV is REAL array, dimension (LDVT,NMAX)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension
!>                      (max(min(MM,NN)))
!> \endverbatim
!>
!> \param[out] SSAV
!> \verbatim
!>          SSAV is REAL array, dimension
!>                      (max(min(MM,NN)))
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension
!>                      (max(min(MM,NN)))
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          max(3*MN+MX,5*MN-4)+2*MN**2 for all pairs
!>          pairs  (MN,MX)=( min(MM(j),NN(j), max(MM(j),NN(j)) )
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension at least 8*min(M,N)
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some MM(j) < 0
!>           -3: Some NN(j) < 0
!>           -4: NTYPES < 0
!>           -7: THRESH < 0
!>          -10: LDA < 1 or LDA < MMAX, where MMAX is max( MM(j) ).
!>          -12: LDU < 1 or LDU < MMAX.
!>          -14: LDVT < 1 or LDVT < NMAX, where NMAX is max( NN(j) ).
!>          -21: LWORK too small.
!>          If  SLATMS, or SGESVD returns an error code, the
!>              absolute value of it is returned.
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
!> \ingroup single_eig
!
!  =====================================================================
   SUBROUTINE SDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, &
                      A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, &
                      SSAV, E, WORK, LWORK, IWORK, NOUT, INFO )
!
   IMPLICIT NONE
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUT, NSIZES, &
                      NTYPES
   REAL               THRESH
!     ..
!     .. Array Arguments ..
   LOGICAL            DOTYPE( * )
   INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
   REAL               A( LDA, * ), ASAV( LDA, * ), E( * ), S( * ), &
                      SSAV( * ), U( LDU, * ), USAV( LDU, * ), &
                      VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXTYP
   PARAMETER          ( MAXTYP = 5 )
!     ..
!     .. Local Scalars ..
   LOGICAL            BADMM, BADNN
   CHARACTER          JOBQ, JOBU, JOBVT, RANGE
   CHARACTER*3        PATH
   INTEGER            I, IINFO, IJQ, IJU, IJVT, IL,IU, IWS, IWTMP, &
                      ITEMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, &
                      MMAX, MNMAX, MNMIN, MTYPES, N, NFAIL, &
                      NMAX, NS, NSI, NSV, NTEST
   REAL               ANORM, DIF, DIV, OVFL, RTUNFL, ULP, &
                      ULPINV, UNFL, VL, VU
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Scalars for DGESVDQ ..
   INTEGER            LIWORK, LRWORK, NUMRANK
!     ..
!     .. Local Arrays for DGESVDQ ..
   REAL               RWORK( 2 )
!     ..
!     .. Local Arrays ..
   CHARACTER          CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 )
   INTEGER            IOLDSD( 4 ), ISEED2( 4 )
   REAL               RESULT( 39 )
!     ..
!     .. External Functions ..
   REAL               SLAMCH, SLARND
   EXTERNAL           SLAMCH, SLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALASVM, SBDT01, SGEJSV, SGESDD, SGESVD, &
                      SGESVDQ, SGESVDX, SGESVJ, SLACPY, SLASET, &
                      SLATMS, SORT01, SORT03, XERBLA
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, NUNIT
!     ..
!     .. Common blocks ..
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               CJOB / 'N', 'O', 'S', 'A' /
   DATA               CJOBR / 'A', 'V', 'I' /
   DATA               CJOBV / 'N', 'V' /
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
   INFO = 0
   NMAX = MAXVAL(NN(1:NSIZES))
   BADNN = any(NN(1:NSIZES) < 0 )
   MMAX = MAXVAL(MM(1:NSIZES))
   BADMM = any(MM(1:NSIZES) < 0 )
   MNMAX = 1
   MINWRK = 1
   DO J = 1, NSIZES
      MNMAX = MAX( MNMAX, MIN( MM( J ), NN( J ) ) )
      MINWRK = MAX( MINWRK, MAX( 3*MIN( MM( J ), &
               NN( J ) )+MAX( MM( J ), NN( J ) ), 5*MIN( MM( J ), &
               NN( J )-4 ) )+2*MIN( MM( J ), NN( J ) )**2 )
   ENDDO
!
!     Check for errors
!
   IF( NSIZES < 0 ) THEN
      INFO = -1
   ELSE IF( BADMM ) THEN
      INFO = -2
   ELSE IF( BADNN ) THEN
      INFO = -3
   ELSE IF( NTYPES < 0 ) THEN
      INFO = -4
   ELSE IF( LDA < MAX( 1, MMAX ) ) THEN
      INFO = -10
   ELSE IF( LDU < MAX( 1, MMAX ) ) THEN
      INFO = -12
   ELSE IF( LDVT < MAX( 1, NMAX ) ) THEN
      INFO = -14
   ELSE IF( MINWRK > LWORK ) THEN
      INFO = -21
   END IF
!
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'SDRVBD', -INFO )
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
!     Initialize constants
!
   PATH( 1: 1 ) = 'Single precision'
   PATH( 2: 3 ) = 'BD'
   NFAIL = 0
   NTEST = 0
   UNFL = SLAMCH( 'Safe minimum' )
   OVFL = 1.0E+0 / UNFL
   ULP = SLAMCH( 'Precision' )
   RTUNFL = SQRT( UNFL )
   ULPINV = 1.0E+0 / ULP
   INFOT = 0
!
!     Loop over sizes, types
!
   DO JSIZE = 1, NSIZES
      M = MM( JSIZE )
      N = NN( JSIZE )
      MNMIN = MIN( M, N )
!
      IF( NSIZES /= 1 ) THEN
         MTYPES = MIN( MAXTYP, NTYPES )
      ELSE
         MTYPES = MIN( MAXTYP+1, NTYPES )
      END IF
!
      DO JTYPE = 1, MTYPES
         IF( .NOT.DOTYPE( JTYPE ) ) GO TO 230
!
         IOLDSD(1:4) = ISEED(1:4)
!
!           Compute "A"
!
         IF( MTYPES > MAXTYP ) GO TO 30
!
         IF( JTYPE == 1 ) THEN
!
!              Zero matrix
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLASET( 'Full', M, N, 0.0E+0, 0.0E+0, A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
         ELSE IF( JTYPE == 2 ) THEN
!
!              Identity matrix
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLASET( 'Full', M, N, 0.0E+0, 1.0E+0, A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLASET : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
!
         ELSE
!
!              (Scaled) random matrix
!
            IF( JTYPE == 3 ) ANORM = 1.0E+0
            IF( JTYPE == 4 ) ANORM = UNFL / ULP
            IF( JTYPE == 5 ) ANORM = OVFL*ULP
            CALL SLATMS( M, N, 'U', ISEED, 'N', S, 4, REAL( MNMIN ), &
                         ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO )
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9996 )'Generator', IINFO, M, N, &
                  JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
         END IF
!
30       CONTINUE
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SLACPY( 'F', M, N, A, LDA, ASAV, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Do for minimal and adequate (for blocking) workspace
!
         DO IWS = 1, 4
!
            RESULT(1:32) = -1.0E0
!
!              Test SGESVD: Factorize A
!
            IWTMP = MAX( 3*MIN( M, N )+MAX( M, N ), 5*MIN( M, N ) )
            LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
            LSWORK = MIN( LSWORK, LWORK )
            LSWORK = MAX( LSWORK, 1 )
            IF( IWS == 4 ) LSWORK = LWORK
!
            IF( IWS > 1 ) CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
            SRNAMT = 'SGESVD'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGESVD( 'A', 'A', M, N, A, LDA, SSAV, USAV, LDU, &
                         VTSAV, LDVT, WORK, LSWORK, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGESVD : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, &
                  LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Do tests 1--4
!
            CALL SBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                         VTSAV, LDVT, WORK, RESULT( 1 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL SORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, &
                            RESULT( 2 ) )
               CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, &
                            RESULT( 3 ) )
            END IF
            RESULT( 4 ) = 0.0E+0
            DO I = 1, MNMIN - 1
               IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 4 ) = ULPINV
               IF( SSAV( I ) < 0.0E+0 ) RESULT( 4 ) = ULPINV
            ENDDO
            IF( MNMIN >= 1 ) THEN
               IF( SSAV( MNMIN ) < 0.0E+0 ) RESULT( 4 ) = ULPINV
            END IF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
            RESULT( 5:7 ) = 0.0E0
            DO IJU = 0, 3
               DO IJVT = 0, 3
                  IF( ( IJU == 3 .AND. IJVT == 3 ) .OR. &
                      ( IJU == 1 .AND. IJVT == 1 ) )GO TO 70
                  JOBU = CJOB( IJU+1 )
                  JOBVT = CJOB( IJVT+1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
                  SRNAMT = 'SGESVD'
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
                               VT, LDVT, WORK, LSWORK, IINFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGESVD : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compare U
!
                  DIF = 0.0E+0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJU == 1 ) THEN
                        CALL SORT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, A, LDA, WORK, LWORK, DIF, &
                                     IINFO )
                     ELSE IF( IJU == 2 ) THEN
                        CALL SORT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, U, LDU, WORK, LWORK, DIF, &
                                     IINFO )
                     ELSE IF( IJU == 3 ) THEN
                        CALL SORT03( 'C', M, M, M, MNMIN, USAV, LDU, &
                                     U, LDU, WORK, LWORK, DIF, &
                                     IINFO )
                     END IF
                  END IF
                  RESULT( 5 ) = MAX( RESULT( 5 ), DIF )
!
!                    Compare VT
!
                  DIF = 0.0E+0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJVT == 1 ) THEN
                        CALL SORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, A, LDA, WORK, LWORK, DIF, &
                                     IINFO )
                     ELSE IF( IJVT == 2 ) THEN
                        CALL SORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     DIF, IINFO )
                     ELSE IF( IJVT == 3 ) THEN
                        CALL SORT03( 'R', N, N, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     DIF, IINFO )
                     END IF
                  END IF
                  RESULT( 6 ) = MAX( RESULT( 6 ), DIF )
!
!                    Compare S
!
                  DIF = 0.0E+0
                  DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                  DO I = 1, MNMIN - 1
                     IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV
                     IF( SSAV( I ) < 0.0E+0 ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                  ENDDO
                  RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
70             CONTINUE
               ENDDO
            ENDDO
!
!              Test SGESDD: Factorize A
!
            IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
            LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
            LSWORK = MIN( LSWORK, LWORK )
            LSWORK = MAX( LSWORK, 1 )
            IF( IWS == 4 ) &
               LSWORK = LWORK
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            SRNAMT = 'SGESDD'
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGESDD( 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, &
                         LDVT, WORK, LSWORK, IWORK, IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGESDD : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, &
                  LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Do tests 8--11
!
            CALL SBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                         VTSAV, LDVT, WORK, RESULT( 8 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL SORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, &
                            RESULT( 9 ) )
               CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, &
                            RESULT( 10 ) )
            END IF
            RESULT( 11 ) = 0.0E+0
            DO I = 1, MNMIN - 1
               IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 11 ) = ULPINV
               IF( SSAV( I ) < 0.0E+0 ) RESULT( 11 ) = ULPINV
            ENDDO
            IF( MNMIN >= 1 ) THEN
               IF( SSAV( MNMIN ) < 0.0E+0 ) RESULT( 11 ) = ULPINV
            END IF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
            RESULT(12:14) = 0.0E+0
            DO IJQ = 0, 2
               JOBQ = CJOB( IJQ+1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'SGESDD'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGESDD( JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, &
                            WORK, LSWORK, IWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGESDD : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 Compare U
!
               DIF = 0.0E+0
               IF( M > 0 .AND. N > 0 ) THEN
                  IF( IJQ == 1 ) THEN
                     IF( M >= N ) THEN
                        CALL SORT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, A, LDA, WORK, LWORK, DIF, &
                                     INFO )
                     ELSE
                        CALL SORT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, U, LDU, WORK, LWORK, DIF, &
                                     INFO )
                     END IF
                  ELSE IF( IJQ == 2 ) THEN
                     CALL SORT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, &
                                  U, LDU, WORK, LWORK, DIF, INFO )
                  END IF
               END IF
               RESULT( 12 ) = MAX( RESULT( 12 ), DIF )
!
!                 Compare VT
!
               DIF = 0.0E+0
               IF( M > 0 .AND. N > 0 ) THEN
                  IF( IJQ == 1 ) THEN
                     IF( M >= N ) THEN
                        CALL SORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     DIF, INFO )
                     ELSE
                        CALL SORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, A, LDA, WORK, LWORK, DIF, &
                                     INFO )
                     END IF
                  ELSE IF( IJQ == 2 ) THEN
                     CALL SORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                  LDVT, VT, LDVT, WORK, LWORK, DIF, &
                                  INFO )
                  END IF
               END IF
               RESULT( 13 ) = MAX( RESULT( 13 ), DIF )
!
!                 Compare S
!
               DIF = 0.0E+0
               DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) &
                     DIF = ULPINV
                  IF( SSAV( I ) < 0.0E+0 ) &
                     DIF = ULPINV
                  DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                  ENDDO
               RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
               ENDDO
!
!              Test SGESVDQ
!              Note: SGESVDQ only works for M >= N
!
            RESULT( 36:39 ) = 0.0E0
!
            IF( M >= N ) THEN
               IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWS == 4 ) LSWORK = LWORK
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'SGESVDQ'
!
               LRWORK = 2
               LIWORK = MAX( N, 1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGESVDQ( 'H', 'N', 'N', 'A', 'A', &
                             M, N, A, LDA, SSAV, USAV, LDU, &
                             VTSAV, LDVT, NUMRANK, IWORK, LIWORK, &
                             WORK, LWORK, RWORK, LRWORK, IINFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGESVDQ : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'SGESVDQ', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
!
!                 Do tests 36--39
!
               CALL SBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                            VTSAV, LDVT, WORK, RESULT( 36 ) )
               IF( M /= 0 .AND. N /= 0 ) THEN
                  CALL SORT01( 'Columns', M, M, USAV, LDU, WORK, &
                               LWORK, RESULT( 37 ) )
                  CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, &
                               LWORK, RESULT( 38 ) )
               END IF
               RESULT( 39 ) = 0.0E+0
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 39 ) = ULPINV
                  IF( SSAV( I ) < 0.0E+0 ) RESULT( 39 ) = ULPINV
                  ENDDO
               IF( MNMIN >= 1 ) THEN
                  IF( SSAV( MNMIN ) < 0.0E+0 ) RESULT( 39 ) = ULPINV
               END IF
            END IF
!
!              Test SGESVJ
!              Note: SGESVJ only works for M >= N
!
            RESULT( 15:18 ) = 0.0E0
!
            IF( M >= N ) THEN
               IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWS == 4 ) &
                  LSWORK = LWORK
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( 'F', M, N, ASAV, LDA, USAV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'SGESVJ'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGESVJ( 'G', 'U', 'V', M, N, USAV, LDA, SSAV, &
                           0, A, LDVT, WORK, LWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGESVJ : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 SGESVJ returns V not VT
!
               DO J=1,N
                  VTSAV(J,1:N) = A(1:N,J)
               END DO
!
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'GESVJ', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
!
!                 Do tests 15--18
!
               CALL SBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                            VTSAV, LDVT, WORK, RESULT( 15 ) )
               IF( M /= 0 .AND. N /= 0 ) THEN
                  CALL SORT01( 'Columns', M, M, USAV, LDU, WORK, &
                               LWORK, RESULT( 16 ) )
                  CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, &
                               LWORK, RESULT( 17 ) )
               END IF
               RESULT( 18 ) = 0.0E+0
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) &
                     RESULT( 18 ) = ULPINV
                  IF( SSAV( I ) < 0.0E+0 ) &
                     RESULT( 18 ) = ULPINV
                  ENDDO
               IF( MNMIN >= 1 ) THEN
                  IF( SSAV( MNMIN ) < 0.0E+0 ) &
                     RESULT( 18 ) = ULPINV
               END IF
            END IF
!
!              Test SGEJSV
!              Note: SGEJSV only works for M >= N
!
            RESULT( 19:22 ) = 0.0E0
            IF( M >= N ) THEN
               IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWS == 4 ) &
                  LSWORK = LWORK
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SLACPY( 'F', M, N, ASAV, LDA, VTSAV, LDA )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               SRNAMT = 'SGEJSV'
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGEJSV( 'G', 'U', 'V', 'R', 'N', 'N', &
                      M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, &
                      WORK, LWORK, IWORK, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGEJSV : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
!                 SGEJSV returns V not VT
!
               DO J=1,N
                  VTSAV(J,1:N) = A(1:N,J)
               END DO
!
               IF( IINFO /= 0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'GEJSV', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
!
!                 Do tests 19--22
!
               CALL SBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                            VTSAV, LDVT, WORK, RESULT( 19 ) )
               IF( M /= 0 .AND. N /= 0 ) THEN
                  CALL SORT01( 'Columns', M, M, USAV, LDU, WORK, &
                               LWORK, RESULT( 20 ) )
                  CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, &
                               LWORK, RESULT( 21 ) )
               END IF
               RESULT( 22 ) = 0.0E+0
               DO I = 1, MNMIN - 1
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 22 ) = ULPINV
                  IF( SSAV( I ) < 0.0E+0 ) RESULT( 22 ) = ULPINV
                  ENDDO
               IF( MNMIN >= 1 ) THEN
                  IF( SSAV( MNMIN ) < 0.0E+0 ) RESULT( 22 ) = ULPINV
               END IF
            END IF
!
!              Test SGESVDX
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGESVDX( 'V', 'V', 'A', M, N, A, LDA, &
                          VL, VU, IL, IU, NS, SSAV, USAV, LDU, &
                          VTSAV, LDVT, WORK, LWORK, IWORK, &
                          IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGESVDX : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
!              Do tests 23--29
!
            RESULT( 23:25 ) = 0.0E0
            CALL SBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, &
                         VTSAV, LDVT, WORK, RESULT( 23 ) )
            IF( M /= 0 .AND. N /= 0 ) THEN
               CALL SORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, &
                            RESULT( 24 ) )
               CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, &
                            RESULT( 25 ) )
            END IF
            RESULT( 26 ) = 0.0E+0
            DO I = 1, MNMIN - 1
               IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 26 ) = ULPINV
               IF( SSAV( I ) < 0.0E+0 ) RESULT( 26 ) = ULPINV
               ENDDO
            IF( MNMIN >= 1 ) THEN
               IF( SSAV( MNMIN ) < 0.0E+0 ) RESULT( 26 ) = ULPINV
            END IF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
            RESULT( 27:29 ) = 0.0E0
            DO IJU = 0, 1
               DO IJVT = 0, 1
                  IF( ( IJU == 0 .AND. IJVT == 0 ) .OR. &
                      ( IJU == 1 .AND. IJVT == 1 ) )GO TO 170
                  JOBU = CJOBV( IJU+1 )
                  JOBVT = CJOBV( IJVT+1 )
                  RANGE = CJOBR( 1 )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
                  CALL SGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, &
                                VL, VU, IL, IU, NS, S, U, LDU, &
                                VT, LDVT, WORK, LWORK, IWORK, &
                                IINFO )
#ifdef _TIMER
                  call system_clock(count_rate=nb_periods_sec,count=S2_time)
                  open(file='results.out', unit=10, position = 'append')
                  write(10,'(A,F16.10,A)') 'Total time : SGESVDX : ',&
                        real(S2_time-S1_time)/real(nb_periods_sec), ' s'
                  close(10)
#endif
!
!                    Compare U
!
                  DIF = 0.0E+0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJU == 1 ) THEN
                        CALL SORT03( 'C', M, MNMIN, M, MNMIN, USAV, &
                                     LDU, U, LDU, WORK, LWORK, DIF, &
                                     IINFO )
                     END IF
                  END IF
                  RESULT( 27 ) = MAX( RESULT( 27 ), DIF )
!
!                    Compare VT
!
                  DIF = 0.0E+0
                  IF( M > 0 .AND. N > 0 ) THEN
                     IF( IJVT == 1 ) THEN
                        CALL SORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, &
                                     LDVT, VT, LDVT, WORK, LWORK, &
                                     DIF, IINFO )
                     END IF
                  END IF
                  RESULT( 28 ) = MAX( RESULT( 28 ), DIF )
!
!                    Compare S
!
                  DIF = 0.0E+0
                  DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                  DO I = 1, MNMIN - 1
                     IF( SSAV( I ) < SSAV( I+1 ) ) &
                        DIF = ULPINV
                     IF( SSAV( I ) < 0.0E+0 ) &
                        DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                     ENDDO
                  RESULT( 29 ) = MAX( RESULT( 29 ), DIF )
  170             CONTINUE
                  ENDDO
               ENDDO
!
!              Do tests 30--32: SGESVDX( 'V', 'V', 'I' )
!
            ISEED2(1:4) = ISEED(1:4)
            IF( MNMIN <= 1 ) THEN
               IL = 1
               IU = MAX( 1, MNMIN )
            ELSE
               IL = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) )
               IF( IU < IL ) THEN
                  ITEMP = IU
                  IU = IL
                  IL = ITEMP
               END IF
            END IF
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGESVDX( 'V', 'V', 'I', M, N, A, LDA, &
                          VL, VU, IL, IU, NSI, S, U, LDU, &
                          VT, LDVT, WORK, LWORK, IWORK, &
                          IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGESVDX : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
            RESULT( 30:32 ) = 0.0E0
            CALL SBDT05( M, N, ASAV, LDA, S, NSI, U, LDU, &
                         VT, LDVT, WORK, RESULT( 30 ) )
            CALL SORT01( 'Columns', M, NSI, U, LDU, WORK, LWORK, &
                         RESULT( 31 ) )
            CALL SORT01( 'Rows', NSI, N, VT, LDVT, WORK, LWORK, &
                         RESULT( 32 ) )
!
!              Do tests 33--35: SGESVDX( 'V', 'V', 'V' )
!
            IF( MNMIN > 0 .AND. NSI > 1 ) THEN
               IF( IL /= 1 ) THEN
                  VU = SSAV( IL ) + &
                       MAX( 0.5E+0*ABS( SSAV( IL )-SSAV( IL-1 ) ), &
                       ULP*ANORM, 2.0E+0*RTUNFL )
               ELSE
                  VU = SSAV( 1 ) + &
                       MAX( 0.5E+0*ABS( SSAV( NS )-SSAV( 1 ) ), &
                       ULP*ANORM, 2.0E+0*RTUNFL )
               END IF
               IF( IU /= NS ) THEN
                  VL = SSAV( IU ) - MAX( ULP*ANORM, 2.0E+0*RTUNFL, &
                       0.5E+0*ABS( SSAV( IU+1 )-SSAV( IU ) ) )
               ELSE
                  VL = SSAV( NS ) - MAX( ULP*ANORM, 2.0E+0*RTUNFL, &
                       0.5E+0*ABS( SSAV( NS )-SSAV( 1 ) ) )
               END IF
               VL = MAX( VL,0.0E+0 )
               VU = MAX( VU,0.0E+0 )
               IF( VL >= VU ) VU = MAX( VU*2, VU+VL+0.5E+0 )
            ELSE
               VL = 0.0E+0
               VU = 1.0E+0
            END IF
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SLACPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGESVDX( 'V', 'V', 'V', M, N, A, LDA, &
                          VL, VU, IL, IU, NSV, S, U, LDU, &
                          VT, LDVT, WORK, LWORK, IWORK, &
                          IINFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGESVDX : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( IINFO /= 0 ) THEN
               WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, &
                  JTYPE, LSWORK, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
!
            RESULT( 33:35 ) = 0.0E0
            CALL SBDT05( M, N, ASAV, LDA, S, NSV, U, LDU, &
                         VT, LDVT, WORK, RESULT( 33 ) )
            CALL SORT01( 'Columns', M, NSV, U, LDU, WORK, LWORK, &
                         RESULT( 34 ) )
            CALL SORT01( 'Rows', NSV, N, VT, LDVT, WORK, LWORK, &
                         RESULT( 35 ) )
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
            DO J = 1, 39
               IF( RESULT( J ) >= THRESH ) THEN
                  IF( NFAIL == 0 ) THEN
                     WRITE( NOUT, FMT = 9999 )
                     WRITE( NOUT, FMT = 9998 )
                  END IF
                  WRITE( NOUT, FMT = 9997 )M, N, JTYPE, IWS, IOLDSD, &
                     J, RESULT( J )
                  NFAIL = NFAIL + 1
               END IF
               ENDDO
            NTEST = NTEST + 39
            ENDDO
  230    CONTINUE
         ENDDO
      ENDDO
!
!     Summary
!
   CALL ALASVM( PATH, NOUT, NFAIL, NTEST, 0 )
!
 9999 FORMAT( ' SVD -- Real Singular Value Decomposition Driver ', &
         / ' Matrix types (see SDRVBD for details):', &
         / / ' 1 = Zero matrix', / ' 2 = Identity matrix', &
         / ' 3 = Evenly spaced singular values near 1', &
         / ' 4 = Evenly spaced singular values near underflow', &
         / ' 5 = Evenly spaced singular values near overflow', / / &
         ' Tests performed: ( A is dense, U and V are orthogonal,', &
         / 19X, ' S is an array, and Upartial, VTpartial, and', &
         / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / ' 2 = | I - U**T U | / ( M ulp ) ', &
         / ' 3 = | I - VT VT**T | / ( N ulp ) ', &
         / ' 4 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / ' 5 = | U - Upartial | / ( M ulp )', &
         / ' 6 = | VT - VTpartial | / ( N ulp )', &
         / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )', &
         / ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / ' 9 = | I - U**T U | / ( M ulp ) ', &
         / '10 = | I - VT VT**T | / ( N ulp ) ', &
         / '11 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / '12 = | U - Upartial | / ( M ulp )', &
         / '13 = | VT - VTpartial | / ( N ulp )', &
         / '14 = | S - Spartial | / ( min(M,N) ulp |S| )', &
         / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / '16 = | I - U**T U | / ( M ulp ) ', &
         / '17 = | I - VT VT**T | / ( N ulp ) ', &
         / '18 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / '19 = | U - Upartial | / ( M ulp )', &
         / '20 = | VT - VTpartial | / ( N ulp )', &
         / '21 = | S - Spartial | / ( min(M,N) ulp |S| )', &
         / '22 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         ' SGESVDX(V,V,A) ', &
         / '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),' &
         / '24 = | I - U**T U | / ( M ulp ) ', &
         / '25 = | I - VT VT**T | / ( N ulp ) ', &
         / '26 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / '27 = | U - Upartial | / ( M ulp )', &
         / '28 = | VT - VTpartial | / ( N ulp )', &
         / '29 = | S - Spartial | / ( min(M,N) ulp |S| )', &
         / '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),', &
         ' SGESVDX(V,V,I) ', &
         / '31 = | I - U**T U | / ( M ulp ) ', &
         / '32 = | I - VT VT**T | / ( N ulp ) ', &
         / '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),', &
         ' SGESVDX(V,V,V) ', &
         / '34 = | I - U**T U | / ( M ulp ) ', &
         / '35 = | I - VT VT**T | / ( N ulp ) ', &
         ' SGESVDQ(H,N,N,A,A', &
         / '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', &
         / '37 = | I - U**T U | / ( M ulp ) ', &
         / '38 = | I - VT VT**T | / ( N ulp ) ', &
         / '39 = 0 if S contains min(M,N) nonnegative values in', &
         ' decreasing order, else 1/ulp', &
         / / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1, &
         ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' SDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', &
         I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), &
         I5, ')' )
 9995 FORMAT( ' SDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', &
         I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X, &
         'ISEED=(', 3( I5, ',' ), I5, ')' )
!
   RETURN
!
!     End of SDRVBD
!
END



