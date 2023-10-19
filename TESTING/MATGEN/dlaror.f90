!> \brief \b DLAROR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAROR( SIDE, INIT, M, N, A, LDA, ISEED, X, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          INIT, SIDE
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAROR pre- or post-multiplies an M by N matrix A by a random
!> orthogonal matrix U, overwriting A.  A may optionally be initialized
!> to the identity matrix before multiplying by U.  U is generated using
!> the method of G.W. Stewart (SIAM J. Numer. Anal. 17, 1980, 403-409).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          Specifies whether A is multiplied on the left or right by U.
!>          = 'L':         Multiply A on the left (premultiply) by U
!>          = 'R':         Multiply A on the right (postmultiply) by U'
!>          = 'C' or 'T':  Multiply A on the left by U and the right
!>                          by U' (Here, U' means U-transpose.)
!> \endverbatim
!>
!> \param[in] INIT
!> \verbatim
!>          INIT is CHARACTER*1
!>          Specifies whether or not A should be initialized to the
!>          identity matrix.
!>          = 'I':  Initialize A to (a section of) the identity matrix
!>                   before applying U.
!>          = 'N':  No initialization.  Apply U to the input matrix A.
!>
!>          INIT = 'I' may be used to generate square or rectangular
!>          orthogonal matrices:
!>
!>          For M = N and SIDE = 'L' or 'R', the rows will be orthogonal
!>          to each other, as will the columns.
!>
!>          If M < N, SIDE = 'R' produces a dense matrix whose rows are
!>          orthogonal and whose columns are not, while SIDE = 'L'
!>          produces a matrix whose rows are orthogonal, and whose first
!>          M columns are orthogonal, and whose remaining columns are
!>          zero.
!>
!>          If M > N, SIDE = 'L' produces a dense matrix whose columns
!>          are orthogonal and whose rows are not, while SIDE = 'R'
!>          produces a matrix whose columns are orthogonal, and whose
!>          first M rows are orthogonal, and whose remaining rows are
!>          zero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          On entry, the array A.
!>          On exit, overwritten by U A ( if SIDE = 'L' ),
!>           or by A U ( if SIDE = 'R' ),
!>           or by U A U' ( if SIDE = 'C' or 'T').
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
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
!>          next call to DLAROR to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (3*MAX( M, N ))
!>          Workspace of length
!>              2*M + N if SIDE = 'L',
!>              2*N + M if SIDE = 'R',
!>              3*N     if SIDE = 'C' or 'T'.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          An error flag.  It is set to:
!>          = 0:  normal return
!>          < 0:  if INFO = -k, the k-th argument had an illegal value
!>          = 1:  if the random numbers generated by DLARND are bad.
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
   SUBROUTINE DLAROR( SIDE, INIT, M, N, A, LDA, ISEED, X, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          INIT, SIDE
   INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   DOUBLE PRECISION   A( LDA, * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE, TOOSML
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, &
                      TOOSML = 1.0D-20 )
!     ..
!     .. Local Scalars ..
   INTEGER            IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM
   DOUBLE PRECISION   FACTOR, XNORM, XNORMS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLARND, DNRM2
   EXTERNAL           LSAME, DLARND, DNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMV, DGER, DLASET, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, SIGN
!     ..
!     .. Executable Statements ..
!
   INFO = 0
   IF( N == 0 .OR. M == 0 ) &
      RETURN
!
   ITYPE = 0
   IF( LSAME( SIDE, 'L' ) ) THEN
      ITYPE = 1
   ELSE IF( LSAME( SIDE, 'R' ) ) THEN
      ITYPE = 2
   ELSE IF( LSAME( SIDE, 'C' ) .OR. LSAME( SIDE, 'T' ) ) THEN
      ITYPE = 3
   END IF
!
!     Check for argument errors.
!
   IF( ITYPE == 0 ) THEN
      INFO = -1
   ELSE IF( M < 0 ) THEN
      INFO = -3
   ELSE IF( N < 0 .OR. ( ITYPE == 3 .AND. N /= M ) ) THEN
      INFO = -4
   ELSE IF( LDA < M ) THEN
      INFO = -6
   END IF
   IF( INFO /= 0 ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DLAROR', -INFO )
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
   IF( ITYPE == 1 ) THEN
      NXFRM = M
   ELSE
      NXFRM = N
   END IF
!
!     Initialize A to the identity matrix if desired
!
   IF( LSAME( INIT, 'I' ) )  THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLASET( 'Full', M, N, ZERO, ONE, A, LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDIF
!
!     If no rotation possible, multiply by random +/-1
!
!     Compute rotation by computing Householder transformations
!     H(2), H(3), ..., H(nhouse)
!
   DO J = 1, NXFRM
      X( J ) = ZERO
   ENDDO
!
   DO IXFRM = 2, NXFRM
      KBEG = NXFRM - IXFRM + 1
!
!        Generate independent normal( 0, 1 ) random numbers
!
      DO J = KBEG, NXFRM
         X( J ) = DLARND( 3, ISEED )
      ENDDO
!
!        Generate a Householder transformation from the random vector X
!
      XNORM = DNRM2( IXFRM, X( KBEG ), 1 )
      XNORMS = SIGN( XNORM, X( KBEG ) )
      X( KBEG+NXFRM ) = SIGN( ONE, -X( KBEG ) )
      FACTOR = XNORMS*( XNORMS+X( KBEG ) )
      IF( ABS( FACTOR ) < TOOSML ) THEN
         INFO = 1
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL XERBLA( 'DLAROR', INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         RETURN
      ELSE
         FACTOR = ONE / FACTOR
      END IF
      X( KBEG ) = X( KBEG ) + XNORMS
!
!        Apply Householder transformation to A
!
      IF( ITYPE == 1 .OR. ITYPE == 3 ) THEN
!
!           Apply H(k) from the left.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEMV( 'T', IXFRM, N, ONE, A( KBEG, 1 ), LDA, &
                     X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 )
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
         CALL DGER( IXFRM, N, -FACTOR, X( KBEG ), 1, X( 2*NXFRM+1 ), &
                    1, A( KBEG, 1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGER : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
      END IF
!
      IF( ITYPE == 2 .OR. ITYPE == 3 ) THEN
!
!           Apply H(k) from the right.
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DGEMV( 'N', M, IXFRM, ONE, A( 1, KBEG ), LDA, &
                     X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 )
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
         CALL DGER( M, IXFRM, -FACTOR, X( 2*NXFRM+1 ), 1, X( KBEG ), &
                    1, A( 1, KBEG ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DGER : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
      END IF
   ENDDO
!
   X( 2*NXFRM ) = SIGN( ONE, DLARND( 3, ISEED ) )
!
!     Scale the matrix A by D.
!
   IF( ITYPE == 1 .OR. ITYPE == 3 ) THEN
      DO IROW = 1, M
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSCAL( N, X( NXFRM+IROW ), A( IROW, 1 ), LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
!
   IF( ITYPE == 2 .OR. ITYPE == 3 ) THEN
      DO JCOL = 1, N
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSCAL( M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
   END IF
   RETURN
!
!     End of DLAROR
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        


