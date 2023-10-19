!> \brief \b DQRT15
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQRT15( SCALE, RKSEL, M, N, NRHS, A, LDA, B, LDB, S,
!                          RANK, NORMA, NORMB, ISEED, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE
!       DOUBLE PRECISION   NORMA, NORMB
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQRT15 generates a matrix with full or deficient rank and of various
!> norms.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SCALE
!> \verbatim
!>          SCALE is INTEGER
!>          SCALE = 1: normally scaled matrix
!>          SCALE = 2: matrix scaled up
!>          SCALE = 3: matrix scaled down
!> \endverbatim
!>
!> \param[in] RKSEL
!> \verbatim
!>          RKSEL is INTEGER
!>          RKSEL = 1: full rank matrix
!>          RKSEL = 2: rank-deficient matrix
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of A.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB, NRHS)
!>          A matrix that is in the range space of matrix A.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension MIN(M,N)
!>          Singular values of A.
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>          number of nonzero singular values of A.
!> \endverbatim
!>
!> \param[out] NORMA
!> \verbatim
!>          NORMA is DOUBLE PRECISION
!>          one-norm of A.
!> \endverbatim
!>
!> \param[out] NORMB
!> \verbatim
!>          NORMB is DOUBLE PRECISION
!>          one-norm of B.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is integer array, dimension (4)
!>          seed for random number generator.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          length of work space required.
!>          LWORK >= MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)
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
!> \ingroup double_lin
!
!  =====================================================================
   SUBROUTINE DQRT15( SCALE, RKSEL, M, N, NRHS, A, LDA, B, LDB, S, &
                      RANK, NORMA, NORMB, ISEED, WORK, LWORK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE
   DOUBLE PRECISION   NORMA, NORMB
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( LWORK )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE, TWO, SVMIN
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                      SVMIN = 0.1D0 )
!     ..
!     .. Local Scalars ..
   INTEGER            INFO, J, MN
   DOUBLE PRECISION   BIGNUM, EPS, SMLNUM, TEMP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   DUMMY( 1 )
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DASUM, DLAMCH, DLANGE, DLARND, DNRM2
   EXTERNAL           DASUM, DLAMCH, DLANGE, DLARND, DNRM2
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMM, DLAORD, DLARF, DLARNV, DLAROR, DLASCL, &
                      DLASET, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
   MN = MIN( M, N )
   IF( LWORK < MAX( M+MN, MN*NRHS, 2*N+M ) ) THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DQRT15', 16 )
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
   SMLNUM = DLAMCH( 'Safe minimum' )
   BIGNUM = ONE / SMLNUM
   EPS = DLAMCH( 'Epsilon' )
   SMLNUM = ( SMLNUM / EPS ) / EPS
   BIGNUM = ONE / SMLNUM
!
!     Determine rank and (unscaled) singular values
!
   IF( RKSEL == 1 ) THEN
      RANK = MN
   ELSE IF( RKSEL == 2 ) THEN
      RANK = ( 3*MN ) / 4
      DO J = RANK + 1, MN
         S( J ) = ZERO
      ENDDO
   ELSE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'DQRT15', 2 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   END IF
!
   IF( RANK > 0 ) THEN
!
!        Nontrivial case
!
      S( 1 ) = ONE
      DO J = 2, RANK
20       CONTINUE
         TEMP = DLARND( 1, ISEED )
         IF( TEMP > SVMIN ) THEN
            S( J ) = ABS( TEMP )
         ELSE
            GO TO 20
         END IF
      ENDDO
      CALL DLAORD( 'Decreasing', RANK, S, 1 )
!
!        Generate 'rank' columns of a random orthogonal matrix in A
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, M, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DSCAL( M, ONE / DNRM2( M, WORK, 1 ), WORK, 1 )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLASET( 'Full', M, RANK, ZERO, ONE, A, LDA )
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
      CALL DLARF( 'Left', M, RANK, WORK, 1, TWO, A, LDA, &
                  WORK( M+1 ) )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARF : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        workspace used: m+mn
!
!        Generate consistent rhs in the range space of A
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLARNV( 2, ISEED, RANK*NRHS, WORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DGEMM( 'No transpose', 'No transpose', M, NRHS, RANK, ONE, &
                  A, LDA, WORK, RANK, ZERO, B, LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        work space used: <= mn *nrhs
!
!        generate (unscaled) matrix A
!
      DO J = 1, RANK
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DSCAL( M, S( J ), A( 1, J ), 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DSCAL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
      IF( RANK < N )  THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLASET( 'Full', M, N-RANK, ZERO, ZERO, A( 1, RANK+1 ), &
                      LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDIF
      CALL DLAROR( 'Right', 'No initialization', M, N, A, LDA, ISEED, &
                   WORK, INFO )
!
   ELSE
!
!        work space used 2*n+m
!
!        Generate null matrix and rhs
!
      DO J = 1, MN
         S( J ) = ZERO
      ENDDO
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLASET( 'Full', M, N, ZERO, ZERO, A, LDA )
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
      CALL DLASET( 'Full', M, NRHS, ZERO, ZERO, B, LDB )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
   END IF
!
!     Scale the matrix
!
   IF( SCALE /= 1 ) THEN
      NORMA = DLANGE( 'Max', M, N, A, LDA, DUMMY )
      IF( NORMA /= ZERO ) THEN
         IF( SCALE == 2 ) THEN
!
!              matrix scaled up
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, M, N, A, &
                         LDA, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, MN, 1, S, &
                         MN, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, M, NRHS, B, &
                         LDB, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ELSE IF( SCALE == 3 ) THEN
!
!              matrix scaled down
!
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, M, N, A, &
                         LDA, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, MN, 1, S, &
                         MN, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, M, NRHS, B, &
                         LDB, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DLASCL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ELSE
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL XERBLA( 'DQRT15', 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            RETURN
         END IF
      END IF
   END IF
!
   NORMA = DASUM( MN, S, 1 )
   NORMB = DLANGE( 'One-norm', M, NRHS, B, LDB, DUMMY )
!
   RETURN
!
!     End of DQRT15
!
END
                                                                                                                                                                                                                                                                                                            




