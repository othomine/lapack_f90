!> \brief \b ZQRT14
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZQRT14( TRANS, M, N, NRHS, A, LDA, X,
!                        LDX, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDA, LDX, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), WORK( LWORK ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZQRT14 checks whether X is in the row space of A or A'.  It does so
!> by scaling both X and A such that their norms are in the range
!> [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
!> (if TRANS = 'C') or an LQ factorization of [A',X]' (if TRANS = 'N'),
!> and returning the norm of the trailing triangle, scaled by
!> MAX(M,N,NRHS)*eps.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, check for X in the row space of A
!>          = 'C':  Conjugate transpose, check for X in row space of A'.
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
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of X.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>          If TRANS = 'N', the N-by-NRHS matrix X.
!>          IF TRANS = 'C', the M-by-NRHS matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          length of workspace array required
!>          If TRANS = 'N', LWORK >= (M+NRHS)*(N+2);
!>          if TRANS = 'C', LWORK >= (N+NRHS)*(M+2).
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
!> \ingroup complex16_lin
!
!  =====================================================================
   DOUBLE PRECISION FUNCTION ZQRT14( TRANS, M, N, NRHS, A, LDA, X, &
                    LDX, WORK, LWORK )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          TRANS
   INTEGER            LDA, LDX, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
   COMPLEX*16         A( LDA, * ), WORK( LWORK ), X( LDX, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            TPSD
   INTEGER            I, INFO, J, LDWORK
   DOUBLE PRECISION   ANRM, ERR, XNRM
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   RWORK( 1 )
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, ZLANGE
   EXTERNAL           LSAME, DLAMCH, ZLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           XERBLA, ZGELQ2, ZGEQR2, ZLACPY, ZLASCL
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, DBLE, DCONJG, MAX, MIN
!     ..
!     .. Executable Statements ..
!
   ZQRT14 = ZERO
   IF( LSAME( TRANS, 'N' ) ) THEN
      LDWORK = M + NRHS
      TPSD = .FALSE.
      IF( LWORK < ( M+NRHS )*( N+2 ) ) THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL XERBLA( 'ZQRT14', 10 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         RETURN
      ELSE IF( N <= 0 .OR. NRHS <= 0 ) THEN
         RETURN
      END IF
   ELSE IF( LSAME( TRANS, 'C' ) ) THEN
      LDWORK = M
      TPSD = .TRUE.
      IF( LWORK < ( N+NRHS )*( M+2 ) ) THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL XERBLA( 'ZQRT14', 10 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : XERBLA : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         RETURN
      ELSE IF( M <= 0 .OR. NRHS <= 0 ) THEN
         RETURN
      END IF
   ELSE
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL XERBLA( 'ZQRT14', 1 )
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
!     Copy and scale A
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLACPY( 'All', M, N, A, LDA, WORK, LDWORK )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   ANRM = ZLANGE( 'M', M, N, WORK, LDWORK, RWORK )
   IF( ANRM /= ZERO )  THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLASCL( 'G', 0, 0, ANRM, ONE, M, N, WORK, LDWORK, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLASCL : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDIF
!
!     Copy X or X' into the right place and scale it
!
   IF( TPSD ) THEN
!
!        Copy X into columns n+1:n+nrhs of work
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZLACPY( 'All', M, NRHS, X, LDX, WORK( N*LDWORK+1 ), &
                   LDWORK )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
      XNRM = ZLANGE( 'M', M, NRHS, WORK( N*LDWORK+1 ), LDWORK, &
             RWORK )
      IF( XNRM /= ZERO )  THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLASCL( 'G', 0, 0, XNRM, ONE, M, NRHS, &
                      WORK( N*LDWORK+1 ), LDWORK, INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLASCL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDIF
!
!        Compute QR factorization of X
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGEQR2( M, N+NRHS, WORK, LDWORK, &
                   WORK( LDWORK*( N+NRHS )+1 ), &
                   WORK( LDWORK*( N+NRHS )+MIN( M, N+NRHS )+1 ), &
                   INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGEQR2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        Compute largest entry in upper triangle of
!        work(n+1:m,n+1:n+nrhs)
!
      ERR = ZERO
      DO J = N + 1, N + NRHS
         DO I = N + 1, MIN( M, J )
            ERR = MAX( ERR, ABS( WORK( I+( J-1 )*M ) ) )
         ENDDO
      ENDDO
!
   ELSE
!
!        Copy X' into rows m+1:m+nrhs of work
!
      DO I = 1, N
         DO J = 1, NRHS
            WORK( M+J+( I-1 )*LDWORK ) = DCONJG( X( I, J ) )
         ENDDO
      ENDDO
!
      XNRM = ZLANGE( 'M', NRHS, N, WORK( M+1 ), LDWORK, RWORK )
      IF( XNRM /= ZERO )  THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL ZLASCL( 'G', 0, 0, XNRM, ONE, NRHS, N, WORK( M+1 ), &
                      LDWORK, INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZLASCL : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDIF
!
!        Compute LQ factorization of work
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL ZGELQ2( LDWORK, N, WORK, LDWORK, WORK( LDWORK*N+1 ), &
                   WORK( LDWORK*( N+1 )+1 ), INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZGELQ2 : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
!        Compute largest entry in lower triangle in
!        work(m+1:m+nrhs,m+1:n)
!
      ERR = ZERO
      DO J = M + 1, N
         DO I = J, LDWORK
            ERR = MAX( ERR, ABS( WORK( I+( J-1 )*LDWORK ) ) )
         ENDDO
      ENDDO
!
   END IF
!
   ZQRT14 = ERR / ( DBLE( MAX( M, N, NRHS ) )*DLAMCH( 'Epsilon' ) )
!
   RETURN
!
!     End of ZQRT14
!
END
                                                                                                                                                                                                                                                                                                            




