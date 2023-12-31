!> \brief \b DQLT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQLT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AF( LDA, * ), C( LDA, * ), CC( LDA, * ),
!      $                   Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQLT03 tests DORMQL, which computes Q*C, Q'*C, C*Q or C*Q'.
!>
!> DQLT03 compares the results of a call to DORMQL with the results of
!> forming Q explicitly by a call to DORGQL and then performing matrix
!> multiplication by a call to DGEMM.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the orthogonal matrix Q.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows or columns of the matrix C; C is m-by-n if
!>          Q is applied from the left, or n-by-m if Q is applied from
!>          the right.  N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          orthogonal matrix Q.  M >= K >= 0.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the QL factorization of an m-by-n matrix, as
!>          returned by DGEQLF. See SGEQLF for further details.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] CC
!> \verbatim
!>          CC is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,M)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays AF, C, CC, and Q.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QL factorization in AF.
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
!>          The length of WORK.  LWORK must be at least M, and should be
!>          M*NB, where NB is the blocksize for this environment.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (4)
!>          The test ratios compare two techniques for multiplying a
!>          random matrix C by an m-by-m orthogonal matrix Q.
!>          RESULT(1) = norm( Q*C - Q*C )  / ( M * norm(C) * EPS )
!>          RESULT(2) = norm( C*Q - C*Q )  / ( M * norm(C) * EPS )
!>          RESULT(3) = norm( Q'*C - Q'*C )/ ( M * norm(C) * EPS )
!>          RESULT(4) = norm( C*Q' - C*Q' )/ ( M * norm(C) * EPS )
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
   SUBROUTINE DQLT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, &
                      RWORK, RESULT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AF( LDA, * ), C( LDA, * ), CC( LDA, * ), &
                      Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), &
                      WORK( LWORK )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
   DOUBLE PRECISION   ROGUE
   PARAMETER          ( ROGUE = -1.0D+10 )
!     ..
!     .. Local Scalars ..
   CHARACTER          SIDE, TRANS
   INTEGER            INFO, ISIDE, ITRANS, J, MC, MINMN, NC
   DOUBLE PRECISION   CNORM, EPS, RESID
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   DOUBLE PRECISION   DLAMCH, DLANGE
   EXTERNAL           LSAME, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           DGEMM, DLACPY, DLARNV, DLASET, DORGQL, DORMQL
!     ..
!     .. Local Arrays ..
   INTEGER            ISEED( 4 )
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Scalars in Common ..
   CHARACTER*32       SRNAMT
!     ..
!     .. Common blocks ..
   COMMON             / SRNAMC / SRNAMT
!     ..
!     .. Data statements ..
   DATA               ISEED / 1988, 1989, 1990, 1991 /
!     ..
!     .. Executable Statements ..
!
   EPS = DLAMCH( 'Epsilon' )
   MINMN = MIN( M, N )
!
!     Quick return if possible
!
   IF( MINMN == 0 ) THEN
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      RESULT( 3 ) = ZERO
      RESULT( 4 ) = ZERO
      RETURN
   END IF
!
!     Copy the last k columns of the factorization to the array Q
!
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( K > 0 .AND. M > K )  THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLACPY( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, &
                   Q( 1, M-K+1 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDIF
   IF( K > 1 )  THEN
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL DLACPY( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, &
                   Q( M-K+1, M-K+2 ), LDA )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
   ENDIF
!
!     Generate the m-by-m matrix Q
!
   SRNAMT = 'DORGQL'
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL DORGQL( M, M, K, Q, LDA, TAU( MINMN-K+1 ), WORK, LWORK, &
                INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : DORGQL : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
   DO ISIDE = 1, 2
      IF( ISIDE == 1 ) THEN
         SIDE = 'L'
         MC = M
         NC = N
      ELSE
         SIDE = 'R'
         MC = N
         NC = M
      END IF
!
!        Generate MC by NC matrix C
!
      DO J = 1, NC
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLARNV( 2, ISEED, MC, C( 1, J ) )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLARNV : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
      ENDDO
      CNORM = DLANGE( '1', MC, NC, C, LDA, RWORK )
      IF( CNORM == 0.0D0 ) &
         CNORM = ONE
!
      DO ITRANS = 1, 2
         IF( ITRANS == 1 ) THEN
            TRANS = 'N'
         ELSE
            TRANS = 'T'
         END IF
!
!           Copy C
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL DLACPY( 'Full', MC, NC, C, LDA, CC, LDA )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : DLACPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
!           Apply Q or Q' to C
!
         SRNAMT = 'DORMQL'
         IF( K > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DORMQL( SIDE, TRANS, MC, NC, K, AF( 1, N-K+1 ), LDA, &
                         TAU( MINMN-K+1 ), CC, LDA, WORK, LWORK, &
                         INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DORMQL : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
!           Form explicit product and subtract
!
         IF( LSAME( SIDE, 'L' ) ) THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGEMM( TRANS, 'No transpose', MC, NC, MC, -ONE, Q, &
                        LDA, C, LDA, ONE, CC, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ELSE
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL DGEMM( 'No transpose', TRANS, MC, NC, NC, -ONE, C, &
                        LDA, Q, LDA, ONE, CC, LDA )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : DGEMM : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         END IF
!
!           Compute error in the difference
!
         RESID = DLANGE( '1', MC, NC, CC, LDA, RWORK )
         RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID / &
            ( DBLE( MAX( 1, M ) )*CNORM*EPS )
!
      ENDDO
   ENDDO
!
   RETURN
!
!     End of DQLT03
!
END
                                                                                                                                                                                                                                                                                                            




