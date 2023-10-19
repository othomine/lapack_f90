!> \brief \b SBDT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            KD, LDA, LDPT, LDQ, M, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), D( * ), E( * ), PT( LDPT, * ),
!      $                   Q( LDQ, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SBDT01 reconstructs a general matrix A from its bidiagonal form
!>    A = Q * B * P**T
!> where Q (m by min(m,n)) and P**T (min(m,n) by n) are orthogonal
!> matrices and B is bidiagonal.
!>
!> The test ratio to test the reduction is
!>    RESID = norm(A - Q * B * P**T) / ( n * norm(A) * EPS )
!> where EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and P**T.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          If KD = 0, B is diagonal and the array E is not referenced.
!>          If KD = 1, the reduction was performed by xGEBRD; B is upper
!>          bidiagonal if M >= N, and lower bidiagonal if M < N.
!>          If KD = -1, the reduction was performed by xGBBRD; B is
!>          always upper bidiagonal.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,N)
!>          The m by min(m,n) orthogonal matrix Q in the reduction
!>          A = Q * B * P**T.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,M).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (min(M,N))
!>          The diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (min(M,N)-1)
!>          The superdiagonal elements of the bidiagonal matrix B if
!>          m >= n, or the subdiagonal elements of B if m < n.
!> \endverbatim
!>
!> \param[in] PT
!> \verbatim
!>          PT is REAL array, dimension (LDPT,N)
!>          The min(m,n) by n orthogonal matrix P**T in the reduction
!>          A = Q * B * P**T.
!> \endverbatim
!>
!> \param[in] LDPT
!> \verbatim
!>          LDPT is INTEGER
!>          The leading dimension of the array PT.
!>          LDPT >= max(1,min(M,N)).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (M+N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          The test ratio:
!>          norm(A - Q * B * P**T) / ( n * norm(A) * EPS )
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
   SUBROUTINE SBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KD, LDA, LDPT, LDQ, M, N
   REAL               RESID
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), D( * ), E( * ), PT( LDPT, * ), &
                      Q( LDQ, * ), WORK( * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, J
   REAL               ANORM, EPS
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   REAL               SASUM, SLAMCH, SLANGE
   EXTERNAL           SASUM, SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
   EXTERNAL           SCOPY, SGEMV
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
   IF( M <= 0 .OR. N <= 0 ) THEN
      RESID = 0.0E+0
      RETURN
   END IF
!
!     Compute A - Q * B * P**T one column at a time.
!
   RESID = 0.0E+0
   IF( KD /= 0 ) THEN
!
!        B is bidiagonal.
!
      IF( KD /= 0 .AND. M >= N ) THEN
!
!           B is upper bidiagonal and M >= N.
!
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            DO I = 1, N - 1
               WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
            ENDDO
            WORK( M+N ) = D( N )*PT( N, J )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'No transpose', M, N, -1.0E+0, Q, LDQ, &
                        WORK( M+1 ), 1, 1.0E+0, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
         ENDDO
      ELSE IF( KD < 0 ) THEN
!
!           B is upper bidiagonal and M < N.
!
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            DO I = 1, M - 1
               WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
            ENDDO
            WORK( M+M ) = D( M )*PT( M, J )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'No transpose', M, M, -1.0E+0, Q, LDQ, &
                        WORK( M+1 ), 1, 1.0E+0, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
         ENDDO
      ELSE
!
!           B is lower bidiagonal.
!
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            WORK( M+1 ) = D( 1 )*PT( 1, J )
            DO I = 2, M
               WORK( M+I ) = E( I-1 )*PT( I-1, J ) + &
                             D( I )*PT( I, J )
            ENDDO
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'No transpose', M, M, -1.0E+0, Q, LDQ, &
                        WORK( M+1 ), 1, 1.0E+0, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
         ENDDO
      END IF
   ELSE
!
!        B is diagonal.
!
      IF( M >= N ) THEN
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            DO I = 1, N
               WORK( M+I ) = D( I )*PT( I, J )
            ENDDO
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'No transpose', M, N, -1.0E+0, Q, LDQ, &
                        WORK( M+1 ), 1, 1.0E+0, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
         ENDDO
      ELSE
         DO J = 1, N
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            DO I = 1, M
               WORK( M+I ) = D( I )*PT( I, J )
            ENDDO
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SGEMV( 'No transpose', M, M, -1.0E+0, Q, LDQ, &
                        WORK( M+1 ), 1, 1.0E+0, WORK, 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SGEMV : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
            ENDDO
      END IF
   END IF
!
!     Compute norm(A - Q * B * P**T) / ( n * norm(A) * EPS )
!
   ANORM = SLANGE( '1', M, N, A, LDA, WORK )
   EPS = SLAMCH( 'Precision' )
!
   IF( ANORM <= 0.0E+0 ) THEN
      IF( RESID /= 0.0E+0 ) &
         RESID = 1.0E+0 / EPS
   ELSE
      IF( ANORM >= RESID ) THEN
         RESID = ( RESID / ANORM ) / ( REAL( N )*EPS )
      ELSE
         IF( ANORM < 1.0E+0 ) THEN
            RESID = ( MIN( RESID, REAL( N )*ANORM ) / ANORM ) / &
                    ( REAL( N )*EPS )
         ELSE
            RESID = MIN( RESID / ANORM, REAL( N ) ) / &
                    ( REAL( N )*EPS )
         END IF
      END IF
   END IF
!
   RETURN
!
!     End of SBDT01
!
END




