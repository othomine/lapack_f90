!> \brief \b SCHKEQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKEQ( THRESH, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NOUT
!       REAL               THRESH
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKEQ tests SGEEQU, SGBEQU, SPOEQU, SPPEQU and SPBEQU
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          Threshold for testing routines. Should be between 2 and 10.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
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
!> \ingroup single_lin
!
!  =====================================================================
   SUBROUTINE SCHKEQ( THRESH, NOUT )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            NOUT
   REAL               THRESH
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE, TEN
   PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E+0, TEN = 1.0E1 )
   INTEGER            NSZ, NSZB
   PARAMETER          ( NSZ = 5, NSZB = 3*NSZ-2 )
   INTEGER            NSZP, NPOW
   PARAMETER          ( NSZP = ( NSZ*( NSZ+1 ) ) / 2, &
                      NPOW = 2*NSZ+1 )
!     ..
!     .. Local Scalars ..
   LOGICAL            OK
   CHARACTER*3        PATH
   INTEGER            I, INFO, J, KL, KU, M, N
   REAL               CCOND, EPS, NORM, RATIO, RCMAX, RCMIN, RCOND
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. Local Arrays ..
   REAL               A( NSZ, NSZ ), AB( NSZB, NSZ ), AP( NSZP ), &
                      C( NSZ ), POW( NPOW ), R( NSZ ), RESLTS( 5 ), &
                      RPOW( NPOW )
!     ..
!     .. External Functions ..
   REAL               SLAMCH
   EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SGBEQU, SGEEQU, SPBEQU, SPOEQU, SPPEQU
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
   PATH( 1:1 ) = 'Single precision'
   PATH( 2:3 ) = 'EQ'
!
   EPS = SLAMCH( 'P' )
   DO I = 1, 5
      RESLTS( I ) = ZERO
   ENDDO
   DO I = 1, NPOW
      POW( I ) = TEN**( I-1 )
      RPOW( I ) = ONE / POW( I )
   ENDDO
!
!     Test SGEEQU
!
   DO N = 0, NSZ
      DO M = 0, NSZ
!
         DO J = 1, NSZ
            DO I = 1, NSZ
               IF( I <= M .AND. J <= N ) THEN
                  A( I, J ) = POW( I+J+1 )*( -1 )**( I+J )
               ELSE
                  A( I, J ) = ZERO
               END IF
            ENDDO
         ENDDO
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SGEEQU( M, N, A, NSZ, R, C, RCOND, CCOND, NORM, INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SGEEQU : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         IF( INFO /= 0 ) THEN
            RESLTS( 1 ) = ONE
         ELSE
            IF( N /= 0 .AND. M /= 0 ) THEN
               RESLTS( 1 ) = MAX( RESLTS( 1 ), &
                             ABS( ( RCOND-RPOW( M ) ) / RPOW( M ) ) )
               RESLTS( 1 ) = MAX( RESLTS( 1 ), &
                             ABS( ( CCOND-RPOW( N ) ) / RPOW( N ) ) )
               RESLTS( 1 ) = MAX( RESLTS( 1 ), &
                             ABS( ( NORM-POW( N+M+1 ) ) / POW( N+M+ &
                             1 ) ) )
               DO I = 1, M
                  RESLTS( 1 ) = MAX( RESLTS( 1 ), &
                                ABS( ( R( I )-RPOW( I+N+1 ) ) / &
                                RPOW( I+N+1 ) ) )
               ENDDO
               DO J = 1, N
                  RESLTS( 1 ) = MAX( RESLTS( 1 ), &
                                ABS( ( C( J )-POW( N-J+1 ) ) / &
                                POW( N-J+1 ) ) )
               ENDDO
            END IF
         END IF
!
      ENDDO
   ENDDO
!
!     Test with zero rows and columns
!
   DO J = 1, NSZ
      A( MAX( NSZ-1, 1 ), J ) = ZERO
   ENDDO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEEQU( NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEEQU : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( INFO /= MAX( NSZ-1, 1 ) ) &
      RESLTS( 1 ) = ONE
!
   DO J = 1, NSZ
      A( MAX( NSZ-1, 1 ), J ) = ONE
      ENDDO
   DO I = 1, NSZ
      A( I, MAX( NSZ-1, 1 ) ) = ZERO
      ENDDO
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SGEEQU( NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SGEEQU : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( INFO /= NSZ+MAX( NSZ-1, 1 ) ) &
      RESLTS( 1 ) = ONE
   RESLTS( 1 ) = RESLTS( 1 ) / EPS
!
!     Test SGBEQU
!
   DO N = 0, NSZ
      DO M = 0, NSZ
         DO KL = 0, MAX( M-1, 0 )
            DO KU = 0, MAX( N-1, 0 )
!
               DO J = 1, NSZ
                  DO I = 1, NSZB
                     AB( I, J ) = ZERO
                     ENDDO
                  ENDDO
               DO J = 1, N
                  DO I = 1, M
                     IF( I <= MIN( M, J+KL ) .AND. I >= &
                         MAX( 1, J-KU ) .AND. J <= N ) THEN
                        AB( KU+1+I-J, J ) = POW( I+J+1 )* &
                                            ( -1 )**( I+J )
                     END IF
                     ENDDO
                  ENDDO
!
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SGBEQU( M, N, KL, KU, AB, NSZB, R, C, RCOND, &
                            CCOND, NORM, INFO )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SGBEQU : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
!
               IF( INFO /= 0 ) THEN
                  IF( .NOT.( ( N+KL < M .AND. INFO == N+KL+1 ) .OR. &
                      ( M+KU < N .AND. INFO == 2*M+KU+1 ) ) ) THEN
                     RESLTS( 2 ) = ONE
                  END IF
               ELSE
                  IF( N /= 0 .AND. M /= 0 ) THEN
!
                     RCMIN = R( 1 )
                     RCMAX = R( 1 )
                     DO I = 1, M
                        RCMIN = MIN( RCMIN, R( I ) )
                        RCMAX = MAX( RCMAX, R( I ) )
                        ENDDO
                     RATIO = RCMIN / RCMAX
                     RESLTS( 2 ) = MAX( RESLTS( 2 ), &
                                   ABS( ( RCOND-RATIO ) / RATIO ) )
!
                     RCMIN = C( 1 )
                     RCMAX = C( 1 )
                     DO J = 1, N
                        RCMIN = MIN( RCMIN, C( J ) )
                        RCMAX = MAX( RCMAX, C( J ) )
                        ENDDO
                     RATIO = RCMIN / RCMAX
                     RESLTS( 2 ) = MAX( RESLTS( 2 ), &
                                   ABS( ( CCOND-RATIO ) / RATIO ) )
!
                     RESLTS( 2 ) = MAX( RESLTS( 2 ), &
                                   ABS( ( NORM-POW( N+M+1 ) ) / &
                                   POW( N+M+1 ) ) )
                     DO I = 1, M
                        RCMAX = ZERO
                        DO J = 1, N
                           IF( I <= J+KL .AND. I >= J-KU ) THEN
                              RATIO = ABS( R( I )*POW( I+J+1 )* &
                                      C( J ) )
                              RCMAX = MAX( RCMAX, RATIO )
                           END IF
                           ENDDO
                        RESLTS( 2 ) = MAX( RESLTS( 2 ), &
                                      ABS( ONE-RCMAX ) )
                        ENDDO
!
                     DO J = 1, N
                        RCMAX = ZERO
                        DO I = 1, M
                           IF( I <= J+KL .AND. I >= J-KU ) THEN
                              RATIO = ABS( R( I )*POW( I+J+1 )* &
                                      C( J ) )
                              RCMAX = MAX( RCMAX, RATIO )
                           END IF
                           ENDDO
                        RESLTS( 2 ) = MAX( RESLTS( 2 ), &
                                      ABS( ONE-RCMAX ) )
                        ENDDO
                  END IF
               END IF
!
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   RESLTS( 2 ) = RESLTS( 2 ) / EPS
!
!     Test SPOEQU
!
   DO N = 0, NSZ
!
      DO I = 1, NSZ
         DO J = 1, NSZ
            IF( I <= N .AND. J == I ) THEN
               A( I, J ) = POW( I+J+1 )*( -1 )**( I+J )
            ELSE
               A( I, J ) = ZERO
            END IF
            ENDDO
         ENDDO
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SPOEQU( N, A, NSZ, R, RCOND, NORM, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SPOEQU : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      IF( INFO /= 0 ) THEN
         RESLTS( 3 ) = ONE
      ELSE
         IF( N /= 0 ) THEN
            RESLTS( 3 ) = MAX( RESLTS( 3 ), &
                          ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )
            RESLTS( 3 ) = MAX( RESLTS( 3 ), &
                          ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ &
                          1 ) ) )
            DO I = 1, N
               RESLTS( 3 ) = MAX( RESLTS( 3 ), &
                             ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ &
                             1 ) ) )
               ENDDO
         END IF
      END IF
      ENDDO
   A( MAX( NSZ-1, 1 ), MAX( NSZ-1, 1 ) ) = -ONE
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SPOEQU( NSZ, A, NSZ, R, RCOND, NORM, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SPOEQU : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( INFO /= MAX( NSZ-1, 1 ) ) &
      RESLTS( 3 ) = ONE
   RESLTS( 3 ) = RESLTS( 3 ) / EPS
!
!     Test SPPEQU
!
   DO N = 0, NSZ
!
!        Upper triangular packed storage
!
      DO I = 1, ( N*( N+1 ) ) / 2
         AP( I ) = ZERO
         ENDDO
      DO I = 1, N
         AP( ( I*( I+1 ) ) / 2 ) = POW( 2*I+1 )
         ENDDO
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SPPEQU( 'U', N, AP, R, RCOND, NORM, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SPPEQU : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      IF( INFO /= 0 ) THEN
         RESLTS( 4 ) = ONE
      ELSE
         IF( N /= 0 ) THEN
            RESLTS( 4 ) = MAX( RESLTS( 4 ), &
                          ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )
            RESLTS( 4 ) = MAX( RESLTS( 4 ), &
                          ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ &
                          1 ) ) )
            DO I = 1, N
               RESLTS( 4 ) = MAX( RESLTS( 4 ), &
                             ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ &
                             1 ) ) )
               ENDDO
         END IF
      END IF
!
!        Lower triangular packed storage
!
      DO I = 1, ( N*( N+1 ) ) / 2
         AP( I ) = ZERO
         ENDDO
      J = 1
      DO I = 1, N
         AP( J ) = POW( 2*I+1 )
         J = J + ( N-I+1 )
         ENDDO
!
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
      CALL SPPEQU( 'L', N, AP, R, RCOND, NORM, INFO )
#ifdef _TIMER
      call system_clock(count_rate=nb_periods_sec,count=S2_time)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : SPPEQU : ',&
            real(S2_time-S1_time)/real(nb_periods_sec), ' s'
      close(10)
#endif
!
      IF( INFO /= 0 ) THEN
         RESLTS( 4 ) = ONE
      ELSE
         IF( N /= 0 ) THEN
            RESLTS( 4 ) = MAX( RESLTS( 4 ), &
                          ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )
            RESLTS( 4 ) = MAX( RESLTS( 4 ), &
                          ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ &
                          1 ) ) )
            DO I = 1, N
               RESLTS( 4 ) = MAX( RESLTS( 4 ), &
                             ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ &
                             1 ) ) )
               ENDDO
         END IF
      END IF
!
      ENDDO
   I = ( NSZ*( NSZ+1 ) ) / 2 - 2
   AP( I ) = -ONE
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL SPPEQU( 'L', NSZ, AP, R, RCOND, NORM, INFO )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : SPPEQU : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
   IF( INFO /= MAX( NSZ-1, 1 ) ) &
      RESLTS( 4 ) = ONE
   RESLTS( 4 ) = RESLTS( 4 ) / EPS
!
!     Test SPBEQU
!
   DO N = 0, NSZ
      DO KL = 0, MAX( N-1, 0 )
!
!           Test upper triangular storage
!
         DO J = 1, NSZ
            DO I = 1, NSZB
               AB( I, J ) = ZERO
               ENDDO
            ENDDO
         DO J = 1, N
            AB( KL+1, J ) = POW( 2*J+1 )
            ENDDO
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SPBEQU( 'U', N, KL, AB, NSZB, R, RCOND, NORM, INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SPBEQU : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         IF( INFO /= 0 ) THEN
            RESLTS( 5 ) = ONE
         ELSE
            IF( N /= 0 ) THEN
               RESLTS( 5 ) = MAX( RESLTS( 5 ), &
                             ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )
               RESLTS( 5 ) = MAX( RESLTS( 5 ), &
                             ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ &
                             1 ) ) )
               DO I = 1, N
                  RESLTS( 5 ) = MAX( RESLTS( 5 ), &
                                ABS( ( R( I )-RPOW( I+1 ) ) / &
                                RPOW( I+1 ) ) )
                  ENDDO
            END IF
         END IF
         IF( N /= 0 ) THEN
            AB( KL+1, MAX( N-1, 1 ) ) = -ONE
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SPBEQU( 'U', N, KL, AB, NSZB, R, RCOND, NORM, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SPBEQU : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( INFO /= MAX( N-1, 1 ) ) &
               RESLTS( 5 ) = ONE
         END IF
!
!           Test lower triangular storage
!
         DO J = 1, NSZ
            DO I = 1, NSZB
               AB( I, J ) = ZERO
               ENDDO
            ENDDO
         DO J = 1, N
            AB( 1, J ) = POW( 2*J+1 )
            ENDDO
!
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SPBEQU( 'L', N, KL, AB, NSZB, R, RCOND, NORM, INFO )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SPBEQU : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
!
         IF( INFO /= 0 ) THEN
            RESLTS( 5 ) = ONE
         ELSE
            IF( N /= 0 ) THEN
               RESLTS( 5 ) = MAX( RESLTS( 5 ), &
                             ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )
               RESLTS( 5 ) = MAX( RESLTS( 5 ), &
                             ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ &
                             1 ) ) )
               DO I = 1, N
                  RESLTS( 5 ) = MAX( RESLTS( 5 ), &
                                ABS( ( R( I )-RPOW( I+1 ) ) / &
                                RPOW( I+1 ) ) )
                  ENDDO
            END IF
         END IF
         IF( N /= 0 ) THEN
            AB( 1, MAX( N-1, 1 ) ) = -ONE
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SPBEQU( 'L', N, KL, AB, NSZB, R, RCOND, NORM, INFO )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SPBEQU : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
            IF( INFO /= MAX( N-1, 1 ) ) &
               RESLTS( 5 ) = ONE
         END IF
         ENDDO
      ENDDO
   RESLTS( 5 ) = RESLTS( 5 ) / EPS
   OK = ( RESLTS( 1 ) <= THRESH ) .AND. &
        ( RESLTS( 2 ) <= THRESH ) .AND. &
        ( RESLTS( 3 ) <= THRESH ) .AND. &
        ( RESLTS( 4 ) <= THRESH ) .AND. ( RESLTS( 5 ) <= THRESH )
   WRITE( NOUT, FMT = * )
   IF( OK ) THEN
      WRITE( NOUT, FMT = 9999 )PATH
   ELSE
      IF( RESLTS( 1 ) > THRESH ) &
         WRITE( NOUT, FMT = 9998 )RESLTS( 1 ), THRESH
      IF( RESLTS( 2 ) > THRESH ) &
         WRITE( NOUT, FMT = 9997 )RESLTS( 2 ), THRESH
      IF( RESLTS( 3 ) > THRESH ) &
         WRITE( NOUT, FMT = 9996 )RESLTS( 3 ), THRESH
      IF( RESLTS( 4 ) > THRESH ) &
         WRITE( NOUT, FMT = 9995 )RESLTS( 4 ), THRESH
      IF( RESLTS( 5 ) > THRESH ) &
         WRITE( NOUT, FMT = 9994 )RESLTS( 5 ), THRESH
   END IF
 9999 FORMAT( 1X, 'All tests for ', A3, &
         ' routines passed the threshold' )
 9998 FORMAT( ' SGEEQU failed test with value ', E10.3, ' exceeding', &
         ' threshold ', E10.3 )
 9997 FORMAT( ' SGBEQU failed test with value ', E10.3, ' exceeding', &
         ' threshold ', E10.3 )
 9996 FORMAT( ' SPOEQU failed test with value ', E10.3, ' exceeding', &
         ' threshold ', E10.3 )
 9995 FORMAT( ' SPPEQU failed test with value ', E10.3, ' exceeding', &
         ' threshold ', E10.3 )
 9994 FORMAT( ' SPBEQU failed test with value ', E10.3, ' exceeding', &
         ' threshold ', E10.3 )
   RETURN
!
!     End of SCHKEQ
!
END
                                                                                                                                                                                                                                                                                                            




