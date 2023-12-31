!> \brief \b SGBT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGBT01( M, N, KL, KU, A, LDA, AFAC, LDAFAC, IPIV, WORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            KL, KU, LDA, LDAFAC, M, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), AFAC( LDAFAC, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGBT01 reconstructs a band matrix A from its L*U factorization and
!> computes the residual:
!>    norm(L*U - A) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon.
!>
!> The expression L*U - A is computed one column at a time, so A and
!> AFAC are not modified.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The original matrix A in band storage, stored in rows 1 to
!>          KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER.
!>          The leading dimension of the array A.  LDA >= max(1,KL+KU+1).
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the banded
!>          factors L and U from the L*U factorization, as computed by
!>          SGBTRF.  U is stored as an upper triangular band matrix with
!>          KL+KU superdiagonals in rows 1 to KL+KU+1, and the
!>          multipliers used during the factorization are stored in rows
!>          KL+KU+2 to 2*KL+KU+1.  See SGBTRF for further details.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.
!>          LDAFAC >= max(1,2*KL*KU+1).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices from SGBTRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*KL+KU+1)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          norm(L*U - A) / ( N * norm(A) * EPS )
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
   SUBROUTINE SGBT01( M, N, KL, KU, A, LDA, AFAC, LDAFAC, IPIV, WORK, &
                      RESID )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   INTEGER            KL, KU, LDA, LDAFAC, M, N
   REAL               RESID
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   REAL               A( LDA, * ), AFAC( LDAFAC, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, I1, I2, IL, IP, IW, J, JL, JU, JUA, KD, LENJ
   REAL               ANORM, EPS, T
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   REAL               SASUM, SLAMCH
   EXTERNAL           SASUM, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           SAXPY, SCOPY
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MAX, MIN, REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if M = 0 or N = 0.
!
   RESID = ZERO
   IF( M <= 0 .OR. N <= 0 ) &
      RETURN
!
!     Determine EPS and the norm of A.
!
   EPS = SLAMCH( 'Epsilon' )
   KD = KU + 1
   ANORM = ZERO
   DO J = 1, N
      I1 = MAX( KD+1-J, 1 )
      I2 = MIN( KD+M-J, KL+KD )
      IF( I2 >= I1 ) &
         ANORM = MAX( ANORM, SASUM( I2-I1+1, A( I1, J ), 1 ) )
   ENDDO
!
!     Compute one column at a time of L*U - A.
!
   KD = KL + KU + 1
   DO J = 1, N
!
!        Copy the J-th column of U to WORK.
!
      JU = MIN( KL+KU, J-1 )
      JL = MIN( KL, M-J )
      LENJ = MIN( M, J ) - J + JU + 1
      IF( LENJ > 0 ) THEN
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
         CALL SCOPY( LENJ, AFAC( KD-JU, J ), 1, WORK, 1 )
#ifdef _TIMER
         call system_clock(count_rate=nb_periods_sec,count=S2_time)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : SCOPY : ',&
               real(S2_time-S1_time)/real(nb_periods_sec), ' s'
         close(10)
#endif
         DO I = LENJ + 1, JU + JL + 1
            WORK( I ) = ZERO
         ENDDO
!
!           Multiply by the unit lower triangular matrix L.  Note that L
!           is stored as a product of transformations and permutations.
!
         DO I = MIN( M-1, J ), J - JU, -1
            IL = MIN( KL, M-I )
            IF( IL > 0 ) THEN
               IW = I - J + JU + 1
               T = WORK( IW )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
               CALL SAXPY( IL, T, AFAC( KD+1, I ), 1, WORK( IW+1 ), &
                           1 )
#ifdef _TIMER
               call system_clock(count_rate=nb_periods_sec,count=S2_time)
               open(file='results.out', unit=10, position = 'append')
               write(10,'(A,F16.10,A)') 'Total time : SAXPY : ',&
                     real(S2_time-S1_time)/real(nb_periods_sec), ' s'
               close(10)
#endif
               IP = IPIV( I )
               IF( I /= IP ) THEN
                  IP = IP - J + JU + 1
                  WORK( IW ) = WORK( IP )
                  WORK( IP ) = T
               END IF
            END IF
         ENDDO
!
!           Subtract the corresponding column of A.
!
         JUA = MIN( JU, KU )
         IF( JUA+JL+1 > 0 )  THEN
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
            CALL SAXPY( JUA+JL+1, -ONE, A( KU+1-JUA, J ), 1, &
                        WORK( JU+1-JUA ), 1 )
#ifdef _TIMER
            call system_clock(count_rate=nb_periods_sec,count=S2_time)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : SAXPY : ',&
                  real(S2_time-S1_time)/real(nb_periods_sec), ' s'
            close(10)
#endif
         ENDIF
!
!           Compute the 1-norm of the column.
!
         RESID = MAX( RESID, SASUM( JU+JL+1, WORK, 1 ) )
      END IF
   ENDDO
!
!     Compute norm(L*U - A) / ( N * norm(A) * EPS )
!
   IF( ANORM <= ZERO ) THEN
      IF( RESID /= ZERO ) &
         RESID = ONE / EPS
   ELSE
      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
   END IF
!
   RETURN
!
!     End of SGBT01
!
END
                                                                                                                                                                                                                                                                                                            




