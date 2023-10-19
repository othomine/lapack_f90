!> \brief \b ZLATM4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATM4( ITYPE, N, NZ1, NZ2, RSIGN, AMAGN, RCOND,
!                          TRIANG, IDIST, ISEED, A, LDA )
!
!       .. Scalar Arguments ..
!       LOGICAL            RSIGN
!       INTEGER            IDIST, ITYPE, LDA, N, NZ1, NZ2
!       DOUBLE PRECISION   AMAGN, RCOND, TRIANG
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLATM4 generates basic square matrices, which may later be
!> multiplied by others in order to produce test matrices.  It is
!> intended mainly to be used to test the generalized eigenvalue
!> routines.
!>
!> It first generates the diagonal and (possibly) subdiagonal,
!> according to the value of ITYPE, NZ1, NZ2, RSIGN, AMAGN, and RCOND.
!> It then fills in the upper triangle with random numbers, if TRIANG is
!> non-zero.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          The "type" of matrix on the diagonal and sub-diagonal.
!>          If ITYPE < 0, then type abs(ITYPE) is generated and then
!>             swapped end for end (A(I,J) := A'(N-J,N-I).)  See also
!>             the description of AMAGN and RSIGN.
!>
!>          Special types:
!>          = 0:  the zero matrix.
!>          = 1:  the identity.
!>          = 2:  a transposed Jordan block.
!>          = 3:  If N is odd, then a k+1 x k+1 transposed Jordan block
!>                followed by a k x k identity block, where k=(N-1)/2.
!>                If N is even, then k=(N-2)/2, and a zero diagonal entry
!>                is tacked onto the end.
!>
!>          Diagonal types.  The diagonal consists of NZ1 zeros, then
!>             k=N-NZ1-NZ2 nonzeros.  The subdiagonal is zero.  ITYPE
!>             specifies the nonzero diagonal entries as follows:
!>          = 4:  1, ..., k
!>          = 5:  1, RCOND, ..., RCOND
!>          = 6:  1, ..., 1, RCOND
!>          = 7:  1, a, a^2, ..., a^(k-1)=RCOND
!>          = 8:  1, 1-d, 1-2*d, ..., 1-(k-1)*d=RCOND
!>          = 9:  random numbers chosen from (RCOND,1)
!>          = 10: random numbers with distribution IDIST (see ZLARND.)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.
!> \endverbatim
!>
!> \param[in] NZ1
!> \verbatim
!>          NZ1 is INTEGER
!>          If abs(ITYPE) > 3, then the first NZ1 diagonal entries will
!>          be zero.
!> \endverbatim
!>
!> \param[in] NZ2
!> \verbatim
!>          NZ2 is INTEGER
!>          If abs(ITYPE) > 3, then the last NZ2 diagonal entries will
!>          be zero.
!> \endverbatim
!>
!> \param[in] RSIGN
!> \verbatim
!>          RSIGN is LOGICAL
!>          = .TRUE.:  The diagonal and subdiagonal entries will be
!>                     multiplied by random numbers of magnitude 1.
!>          = .FALSE.: The diagonal and subdiagonal entries will be
!>                     left as they are (usually non-negative real.)
!> \endverbatim
!>
!> \param[in] AMAGN
!> \verbatim
!>          AMAGN is DOUBLE PRECISION
!>          The diagonal and subdiagonal entries will be multiplied by
!>          AMAGN.
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          If abs(ITYPE) > 4, then the smallest diagonal entry will be
!>          RCOND.  RCOND must be between 0 and 1.
!> \endverbatim
!>
!> \param[in] TRIANG
!> \verbatim
!>          TRIANG is DOUBLE PRECISION
!>          The entries above the diagonal will be random numbers with
!>          magnitude bounded by TRIANG (i.e., random numbers multiplied
!>          by TRIANG.)
!> \endverbatim
!>
!> \param[in] IDIST
!> \verbatim
!>          IDIST is INTEGER
!>          On entry, DIST specifies the type of distribution to be used
!>          to generate a random matrix .
!>          = 1: real and imaginary parts each UNIFORM( 0, 1 )
!>          = 2: real and imaginary parts each UNIFORM( -1, 1 )
!>          = 3: real and imaginary parts each NORMAL( 0, 1 )
!>          = 4: complex number uniform in DISK( 0, 1 )
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator.  The values of ISEED are changed on exit, and can
!>          be used in the next call to ZLATM4 to continue the same
!>          random number sequence.
!>          Note: ISEED(4) should be odd, for the random number generator
!>          used at present.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          Array to be computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          Leading dimension of A.  Must be at least 1 and at least N.
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
   SUBROUTINE ZLATM4( ITYPE, N, NZ1, NZ2, RSIGN, AMAGN, RCOND, &
                      TRIANG, IDIST, ISEED, A, LDA )
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   LOGICAL            RSIGN
   INTEGER            IDIST, ITYPE, LDA, N, NZ1, NZ2
   DOUBLE PRECISION   AMAGN, RCOND, TRIANG
!     ..
!     .. Array Arguments ..
   INTEGER            ISEED( 4 )
   COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   INTEGER            I, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN
   DOUBLE PRECISION   ALPHA
   COMPLEX*16         CTEMP
#ifdef _TIMER
      INTEGER(8)         nb_periods_sec, S1_time, S2_time
#endif
!     ..
!     .. External Functions ..
   DOUBLE PRECISION   DLARAN
   COMPLEX*16         ZLARND
   EXTERNAL           DLARAN, ZLARND
!     ..
!     .. External Subroutines ..
   EXTERNAL           ZLASET
!     ..
!     .. Executable Statements ..
!
   IF( N <= 0 ) RETURN
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S1_time)
#endif
   CALL ZLASET( 'Full', N, N, (0.0D+0,0.0D+0), (0.0D+0,0.0D+0), A, LDA )
#ifdef _TIMER
   call system_clock(count_rate=nb_periods_sec,count=S2_time)
   open(file='results.out', unit=10, position = 'append')
   write(10,'(A,F16.10,A)') 'Total time : ZLASET : ',&
         real(S2_time-S1_time)/real(nb_periods_sec), ' s'
   close(10)
#endif
!
!     Insure a correct ISEED
!
   IF( MOD( ISEED( 4 ), 2 ) /= 1 ) ISEED( 4 ) = ISEED( 4 ) + 1
!
!     Compute diagonal and subdiagonal according to ITYPE, NZ1, NZ2,
!     and RCOND
!
   IF( ITYPE /= 0 ) THEN
      IF( ABS( ITYPE ) >= 4 ) THEN
         KBEG = MAX( 1, MIN( N, NZ1+1 ) )
         KEND = MAX( KBEG, MIN( N, N-NZ2 ) )
         KLEN = KEND + 1 - KBEG
      ELSE
         KBEG = 1
         KEND = N
         KLEN = N
      END IF
      ISDB = 1
      ISDE = 0
      SELECT CASE (ABS( ITYPE ))
       CASE (1)
!
!        abs(ITYPE) = 1: Identity
!
        FORALL (JD = 1:N) A( JD, JD ) = (1.0D0,0.0D0)
       CASE (2)
!
!        abs(ITYPE) = 2: Transposed Jordan block
!
        FORALL (JD = 1:N-1) A( JD+1, JD ) = (1.0D0,0.0D0)
        ISDB = 1
        ISDE = N - 1
       CASE (3)
!
!        abs(ITYPE) = 3: Transposed Jordan block, followed by the
!                        identity.
!
        K = ( N-1 ) / 2
        FORALL (JD = 1:K) A( JD+1, JD ) = (1.0D0,0.0D0)
        ISDB = 1
        ISDE = K
        FORALL (JD = K + 2:2*K + 1) A( JD, JD ) = (1.0D0,0.0D0)
       CASE (4)
!
!        abs(ITYPE) = 4: 1,...,k
!
        FORALL (JD = KBEG:KEND) A( JD, JD ) = DCMPLX( JD-NZ1 )
       CASE (5)
!
!        abs(ITYPE) = 5: One large D value:
!
        FORALL (JD = KBEG+1:KEND) A( JD, JD ) = DCMPLX( RCOND )
        A( KBEG, KBEG ) = (1.0D0,0.0D0)
       CASE (6)
!
!        abs(ITYPE) = 6: One small D value:
!
        FORALL (JD = KBEG:KEND-1) A( JD, JD ) = (1.0D0,0.0D0)
        A( KEND, KEND ) = DCMPLX( RCOND )
       CASE (7)
!
!        abs(ITYPE) = 7: Exponentially distributed D values:
!
      A( KBEG, KBEG ) = (1.0D0,0.0D0)
      IF( KLEN > 1 ) THEN
         ALPHA = RCOND**( 1.0D0 / DBLE( KLEN-1 ) )
         DO I = 2, KLEN
            A( NZ1+I, NZ1+I ) = DCMPLX( ALPHA**DBLE( I-1 ) )
         ENDDO
      END IF
       CASE (8)
!
!        abs(ITYPE) = 8: Arithmetically distributed D values:
!
      A( KBEG, KBEG ) = (1.0D0,0.0D0)
      IF( KLEN > 1 ) THEN
         ALPHA = ( 1.0D0-RCOND ) / DBLE( KLEN-1 )
         DO I = 2, KLEN
            A( NZ1+I, NZ1+I ) = DCMPLX( DBLE( KLEN-I )*ALPHA+RCOND )
         ENDDO
      END IF
       CASE (9)
!
!        abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):
!
      ALPHA = LOG( RCOND )
      DO JD = KBEG, KEND
         A( JD, JD ) = EXP( ALPHA*DLARAN( ISEED ) )
      ENDDO
       CASE (10)
!
!        abs(ITYPE) = 10: Randomly distributed D values from DIST
!
      DO JD = KBEG, KEND
         A( JD, JD ) = ZLARND( IDIST, ISEED )
      ENDDO
!
       END SELECT
!
!        Scale by AMAGN
!
      FORALL (JD = KBEG:KEND) A( JD, JD ) = AMAGN*DBLE( A( JD, JD ) )
      FORALL (JD = ISDB:ISDE) A( JD+1, JD ) = AMAGN*DBLE( A( JD+1, JD ) )
!
!        If RSIGN = .TRUE., assign random signs to diagonal and
!        subdiagonal
!
      IF( RSIGN ) THEN
         DO JD = KBEG, KEND
            IF( DBLE( A( JD, JD ) ) /= 0.0D0 ) THEN
               CTEMP = ZLARND( 3, ISEED )
               CTEMP = CTEMP / ABS( CTEMP )
               A( JD, JD ) = CTEMP*DBLE( A( JD, JD ) )
            END IF
         ENDDO
         DO JD = ISDB, ISDE
            IF( DBLE( A( JD+1, JD ) ) /= 0.0D0 ) THEN
               CTEMP = ZLARND( 3, ISEED )
               CTEMP = CTEMP / ABS( CTEMP )
               A( JD+1, JD ) = CTEMP*DBLE( A( JD+1, JD ) )
            END IF
         ENDDO
      END IF
!
!        Reverse if ITYPE < 0
!
      IF( ITYPE < 0 ) THEN
         DO JD = KBEG, ( KBEG+KEND-1 ) / 2
            CTEMP = A( JD, JD )
            A( JD, JD ) = A( KBEG+KEND-JD, KBEG+KEND-JD )
            A( KBEG+KEND-JD, KBEG+KEND-JD ) = CTEMP
         ENDDO
         DO JD = 1, ( N-1 ) / 2
            CTEMP = A( JD+1, JD )
            A( JD+1, JD ) = A( N+1-JD, N-JD )
            A( N+1-JD, N-JD ) = CTEMP
         ENDDO
      END IF
!
   END IF
!
!     Fill in upper triangle
!
   IF( TRIANG /= 0.0D0 ) THEN
      DO JC = 2, N
         DO JR = 1, JC - 1
            A( JR, JC ) = TRIANG*ZLARND( IDIST, ISEED )
            ENDDO
         ENDDO
   END IF
!
   RETURN
!
!     End of ZLATM4
!
END




