!> \brief \b CSYTRI2X
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYTRI2X + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytri2x.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytri2x.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytri2x.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), WORK( N+NB+1,* )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYTRI2X computes the inverse of a real symmetric indefinite matrix
!> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!> CSYTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the NNB diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by CSYTRF.
!>
!>          On exit, if INFO = 0, the (symmetric) inverse of the original
!>          matrix.  If UPLO = 'U', the upper triangular part of the
!>          inverse is formed and the part of A below the diagonal is not
!>          referenced; if UPLO = 'L' the lower triangular part of the
!>          inverse is formed and the part of A above the diagonal is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the NNB structure of D
!>          as determined by CSYTRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N+NB+1,NB+3)
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Block size
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
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
!> \ingroup hetri2x
!
!  =====================================================================
   SUBROUTINE CSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
   IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, LDA, N, NB
!     ..
!     .. Array Arguments ..
   INTEGER            IPIV( * )
   COMPLEX            A( LDA, * ), WORK( N+NB+1,* )
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I, IINFO, IP, K, CUT, NNB
   INTEGER            J, U11, INVD

   COMPLEX   AK, AKKP1, AKP1, uoD, uoT
   COMPLEX   U01_I_J, U01_IP1_J
   COMPLEX   U11_I_J, U11_IP1_J
!     ..
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
   EXTERNAL           CSYCONV, XERBLA, CTRTRI
   EXTERNAL           CGEMM, CTRMM, CSYSWAPR
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N < 0 ) THEN
      INFO = -2
   ELSE IF( LDA < MAX( 1, N ) ) THEN
      INFO = -4
   END IF
!
!     Quick return if possible
!
!
   IF( INFO /= 0 ) THEN
      CALL XERBLA( 'CSYTRI2X', -INFO )
      RETURN
   END IF
   IF( N == 0 ) RETURN
!
!     Convert A
!     Workspace got Non-diag elements of D
!
   CALL CSYCONV( UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO )
!
!     Check that the diagonal matrix D is nonsingular.
!
   IF( UPPER ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
      DO INFO = N, 1, -1
         IF( IPIV( INFO ) > 0 .AND. A( INFO, INFO ) == (0.0E+0,0.0E+0) ) RETURN
      END DO
   ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
      DO INFO = 1, N
         IF( IPIV( INFO ) > 0 .AND. A( INFO, INFO ) == (0.0E+0,0.0E+0) ) RETURN
      END DO
   END IF
   INFO = 0
!
!  Splitting Workspace
!     U01 is a block (N,NB+1)
!     The first element of U01 is in WORK(1,1)
!     U11 is a block (NB+1,NB+1)
!     The first element of U11 is in WORK(N+1,1)
   U11 = N
!     INVD is a block (N,2)
!     The first element of INVD is in WORK(1,INVD)
   INVD = NB+2

   IF( UPPER ) THEN
!
!        invA = P * inv(U**T)*inv(D)*inv(U)*P**T.
!
     CALL CTRTRI( UPLO, 'U', N, A, LDA, INFO )
!
!       inv(D) and inv(D)*inv(U)
!
     K=1
     DO WHILE ( K  <=  N )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal NNB
          WORK(K,INVD) = (1.0E+0,0.0E+0) /  A( K, K )
          WORK(K,INVD+1) = 0
         K=K+1
      ELSE
!           2 x 2 diagonal NNB
          uoT = (1.0E+0,0.0E+0)/WORK(K+1,1)
          AK = A( K, K ) * uoT
          AKP1 = A( K+1, K+1 ) * uoT
          AKKP1 = WORK(K+1,1)  * uoT
          uoD = uoT/( AK*AKP1-(1.0E+0,0.0E+0) )
          WORK(K,INVD) = AKP1 * uoD
          WORK(K+1,INVD+1) = AK * uoD
          WORK(K,INVD+1) = -AKKP1 * uoD
          WORK(K+1,INVD) = -AKKP1 * uoD
         K=K+2
      END IF
     END DO
!
!       inv(U**T) = (inv(U))**T
!
!       inv(U**T)*inv(D)*inv(U)
!
     CUT=N
     DO WHILE (CUT  >  0)
        NNB=NB
        IF (CUT  <=  NNB) THEN
           NNB=CUT
        ELSE
!             count negative elements,
!             need a even number for a clear cut
           IF (MOD(COUNT(IPIV(CUT+1-NNB:CUT)  <  0),2)  ==  1) NNB=NNB+1
        END IF

        CUT=CUT-NNB
!
!          U01 Block
!
        WORK(1:CUT,1:NNB)=A(1:CUT,CUT+1:CUT+NNB)
!
!          U11 Block
!
        DO I=1,NNB
          WORK(U11+I,1:I-1)=(0.0E+0,0.0E+0)
          WORK(U11+I,I)=(1.0E+0,0.0E+0)
          WORK(U11+I,I+1:NNB)=A(CUT+I,CUT+I+1:CUT+NNB)
        END DO
!
!          invD*U01
!
        I=1
        DO WHILE (I  <=  CUT)
          IF (IPIV(I) > 0) THEN
             WORK(I,1:NNB)=WORK(I,INVD)*WORK(I,1:NNB)
             I=I+1
          ELSE
             DO J=1,NNB
                U01_I_J = WORK(I,J)
                U01_IP1_J = WORK(I+1,J)
                WORK(I,J)=WORK(I,INVD)*U01_I_J+ WORK(I,INVD+1)*U01_IP1_J
                WORK(I+1,J)=WORK(I+1,INVD)*U01_I_J+ WORK(I+1,INVD+1)*U01_IP1_J
             END DO
             I=I+2
          END IF
        END DO
!
!        invD1*U11
!
        I=1
        DO WHILE (I  <=  NNB)
          IF (IPIV(CUT+I) > 0) THEN
             WORK(U11+I,I:NNB)=WORK(CUT+I,INVD)*WORK(U11+I,I:NNB)
             I=I+1
          ELSE
             DO J=I,NNB
                U11_I_J = WORK(U11+I,J)
                U11_IP1_J = WORK(U11+I+1,J)
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*WORK(U11+I+1,J)
                WORK(U11+I+1,J)=WORK(CUT+I+1,INVD)*U11_I_J+ WORK(CUT+I+1,INVD+1)*U11_IP1_J
             END DO
             I=I+2
          END IF
        END DO
!
!       U11**T*invD1*U11->U11
!
     CALL CTRMM('L','U','T','U',NNB, NNB, &
                (1.0E+0,0.0E+0),A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1)
!
      DO I=1,NNB
         A(CUT+I,CUT+I:CUT+NNB)=WORK(U11+I,I:NNB)
      END DO
!
!          U01**T*invD*U01->A(CUT+I,CUT+J)
!
      CALL CGEMM('T','N',NNB,NNB,CUT,(1.0E+0,0.0E+0),A(1,CUT+1),LDA, &
                 WORK,N+NB+1, (0.0E+0,0.0E+0), WORK(U11+1,1), N+NB+1)
!
!        U11 =  U11**T*invD1*U11 + U01**T*invD*U01
!
      DO I=1,NNB
         A(CUT+I,CUT+I:CUT+NNB)=A(CUT+I,CUT+I:CUT+NNB)+WORK(U11+I,I:NNB)
      END DO
!
!        U01 =  U00**T*invD0*U01
!
      CALL CTRMM('L',UPLO,'T','U',CUT, NNB, (1.0E+0,0.0E+0),A,LDA,WORK,N+NB+1)

!
!        Update U01
!
!       A(1:CUT,CUT+1:CUT+NNB)=WORK(1:CUT,1:NNB)
      DO I=1,CUT
        DO J=1,NNB
         A(I,CUT+J)=WORK(I,J)
        END DO
      END DO
!
!      Next Block
!
    END DO
!
!        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T
!
         I=1
         DO WHILE ( I  <=  N )
            IF( IPIV(I)  >  0 ) THEN
               IP=IPIV(I)
              IF (I  <  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, I ,IP )
              ELSEIF (I  >  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, IP ,I )
              ENDIF
            ELSE
              IP=-IPIV(I)
              I=I+1
              IF ( (I-1)  <  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, I-1 ,IP )
              ELSEIF ( (I-1)  >  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, IP ,I-1 )
              ENDIF
           ENDIF
            I=I+1
         END DO
   ELSE
!
!        LOWER...
!
!        invA = P * inv(U**T)*inv(D)*inv(U)*P**T.
!
      CALL CTRTRI( UPLO, 'U', N, A, LDA, INFO )
!
!       inv(D) and inv(D)*inv(U)
!
     K=N
     DO WHILE ( K  >=  1 )
      IF( IPIV( K ) > 0 ) THEN
!           1 x 1 diagonal NNB
          WORK(K,INVD) = (1.0E+0,0.0E+0) /  A( K, K )
          WORK(K,INVD+1) = 0
         K=K-1
      ELSE
!           2 x 2 diagonal NNB
          uoT = (1.0E+0,0.0E+0)/WORK(K-1,1)
          AK = A( K-1, K-1 ) * uoT
          AKP1 = A( K, K ) * uoT
          AKKP1 = WORK(K-1,1) * uoT
          uoD = uoT/( AK*AKP1-(1.0E+0,0.0E+0) )
          WORK(K-1,INVD) = AKP1 * uoD
          WORK(K,INVD) = AK * uoD
          WORK(K,INVD+1) = -AKKP1 * uoD
          WORK(K-1,INVD+1) = -AKKP1 * uoD
         K=K-2
      END IF
     END DO
!
!       inv(U**T) = (inv(U))**T
!
!       inv(U**T)*inv(D)*inv(U)
!
     CUT=0
     DO WHILE (CUT  <  N)
        NNB=NB
        IF (CUT + NNB  >=  N) THEN
           NNB=N-CUT
        ELSE
!             count negative elements,
!             need a even number for a clear cut
           IF (MOD(COUNT(IPIV(CUT+1:CUT+NNB)  <  0),2)  ==  1) NNB=NNB+1
        END IF
!      L21 Block
       WORK(1:N-CUT-NNB,1:NNB)=A(CUT+NNB+1:N,CUT+1:CUT+NNB)
!     L11 Block
        DO I=1,NNB
          WORK(U11+I,1:I-1)=A(CUT+I,CUT+1:CUT+I-1)
          WORK(U11+I,I)=(1.0E+0,0.0E+0)
          WORK(U11+I,I+1:NNB)=(0.0E+0,0.0E+0)
        END DO
!
!          invD*L21
!
        I=N-CUT-NNB
        DO WHILE (I  >=  1)
          IF (IPIV(CUT+NNB+I) > 0) THEN
             WORK(I,1:NNB)=WORK(CUT+NNB+I,INVD)*WORK(I,1:NNB)
             I=I-1
          ELSE
             DO J=1,NNB
                U01_I_J = WORK(I,J)
                U01_IP1_J = WORK(I-1,J)
                WORK(I,J)=WORK(CUT+NNB+I,INVD)*U01_I_J+ &
                           WORK(CUT+NNB+I,INVD+1)*U01_IP1_J
                WORK(I-1,J)=WORK(CUT+NNB+I-1,INVD+1)*U01_I_J+ &
                           WORK(CUT+NNB+I-1,INVD)*U01_IP1_J
             END DO
             I=I-2
          END IF
        END DO
!
!        invD1*L11
!
        I=NNB
        DO WHILE (I  >=  1)
          IF (IPIV(CUT+I) > 0) THEN
             WORK(U11+I,1:NNB)=WORK(CUT+I,INVD)*WORK(U11+I,1:NNB)
             I=I-1
          ELSE
             DO J=1,NNB
                U11_I_J = WORK(U11+I,J)
                U11_IP1_J = WORK(U11+I-1,J)
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*U11_IP1_J
                WORK(U11+I-1,J)=WORK(CUT+I-1,INVD+1)*U11_I_J+ WORK(CUT+I-1,INVD)*U11_IP1_J
             END DO
             I=I-2
          END IF
        END DO
!
!       L11**T*invD1*L11->L11
!
     CALL CTRMM('L',UPLO,'T','U',NNB, NNB, &
                (1.0E+0,0.0E+0),A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1)
!
      DO I=1,NNB
         A(CUT+I,CUT+1:CUT+I)=WORK(U11+I,1:I)
      END DO
!
     IF ( (CUT+NNB)  <  N ) THEN
!
!          L21**T*invD2*L21->A(CUT+I,CUT+J)
!
      CALL CGEMM('T','N',NNB,NNB,N-NNB-CUT,(1.0E+0,0.0E+0),A(CUT+NNB+1,CUT+1) &
                ,LDA,WORK,N+NB+1, (0.0E+0,0.0E+0), WORK(U11+1,1), N+NB+1)

!
!        L11 =  L11**T*invD1*L11 + U01**T*invD*U01
!
      DO I=1,NNB
         A(CUT+I,CUT+1:CUT+I)=A(CUT+I,CUT+1:CUT+I)+WORK(U11+I,1:I)
      END DO
!
!        L01 =  L22**T*invD2*L21
!
      CALL CTRMM('L',UPLO,'T','U', N-NNB-CUT, NNB, &
                (1.0E+0,0.0E+0),A(CUT+NNB+1,CUT+NNB+1),LDA,WORK,N+NB+1)

!      Update L21
      A(CUT+NNB+1:N,CUT+1:CUT+NNB)=WORK(1:N-CUT-NNB,1:NNB)
    ELSE
!
!        L11 =  L11**T*invD1*L11
!
      DO I=1,NNB
         A(CUT+I,CUT+1:CUT+I)=WORK(U11+I,1:I)
      END DO
    END IF
!
!      Next Block
!
        CUT=CUT+NNB
    END DO
!
!        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T
!
         I=N
         DO WHILE ( I  >=  1 )
            IF( IPIV(I)  >  0 ) THEN
               IP=IPIV(I)
              IF (I  <  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, I ,IP  )
              ELSEIF (I  >  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, IP ,I )
              ENDIF
            ELSE
              IP=-IPIV(I)
              IF ( I  <  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, I ,IP )
              ELSEIF ( I  >  IP) THEN
                 CALL CSYSWAPR( UPLO, N, A, LDA, IP ,I )
              ENDIF
              I=I-1
            ENDIF
            I=I-1
         END DO
   END IF
!
   RETURN
!
!     End of CSYTRI2X
!
   END
