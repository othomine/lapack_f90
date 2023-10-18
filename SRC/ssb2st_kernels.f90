!> \brief \b SSB2ST_KERNELS
!
!  @generated from zhb2st_kernels.f, fortran z -> s, Wed Dec  7 08:22:40 2016
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSB2ST_KERNELS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssb2st_kernels.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssb2st_kernels.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssb2st_kernels.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE  SSB2ST_KERNELS( UPLO, WANTZ, TTYPE,
!                                   ST, ED, SWEEP, N, NB, IB,
!                                   A, LDA, V, TAU, LDVT, WORK)
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       LOGICAL            WANTZ
!       INTEGER            TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), V( * ),
!                          TAU( * ), WORK( * )
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSB2ST_KERNELS is an internal routine used by the SSYTRD_SB2ST
!> subroutine.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL which indicate if Eigenvalue are requested or both
!>          Eigenvalue/Eigenvectors.
!> \endverbatim
!>
!> \param[in] TTYPE
!> \verbatim
!>          TTYPE is INTEGER
!> \endverbatim
!>
!> \param[in] ST
!> \verbatim
!>          ST is INTEGER
!>          internal parameter for indices.
!> \endverbatim
!>
!> \param[in] ED
!> \verbatim
!>          ED is INTEGER
!>          internal parameter for indices.
!> \endverbatim
!>
!> \param[in] SWEEP
!> \verbatim
!>          SWEEP is INTEGER
!>          internal parameter for indices.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER. The order of the matrix A.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER. The size of the band.
!> \endverbatim
!>
!> \param[in] IB
!> \verbatim
!>          IB is INTEGER.
!> \endverbatim
!>
!> \param[in, out] A
!> \verbatim
!>          A is REAL array. A pointer to the matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER. The leading dimension of the matrix A.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is REAL array, dimension 2*n if eigenvalues only are
!>          requested or to be queried for vectors.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (2*n).
!>          The scalar factors of the Householder reflectors are stored
!>          in this array.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array. Workspace of size nb.
!> \endverbatim
!>
!>
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Implemented by Azzam Haidar.
!>
!>  All details are available on technical report, SC11, SC13 papers.
!>
!>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
!>  Parallel reduction to condensed forms for symmetric eigenvalue problems
!>  using aggregated fine-grained and memory-aware kernels. In Proceedings
!>  of 2011 International Conference for High Performance Computing,
!>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
!>  Article 8 , 11 pages.
!>  http://doi.acm.org/10.1145/2063384.2063394
!>
!>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
!>  An improved parallel singular value algorithm and its implementation
!>  for multicore hardware, In Proceedings of 2013 International Conference
!>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
!>  Denver, Colorado, USA, 2013.
!>  Article 90, 12 pages.
!>  http://doi.acm.org/10.1145/2503210.2503292
!>
!>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
!>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure
!>  calculations based on fine-grained memory aware tasks.
!>  International Journal of High Performance Computing Applications.
!>  Volume 28 Issue 2, Pages 196-209, May 2014.
!>  http://hpc.sagepub.com/content/28/2/196
!>
!> \endverbatim
!>
!> \ingroup hb2st_kernels
!>
!  =====================================================================
   SUBROUTINE  SSB2ST_KERNELS( UPLO, WANTZ, TTYPE, &
                               ST, ED, SWEEP, N, NB, IB, &
                               A, LDA, V, TAU, LDVT, WORK)
!
   IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   LOGICAL            WANTZ
   INTEGER            TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT
!     ..
!     .. Array Arguments ..
   REAL               A( LDA, * ), V( * ), &
                      TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   REAL               ZERO, ONE
   PARAMETER          ( ZERO = 0.0E+0, &
                      ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I, J1, J2, LM, LN, VPOS, TAUPOS, &
                      DPOS, OFDPOS, AJETER
   REAL               CTMP
!     ..
!     .. External Subroutines ..
   EXTERNAL           SLARFG, SLARFX, SLARFY
!     ..
!     .. Intrinsic Functions ..
   INTRINSIC          MOD
!     .. External Functions ..
   LOGICAL            LSAME
   EXTERNAL           LSAME
!     ..
!     ..
!     .. Executable Statements ..
!
   AJETER = IB + LDVT
   UPPER = LSAME( UPLO, 'U' )

   IF( UPPER ) THEN
       DPOS    = 2 * NB + 1
       OFDPOS  = 2 * NB
   ELSE
       DPOS    = 1
       OFDPOS  = 2
   ENDIF

!
!     Upper case
!
   IF( UPPER ) THEN
!
       IF( WANTZ ) THEN
           VPOS   = MOD( SWEEP-1, 2 ) * N + ST
           TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
       ELSE
           VPOS   = MOD( SWEEP-1, 2 ) * N + ST
           TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
       ENDIF
!
       IF( TTYPE == 1 ) THEN
           LM = ED - ST + 1
!
           V( VPOS ) = ONE
           DO I = 1, LM-1
               V( VPOS+I )         = ( A( OFDPOS-I, ST+I ) )
               A( OFDPOS-I, ST+I ) = ZERO
           ENDDO
           CTMP = ( A( OFDPOS, ST ) )
           CALL SLARFG( LM, CTMP, V( VPOS+1 ), 1, &
                                          TAU( TAUPOS ) )
           A( OFDPOS, ST ) = CTMP
!
           LM = ED - ST + 1
           CALL SLARFY( UPLO, LM, V( VPOS ), 1, &
                        ( TAU( TAUPOS ) ), &
                        A( DPOS, ST ), LDA-1, WORK)
       ENDIF
!
       IF( TTYPE == 3 ) THEN
!
           LM = ED - ST + 1
           CALL SLARFY( UPLO, LM, V( VPOS ), 1, &
                        ( TAU( TAUPOS ) ), &
                        A( DPOS, ST ), LDA-1, WORK)
       ENDIF
!
       IF( TTYPE == 2 ) THEN
           J1 = ED+1
           J2 = MIN( ED+NB, N )
           LN = ED-ST+1
           LM = J2-J1+1
           IF( LM > 0) THEN
               CALL SLARFX( 'Left', LN, LM, V( VPOS ), &
                            ( TAU( TAUPOS ) ), &
                            A( DPOS-NB, J1 ), LDA-1, WORK)
!
               IF( WANTZ ) THEN
                   VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                   TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
               ELSE
                   VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                   TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
               ENDIF
!
               V( VPOS ) = ONE
               DO I = 1, LM-1
                   V( VPOS+I )          = &
                                       ( A( DPOS-NB-I, J1+I ) )
                   A( DPOS-NB-I, J1+I ) = ZERO
               ENDDO
               CTMP = ( A( DPOS-NB, J1 ) )
               CALL SLARFG( LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) )
               A( DPOS-NB, J1 ) = CTMP
!
               CALL SLARFX( 'Right', LN-1, LM, V( VPOS ), &
                            TAU( TAUPOS ), &
                            A( DPOS-NB+1, J1 ), LDA-1, WORK)
           ENDIF
       ENDIF
!
!     Lower case
!
   ELSE
!
       IF( WANTZ ) THEN
           VPOS   = MOD( SWEEP-1, 2 ) * N + ST
           TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
       ELSE
           VPOS   = MOD( SWEEP-1, 2 ) * N + ST
           TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
       ENDIF
!
       IF( TTYPE == 1 ) THEN
           LM = ED - ST + 1
!
           V( VPOS ) = ONE
           DO I = 1, LM-1
               V( VPOS+I )         = A( OFDPOS+I, ST-1 )
               A( OFDPOS+I, ST-1 ) = ZERO
           ENDDO
           CALL SLARFG( LM, A( OFDPOS, ST-1 ), V( VPOS+1 ), 1, &
                                          TAU( TAUPOS ) )
!
           LM = ED - ST + 1
!
           CALL SLARFY( UPLO, LM, V( VPOS ), 1, &
                        ( TAU( TAUPOS ) ), &
                        A( DPOS, ST ), LDA-1, WORK)

       ENDIF
!
       IF( TTYPE == 3 ) THEN
           LM = ED - ST + 1
!
           CALL SLARFY( UPLO, LM, V( VPOS ), 1, &
                        ( TAU( TAUPOS ) ), &
                        A( DPOS, ST ), LDA-1, WORK)

       ENDIF
!
       IF( TTYPE == 2 ) THEN
           J1 = ED+1
           J2 = MIN( ED+NB, N )
           LN = ED-ST+1
           LM = J2-J1+1
!
           IF( LM > 0) THEN
               CALL SLARFX( 'Right', LM, LN, V( VPOS ), &
                            TAU( TAUPOS ), A( DPOS+NB, ST ), &
                            LDA-1, WORK)
!
               IF( WANTZ ) THEN
                   VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                   TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
               ELSE
                   VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                   TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
               ENDIF
!
               V( VPOS ) = ONE
               DO I = 1, LM-1
                   V( VPOS+I )        = A( DPOS+NB+I, ST )
                   A( DPOS+NB+I, ST ) = ZERO
               ENDDO
               CALL SLARFG( LM, A( DPOS+NB, ST ), V( VPOS+1 ), 1, &
                                           TAU( TAUPOS ) )
!
               CALL SLARFX( 'Left', LM, LN-1, V( VPOS ), &
                            ( TAU( TAUPOS ) ), &
                            A( DPOS+NB-1, ST+1 ), LDA-1, WORK)

           ENDIF
       ENDIF
   ENDIF
!
   RETURN
!
!     End of SSB2ST_KERNELS
!
END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
