!> \brief \b ZCHKEE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM ZCHKEE
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKEE tests the COMPLEX*16 LAPACK subroutines for the matrix
!> eigenvalue problem.  The test paths in this version are
!>
!> NEP (Nonsymmetric Eigenvalue Problem):
!>     Test ZGEHRD, ZUNGHR, ZHSEQR, ZTREVC, ZHSEIN, and ZUNMHR
!>
!> SEP (Hermitian Eigenvalue Problem):
!>     Test ZHETRD, ZUNGTR, ZSTEQR, ZSTERF, ZSTEIN, ZSTEDC,
!>     and drivers ZHEEV(X), ZHBEV(X), ZHPEV(X),
!>                 ZHEEVD,   ZHBEVD,   ZHPEVD
!>
!> SVD (Singular Value Decomposition):
!>     Test ZGEBRD, ZUNGBR, and ZBDSQR
!>     and the drivers ZGESVD, ZGESDD
!>
!> ZEV (Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test ZGEEV
!>
!> ZES (Nonsymmetric Schur form Driver):
!>     Test ZGEES
!>
!> ZVX (Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test ZGEEVX
!>
!> ZSX (Nonsymmetric Schur form Expert Driver):
!>     Test ZGEESX
!>
!> ZGG (Generalized Nonsymmetric Eigenvalue Problem):
!>     Test ZGGHD3, ZGGBAL, ZGGBAK, ZHGEQZ, and ZTGEVC
!>
!> ZGS (Generalized Nonsymmetric Schur form Driver):
!>     Test ZGGES
!>
!> ZGV (Generalized Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test ZGGEV
!>
!> ZGX (Generalized Nonsymmetric Schur form Expert Driver):
!>     Test ZGGESX
!>
!> ZXV (Generalized Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test ZGGEVX
!>
!> ZSG (Hermitian Generalized Eigenvalue Problem):
!>     Test ZHEGST, ZHEGV, ZHEGVD, ZHEGVX, ZHPGST, ZHPGV, ZHPGVD,
!>     ZHPGVX, ZHBGST, ZHBGV, ZHBGVD, and ZHBGVX
!>
!> ZHB (Hermitian Band Eigenvalue Problem):
!>     Test ZHBTRD
!>
!> ZBB (Band Singular Value Decomposition):
!>     Test ZGBBRD
!>
!> ZEC (Eigencondition estimation):
!>     Test ZTRSYL, ZTREXC, ZTRSNA, and ZTRSEN
!>
!> ZBL (Balancing a general matrix)
!>     Test ZGEBAL
!>
!> ZBK (Back transformation on a balanced matrix)
!>     Test ZGEBAK
!>
!> ZGL (Balancing a matrix pair)
!>     Test ZGGBAL
!>
!> ZGK (Back transformation on a matrix pair)
!>     Test ZGGBAK
!>
!> GLM (Generalized Linear Regression Model):
!>     Tests ZGGGLM
!>
!> GQR (Generalized QR and RQ factorizations):
!>     Tests ZGGQRF and ZGGRQF
!>
!> GSV (Generalized Singular Value Decomposition):
!>     Tests ZGGSVD, ZGGSVP, ZTGSJA, ZLAGS2, ZLAPLL, and ZLAPMT
!>
!> CSD (CS decomposition):
!>     Tests ZUNCSD
!>
!> LSE (Constrained Linear Least Squares):
!>     Tests ZGGLSE
!>
!> Each test path has a different set of inputs, but the data sets for
!> the driver routines xEV, xES, xVX, and xSX can be concatenated in a
!> single input file.  The first line of input should contain one of the
!> 3-character path names in columns 1-3.  The number of remaining lines
!> depends on what is found on the first line.
!>
!> The number of matrix types used in testing is often controllable from
!> the input file.  The number of matrix types for each path, and the
!> test routine that describes them, is as follows:
!>
!> Path name(s)  Types    Test routine
!>
!> ZHS or NEP      21     ZCHKHS
!> ZST or SEP      21     ZCHKST (routines)
!>                 18     ZDRVST (drivers)
!> ZBD or SVD      16     ZCHKBD (routines)
!>                  5     ZDRVBD (drivers)
!> ZEV             21     ZDRVEV
!> ZES             21     ZDRVES
!> ZVX             21     ZDRVVX
!> ZSX             21     ZDRVSX
!> ZGG             26     ZCHKGG (routines)
!> ZGS             26     ZDRGES
!> ZGX              5     ZDRGSX
!> ZGV             26     ZDRGEV
!> ZXV              2     ZDRGVX
!> ZSG             21     ZDRVSG
!> ZHB             15     ZCHKHB
!> ZBB             15     ZCHKBB
!> ZEC              -     ZCHKEC
!> ZBL              -     ZCHKBL
!> ZBK              -     ZCHKBK
!> ZGL              -     ZCHKGL
!> ZGK              -     ZCHKGK
!> GLM              8     ZCKGLM
!> GQR              8     ZCKGQR
!> GSV              8     ZCKGSV
!> CSD              3     ZCKCSD
!> LSE              8     ZCKLSE
!>
!>-----------------------------------------------------------------------
!>
!> NEP input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, NX, NS, and
!>          MAXB.
!>
!> line 5:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 6:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for the minimum blocksize NBMIN.
!>
!> line 7:  NXVAL, INTEGER array, dimension (NPARMS)
!>          The values for the crossover point NX.
!>
!> line 8:  INMIN, INTEGER array, dimension (NPARMS)
!>          LAHQR vs TTQRE crossover point, >= 11
!>
!> line 9:  INWIN, INTEGER array, dimension (NPARMS)
!>          recommended deflation window size
!>
!> line 10: INIBL, INTEGER array, dimension (NPARMS)
!>          nibble crossover point
!>
!> line 11:  ISHFTS, INTEGER array, dimension (NPARMS)
!>          number of simultaneous shifts)
!>
!> line 12:  IACC22, INTEGER array, dimension (NPARMS)
!>          select structured matrix multiply: 0, 1 or 2)
!>
!> line 13: THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.  To have all of the test
!>          ratios printed, use THRESH = 0.0 .
!>
!> line 14: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 14 was 2:
!>
!> line 15: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 15-EOF:  The remaining lines occur in sets of 1 or 2 and allow
!>          the user to specify the matrix types.  Each line contains
!>          a 3-character path name in columns 1-3, and the number
!>          of matrix types must be the first nonblank item in columns
!>          4-80.  If the number of matrix types is at least 1 but is
!>          less than the maximum number of possible types, a second
!>          line will be read to get the numbers of the matrix types to
!>          be used.  For example,
!> NEP 21
!>          requests all of the matrix types for the nonsymmetric
!>          eigenvalue problem, while
!> NEP  4
!> 9 10 11 12
!>          requests only matrices of type 9, 10, 11, and 12.
!>
!>          The valid 3-character path names are 'NEP' or 'ZHS' for the
!>          nonsymmetric eigenvalue routines.
!>
!>-----------------------------------------------------------------------
!>
!> SEP or ZSG input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, and NX.
!>
!> line 5:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 6:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for the minimum blocksize NBMIN.
!>
!> line 7:  NXVAL, INTEGER array, dimension (NPARMS)
!>          The values for the crossover point NX.
!>
!> line 8:  THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 9:  TSTCHK, LOGICAL
!>          Flag indicating whether or not to test the LAPACK routines.
!>
!> line 10: TSTDRV, LOGICAL
!>          Flag indicating whether or not to test the driver routines.
!>
!> line 11: TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 12: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 12 was 2:
!>
!> line 13: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 13-EOF:  Lines specifying matrix types, as for NEP.
!>          The valid 3-character path names are 'SEP' or 'ZST' for the
!>          Hermitian eigenvalue routines and driver routines, and
!>          'ZSG' for the routines for the Hermitian generalized
!>          eigenvalue problem.
!>
!>-----------------------------------------------------------------------
!>
!> SVD input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of M and N.
!>
!> line 3:  MVAL, INTEGER array, dimension (NN)
!>          The values for the matrix row dimension M.
!>
!> line 4:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix column dimension N.
!>
!> line 5:  NPARMS, INTEGER
!>          Number of values of the parameter NB, NBMIN, NX, and NRHS.
!>
!> line 6:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 7:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for the minimum blocksize NBMIN.
!>
!> line 8:  NXVAL, INTEGER array, dimension (NPARMS)
!>          The values for the crossover point NX.
!>
!> line 9:  NSVAL, INTEGER array, dimension (NPARMS)
!>          The values for the number of right hand sides NRHS.
!>
!> line 10: THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 11: TSTCHK, LOGICAL
!>          Flag indicating whether or not to test the LAPACK routines.
!>
!> line 12: TSTDRV, LOGICAL
!>          Flag indicating whether or not to test the driver routines.
!>
!> line 13: TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 14: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 14 was 2:
!>
!> line 15: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 15-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path names are 'SVD' or 'ZBD' for both the
!>          SVD routines and the SVD driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> ZEV and ZES data files:
!>
!> line 1:  'ZEV' or 'ZES' in columns 1 to 3.
!>
!> line 2:  NSIZES, INTEGER
!>          Number of sizes of matrices to use. Should be at least 0
!>          and at most 20. If NSIZES = 0, no testing is done
!>          (although the remaining  3 lines are still read).
!>
!> line 3:  NN, INTEGER array, dimension(NSIZES)
!>          Dimensions of matrices to be tested.
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHSEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 5:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          If it is 0., all test case data will be printed.
!>
!> line 6:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 6 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 8 and following:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'ZEV' to test CGEEV, or
!>          'ZES' to test CGEES.
!>
!>-----------------------------------------------------------------------
!>
!> The ZVX data has two parts. The first part is identical to ZEV,
!> and the second part consists of test matrices with precomputed
!> solutions.
!>
!> line 1:  'ZVX' in columns 1-3.
!>
!> line 2:  NSIZES, INTEGER
!>          If NSIZES = 0, no testing of randomly generated examples
!>          is done, but any precomputed examples are tested.
!>
!> line 3:  NN, INTEGER array, dimension(NSIZES)
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>
!> line 5:  THRESH, REAL
!>
!> line 6:  NEWSD, INTEGER
!>
!> If line 6 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>
!> lines 8 and following: The first line contains 'ZVX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 1+N+N**2 lines, where N is
!>          its dimension. The first line contains the dimension N and
!>          ISRT (two integers). ISRT indicates whether the last N lines
!>          are sorted by increasing real part of the eigenvalue
!>          (ISRT=0) or by increasing imaginary part (ISRT=1). The next
!>          N**2 lines contain the matrix rowwise, one entry per line.
!>          The last N lines correspond to each eigenvalue. Each of
!>          these last N lines contains 4 real values: the real part of
!>          the eigenvalues, the imaginary part of the eigenvalue, the
!>          reciprocal condition number of the eigenvalues, and the
!>          reciprocal condition number of the vector eigenvector. The
!>          end of data is indicated by dimension N=0. Even if no data
!>          is to be tested, there must be at least one line containing
!>          N=0.
!>
!>-----------------------------------------------------------------------
!>
!> The ZSX data is like ZVX. The first part is identical to ZEV, and the
!> second part consists of test matrices with precomputed solutions.
!>
!> line 1:  'ZSX' in columns 1-3.
!>
!> line 2:  NSIZES, INTEGER
!>          If NSIZES = 0, no testing of randomly generated examples
!>          is done, but any precomputed examples are tested.
!>
!> line 3:  NN, INTEGER array, dimension(NSIZES)
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>
!> line 5:  THRESH, REAL
!>
!> line 6:  NEWSD, INTEGER
!>
!> If line 6 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>
!> lines 8 and following: The first line contains 'ZSX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 3+N**2 lines, where N is
!>          its dimension. The first line contains the dimension N, the
!>          dimension M of an invariant subspace, and ISRT. The second
!>          line contains M integers, identifying the eigenvalues in the
!>          invariant subspace (by their position in a list of
!>          eigenvalues ordered by increasing real part (if ISRT=0) or
!>          by increasing imaginary part (if ISRT=1)). The next N**2
!>          lines contain the matrix rowwise. The last line contains the
!>          reciprocal condition number for the average of the selected
!>          eigenvalues, and the reciprocal condition number for the
!>          corresponding right invariant subspace. The end of data in
!>          indicated by a line containing N=0, M=0, and ISRT = 0.  Even
!>          if no data is to be tested, there must be at least one line
!>          containing N=0, M=0 and ISRT=0.
!>
!>-----------------------------------------------------------------------
!>
!> ZGG input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, NBCOL, NS, and
!>          MAXB.
!>
!> line 5:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 6:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for NBMIN, the minimum row dimension for blocks.
!>
!> line 7:  NSVAL, INTEGER array, dimension (NPARMS)
!>          The values for the number of shifts.
!>
!> line 8:  MXBVAL, INTEGER array, dimension (NPARMS)
!>          The values for MAXB, used in determining minimum blocksize.
!>
!> line 9:  IACC22, INTEGER array, dimension (NPARMS)
!>          select structured matrix multiply: 1 or 2)
!>
!> line 10: NBCOL, INTEGER array, dimension (NPARMS)
!>          The values for NBCOL, the minimum column dimension for
!>          blocks.
!>
!> line 11: THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 12: TSTCHK, LOGICAL
!>          Flag indicating whether or not to test the LAPACK routines.
!>
!> line 13: TSTDRV, LOGICAL
!>          Flag indicating whether or not to test the driver routines.
!>
!> line 14: TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 15: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 15 was 2:
!>
!> line 16: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 17-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'ZGG' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> ZGS and ZGV input files:
!>
!> line 1:  'ZGS' or 'ZGV' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension(NN)
!>          Dimensions of matrices to be tested.
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHGEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 5:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          If it is 0., all test case data will be printed.
!>
!> line 6:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits.
!>
!> line 7:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 17 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 7-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'ZGS' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> ZGX input file:
!> line 1:  'ZGX' in columns 1 to 3.
!>
!> line 2:  N, INTEGER
!>          Value of N.
!>
!> line 3:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHGEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 4:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          Information will be printed about each test for which the
!>          test ratio is greater than or equal to the threshold.
!>
!> line 5:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 6:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 6 was 2:
!>
!> line 7: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> If line 2 was 0:
!>
!> line 7-EOF: Precomputed examples are tested.
!>
!> remaining lines : Each example is stored on 3+2*N*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer).  The next line contains an integer k such
!>          that only the last k eigenvalues will be selected and appear
!>          in the leading diagonal blocks of $A$ and $B$. The next N*N
!>          lines contain the matrix A, one element per line. The next N*N
!>          lines contain the matrix B. The last line contains the
!>          reciprocal of the eigenvalue cluster condition number and the
!>          reciprocal of the deflating subspace (associated with the
!>          selected eigencluster) condition number.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> ZXV input files:
!> line 1:  'ZXV' in columns 1 to 3.
!>
!> line 2:  N, INTEGER
!>          Value of N.
!>
!> line 3:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHGEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 4:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          Information will be printed about each test for which the
!>          test ratio is greater than or equal to the threshold.
!>
!> line 5:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 6:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 6 was 2:
!>
!> line 7: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> If line 2 was 0:
!>
!> line 7-EOF: Precomputed examples are tested.
!>
!> remaining lines : Each example is stored on 3+2*N*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer). The next N*N lines contain the matrix A, one
!>          element per line. The next N*N lines contain the matrix B.
!>          The next line contains the reciprocals of the eigenvalue
!>          condition numbers.  The last line contains the reciprocals of
!>          the eigenvector condition numbers.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> ZHB input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NK, INTEGER
!>          Number of values of K.
!>
!> line 5:  KVAL, INTEGER array, dimension (NK)
!>          The values for the matrix dimension K.
!>
!> line 6:  THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 8-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'ZHB'.
!>
!>-----------------------------------------------------------------------
!>
!> ZBB input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of M and N.
!>
!> line 3:  MVAL, INTEGER array, dimension (NN)
!>          The values for the matrix row dimension M.
!>
!> line 4:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix column dimension N.
!>
!> line 4:  NK, INTEGER
!>          Number of values of K.
!>
!> line 5:  KVAL, INTEGER array, dimension (NK)
!>          The values for the matrix bandwidth K.
!>
!> line 6:  NPARMS, INTEGER
!>          Number of values of the parameter NRHS
!>
!> line 7:  NSVAL, INTEGER array, dimension (NPARMS)
!>          The values for the number of right hand sides NRHS.
!>
!> line 8:  THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 9:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 9 was 2:
!>
!> line 10: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 10-EOF:  Lines specifying matrix types, as for SVD.
!>          The 3-character path name is 'ZBB'.
!>
!>-----------------------------------------------------------------------
!>
!> ZEC input file:
!>
!> line  2: THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> lines  3-EOF:
!>
!> Input for testing the eigencondition routines consists of a set of
!> specially constructed test cases and their solutions.  The data
!> format is not intended to be modified by the user.
!>
!>-----------------------------------------------------------------------
!>
!> ZBL and ZBK input files:
!>
!> line 1:  'ZBL' in columns 1-3 to test CGEBAL, or 'ZBK' in
!>          columns 1-3 to test CGEBAK.
!>
!> The remaining lines consist of specially constructed test cases.
!>
!>-----------------------------------------------------------------------
!>
!> ZGL and ZGK input files:
!>
!> line 1:  'ZGL' in columns 1-3 to test ZGGBAL, or 'ZGK' in
!>          columns 1-3 to test ZGGBAK.
!>
!> The remaining lines consist of specially constructed test cases.
!>
!>-----------------------------------------------------------------------
!>
!> GLM data file:
!>
!> line 1:  'GLM' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M (row dimension).
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P (row dimension).
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N (column dimension), note M <= N <= M+P.
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GLM' for the generalized
!>          linear regression model routines.
!>
!>-----------------------------------------------------------------------
!>
!> GQR data file:
!>
!> line 1:  'GQR' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M.
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P.
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N.
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GQR' for the generalized
!>          QR and RQ routines.
!>
!>-----------------------------------------------------------------------
!>
!> GSV data file:
!>
!> line 1:  'GSV' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M (row dimension).
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P (row dimension).
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N (column dimension).
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GSV' for the generalized
!>          SVD routines.
!>
!>-----------------------------------------------------------------------
!>
!> CSD data file:
!>
!> line 1:  'CSD' in columns 1 to 3.
!>
!> line 2:  NM, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NM)
!>          Values of M (row and column dimension of orthogonal matrix).
!>
!> line 4:  PVAL, INTEGER array, dimension(NM)
!>          Values of P (row dimension of top-left block).
!>
!> line 5:  NVAL, INTEGER array, dimension(NM)
!>          Values of N (column dimension of top-left block).
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'CSD' for the CSD routine.
!>
!>-----------------------------------------------------------------------
!>
!> LSE data file:
!>
!> line 1:  'LSE' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M.
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P.
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N, note P <= N <= P+M.
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GSV' for the generalized
!>          SVD routines.
!>
!>-----------------------------------------------------------------------
!>
!> NMAX is currently set to 132 and must be at least 12 for some of the
!> precomputed examples, and LWORK = NMAX*(5*NMAX+20) in the parameter
!> statements below.  For SVD, we assume NRHS may be as big as N.  The
!> parameter NEED is set to 14 to allow for 14 N-by-N matrices for ZGG.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_eig
!
!  =====================================================================
   PROGRAM ZCHKEE
!
#if defined(_OPENMP)
   use omp_lib
#endif
!
!  -- LAPACK test routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            NMAX
   PARAMETER          ( NMAX = 132 )
   INTEGER            NCMAX
   PARAMETER          ( NCMAX = 20 )
   INTEGER            NEED
   PARAMETER          ( NEED = 14 )
   INTEGER            LWORK
   PARAMETER          ( LWORK = NMAX*( 5*NMAX+20 ) )
   INTEGER            LIWORK
   PARAMETER          ( LIWORK = NMAX*( NMAX+20 ) )
   INTEGER            MAXIN
   PARAMETER          ( MAXIN = 20 )
   INTEGER            MAXT
   PARAMETER          ( MAXT = 30 )
   INTEGER            NIN, NOUT
   PARAMETER          ( NIN = 5, NOUT = 6 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ZBK, ZBL, ZES, ZEV, ZGK, ZGL, ZGS, ZGV, ZGX, &
                      ZSX, ZVX, ZXV, CSD, FATAL, GLM, GQR, GSV, LSE, &
                      NEP, SEP, SVD, TSTCHK, TSTDIF, TSTDRV, TSTERR, &
                      ZBB, ZGG, ZHB
   CHARACTER          C1
   CHARACTER*3        C3, PATH
   CHARACTER*32       VNAME
   CHARACTER*10       INTSTR
   CHARACTER*80       LINE
   INTEGER            I, I1, IC, INFO, ITMP, K, LENP, MAXTYP, NEWSD, &
                      NK, NN, NPARMS, NRHS, NTYPES, &
                      VERS_MAJOR, VERS_MINOR, VERS_PATCH
   INTEGER*4          N_THREADS, ONE_THREAD
   DOUBLE PRECISION   EPS, THRESH, THRSHN
   INTEGER(8)         nb_periods_sec, S1, S2, S1T, S2T
   REAL               STOT
!     ..
!     .. Local Arrays ..
   LOGICAL            DOTYPE( MAXT ), LOGWRK( NMAX )
   INTEGER            IOLDSD( 4 ), ISEED( 4 ), IWORK( LIWORK ), &
                      KVAL( MAXIN ), MVAL( MAXIN ), MXBVAL( MAXIN ), &
                      NBCOL( MAXIN ), NBMIN( MAXIN ), NBVAL( MAXIN ), &
                      NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), &
                      PVAL( MAXIN )
   INTEGER            INMIN( MAXIN ), INWIN( MAXIN ), INIBL( MAXIN ), &
                      ISHFTS( MAXIN ), IACC22( MAXIN )
   DOUBLE PRECISION   ALPHA( NMAX ), BETA( NMAX ), DR( NMAX, 12 ), &
                      RESULT( 500 )
   COMPLEX*16         DC( NMAX, 6 ), TAUA( NMAX ), TAUB( NMAX ), &
                      X( 5*NMAX )
!     ..
!     .. Allocatable Arrays ..
   INTEGER AllocateStatus
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK, S
   COMPLEX*16, DIMENSION(:), ALLOCATABLE :: WORK
   COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: A, B, C
!     ..
!     .. External Functions ..
   LOGICAL            LSAMEN
   DOUBLE PRECISION   DLAMCH, DSECND
   EXTERNAL           LSAMEN, DLAMCH, DSECND
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAREQ, XLAENV, ZCHKBB, ZCHKBD, ZCHKBK, ZCHKBL, &
                      ZCHKEC, ZCHKGG, ZCHKGK, ZCHKGL, ZCHKHB, ZCHKHS, &
                      ZCHKST, ZCKCSD, ZCKGLM, ZCKGQR, ZCKGSV, ZCKLSE, &
                      ZDRGES, ZDRGEV, ZDRGSX, ZDRGVX, ZDRVBD, ZDRVES, &
                      ZDRVEV, ZDRVSG, ZDRVST, ZDRVSX, ZDRVVX, &
                      ZERRBD, ZERRED, ZERRGG, ZERRHS, ZERRST, ILAVER, &
                      ZDRGES3, ZDRGEV3, &
                      ZCHKST2STG, ZDRVST2STG, ZCHKHB2STG
!     ..
!     .. Scalars in Common ..
   LOGICAL            LERR, OK
   CHARACTER*32       SRNAMT
   INTEGER            INFOT, MAXB, NPROC, NSHIFT, NUNIT, SELDIM, &
                      SELOPT
!     ..
!     .. Arrays in Common ..
   LOGICAL            SELVAL( 20 )
   INTEGER            IPARMS( 100 )
   DOUBLE PRECISION   SELWI( 20 ), SELWR( 20 )
!     ..
!     .. Common blocks ..
   COMMON             / CENVIR / NPROC, NSHIFT, MAXB
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
   COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
   COMMON             / CLAENV / IPARMS
!     ..
!     .. Data statements ..
   DATA               INTSTR / '0123456789' /
   DATA               IOLDSD / 0, 0, 0, 1 /
!     ..
!     .. Allocate memory dynamically ..
!
   ALLOCATE ( S(NMAX*NMAX), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( A(NMAX*NMAX,NEED), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( B(NMAX*NMAX,5), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( C(NCMAX*NCMAX,NCMAX*NCMAX), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( RWORK(LWORK), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( WORK(LWORK), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!     ..
!     .. Executable Statements ..
!
   A = 0.0
   B = 0.0
   C = 0.0
   DC = 0.0
   call system_clock(count_rate=nb_periods_sec,count=S1T)
   FATAL = .FALSE.
   NUNIT = NOUT
!
!     Return to here to read multiple sets of data
!
   DO
!
!     Read the first line and set the 3-character test path
!
   READ( NIN, FMT = '(A80)', END = 380 )LINE
   PATH = LINE( 1: 3 )
   NEP = LSAMEN( 3, PATH, 'NEP' ) .OR. LSAMEN( 3, PATH, 'ZHS' )
   SEP = LSAMEN( 3, PATH, 'SEP' ) .OR. LSAMEN( 3, PATH, 'ZST' ) .OR. &
         LSAMEN( 3, PATH, 'ZSG' ) .OR. LSAMEN( 3, PATH, 'SE2' )
   SVD = LSAMEN( 3, PATH, 'SVD' ) .OR. LSAMEN( 3, PATH, 'ZBD' )
   ZEV = LSAMEN( 3, PATH, 'ZEV' )
   ZES = LSAMEN( 3, PATH, 'ZES' )
   ZVX = LSAMEN( 3, PATH, 'ZVX' )
   ZSX = LSAMEN( 3, PATH, 'ZSX' )
   ZGG = LSAMEN( 3, PATH, 'ZGG' )
   ZGS = LSAMEN( 3, PATH, 'ZGS' )
   ZGX = LSAMEN( 3, PATH, 'ZGX' )
   ZGV = LSAMEN( 3, PATH, 'ZGV' )
   ZXV = LSAMEN( 3, PATH, 'ZXV' )
   ZHB = LSAMEN( 3, PATH, 'ZHB' )
   ZBB = LSAMEN( 3, PATH, 'ZBB' )
   GLM = LSAMEN( 3, PATH, 'GLM' )
   GQR = LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GRQ' )
   GSV = LSAMEN( 3, PATH, 'GSV' )
   CSD = LSAMEN( 3, PATH, 'CSD' )
   LSE = LSAMEN( 3, PATH, 'LSE' )
   ZBL = LSAMEN( 3, PATH, 'ZBL' )
   ZBK = LSAMEN( 3, PATH, 'ZBK' )
   ZGL = LSAMEN( 3, PATH, 'ZGL' )
   ZGK = LSAMEN( 3, PATH, 'ZGK' )
!
!     Report values of parameters.
!
   IF( PATH == '   ' ) CYCLE
   IF( NEP ) THEN
      WRITE( NOUT, FMT = 9987 )
   ELSE IF( SEP ) THEN
      WRITE( NOUT, FMT = 9986 )
   ELSE IF( SVD ) THEN
      WRITE( NOUT, FMT = 9985 )
   ELSE IF( ZEV ) THEN
      WRITE( NOUT, FMT = 9979 )
   ELSE IF( ZES ) THEN
      WRITE( NOUT, FMT = 9978 )
   ELSE IF( ZVX ) THEN
      WRITE( NOUT, FMT = 9977 )
   ELSE IF( ZSX ) THEN
      WRITE( NOUT, FMT = 9976 )
   ELSE IF( ZGG ) THEN
      WRITE( NOUT, FMT = 9975 )
   ELSE IF( ZGS ) THEN
      WRITE( NOUT, FMT = 9964 )
   ELSE IF( ZGX ) THEN
      WRITE( NOUT, FMT = 9965 )
   ELSE IF( ZGV ) THEN
      WRITE( NOUT, FMT = 9963 )
   ELSE IF( ZXV ) THEN
      WRITE( NOUT, FMT = 9962 )
   ELSE IF( ZHB ) THEN
      WRITE( NOUT, FMT = 9974 )
   ELSE IF( ZBB ) THEN
      WRITE( NOUT, FMT = 9967 )
   ELSE IF( GLM ) THEN
      WRITE( NOUT, FMT = 9971 )
   ELSE IF( GQR ) THEN
      WRITE( NOUT, FMT = 9970 )
   ELSE IF( GSV ) THEN
      WRITE( NOUT, FMT = 9969 )
   ELSE IF( CSD ) THEN
      WRITE( NOUT, FMT = 9960 )
   ELSE IF( LSE ) THEN
      WRITE( NOUT, FMT = 9968 )
   ELSE IF( ZBL ) THEN
!
!        ZGEBAL:  Balancing
!
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCHKBL( NIN, NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCHKBL : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      CYCLE
   ELSE IF( ZBK ) THEN
!
!        ZGEBAK:  Back transformation
!
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCHKBK( NIN, NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCHKBK : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      CYCLE
   ELSE IF( ZGL ) THEN
!
!        ZGGBAL:  Balancing
!
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCHKGL( NIN, NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCHKGL : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      CYCLE
   ELSE IF( ZGK ) THEN
!
!        ZGGBAK:  Back transformation
!
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCHKGK( NIN, NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCHKGK : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      CYCLE
   ELSE IF( LSAMEN( 3, PATH, 'ZEC' ) ) THEN
!
!        ZEC:  Eigencondition estimation
!
      READ( NIN, FMT = * )THRESH
      CALL XLAENV( 1, 1 )
      CALL XLAENV( 12, 1 )
      TSTERR = .TRUE.
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCHKEC( THRESH, TSTERR, NIN, NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCHKEC : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      CYCLE
   ELSE
      WRITE( NOUT, FMT = 9992 )PATH
      CYCLE
   END IF
   CALL ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
   WRITE( NOUT, FMT = 9972 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH
   WRITE( NOUT, FMT = 9984 )
!
!     Read the number of values of M, P, and N.
!
   READ(NIN,*)NN
   IF( NN < 0 ) THEN
      WRITE( NOUT, FMT = 9989 )'   NN ', NN, 1
      NN = 0
      FATAL = .TRUE.
   ELSE IF( NN > MAXIN ) THEN
      WRITE( NOUT, FMT = 9988 )'   NN ', NN, MAXIN
      NN = 0
      FATAL = .TRUE.
   END IF
!
!     Read the values of M
!
   IF( .NOT.( ZGX .OR. ZXV ) ) THEN
      READ(NIN,*) MVAL(1:NN)
      IF( SVD ) THEN
         VNAME = '    M '
      ELSE
         VNAME = '    N '
      END IF
      DO I = 1, NN
         IF( MVAL( I ) < 0 ) THEN
            WRITE( NOUT, FMT = 9989 )VNAME, MVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( MVAL( I ) > NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )VNAME, MVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
      ENDDO
      WRITE( NOUT, FMT = 9983 )'M:    ', MVAL(1:NN)
   END IF
!
!     Read the values of P
!
   IF( GLM .OR. GQR .OR. GSV .OR. CSD .OR. LSE ) THEN
      READ(NIN,*) PVAL(1:NN)
      DO I = 1, NN
         IF( PVAL( I ) < 0 ) THEN
            WRITE( NOUT, FMT = 9989 )' P  ', PVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( PVAL( I ) > NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )' P  ', PVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
      ENDDO
      WRITE( NOUT, FMT = 9983 )'P:    ', PVAL(1:NN)
   END IF
!
!     Read the values of N
!
   IF( SVD .OR. ZBB .OR. GLM .OR. GQR .OR. GSV .OR. CSD .OR. &
       LSE ) THEN
      READ(NIN,*) NVAL(1:NN)
      DO I = 1, NN
         IF( NVAL( I ) < 0 ) THEN
            WRITE( NOUT, FMT = 9989 )'    N ', NVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( NVAL( I ) > NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )'    N ', NVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
      ENDDO
   ELSE
      NVAL(1:NN) = MVAL(1:NN)
   END IF
   IF( .NOT.( ZGX .OR. ZXV ) ) THEN
      WRITE( NOUT, FMT = 9983 )'N:    ', NVAL(1:NN)
   ELSE
      WRITE( NOUT, FMT = 9983 )'N:    ', NN
   END IF
!
!     Read the number of values of K, followed by the values of K
!
   IF( ZHB .OR. ZBB ) THEN
      READ(NIN,*) NK
      READ(NIN,*) KVAL(1:NK)
      DO I = 1, NK
         IF( KVAL( I ) < 0 ) THEN
            WRITE( NOUT, FMT = 9989 )'    K ', KVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( KVAL( I ) > NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )'    K ', KVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
      ENDDO
      WRITE( NOUT, FMT = 9983 )'K:    ', KVAL(1:NK)
   END IF
!
   IF( ZEV .OR. ZES .OR. ZVX .OR. ZSX ) THEN
!
!        For the nonsymmetric QR driver routines, only one set of
!        parameters is allowed.
!
      READ(NIN,*) NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), &
         INMIN( 1 ), INWIN( 1 ), INIBL(1), ISHFTS(1), IACC22(1)
      IF( NBVAL( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( NBMIN( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( NXVAL( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( INMIN( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   INMIN ', INMIN( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( INWIN( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   INWIN ', INWIN( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( INIBL( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   INIBL ', INIBL( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( ISHFTS( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   ISHFTS ', ISHFTS( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( IACC22( 1 ) < 0 ) THEN
         WRITE( NOUT, FMT = 9989 )'   IACC22 ', IACC22( 1 ), 0
         FATAL = .TRUE.
      END IF
      CALL XLAENV( 1, NBVAL( 1 ) )
      CALL XLAENV( 2, NBMIN( 1 ) )
      CALL XLAENV( 3, NXVAL( 1 ) )
      CALL XLAENV(12, MAX( 11, INMIN( 1 ) ) )
      CALL XLAENV(13, INWIN( 1 ) )
      CALL XLAENV(14, INIBL( 1 ) )
      CALL XLAENV(15, ISHFTS( 1 ) )
      CALL XLAENV(16, IACC22( 1 ) )
      WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 )
      WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 )
      WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 )
      WRITE( NOUT, FMT = 9983 )'INMIN:   ', INMIN( 1 )
      WRITE( NOUT, FMT = 9983 )'INWIN: ', INWIN( 1 )
      WRITE( NOUT, FMT = 9983 )'INIBL: ', INIBL( 1 )
      WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ISHFTS( 1 )
      WRITE( NOUT, FMT = 9983 )'IACC22: ', IACC22( 1 )
!
   ELSE IF( ZGS .OR. ZGX .OR. ZGV .OR. ZXV ) THEN
!
!        For the nonsymmetric generalized driver routines, only one set of
!        parameters is allowed.
!
      READ(NIN,*) NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), NSVAL( 1 ), MXBVAL( 1 )
      IF( NBVAL( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( NBMIN( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( NXVAL( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1
         FATAL = .TRUE.
      ELSE IF( NSVAL( 1 ) < 2 ) THEN
         WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( 1 ), 2
         FATAL = .TRUE.
      ELSE IF( MXBVAL( 1 ) < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( 1 ), 1
         FATAL = .TRUE.
      END IF
      CALL XLAENV( 1, NBVAL( 1 ) )
      CALL XLAENV( 2, NBMIN( 1 ) )
      CALL XLAENV( 3, NXVAL( 1 ) )
      CALL XLAENV( 4, NSVAL( 1 ) )
      CALL XLAENV( 8, MXBVAL( 1 ) )
      WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 )
      WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 )
      WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 )
      WRITE( NOUT, FMT = 9983 )'NS:   ', NSVAL( 1 )
      WRITE( NOUT, FMT = 9983 )'MAXB: ', MXBVAL( 1 )
   ELSE IF( .NOT.ZHB .AND. .NOT.GLM .AND. .NOT.GQR .AND. .NOT. &
            GSV .AND. .NOT.CSD .AND. .NOT.LSE ) THEN
!
!        For the other paths, the number of parameters can be varied
!        from the input file.  Read the number of parameter values.
!
      READ(NIN,*) NPARMS
      IF( NPARMS < 1 ) THEN
         WRITE( NOUT, FMT = 9989 )'NPARMS', NPARMS, 1
         NPARMS = 0
         FATAL = .TRUE.
      ELSE IF( NPARMS > MAXIN ) THEN
         WRITE( NOUT, FMT = 9988 )'NPARMS', NPARMS, MAXIN
         NPARMS = 0
         FATAL = .TRUE.
      END IF
!
!        Read the values of NB
!
      IF( .NOT.ZBB ) THEN
         READ(NIN,*) NBVAL(1:NPARMS)
         DO I = 1, NPARMS
            IF( NBVAL( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( NBVAL( I ) > NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )'   NB ', NBVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
         ENDDO
         WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL(1:NPARMS)
      END IF
!
!        Read the values of NBMIN
!
      IF( NEP .OR. SEP .OR. SVD .OR. ZGG ) THEN
         READ(NIN,*) NBMIN(1:NPARMS)
         DO I = 1, NPARMS
            IF( NBMIN( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( I ), 0
               FATAL = .TRUE.
            ELSE IF( NBMIN( I ) > NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )'NBMIN ', NBMIN( I ), NMAX
               FATAL = .TRUE.
            END IF
         ENDDO
         WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN(1:NPARMS)
      ELSE
         NBMIN(1:NPARMS) = 1
      END IF
!
!        Read the values of NX
!
      IF( NEP .OR. SEP .OR. SVD ) THEN
         READ(NIN,*) NXVAL(1:NPARMS)
         DO I = 1, NPARMS
            IF( NXVAL( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( NXVAL( I ) > NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )'   NX ', NXVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL(1:NPARMS)
      ELSE
         NXVAL(1:NPARMS) = 1
      END IF
!
!        Read the values of NSHIFT (if ZGG) or NRHS (if SVD
!        or ZBB).
!
      IF( SVD .OR. ZBB .OR. ZGG ) THEN
         READ(NIN,*) NSVAL(1:NPARMS)
         DO I = 1, NPARMS
            IF( NSVAL( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( NSVAL( I ) > NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )'   NS ', NSVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'NS:   ', NSVAL(1:NPARMS)
      ELSE
         NSVAL(1:NPARMS) = 1
      END IF
!
!        Read the values for MAXB.
!
      IF( ZGG ) THEN
         READ(NIN,*) MXBVAL(1:NPARMS)
         DO I = 1, NPARMS
            IF( MXBVAL( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( MXBVAL( I ) > NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )' MAXB ', MXBVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'MAXB: ', MXBVAL(1:NPARMS)
      ELSE
         MXBVAL(1:NPARMS) = 1
      END IF
!
!        Read the values for INMIN.
!
      IF( NEP ) THEN
         READ(NIN,*) INMIN(1:NPARMS)
         DO I = 1, NPARMS
            IF( INMIN( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )' INMIN ', INMIN( I ), 0
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'INMIN: ', INMIN(1:NPARMS)
      ELSE
         INMIN(1:NPARMS) = 1
      END IF
!
!        Read the values for INWIN.
!
      IF( NEP ) THEN
         READ(NIN,*) INWIN(1:NPARMS)
         DO I = 1, NPARMS
            IF( INWIN( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )' INWIN ', INWIN( I ), 0
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'INWIN: ', INWIN(1:NPARMS)
      ELSE
         INWIN(1:NPARMS) = 1
      END IF
!
!        Read the values for INIBL.
!
      IF( NEP ) THEN
         READ(NIN,*) INIBL(1:NPARMS)
         DO I = 1, NPARMS
            IF( INIBL( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )' INIBL ', INIBL( I ), 0
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'INIBL: ', INIBL(1:NPARMS)
      ELSE
         INIBL(1:NPARMS) = 1
      END IF
!
!        Read the values for ISHFTS.
!
      IF( NEP ) THEN
         READ(NIN,*) ISHFTS(1:NPARMS)
         DO I = 1, NPARMS
            IF( ISHFTS( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )' ISHFTS ', ISHFTS( I ), 0
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ISHFTS(1:NPARMS)
      ELSE
         ISHFTS(1:NPARMS) = 1
      END IF
!
!        Read the values for IACC22.
!
      IF( NEP .OR. ZGG ) THEN
         READ(NIN,*) IACC22(1:NPARMS)
         DO I = 1, NPARMS
            IF( IACC22( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )' IACC22 ', IACC22( I ), 0
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'IACC22: ', IACC22(1:NPARMS)
      ELSE
         IACC22(1:NPARMS) = 1
      END IF
!
!        Read the values for NBCOL.
!
      IF( ZGG ) THEN
         READ(NIN,*) NBCOL(1:NPARMS)
         DO I = 1, NPARMS
            IF( NBCOL( I ) < 0 ) THEN
               WRITE( NOUT, FMT = 9989 )'NBCOL ', NBCOL( I ), 0
               FATAL = .TRUE.
            ELSE IF( NBCOL( I ) > NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )'NBCOL ', NBCOL( I ), NMAX
               FATAL = .TRUE.
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9983 )'NBCOL:', NBCOL(1:NPARMS)
      ELSE
         NBCOL(1:NPARMS) = 1
      END IF
   END IF
!
!     Calculate and print the machine dependent constants.
!
   WRITE( NOUT, FMT = * )
   EPS = DLAMCH( 'Underflow threshold' )
   WRITE( NOUT, FMT = 9981 )'underflow', EPS
   EPS = DLAMCH( 'Overflow threshold' )
   WRITE( NOUT, FMT = 9981 )'overflow ', EPS
   EPS = DLAMCH( 'Epsilon' )
   WRITE( NOUT, FMT = 9981 )'precision', EPS
!
!     Read the threshold value for the test ratios.
!
   READ(NIN,*) THRESH
   WRITE( NOUT, FMT = 9982 )THRESH
   IF( SEP .OR. SVD .OR. ZGG ) THEN
!
!        Read the flag that indicates whether to test LAPACK routines.
!
      READ(NIN,*) TSTCHK
!
!        Read the flag that indicates whether to test driver routines.
!
      READ(NIN,*) TSTDRV
   END IF
!
!     Read the flag that indicates whether to test the error exits.
!
   READ(NIN,*) TSTERR
!
!     Read the code describing how to set the random number seed.
!
   READ(NIN,*) NEWSD
!
!     If NEWSD = 2, read another line with 4 integers for the seed.
!
   IF( NEWSD == 2 ) READ(NIN,*) IOLDSD(1:4)
!
   ISEED(1:4) = IOLDSD(1:4)
!
   IF( FATAL ) THEN
      WRITE( NOUT, FMT = 9999 )
      STOP
   END IF
!
!     Read the input lines indicating the test path and its parameters.
!     The first three characters indicate the test path, and the number
!     of test matrix types must be the first nonblank item in columns
!     4-80.
!
  190 CONTINUE
!
   IF( .NOT.( ZGX .OR. ZXV ) ) THEN
!
  200    CONTINUE
      READ( NIN, FMT = '(A80)', END = 380 )LINE
      C3 = LINE( 1: 3 )
      LENP = LEN( LINE )
      I = 3
      ITMP = 0
      I1 = 0
  210    CONTINUE
      I = I + 1
      IF( I > LENP ) THEN
         IF( I1 > 0 ) THEN
            GO TO 240
         ELSE
            NTYPES = MAXT
            GO TO 240
         END IF
      END IF
      IF( LINE( I: I ) /= ' ' .AND. LINE( I: I ) /= ',' ) THEN
         I1 = I
         C1 = LINE( I1: I1 )
!
!        Check that a valid integer was read
!
         DO K = 1, 10
            IF( C1 == INTSTR( K: K ) ) THEN
               IC = K - 1
               GO TO 230
            END IF
            ENDDO
         WRITE( NOUT, FMT = 9991 )I, LINE
         GO TO 200
  230       CONTINUE
         ITMP = 10*ITMP + IC
         GO TO 210
      ELSE IF( I1 > 0 ) THEN
         GO TO 240
      ELSE
         GO TO 210
      END IF
  240    CONTINUE
      NTYPES = ITMP
!
!     Skip the tests if NTYPES is <= 0.
!
      IF( .NOT.( ZEV .OR. ZES .OR. ZVX .OR. ZSX .OR. ZGV .OR. &
          ZGS ) .AND. NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
         GO TO 200
      END IF
!
   ELSE
      IF( ZGX ) C3 = 'ZGX'
      IF( ZXV ) C3 = 'ZXV'
   END IF
!
!     Reset the random number seed.
!
   IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
!
   IF( LSAMEN( 3, C3, 'ZHS' ) .OR. LSAMEN( 3, C3, 'NEP' ) ) THEN
!
!        -------------------------------------
!        NEP:  Nonsymmetric Eigenvalue Problem
!        -------------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!           NS    = number of shifts
!           MAXB  = minimum submatrix size
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) THEN
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZERRHS( 'ZHSEQR', NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZERRHS : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
    ENDIF
      DO I = 1, NPARMS
         CALL XLAENV( 1, NBVAL( I ) )
         CALL XLAENV( 2, NBMIN( I ) )
         CALL XLAENV( 3, NXVAL( I ) )
         CALL XLAENV(12, MAX( 11, INMIN( I ) ) )
         CALL XLAENV(13, INWIN( I ) )
         CALL XLAENV(14, INIBL( I ) )
         CALL XLAENV(15, ISHFTS( I ) )
         CALL XLAENV(16, IACC22( I ) )
!
         IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
         WRITE( NOUT, FMT = 9961 )C3, NBVAL( I ), NBMIN( I ), &
            NXVAL( I ), MAX( 11, INMIN(I)), &
            INWIN( I ), INIBL( I ), ISHFTS( I ), IACC22( I )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKHS( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), &
                      A( 1, 7 ), DC( 1, 1 ), DC( 1, 2 ), A( 1, 8 ), &
                      A( 1, 9 ), A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), &
                      DC( 1, 3 ), WORK, LWORK, RWORK, IWORK, LOGWRK, &
                      RESULT, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKHS : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCHKHS', INFO
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'ZST' ) .OR. LSAMEN( 3, C3, 'SEP' ) &
                                   .OR. LSAMEN( 3, C3, 'SE2' ) ) THEN
!
!        ----------------------------------
!        SEP:  Symmetric Eigenvalue Problem
!        ----------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      CALL XLAENV( 1, 1 )
      CALL XLAENV( 9, 25 )
      IF( TSTERR ) THEN
#if defined(_OPENMP)
         N_THREADS = OMP_GET_MAX_THREADS()
         ONE_THREAD = 1
         CALL OMP_SET_NUM_THREADS(ONE_THREAD)
#endif
         CALL ZERRST( 'ZST', NOUT )
#if defined(_OPENMP)
         CALL OMP_SET_NUM_THREADS(N_THREADS)
#endif
      END IF
      DO I = 1, NPARMS
         CALL XLAENV( 1, NBVAL( I ) )
         CALL XLAENV( 2, NBMIN( I ) )
         CALL XLAENV( 3, NXVAL( I ) )
!
         IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
         WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I )
         IF( TSTCHK ) THEN
            IF( LSAMEN( 3, C3, 'SE2' ) ) THEN
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZCHKST2STG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                         NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), &
                         DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), &
                         DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), &
                         DR( 1, 7 ), DR( 1, 8 ), DR( 1, 9 ), &
                         DR( 1, 10 ), DR( 1, 11 ), A( 1, 3 ), NMAX, &
                         A( 1, 4 ), A( 1, 5 ), DC( 1, 1 ), A( 1, 6 ), &
                         WORK, LWORK, RWORK, LWORK, IWORK, LIWORK, &
                         RESULT, INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCHKST2STG : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
            ELSE
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZCHKST( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                         NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), &
                         DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), &
                         DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), &
                         DR( 1, 7 ), DR( 1, 8 ), DR( 1, 9 ), &
                         DR( 1, 10 ), DR( 1, 11 ), A( 1, 3 ), NMAX, &
                         A( 1, 4 ), A( 1, 5 ), DC( 1, 1 ), A( 1, 6 ), &
                         WORK, LWORK, RWORK, LWORK, IWORK, LIWORK, &
                         RESULT, INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCHKST : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
            ENDIF
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCHKST', INFO
         END IF
         IF( TSTDRV ) THEN
            IF( LSAMEN( 3, C3, 'SE2' ) ) THEN
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZDRVST2STG( NN, NVAL, 18, DOTYPE, ISEED, THRESH, &
                       NOUT, A( 1, 1 ), NMAX, DR( 1, 3 ), DR( 1, 4 ), &
                       DR( 1, 5 ), DR( 1, 8 ), DR( 1, 9 ), &
                       DR( 1, 10 ), A( 1, 2 ), NMAX, A( 1, 3 ), &
                       DC( 1, 1 ), A( 1, 4 ), WORK, LWORK, RWORK, &
                       LWORK, IWORK, LIWORK, RESULT, INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZDRVST2STG : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
        ELSE
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZDRVST( NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, &
                       A( 1, 1 ), NMAX, DR( 1, 3 ), DR( 1, 4 ), &
                       DR( 1, 5 ), DR( 1, 8 ), DR( 1, 9 ), &
                       DR( 1, 10 ), A( 1, 2 ), NMAX, A( 1, 3 ), &
                       DC( 1, 1 ), A( 1, 4 ), WORK, LWORK, RWORK, &
                       LWORK, IWORK, LIWORK, RESULT, INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZDRVST : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
            ENDIF
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZDRVST', INFO
         END IF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'ZSG' ) ) THEN
!
!        ----------------------------------------------
!        ZSG:  Hermitian Generalized Eigenvalue Problem
!        ----------------------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      CALL XLAENV( 9, 25 )
      DO I = 1, NPARMS
         CALL XLAENV( 1, NBVAL( I ) )
         CALL XLAENV( 2, NBMIN( I ) )
         CALL XLAENV( 3, NXVAL( I ) )
!
         IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
         WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I )
         IF( TSTCHK ) THEN
!               CALL ZDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
!     $                      DR( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
!     $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
!     $                      LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT,
!     $                      INFO )
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZDRVSG2STG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                             NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, &
                             DR( 1, 3 ), DR( 1, 4 ), A( 1, 3 ), NMAX, &
                             A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), &
                             A( 1, 7 ), WORK, LWORK, RWORK, LWORK, &
                             IWORK, LIWORK, RESULT, INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZDRVSG2STG : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
            IF( INFO /= 0 ) &
               WRITE( NOUT, FMT = 9980 )'ZDRVSG', INFO
         END IF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'ZBD' ) .OR. LSAMEN( 3, C3, 'SVD' ) ) THEN
!
!        ----------------------------------
!        SVD:  Singular Value Decomposition
!        ----------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!           NRHS  = number of right hand sides
!
      MAXTYP = 16
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      CALL XLAENV( 9, 25 )
!
!        Test the error exits
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR .AND. TSTCHK ) CALL ZERRBD( 'ZBD', NOUT )
      IF( TSTERR .AND. TSTDRV ) CALL ZERRED( 'ZBD', NOUT )
!
      DO I = 1, NPARMS
         NRHS = NSVAL( I )
         CALL XLAENV( 1, NBVAL( I ) )
         CALL XLAENV( 2, NBMIN( I ) )
         CALL XLAENV( 3, NXVAL( I ) )
         IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
         WRITE( NOUT, FMT = 9995 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), NRHS
         IF( TSTCHK ) THEN
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZCHKBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, NRHS, ISEED, &
                         THRESH, A( 1, 1 ), NMAX, DR( 1, 1 ), &
                         DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), &
                         A( 1, 2 ), NMAX, A( 1, 3 ), A( 1, 4 ), &
                         A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), &
                         A( 1, 8 ), WORK, LWORK, RWORK, NOUT, INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCHKBD : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCHKBD', INFO
         END IF
         IF( TSTDRV ) THEN
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZDRVBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, ISEED, &
                         THRESH, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, &
                         A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), &
                         A( 1, 6 ), DR( 1, 1 ), DR( 1, 2 ), &
                         DR( 1, 3 ), WORK, LWORK, RWORK, IWORK, NOUT, &
                         INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZDRVBD : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
    ENDIF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'ZEV' ) ) THEN
!
!        --------------------------------------------
!        ZEV:  Nonsymmetric Eigenvalue Problem Driver
!              ZGEEV (eigenvalues and eigenvectors)
!        --------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL ZERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVEV( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), DC( 1, 1 ), &
                      DC( 1, 2 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, &
                      A( 1, 5 ), NMAX, RESULT, WORK, LWORK, RWORK, &
                      IWORK, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVEV : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZGEEV', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'ZES' ) ) THEN
!
!        --------------------------------------------
!        ZES:  Nonsymmetric Eigenvalue Problem Driver
!              ZGEES (Schur form)
!        --------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL ZERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVES( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      DC( 1, 1 ), DC( 1, 2 ), A( 1, 4 ), NMAX, &
                      RESULT, WORK, LWORK, RWORK, IWORK, LOGWRK, &
                      INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVES : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZGEES', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'ZVX' ) ) THEN
!
!        --------------------------------------------------------------
!        ZVX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              ZGEEVX (eigenvalues, eigenvectors and condition numbers)
!        --------------------------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL ZERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVVX( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, &
                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), DC( 1, 1 ), &
                      DC( 1, 2 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, &
                      A( 1, 5 ), NMAX, DR( 1, 1 ), DR( 1, 2 ), &
                      DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), &
                      DR( 1, 7 ), DR( 1, 8 ), RESULT, WORK, LWORK, &
                      RWORK, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVVX : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZGEEVX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'ZSX' ) ) THEN
!
!        ---------------------------------------------------
!        ZSX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              ZGEESX (Schur form and condition numbers)
!        ---------------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL ZERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRVSX( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, &
                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      DC( 1, 1 ), DC( 1, 2 ), DC( 1, 3 ), A( 1, 4 ), &
                      NMAX, A( 1, 5 ), RESULT, WORK, LWORK, RWORK, &
                      LOGWRK, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRVSX : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZGEESX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'ZGG' ) ) THEN
!
!        -------------------------------------------------
!        ZGG:  Generalized Nonsymmetric Eigenvalue Problem
!        -------------------------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NS    = number of shifts
!           MAXB  = minimum submatrix size
!           IACC22: structured matrix multiply
!           NBCOL = minimum column dimension for blocks
!
      MAXTYP = 26
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      CALL XLAENV(1,1)
      IF( TSTCHK .AND. TSTERR ) CALL ZERRGG( C3, NOUT )
      DO I = 1, NPARMS
         CALL XLAENV( 1, NBVAL( I ) )
         CALL XLAENV( 2, NBMIN( I ) )
         CALL XLAENV( 4, NSVAL( I ) )
         CALL XLAENV( 8, MXBVAL( I ) )
         CALL XLAENV( 16, IACC22( I ) )
         CALL XLAENV( 5, NBCOL( I ) )
!
         IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
         WRITE( NOUT, FMT = 9996 )C3, NBVAL( I ), NBMIN( I ), &
            NSVAL( I ), MXBVAL( I ), IACC22( I ), NBCOL( I )
         TSTDIF = .FALSE.
         THRSHN = 10.D0
         IF( TSTCHK ) THEN
            call system_clock(count_rate=nb_periods_sec,count=S1)
            CALL ZCHKGG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                         TSTDIF, THRSHN, NOUT, A( 1, 1 ), NMAX, &
                         A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                         A( 1, 6 ), A( 1, 7 ), A( 1, 8 ), A( 1, 9 ), &
                         NMAX, A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), &
                         DC( 1, 1 ), DC( 1, 2 ), DC( 1, 3 ), &
                         DC( 1, 4 ), A( 1, 13 ), A( 1, 14 ), WORK, &
                         LWORK, RWORK, LOGWRK, RESULT, INFO )
            call system_clock(count_rate=nb_periods_sec,count=S2)
            open(file='results.out', unit=10, position = 'append')
            write(10,'(A,F16.10,A)') 'Total time : ZCHKGG : ', &
                  real(S2-S1)/real(nb_periods_sec), ' s'
            close(10)
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCHKGG', INFO
         END IF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'ZGS' ) ) THEN
!
!        -------------------------------------------------
!        ZGS:  Generalized Nonsymmetric Eigenvalue Problem
!              ZGGES (Schur form)
!        -------------------------------------------------
!
      MAXTYP = 26
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL ZERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRGES( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                      DC( 1, 1 ), DC( 1, 2 ), WORK, LWORK, RWORK, &
                      RESULT, LOGWRK, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRGES : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZDRGES', INFO
!
!     Blocked version
!
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRGES3( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                       A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                       A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                       DC( 1, 1 ), DC( 1, 2 ), WORK, LWORK, RWORK, &
                       RESULT, LOGWRK, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRGES3 : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZDRGES3', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( ZGX ) THEN
!
!        -------------------------------------------------
!        ZGX  Generalized Nonsymmetric Eigenvalue Problem
!              ZGGESX (Schur form and condition numbers)
!        -------------------------------------------------
!
      MAXTYP = 5
      NTYPES = MAXTYP
      IF( NN < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL ZERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL XLAENV( 5, 2 )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRGSX( NN, NCMAX, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      A( 1, 6 ), DC( 1, 1 ), DC( 1, 2 ), C, &
                      NCMAX*NCMAX, S, WORK, LWORK, RWORK, IWORK, &
                      LIWORK, LOGWRK, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRGSX : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZDRGSX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'ZGV' ) ) THEN
!
!        -------------------------------------------------
!        ZGV:  Generalized Nonsymmetric Eigenvalue Problem
!              ZGGEV (Eigenvalue/vector form)
!        -------------------------------------------------
!
      MAXTYP = 26
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL ZERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRGEV( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                      A( 1, 9 ), NMAX, DC( 1, 1 ), DC( 1, 2 ), &
                      DC( 1, 3 ), DC( 1, 4 ), WORK, LWORK, RWORK, &
                      RESULT, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRGEV : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZDRGEV', INFO
!
! Blocked version
!
         CALL XLAENV(16,2)
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRGEV3( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                       A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                       A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                       A( 1, 9 ), NMAX, DC( 1, 1 ), DC( 1, 2 ), &
                       DC( 1, 3 ), DC( 1, 4 ), WORK, LWORK, RWORK, &
                       RESULT, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRGEV3 : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) &
            WRITE( NOUT, FMT = 9980 )'ZDRGEV3', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( ZXV ) THEN
!
!        -------------------------------------------------
!        ZXV:  Generalized Nonsymmetric Eigenvalue Problem
!              ZGGEVX (eigenvalue/vector with condition numbers)
!        -------------------------------------------------
!
      MAXTYP = 2
      NTYPES = MAXTYP
      IF( NN < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) &
            CALL ZERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZDRGVX( NN, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), DC( 1, 1 ), &
                      DC( 1, 2 ), A( 1, 5 ), A( 1, 6 ), IWORK( 1 ), &
                      IWORK( 2 ), DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), &
                      DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), WORK, &
                      LWORK, RWORK, IWORK( 3 ), LIWORK-2, RESULT, &
                      LOGWRK, INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZDRGVX : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
!
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZDRGVX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'ZHB' ) ) THEN
!
!        ------------------------------
!        ZHB:  Hermitian Band Reduction
!        ------------------------------
!
      MAXTYP = 15
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      IF( TSTERR ) THEN
#if defined(_OPENMP)
         N_THREADS = OMP_GET_MAX_THREADS()
         ONE_THREAD = 1
         CALL OMP_SET_NUM_THREADS(ONE_THREAD)
#endif
         CALL ZERRST( 'ZHB', NOUT )
#if defined(_OPENMP)
         CALL OMP_SET_NUM_THREADS(N_THREADS)
#endif
      END IF
!         CALL ZCHKHB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                NOUT, A( 1, 1 ), NMAX, DR( 1, 1 ), DR( 1, 2 ),
!     $                A( 1, 2 ), NMAX, WORK, LWORK, RWORK, RESULT,
!     $                INFO )
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCHKHB2STG( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, &
                    THRESH, NOUT, A( 1, 1 ), NMAX, DR( 1, 1 ), &
                    DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), &
                    A( 1, 2 ), NMAX, WORK, LWORK, RWORK, RESULT, &
                    INFO )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCHKHB2STG : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCHKHB', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'ZBB' ) ) THEN
!
!        ------------------------------
!        ZBB:  General Band Reduction
!        ------------------------------
!
      MAXTYP = 15
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      DO I = 1, NPARMS
         NRHS = NSVAL( I )
!
         IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
         WRITE( NOUT, FMT = 9966 )C3, NRHS
         call system_clock(count_rate=nb_periods_sec,count=S1)
         CALL ZCHKBB( NN, MVAL, NVAL, NK, KVAL, MAXTYP, DOTYPE, NRHS, &
                      ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, &
                      A( 1, 2 ), 2*NMAX, DR( 1, 1 ), DR( 1, 2 ), &
                      A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, A( 1, 6 ), &
                      NMAX, A( 1, 7 ), WORK, LWORK, RWORK, RESULT, &
                      INFO )
         call system_clock(count_rate=nb_periods_sec,count=S2)
         open(file='results.out', unit=10, position = 'append')
         write(10,'(A,F16.10,A)') 'Total time : ZCHKBB : ', &
               real(S2-S1)/real(nb_periods_sec), ' s'
         close(10)
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCHKBB', INFO
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'GLM' ) ) THEN
!
!        -----------------------------------------
!        GLM:  Generalized Linear Regression Model
!        -----------------------------------------
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL ZERRGG( 'GLM', NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCKGLM( NN, NVAL, MVAL, PVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, &
                   WORK, DR( 1, 1 ), NIN, NOUT, INFO )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCKGLM : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCKGLM', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'GQR' ) ) THEN
!
!        ------------------------------------------
!        GQR:  Generalized QR and RQ factorizations
!        ------------------------------------------
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL ZERRGG( 'GQR', NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCKGQR( NN, MVAL, NN, PVAL, NN, NVAL, NTYPES, ISEED, &
                   THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                   A( 1, 4 ), TAUA, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                   B( 1, 4 ), B( 1, 5 ), TAUB, WORK, DR( 1, 1 ), NIN, &
                   NOUT, INFO )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCKGQR : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCKGQR', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'GSV' ) ) THEN
!
!        ----------------------------------------------
!        GSV:  Generalized Singular Value Decomposition
!        ----------------------------------------------
!
      CALL XLAENV(1,1)
      IF( TSTERR ) CALL ZERRGG( 'GSV', NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCKGSV( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                   A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), ALPHA, BETA, &
                   B( 1, 4 ), IWORK, WORK, DR( 1, 1 ), NIN, NOUT, &
                   INFO )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCKGSV : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCKGSV', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'CSD' ) ) THEN
!
!        ----------------------------------------------
!        CSD:  CS Decomposition
!        ----------------------------------------------
!
      CALL XLAENV(1,1)
      IF( TSTERR ) CALL ZERRGG( 'CSD', NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCKCSD( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), &
                   A( 1, 5 ), A( 1, 6 ), RWORK, IWORK, WORK, &
                   DR( 1, 1 ), NIN, NOUT, INFO )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCKCSD : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCKCSD', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'LSE' ) ) THEN
!
!        --------------------------------------
!        LSE:  Constrained Linear Least Squares
!        --------------------------------------
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL ZERRGG( 'LSE', NOUT )
      call system_clock(count_rate=nb_periods_sec,count=S1)
      CALL ZCKLSE( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, &
                   WORK, DR( 1, 1 ), NIN, NOUT, INFO )
      call system_clock(count_rate=nb_periods_sec,count=S2)
      open(file='results.out', unit=10, position = 'append')
      write(10,'(A,F16.10,A)') 'Total time : ZCKLSE : ', &
            real(S2-S1)/real(nb_periods_sec), ' s'
      close(10)
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'ZCKLSE', INFO
   ELSE
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9992 )C3
   END IF
   ENDDO
   IF( .NOT.( ZGX .OR. ZXV ) ) GO TO 190
  380 CONTINUE
   WRITE( NOUT, FMT = 9994 )
   call system_clock(count_rate=nb_periods_sec,count=S2T)
   WRITE( NOUT, FMT = 9993 ) real(S2T - S1T)/real(nb_periods_sec)
!
   DEALLOCATE (S, STAT = AllocateStatus)
   DEALLOCATE (A, STAT = AllocateStatus)
   DEALLOCATE (B, STAT = AllocateStatus)
   DEALLOCATE (C, STAT = AllocateStatus)
   DEALLOCATE (RWORK, STAT = AllocateStatus)
   DEALLOCATE (WORK,  STAT = AllocateStatus)
!
 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9997 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4 )
 9996 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NS =', I4, &
         ', MAXB =', I4, ', IACC22 =', I4, ', NBCOL =', I4 )
 9995 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, &
         ', NRHS =', I4 )
 9994 FORMAT( / / ' End of tests' )
 9993 FORMAT( ' Total time used = ', F16.8, ' seconds', / )
 9992 FORMAT( 1X, A3, ':  Unrecognized path name' )
 9991 FORMAT( / / ' *** Invalid integer value in column ', I2, &
         ' of input', ' line:', / A79 )
 9990 FORMAT( / / 1X, A3, ' routines were not tested' )
 9989 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be >=', &
         I6 )
 9988 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be <=', &
         I6 )
 9987 FORMAT( ' Tests of the Nonsymmetric Eigenvalue Problem routines' )
 9986 FORMAT( ' Tests of the Hermitian Eigenvalue Problem routines' )
 9985 FORMAT( ' Tests of the Singular Value Decomposition routines' )
 9984 FORMAT( / ' The following parameter values will be used:' )
 9983 FORMAT( 4X, A, 10I6, / 10X, 10I6 )
 9982 FORMAT( / ' Routines pass computational tests if test ratio is ', &
         'less than', F8.2, / )
 9981 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )
 9980 FORMAT( ' *** Error code from ', A, ' = ', I4 )
 9979 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', &
         / '    ZGEEV (eigenvalues and eigevectors)' )
 9978 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', &
         / '    ZGEES (Schur form)' )
 9977 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', &
         ' Driver', / '    ZGEEVX (eigenvalues, eigenvectors and', &
         ' condition numbers)' )
 9976 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', &
         ' Driver', / '    ZGEESX (Schur form and condition', &
         ' numbers)' )
 9975 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem routines' )
 9974 FORMAT( ' Tests of ZHBTRD', / ' (reduction of a Hermitian band ', &
         'matrix to real tridiagonal form)' )
 9973 FORMAT( / 1X, 71( '-' ) )
 9972 FORMAT( / ' LAPACK VERSION ', I1, '.', I1, '.', I1 )
 9971 FORMAT( / ' Tests of the Generalized Linear Regression Model ', &
         'routines' )
 9970 FORMAT( / ' Tests of the Generalized QR and RQ routines' )
 9969 FORMAT( / ' Tests of the Generalized Singular Value', &
         ' Decomposition routines' )
 9968 FORMAT( / ' Tests of the Linear Least Squares routines' )
 9967 FORMAT( ' Tests of ZGBBRD', / ' (reduction of a general band ', &
         'matrix to real bidiagonal form)' )
 9966 FORMAT( / / 1X, A3, ':  NRHS =', I4 )
 9965 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Expert Driver ZGGESX' )
 9964 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Driver ZGGES' )
 9963 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Driver ZGGEV' )
 9962 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Expert Driver ZGGEVX' )
 9961 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, &
         ', INMIN=', I4, &
         ', INWIN =', I4, ', INIBL =', I4, ', ISHFTS =', I4, &
         ', IACC22 =', I4)
 9960 FORMAT( / ' Tests of the CS Decomposition routines' )
!
!     End of ZCHKEE
!
   END
