!> \brief \b SCHKEE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM SCHKEE
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKEE tests the REAL LAPACK subroutines for the matrix
!> eigenvalue problem.  The test paths in this version are
!>
!> NEP (Nonsymmetric Eigenvalue Problem):
!>     Test SGEHRD, SORGHR, SHSEQR, STREVC, SHSEIN, and SORMHR
!>
!> SEP (Symmetric Eigenvalue Problem):
!>     Test SSYTRD, SORGTR, SSTEQR, SSTERF, SSTEIN, SSTEDC,
!>     and drivers SSYEV(X), SSBEV(X), SSPEV(X), SSTEV(X),
!>                 SSYEVD,   SSBEVD,   SSPEVD,   SSTEVD
!>
!> SVD (Singular Value Decomposition):
!>     Test SGEBRD, SORGBR, SBDSQR, SBDSDC
!>     and the drivers SGESVD, SGESDD
!>
!> SEV (Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test SGEEV
!>
!> SES (Nonsymmetric Schur form Driver):
!>     Test SGEES
!>
!> SVX (Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test SGEEVX
!>
!> SSX (Nonsymmetric Schur form Expert Driver):
!>     Test SGEESX
!>
!> SGG (Generalized Nonsymmetric Eigenvalue Problem):
!>     Test SGGHD3, SGGBAL, SGGBAK, SHGEQZ, and STGEVC
!>
!> SGS (Generalized Nonsymmetric Schur form Driver):
!>     Test SGGES
!>
!> SGV (Generalized Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test SGGEV
!>
!> SGX (Generalized Nonsymmetric Schur form Expert Driver):
!>     Test SGGESX
!>
!> SXV (Generalized Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test SGGEVX
!>
!> SSG (Symmetric Generalized Eigenvalue Problem):
!>     Test SSYGST, SSYGV, SSYGVD, SSYGVX, SSPGST, SSPGV, SSPGVD,
!>     SSPGVX, SSBGST, SSBGV, SSBGVD, and SSBGVX
!>
!> SSB (Symmetric Band Eigenvalue Problem):
!>     Test SSBTRD
!>
!> SBB (Band Singular Value Decomposition):
!>     Test SGBBRD
!>
!> SEC (Eigencondition estimation):
!>     Test SLALN2, SLASY2, SLAEQU, SLAEXC, STRSYL, STREXC, STRSNA,
!>     STRSEN, and SLAQTR
!>
!> SBL (Balancing a general matrix)
!>     Test SGEBAL
!>
!> SBK (Back transformation on a balanced matrix)
!>     Test SGEBAK
!>
!> SGL (Balancing a matrix pair)
!>     Test SGGBAL
!>
!> SGK (Back transformation on a matrix pair)
!>     Test SGGBAK
!>
!> GLM (Generalized Linear Regression Model):
!>     Tests SGGGLM
!>
!> GQR (Generalized QR and RQ factorizations):
!>     Tests SGGQRF and SGGRQF
!>
!> GSV (Generalized Singular Value Decomposition):
!>     Tests SGGSVD, SGGSVP, STGSJA, SLAGS2, SLAPLL, and SLAPMT
!>
!> CSD (CS decomposition):
!>     Tests SORCSD
!>
!> LSE (Constrained Linear Least Squares):
!>     Tests SGGLSE
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
!> SHS or NEP      21     SCHKHS
!> SST or SEP      21     SCHKST (routines)
!>                 18     SDRVST (drivers)
!> SBD or SVD      16     SCHKBD (routines)
!>                  5     SDRVBD (drivers)
!> SEV             21     SDRVEV
!> SES             21     SDRVES
!> SVX             21     SDRVVX
!> SSX             21     SDRVSX
!> SGG             26     SCHKGG (routines)
!> SGS             26     SDRGES
!> SGX              5     SDRGSX
!> SGV             26     SDRGEV
!> SXV              2     SDRGVX
!> SSG             21     SDRVSG
!> SSB             15     SCHKSB
!> SBB             15     SCHKBB
!> SEC              -     SCHKEC
!> SBL              -     SCHKBL
!> SBK              -     SCHKBK
!> SGL              -     SCHKGL
!> SGK              -     SCHKGK
!> GLM              8     SCKGLM
!> GQR              8     SCKGQR
!> GSV              8     SCKGSV
!> CSD              3     SCKCSD
!> LSE              8     SCKLSE
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
!>          The valid 3-character path names are 'NEP' or 'SHS' for the
!>          nonsymmetric eigenvalue routines.
!>
!>-----------------------------------------------------------------------
!>
!> SEP or SSG input file:
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
!>          The 3-character path names are 'SEP' or 'SST' for the
!>          symmetric eigenvalue routines and driver routines, and
!>          'SSG' for the routines for the symmetric generalized
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
!>          The 3-character path names are 'SVD' or 'SBD' for both the
!>          SVD routines and the SVD driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> SEV and SES data files:
!>
!> line 1:  'SEV' or 'SES' in columns 1 to 3.
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
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9 and following:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'SEV' to test SGEEV, or
!>          'SES' to test SGEES.
!>
!>-----------------------------------------------------------------------
!>
!> The SVX data has two parts. The first part is identical to SEV,
!> and the second part consists of test matrices with precomputed
!> solutions.
!>
!> line 1:  'SVX' in columns 1-3.
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
!> line 6:  TSTERR, LOGICAL
!>
!> line 7:  NEWSD, INTEGER
!>
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>
!> lines 9 and following: The first line contains 'SVX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 1+2*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer). The next N lines contain the matrix, one
!>          row per line. The last N lines correspond to each
!>          eigenvalue. Each of these last N lines contains 4 real
!>          values: the real part of the eigenvalue, the imaginary
!>          part of the eigenvalue, the reciprocal condition number of
!>          the eigenvalues, and the reciprocal condition number of the
!>          eigenvector.  The end of data is indicated by dimension N=0.
!>          Even if no data is to be tested, there must be at least one
!>          line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> The SSX data is like SVX. The first part is identical to SEV, and the
!> second part consists of test matrices with precomputed solutions.
!>
!> line 1:  'SSX' in columns 1-3.
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
!> line 6:  TSTERR, LOGICAL
!>
!> line 7:  NEWSD, INTEGER
!>
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>
!> lines 9 and following: The first line contains 'SSX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 3+N lines, where N is its
!>          dimension. The first line contains the dimension N and the
!>          dimension M of an invariant subspace. The second line
!>          contains M integers, identifying the eigenvalues in the
!>          invariant subspace (by their position in a list of
!>          eigenvalues ordered by increasing real part). The next N
!>          lines contain the matrix. The last line contains the
!>          reciprocal condition number for the average of the selected
!>          eigenvalues, and the reciprocal condition number for the
!>          corresponding right invariant subspace. The end of data is
!>          indicated by a line containing N=0 and M=0. Even if no data
!>          is to be tested, there must be at least one line containing
!>          N=0 and M=0.
!>
!>-----------------------------------------------------------------------
!>
!> SGG input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, NS, MAXB, and
!>          NBCOL.
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
!>          The 3-character path name is 'SGG' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> SGS and SGV input files:
!>
!> line 1:  'SGS' or 'SGV' in columns 1 to 3.
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
!>          The 3-character path name is 'SGS' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> SXV input files:
!>
!> line 1:  'SXV' in columns 1 to 3.
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
!> remaining lines : Each example is stored on 3+2*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer). The next N lines contain the matrix A, one
!>          row per line. The next N lines contain the matrix B.  The
!>          next line contains the reciprocals of the eigenvalue
!>          condition numbers.  The last line contains the reciprocals of
!>          the eigenvector condition numbers.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> SGX input files:
!>
!> line 1:  'SGX' in columns 1 to 3.
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
!> remaining lines : Each example is stored on 3+2*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer).  The next line contains an integer k such
!>          that only the last k eigenvalues will be selected and appear
!>          in the leading diagonal blocks of $A$ and $B$. The next N
!>          lines contain the matrix A, one row per line.  The next N
!>          lines contain the matrix B.  The last line contains the
!>          reciprocal of the eigenvalue cluster condition number and the
!>          reciprocal of the deflating subspace (associated with the
!>          selected eigencluster) condition number.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> SSB input file:
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
!>          The 3-character path name is 'SSB'.
!>
!>-----------------------------------------------------------------------
!>
!> SBB input file:
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
!>          The 3-character path name is 'SBB'.
!>
!>-----------------------------------------------------------------------
!>
!> SEC input file:
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
!> SBL and SBK input files:
!>
!> line 1:  'SBL' in columns 1-3 to test SGEBAL, or 'SBK' in
!>          columns 1-3 to test SGEBAK.
!>
!> The remaining lines consist of specially constructed test cases.
!>
!>-----------------------------------------------------------------------
!>
!> SGL and SGK input files:
!>
!> line 1:  'SGL' in columns 1-3 to test SGGBAL, or 'SGK' in
!>          columns 1-3 to test SGGBAK.
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
!> precomputed examples, and LWORK = NMAX*(5*NMAX+5)+1 in the parameter
!> statements below.  For SVD, we assume NRHS may be as big as N.  The
!> parameter NEED is set to 14 to allow for 14 N-by-N matrices for SGG.
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
!> \author Olivier Thomine [F90 conversion, profiling & optimization]
!
!> \ingroup single_eig
!
!  =====================================================================
   PROGRAM SCHKEE
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
   PARAMETER          ( LWORK = NMAX*( 5*NMAX+5 )+1 )
   INTEGER            LIWORK
   PARAMETER          ( LIWORK = NMAX*( 5*NMAX+20 ) )
   INTEGER            MAXIN
   PARAMETER          ( MAXIN = 20 )
   INTEGER            MAXT
   PARAMETER          ( MAXT = 30 )
   INTEGER            NIN, NOUT
   PARAMETER          ( NIN = 5, NOUT = 6 )
!     ..
!     .. Local Scalars ..
   LOGICAL            CSD, FATAL, GLM, GQR, GSV, LSE, NEP, SBB, SBK, &
                      SBL, SEP, SES, SEV, SGG, SGK, SGL, SGS, SGV, &
                      SGX, SSB, SSX, SVD, SVX, SXV, TSTCHK, TSTDIF, &
                      TSTDRV, TSTERR
   CHARACTER          C1
   CHARACTER*3        C3, PATH
   CHARACTER*32       VNAME
   CHARACTER*10       INTSTR
   CHARACTER*80       LINE
   INTEGER            I, I1, IC, INFO, ITMP, K, LENP, MAXTYP, NEWSD, &
                      NK, NN, NPARMS, NRHS, NTYPES, &
                      VERS_MAJOR, VERS_MINOR, VERS_PATCH
   INTEGER*4          N_THREADS, ONE_THREAD
   REAL               EPS, THRESH, THRSHN
   INTEGER(8) nb_periods_sec, S1T_time, S2T_time
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
   REAL               D( NMAX, 12 ), RESULT( 500 ), TAUA( NMAX ), &
                      TAUB( NMAX ), X( 5*NMAX )
!     ..
!     .. Allocatable Arrays ..
   INTEGER AllocateStatus
   REAL, DIMENSION(:), ALLOCATABLE :: WORK
   REAL, DIMENSION(:,:), ALLOCATABLE :: A, B, C
!     ..
!     .. External Functions ..
   LOGICAL            LSAMEN
   REAL               SECOND, SLAMCH
   EXTERNAL           LSAMEN, SECOND, SLAMCH
!     ..
!     .. External Subroutines ..
   EXTERNAL           ALAREQ, SCHKBB, SCHKBD, SCHKBK, SCHKBL, SCHKEC, &
                      SCHKGG, SCHKGK, SCHKGL, SCHKHS, SCHKSB, SCHKST, &
                      SCKCSD, SCKGLM, SCKGQR, SCKGSV, SCKLSE, SDRGES, &
                      SDRGEV, SDRGSX, SDRGVX, SDRVBD, SDRVES, SDRVEV, &
                      SDRVSG, SDRVST, SDRVSX, SDRVVX, SERRBD, &
                      SERRED, SERRGG, SERRHS, SERRST, ILAVER, XLAENV, &
                      SDRGES3, SDRGEV3, &
                      SCHKST2STG, SDRVST2STG, SCHKSB2STG, SDRVSG2STG
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
   REAL               SELWI( 20 ), SELWR( 20 )
!     ..
!     .. Common blocks ..
   COMMON             / CENVIR / NPROC, NSHIFT, MAXB
   COMMON             / CLAENV / IPARMS
   COMMON             / INFOC / INFOT, NUNIT, OK, LERR
   COMMON             / SRNAMC / SRNAMT
   COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
!     ..
!     .. Data statements ..
   DATA               INTSTR / '0123456789' /
   DATA               IOLDSD / 0, 0, 0, 1 /
!     ..
!     .. Allocate memory dynamically ..
!
   ALLOCATE ( A(NMAX*NMAX,NEED), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( B(NMAX*NMAX,5), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( C(NCMAX*NCMAX,NCMAX*NCMAX), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   ALLOCATE ( WORK(LWORK), STAT = AllocateStatus )
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!     ..
!     .. Executable Statements ..
!
   A = 0.0
   B = 0.0
   C = 0.0
   D = 0.0
    call system_clock(count_rate=nb_periods_sec,count=S1T_time)
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
   PATH = LINE(1: 3)
   NEP = PATH == 'NEP' .OR. PATH == 'SHS'
   SEP = PATH == 'SEP' .OR. PATH == 'SST' .OR. &
         PATH == 'SSG' .OR. PATH == 'SE2'
   SVD = PATH == 'SVD' .OR. PATH == 'SBD'
   SEV = PATH == 'SEV'
   SES = PATH == 'SES'
   SVX = PATH == 'SVX'
   SSX = PATH == 'SSX'
   SGG = PATH == 'SGG'
   SGS = PATH == 'SGS'
   SGX = PATH == 'SGX'
   SGV = PATH == 'SGV'
   SXV = PATH == 'SXV'
   SSB = PATH == 'SSB'
   SBB = PATH == 'SBB'
   GLM = PATH == 'GLM'
   GQR = PATH == 'GQR' .OR. PATH == 'GRQ'
   GSV = PATH == 'GSV'
   CSD = PATH == 'CSD'
   LSE = PATH == 'LSE'
   SBL = PATH == 'SBL'
   SBK = PATH == 'SBK'
   SGL = PATH == 'SGL'
   SGK = PATH == 'SGK'
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
   ELSE IF( SEV ) THEN
      WRITE( NOUT, FMT = 9979 )
   ELSE IF( SES ) THEN
      WRITE( NOUT, FMT = 9978 )
   ELSE IF( SVX ) THEN
      WRITE( NOUT, FMT = 9977 )
   ELSE IF( SSX ) THEN
      WRITE( NOUT, FMT = 9976 )
   ELSE IF( SGG ) THEN
      WRITE( NOUT, FMT = 9975 )
   ELSE IF( SGS ) THEN
      WRITE( NOUT, FMT = 9964 )
   ELSE IF( SGX ) THEN
      WRITE( NOUT, FMT = 9965 )
   ELSE IF( SGV ) THEN
      WRITE( NOUT, FMT = 9963 )
   ELSE IF( SXV ) THEN
      WRITE( NOUT, FMT = 9962 )
   ELSE IF( SSB ) THEN
      WRITE( NOUT, FMT = 9974 )
   ELSE IF( SBB ) THEN
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
   ELSE IF( SBL ) THEN
!
!        SGEBAL:  Balancing
!
      CALL SCHKBL( NIN, NOUT )
      CYCLE
   ELSE IF( SBK ) THEN
!
!        SGEBAK:  Back transformation
!
      CALL SCHKBK( NIN, NOUT )
      CYCLE
   ELSE IF( SGL ) THEN
!
!        SGGBAL:  Balancing
!
      CALL SCHKGL( NIN, NOUT )
      CYCLE
   ELSE IF( SGK ) THEN
!
!        SGGBAK:  Back transformation
!
      CALL SCHKGK( NIN, NOUT )
      CYCLE
   ELSE IF( PATH == 'SEC' ) THEN
!
!        SEC:  Eigencondition estimation
!
      READ( NIN,*)THRESH
      CALL XLAENV( 1, 1 )
      CALL XLAENV( 12, 11 )
      CALL XLAENV( 13, 2 )
      CALL XLAENV( 14, 0 )
      CALL XLAENV( 15, 2 )
      CALL XLAENV( 16, 2 )
      TSTERR = .TRUE.
      CALL SCHKEC( THRESH, TSTERR, NIN, NOUT )
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
   IF( .NOT.( SGX .OR. SXV ) ) THEN
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
      WRITE( NOUT, FMT = 9983 )'P:    ', ( PVAL( I ), I = 1, NN )
   END IF
!
!     Read the values of N
!
   IF( SVD .OR. SBB .OR. GLM .OR. GQR .OR. GSV .OR. CSD .OR. &
       LSE ) THEN
      READ( NIN,*)( NVAL( I ), I = 1, NN )
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
      DO I = 1, NN
         NVAL( I ) = MVAL( I )
      ENDDO
   END IF
   IF( .NOT.( SGX .OR. SXV ) ) THEN
      WRITE( NOUT, FMT = 9983 )'N:    ', ( NVAL( I ), I = 1, NN )
   ELSE
      WRITE( NOUT, FMT = 9983 )'N:    ', NN
   END IF
!
!     Read the number of values of K, followed by the values of K
!
   IF( SSB .OR. SBB ) THEN
      READ( NIN,*)NK
      READ( NIN,*)( KVAL( I ), I = 1, NK )
      DO I = 1, NK
         IF( KVAL( I ) < 0 ) THEN
            WRITE( NOUT, FMT = 9989 )'    K ', KVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( KVAL( I ) > NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )'    K ', KVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
      ENDDO
      WRITE( NOUT, FMT = 9983 )'K:    ', ( KVAL( I ), I = 1, NK )
   END IF
!
   IF( SEV .OR. SES .OR. SVX .OR. SSX ) THEN
!
!        For the nonsymmetric QR driver routines, only one set of
!        parameters is allowed.
!
      READ( NIN,*)NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), &
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
   ELSE IF( SGS .OR. SGX .OR. SGV .OR. SXV ) THEN
!
!        For the nonsymmetric generalized driver routines, only one set
!        of parameters is allowed.
!
      READ( NIN,*)NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), NSVAL( 1 ), MXBVAL( 1 )
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
   ELSE IF( .NOT.SSB .AND. .NOT.GLM .AND. .NOT.GQR .AND. .NOT. &
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
      IF( .NOT.SBB ) THEN
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
      IF( NEP .OR. SEP .OR. SVD .OR. SGG ) THEN
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
!        Read the values of NSHIFT (if SGG) or NRHS (if SVD
!        or SBB).
!
      IF( SVD .OR. SBB .OR. SGG ) THEN
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
      IF( SGG ) THEN
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
      IF( NEP .OR. SGG ) THEN
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
      IF( SGG ) THEN
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
   EPS = SLAMCH( 'Underflow threshold' )
   WRITE( NOUT, FMT = 9981 )'underflow', EPS
   EPS = SLAMCH( 'Overflow threshold' )
   WRITE( NOUT, FMT = 9981 )'overflow ', EPS
   EPS = SLAMCH( 'Epsilon' )
   WRITE( NOUT, FMT = 9981 )'precision', EPS
!
!     Read the threshold value for the test ratios.
!
   READ(NIN,*) THRESH
   WRITE( NOUT, FMT = 9982 )THRESH
   IF( SEP .OR. SVD .OR. SGG ) THEN
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
   IF( .NOT.( SGX .OR. SXV ) ) THEN
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
      IF( .NOT.( SEV .OR. SES .OR. SVX .OR. SSX .OR. SGV .OR. &
          SGS ) .AND. NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
         GO TO 200
      END IF
!
   ELSE
      IF( SXV ) C3 = 'SXV'
      IF( SGX ) C3 = 'SGX'
   END IF
!
!     Reset the random number seed.
!
   IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
!
   IF( C3 == 'SHS' .OR. C3 == 'NEP' ) THEN
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
         CALL SERRHS( 'SHSEQR', NOUT )
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
         CALL SCHKHS( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), &
                      A( 1, 7 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), &
                      D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 8 ), &
                      A( 1, 9 ), A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), &
                      D( 1, 7 ), WORK, LWORK, IWORK, LOGWRK, RESULT, &
                      INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCHKHS', INFO
         ENDDO
!
   ELSE IF( C3 == 'SST' .OR. C3 == 'SEP' .OR. C3 == 'SE2' ) THEN
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
         CALL SERRST( 'SST', NOUT )
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
            CALL SCHKST2STG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                         NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), &
                         D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), &
                         D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), &
                         D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, &
                         A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), &
                         WORK, LWORK, IWORK, LIWORK, RESULT, INFO )
            ELSE
            CALL SCHKST( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                         NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), &
                         D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), &
                         D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), &
                         D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, &
                         A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), &
                         WORK, LWORK, IWORK, LIWORK, RESULT, INFO )
            ENDIF
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCHKST', INFO
         END IF
         IF( TSTDRV ) THEN
            IF( LSAMEN( 3, C3, 'SE2' ) ) THEN
            CALL SDRVST2STG( NN, NVAL, 18, DOTYPE, ISEED, THRESH, &
                         NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), &
                         D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), &
                         D( 1, 10 ), D( 1, 11), A( 1, 2 ), NMAX, &
                         A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, &
                         LWORK, IWORK, LIWORK, RESULT, INFO )
            ELSE
            CALL SDRVST( NN, NVAL, 18, DOTYPE, ISEED, THRESH, &
                         NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), &
                         D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), &
                         D( 1, 10 ), D( 1, 11), A( 1, 2 ), NMAX, &
                         A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, &
                         LWORK, IWORK, LIWORK, RESULT, INFO )
            ENDIF
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRVST', INFO
         END IF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'SSG' ) ) THEN
!
!        ----------------------------------------------
!        SSG:  Symmetric Generalized Eigenvalue Problem
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
!               CALL SDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
!     $                      D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
!     $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
!     $                      LWORK, IWORK, LIWORK, RESULT, INFO )
            CALL SDRVSG2STG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                             NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, &
                             D( 1, 3 ), D( 1, 3 ), A( 1, 3 ), NMAX, &
                             A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), &
                             A( 1, 7 ), WORK, LWORK, IWORK, LIWORK, &
                             RESULT, INFO )
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRVSG', INFO
         END IF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'SBD' ) .OR. LSAMEN( 3, C3, 'SVD' ) ) THEN
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
      CALL XLAENV( 1, 1 )
      CALL XLAENV( 9, 25 )
!
!        Test the error exits
!
      IF( TSTERR .AND. TSTCHK ) CALL SERRBD( 'SBD', NOUT )
      IF( TSTERR .AND. TSTDRV ) CALL SERRED( 'SBD', NOUT )
!
      DO I = 1, NPARMS
         NRHS = NSVAL( I )
         CALL XLAENV( 1, NBVAL( I ) )
         CALL XLAENV( 2, NBMIN( I ) )
         CALL XLAENV( 3, NXVAL( I ) )
         IF( NEWSD == 0 ) ISEED(1:4) = IOLDSD(1:4)
         WRITE( NOUT, FMT = 9995 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), NRHS
         IF( TSTCHK ) THEN
            CALL SCHKBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, NRHS, ISEED, &
                         THRESH, A( 1, 1 ), NMAX, D( 1, 1 ), &
                         D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 2 ), &
                         NMAX, A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, &
                         A( 1, 6 ), NMAX, A( 1, 7 ), A( 1, 8 ), WORK, &
                         LWORK, IWORK, NOUT, INFO )
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCHKBD', INFO
         END IF
         IF( TSTDRV ) THEN
            CALL SDRVBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, ISEED, &
                         THRESH, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, &
                         A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), &
                         A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), &
                         WORK, LWORK, IWORK, NOUT, INFO )
    ENDIF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'SEV' ) ) THEN
!
!        --------------------------------------------
!        SEV:  Nonsymmetric Eigenvalue Problem Driver
!              SGEEV (eigenvalues and eigenvectors)
!        --------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL SDRVEV( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), &
                      D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), &
                      NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, RESULT, &
                      WORK, LWORK, IWORK, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SGEEV', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'SES' ) ) THEN
!
!        --------------------------------------------
!        SES:  Nonsymmetric Eigenvalue Problem Driver
!              SGEES (Schur form)
!        --------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL SDRVES( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), &
                      A( 1, 4 ), NMAX, RESULT, WORK, LWORK, IWORK, &
                      LOGWRK, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SGEES', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( C3 == 'SVX' ) THEN
!
!        --------------------------------------------------------------
!        SVX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              SGEEVX (eigenvalues, eigenvectors and condition numbers)
!        --------------------------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL SDRVVX( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, &
                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), &
                      D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), &
                      NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, &
                      D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), &
                      D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), D( 1, 12 ), &
                      RESULT, WORK, LWORK, IWORK, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SGEEVX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'SSX' ) ) THEN
!
!        ---------------------------------------------------
!        SSX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              SGEESX (Schur form and condition numbers)
!        ---------------------------------------------------
!
      MAXTYP = 21
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRED( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL SDRVSX( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, &
                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), &
                      D( 1, 5 ), D( 1, 6 ), A( 1, 4 ), NMAX, &
                      A( 1, 5 ), RESULT, WORK, LWORK, IWORK, LOGWRK, &
                      INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SGEESX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'SGG' ) ) THEN
!
!        -------------------------------------------------
!        SGG:  Generalized Nonsymmetric Eigenvalue Problem
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
      IF( TSTCHK .AND. TSTERR ) CALL SERRGG( C3, NOUT )
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
         THRSHN = 10.
         IF( TSTCHK ) THEN
            CALL SCHKGG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, &
                         TSTDIF, THRSHN, NOUT, A( 1, 1 ), NMAX, &
                         A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                         A( 1, 6 ), A( 1, 7 ), A( 1, 8 ), A( 1, 9 ), &
                         NMAX, A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), &
                         D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), &
                         D( 1, 5 ), D( 1, 6 ), A( 1, 13 ), &
                         A( 1, 14 ), WORK, LWORK, LOGWRK, RESULT, &
                         INFO )
            IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCHKGG', INFO
         END IF
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'SGS' ) ) THEN
!
!        -------------------------------------------------
!        SGS:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGES (Schur form)
!        -------------------------------------------------
!
      MAXTYP = 26
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL SDRGES( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                      D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, &
                      RESULT, LOGWRK, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRGES', INFO
!
!     Blocked version
!
         CALL XLAENV(16,1)
         CALL SDRGES3( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                       A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                       A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                       D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, &
                       RESULT, LOGWRK, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRGES3', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( SGX ) THEN
!
!        -------------------------------------------------
!        SGX:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGESX (Schur form and condition numbers)
!        -------------------------------------------------
!
      MAXTYP = 5
      NTYPES = MAXTYP
      IF( NN < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL XLAENV( 5, 2 )
         CALL SDRGSX( NN, NCMAX, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), &
                      A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), &
                      C( 1, 1 ), NCMAX*NCMAX, A( 1, 12 ), WORK, &
                      LWORK, IWORK, LIWORK, LOGWRK, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRGSX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'SGV' ) ) THEN
!
!        -------------------------------------------------
!        SGV:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGEV (Eigenvalue/vector form)
!        -------------------------------------------------
!
      MAXTYP = 26
      NTYPES = MIN( MAXTYP, NTYPES )
      IF( NTYPES <= 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL SDRGEV( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                      A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                      A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                      A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), &
                      D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), &
                      WORK, LWORK, RESULT, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRGEV', INFO
!
! Blocked version
!
         CALL SDRGEV3( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, &
                       A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), &
                       A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), &
                       A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), &
                       D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), &
                       WORK, LWORK, RESULT, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRGEV3', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( SXV ) THEN
!
!        -------------------------------------------------
!        SXV:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGEVX (eigenvalue/vector with condition numbers)
!        -------------------------------------------------
!
      MAXTYP = 2
      NTYPES = MAXTYP
      IF( NN < 0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
      ELSE
         IF( TSTERR ) CALL SERRGG( C3, NOUT )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL SDRGVX( NN, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, &
                      A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), D( 1, 1 ), &
                      D( 1, 2 ), D( 1, 3 ), A( 1, 5 ), A( 1, 6 ), &
                      IWORK( 1 ), IWORK( 2 ), D( 1, 4 ), D( 1, 5 ), &
                      D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), &
                      WORK, LWORK, IWORK( 3 ), LIWORK-2, RESULT, &
                      LOGWRK, INFO )
!
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SDRGVX', INFO
      END IF
      WRITE( NOUT, FMT = 9973 )
      CYCLE
!
   ELSE IF( LSAMEN( 3, C3, 'SSB' ) ) THEN
!
!        ------------------------------
!        SSB:  Symmetric Band Reduction
!        ------------------------------
!
      MAXTYP = 15
      NTYPES = MIN( MAXTYP, NTYPES )
      CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
      IF( TSTERR ) CALL SERRST( 'SSB', NOUT )
!         CALL SCHKSB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ),
!     $                A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
      CALL SCHKSB2STG( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, &
                    THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), &
                    D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), &
                    A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCHKSB', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'SBB' ) ) THEN
!
!        ------------------------------
!        SBB:  General Band Reduction
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
         CALL SCHKBB( NN, MVAL, NVAL, NK, KVAL, MAXTYP, DOTYPE, NRHS, &
                      ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, &
                      A( 1, 2 ), 2*NMAX, D( 1, 1 ), D( 1, 2 ), &
                      A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, A( 1, 6 ), &
                      NMAX, A( 1, 7 ), WORK, LWORK, RESULT, INFO )
         IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCHKBB', INFO
         ENDDO
!
   ELSE IF( LSAMEN( 3, C3, 'GLM' ) ) THEN
!
!        -----------------------------------------
!        GLM:  Generalized Linear Regression Model
!        -----------------------------------------
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL SERRGG( 'GLM', NOUT )
      CALL SCKGLM( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, &
                   WORK, D( 1, 1 ), NIN, NOUT, INFO )
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCKGLM', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'GQR' ) ) THEN
!
!        ------------------------------------------
!        GQR:  Generalized QR and RQ factorizations
!        ------------------------------------------
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL SERRGG( 'GQR', NOUT )
      CALL SCKGQR( NN, MVAL, NN, PVAL, NN, NVAL, NTYPES, ISEED, &
                   THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), &
                   A( 1, 4 ), TAUA, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), &
                   B( 1, 4 ), B( 1, 5 ), TAUB, WORK, D( 1, 1 ), NIN, &
                   NOUT, INFO )
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCKGQR', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'GSV' ) ) THEN
!
!        ----------------------------------------------
!        GSV:  Generalized Singular Value Decomposition
!        ----------------------------------------------
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL SERRGG( 'GSV', NOUT )
      CALL SCKGSV( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), &
                   A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), TAUA, TAUB, &
                   B( 1, 4 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT, &
                   INFO )
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCKGSV', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'CSD' ) ) THEN
!
!        ----------------------------------------------
!        CSD:  CS Decomposition
!        ----------------------------------------------
!
      CALL XLAENV(1,1)
      IF( TSTERR ) CALL SERRGG( 'CSD', NOUT )
      CALL SCKCSD( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), &
                   A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), IWORK, WORK, &
                   D( 1, 1 ), NIN, NOUT, INFO )
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCKCSD', INFO
!
   ELSE IF( LSAMEN( 3, C3, 'LSE' ) ) THEN
!
!        --------------------------------------
!        LSE:  Constrained Linear Least Squares
!        --------------------------------------
!
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL SERRGG( 'LSE', NOUT )
      CALL SCKLSE( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, &
                   A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, &
                   WORK, D( 1, 1 ), NIN, NOUT, INFO )
      IF( INFO /= 0 ) WRITE( NOUT, FMT = 9980 )'SCKLSE', INFO
   ELSE
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9992 )C3
   END IF
   ENDDO
   IF( .NOT.( SGX .OR. SXV ) ) GO TO 190
  380 CONTINUE
   WRITE( NOUT, FMT = 9994 )
    call system_clock(count_rate=nb_periods_sec,count=S2T_time)
   WRITE( NOUT, FMT = 9993 ) real(S2T_time - S1T_time)/real(nb_periods_sec)
!
   DEALLOCATE (A, STAT = AllocateStatus)
   DEALLOCATE (B, STAT = AllocateStatus)
   DEALLOCATE (C, STAT = AllocateStatus)
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
 9986 FORMAT( ' Tests of the Symmetric Eigenvalue Problem routines' )
 9985 FORMAT( ' Tests of the Singular Value Decomposition routines' )
 9984 FORMAT( / ' The following parameter values will be used:' )
 9983 FORMAT( 4X, A, 10I6, / 10X, 10I6 )
 9982 FORMAT( / ' Routines pass computational tests if test ratio is ', &
         'less than', F8.2, / )
 9981 FORMAT( ' Relative machine ', A, ' is taken to be', E16.6 )
 9980 FORMAT( ' *** Error code from ', A, ' = ', I4 )
 9979 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', &
         / '    SGEEV (eigenvalues and eigevectors)' )
 9978 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', &
         / '    SGEES (Schur form)' )
 9977 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', &
         ' Driver', / '    SGEEVX (eigenvalues, eigenvectors and', &
         ' condition numbers)' )
 9976 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', &
         ' Driver', / '    SGEESX (Schur form and condition', &
         ' numbers)' )
 9975 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem routines' )
 9974 FORMAT( ' Tests of SSBTRD', / ' (reduction of a symmetric band ', &
         'matrix to tridiagonal form)' )
 9973 FORMAT( / 1X, 71( '-' ) )
 9972 FORMAT( / ' LAPACK VERSION ', I1, '.', I1, '.', I1 )
 9971 FORMAT( / ' Tests of the Generalized Linear Regression Model ', &
         'routines' )
 9970 FORMAT( / ' Tests of the Generalized QR and RQ routines' )
 9969 FORMAT( / ' Tests of the Generalized Singular Value', &
         ' Decomposition routines' )
 9968 FORMAT( / ' Tests of the Linear Least Squares routines' )
 9967 FORMAT( ' Tests of SGBBRD', / ' (reduction of a general band ', &
         'matrix to real bidiagonal form)' )
 9966 FORMAT( / / 1X, A3, ':  NRHS =', I4 )
 9965 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Expert Driver SGGESX' )
 9964 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Driver SGGES' )
 9963 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Driver SGGEV' )
 9962 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', &
         'Problem Expert Driver SGGEVX' )
 9961 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, &
         ', INMIN=', I4, &
         ', INWIN =', I4, ', INIBL =', I4, ', ISHFTS =', I4, &
         ', IACC22 =', I4)
 9960 FORMAT( / ' Tests of the CS Decomposition routines' )
!
!     End of SCHKEE
!
   END





