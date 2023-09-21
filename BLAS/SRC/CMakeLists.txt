#######################################################################
#  This is the makefile to create a library for the BLAS.
#  The files are grouped as follows:
#
#       SBLAS1 -- Single precision real BLAS routines
#       CBLAS1 -- Single precision complex BLAS routines
#       DBLAS1 -- Double precision real BLAS routines
#       ZBLAS1 -- Double precision complex BLAS routines
#
#       CB1AUX -- Real BLAS routines called by complex routines
#       ZB1AUX -- D.P. real BLAS routines called by d.p. complex
#                 routines
#
#      ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
#
#       SBLAS2 -- Single precision real BLAS2 routines
#       CBLAS2 -- Single precision complex BLAS2 routines
#       DBLAS2 -- Double precision real BLAS2 routines
#       ZBLAS2 -- Double precision complex BLAS2 routines
#
#       SBLAS3 -- Single precision real BLAS3 routines
#       CBLAS3 -- Single precision complex BLAS3 routines
#       DBLAS3 -- Double precision real BLAS3 routines
#       ZBLAS3 -- Double precision complex BLAS3 routines
#
#######################################################################

#---------------------------------------------------------
#  Level 1 BLAS
#---------------------------------------------------------

set(SBLAS1 isamax.f90 sasum.f90 saxpy.f90 scopy.f90 sdot.f90 snrm2.f90
	srot.f90 srotg.f90 sscal.f90 sswap.f90 sdsdot.f90 srotmg.f90 srotm.f90)

set(CBLAS1 scabs1.f90 scasum.f90 scnrm2.f90 icamax.f90 caxpy.f90 ccopy.f90
	cdotc.f90 cdotu.f90 csscal.f90 crotg.f90 cscal.f90 cswap.f90 csrot.f90)

set(DBLAS1 idamax.f90 dasum.f90 daxpy.f90 dcopy.f90 ddot.f90 dnrm2.f90
	drot.f90 drotg.f90 dscal.f90 dsdot.f90 dswap.f90 drotmg.f90 drotm.f90)

set(DB1AUX sscal.f90 isamax.f90)

set(ZBLAS1 dcabs1.f90 dzasum.f90 dznrm2.f90 izamax.f90 zaxpy.f90 zcopy.f90
	zdotc.f90 zdotu.f90 zdscal.f90 zrotg.f90 zscal.f90 zswap.f90 zdrot.f90)

set(CB1AUX
    isamax.f90 idamax.f90
    sasum.f90 saxpy.f90 scopy.f90 sdot.f90 sgemm.f90 sgemv.f90 snrm2.f90 srot.f90 sscal.f90
    sswap.f90)

set(ZB1AUX
    icamax.f90 idamax.f90
    cgemm.f90 cherk.f90 cscal.f90 ctrsm.f90
    dasum.f90 daxpy.f90 dcopy.f90 ddot.f90 dgemm.f90 dgemv.f90 dnrm2.f90 drot.f90 dscal.f90
    dswap.f90
    scabs1.f90)

#---------------------------------------------------------------------
#  Auxiliary routines needed by both the Level 2 and Level 3 BLAS
#---------------------------------------------------------------------
set(ALLBLAS lsame.f90 xerbla.f90 xerbla_array.f90)

#---------------------------------------------------------
#  Level 2 BLAS
#---------------------------------------------------------
set(SBLAS2 sgemv.f90 sgbmv.f90 ssymv.f90 ssbmv.f90 sspmv.f90
	strmv.f90 stbmv.f90 stpmv.f90 strsv.f90 stbsv.f90 stpsv.f90
	sger.f90 ssyr.f90 sspr.f90 ssyr2.f90 sspr2.f90)

set(CBLAS2 cgemv.f90 cgbmv.f90 chemv.f90 chbmv.f90 chpmv.f90
	ctrmv.f90 ctbmv.f90 ctpmv.f90 ctrsv.f90 ctbsv.f90 ctpsv.f90
	cgerc.f90 cgeru.f90 cher.f90 chpr.f90 cher2.f90 chpr2.f90)

set(DBLAS2 dgemv.f90 dgbmv.f90 dsymv.f90 dsbmv.f90 dspmv.f90
	dtrmv.f90 dtbmv.f90 dtpmv.f90 dtrsv.f90 dtbsv.f90 dtpsv.f90
	dger.f90 dsyr.f90 dspr.f90 dsyr2.f90 dspr2.f90)

set(ZBLAS2 zgemv.f90 zgbmv.f90 zhemv.f90 zhbmv.f90 zhpmv.f90
	ztrmv.f90 ztbmv.f90 ztpmv.f90 ztrsv.f90 ztbsv.f90 ztpsv.f90
	zgerc.f90 zgeru.f90 zher.f90 zhpr.f90 zher2.f90 zhpr2.f90)

#---------------------------------------------------------
#  Level 3 BLAS
#---------------------------------------------------------
set(SBLAS3 sgemm.f90 ssymm.f90 ssyrk.f90 ssyr2k.f90 strmm.f90 strsm.f90)

set(CBLAS3 cgemm.f90 csymm.f90 csyrk.f90 csyr2k.f90 ctrmm.f90 ctrsm.f90
	chemm.f90 cherk.f90 cher2k.f90)

set(DBLAS3 dgemm.f90 dsymm.f90 dsyrk.f90 dsyr2k.f90 dtrmm.f90 dtrsm.f90)

set(ZBLAS3 zgemm.f90 zsymm.f90 zsyrk.f90 zsyr2k.f90 ztrmm.f90 ztrsm.f90
	zhemm.f90 zherk.f90 zher2k.f90)


set(SOURCES)
if(BUILD_SINGLE)
  list(APPEND SOURCES ${SBLAS1} ${ALLBLAS} ${SBLAS2} ${SBLAS3})
endif()
if(BUILD_DOUBLE)
  list(APPEND SOURCES
       ${DBLAS1} ${DB1AUX} ${ALLBLAS} ${DBLAS2} ${DBLAS3} ${SBLAS3})
endif()
if(BUILD_COMPLEX)
  list(APPEND SOURCES ${CBLAS1} ${CB1AUX} ${ALLBLAS} ${CBLAS2} ${CBLAS3})
endif()
if(BUILD_COMPLEX16)
  list(APPEND SOURCES ${ZBLAS1} ${ZB1AUX} ${ALLBLAS} ${ZBLAS2} ${ZBLAS3})
endif()
list(REMOVE_DUPLICATES SOURCES)

add_library(${BLASLIB}_obj OBJECT ${SOURCES})
set_target_properties(${BLASLIB}_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(BUILD_INDEX64_EXT_API)
  set(SOURCES_64_F)
  # Copy files so we can set source property specific to /${BLASLIB}_64_obj target
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BLASLIB}_64_obj)
  file(COPY ${SOURCES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${BLASLIB}_64_obj)
  file(GLOB SOURCES_64_F ${CMAKE_CURRENT_BINARY_DIR}/${BLASLIB}_64_obj/*.f*)
  add_library(${BLASLIB}_64_obj OBJECT ${SOURCES_64_F})
  target_compile_options(${BLASLIB}_64_obj PRIVATE ${FOPT_ILP64})
  set_target_properties(${BLASLIB}_64_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
  #Add _64 suffix to all Fortran functions via macros
  foreach(F IN LISTS SOURCES_64_F)
      set(COPT_64_F -cpp)
      file(STRINGS ${F} ${F}.lst)
      list(FILTER ${F}.lst INCLUDE REGEX "subroutine|SUBROUTINE|external|EXTERNAL|function|FUNCTION")
      list(FILTER ${F}.lst EXCLUDE REGEX "^!.*")
      list(FILTER ${F}.lst EXCLUDE REGEX "^[*].*")
      list(FILTER ${F}.lst EXCLUDE REGEX "end|END")
      foreach(FUNC IN LISTS ${F}.lst)
        string(REGEX REPLACE "^[a-zA-Z0-9_ *]*(subroutine|SUBROUTINE|external|EXTERNAL|function|FUNCTION)[ ]*[*]?" "" FUNC ${FUNC})
        string(REGEX REPLACE "[(][a-zA-Z0-9_, )]*$" "" FUNC ${FUNC})
        string(STRIP ${FUNC} FUNC)
        list(APPEND COPT_64_F "-D${FUNC}=${FUNC}_64")
      endforeach()
      list(REMOVE_DUPLICATES COPT_64_F)
      set_source_files_properties(${F} PROPERTIES COMPILE_OPTIONS "${COPT_64_F}")
  endforeach()
endif()

add_library(${BLASLIB}
        $<TARGET_OBJECTS:${BLASLIB}_obj>
        $<$<BOOL:${BUILD_INDEX64_EXT_API}>: $<TARGET_OBJECTS:${BLASLIB}_64_obj>>)

set_target_properties(
  ${BLASLIB} PROPERTIES
  VERSION ${LAPACK_VERSION}
  SOVERSION ${LAPACK_MAJOR_VERSION}
  POSITION_INDEPENDENT_CODE ON
  )
lapack_install_library(${BLASLIB})

if( TEST_FORTRAN_COMPILER )
  add_dependencies( ${BLASLIB} run_test_zcomplexabs run_test_zcomplexdiv run_test_zcomplexmult run_test_zminMax )
endif()
