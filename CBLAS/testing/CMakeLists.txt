#######################################################################
#  This CMakeLists.txt creates the test programs for the CBLAS routines.
#
#######################################################################

macro(add_cblas_test output input target)
  if(NOT "${input}" STREQUAL "")
    set(TEST_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/${input}")
  endif()
  set(TEST_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${output}")
  set(testName "${target}")

  if(DEFINED TEST_INPUT AND EXISTS "${TEST_INPUT}")
    add_test(NAME CBLAS-${testName} COMMAND "${CMAKE_COMMAND}"
      -DTEST=$<TARGET_FILE:${target}>
      -DINPUT=${TEST_INPUT}
      -DOUTPUT=${TEST_OUTPUT}
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
  else()
    add_test(NAME CBLAS-${testName} COMMAND "${CMAKE_COMMAND}"
      -DTEST=$<TARGET_FILE:${target}>
      -DOUTPUT=${TEST_OUTPUT}
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
  endif()
endmacro()


# Object files for single precision real
set(STESTL1O c_sblas1.c)
set(STESTL2O c_sblas2.c c_s2chke.c auxiliary.c c_xerbla.c)
set(STESTL3O c_sblas3.c c_s3chke.c auxiliary.c c_xerbla.c)

# Object files for double precision real
set(DTESTL1O c_dblas1.c)
set(DTESTL2O c_dblas2.c c_d2chke.c auxiliary.c c_xerbla.c)
set(DTESTL3O c_dblas3.c c_d3chke.c auxiliary.c c_xerbla.c)

# Object files for single precision complex
set(CTESTL1O c_cblat1.f c_cblas1.c)
set(CTESTL2O c_cblas2.c c_c2chke.c auxiliary.c c_xerbla.c)
set(CTESTL3O c_cblas3.c c_c3chke.c auxiliary.c c_xerbla.c)

# Object files for double precision complex
set(ZTESTL1O c_zblas1.c)
set(ZTESTL2O c_zblas2.c c_z2chke.c auxiliary.c c_xerbla.c)
set(ZTESTL3O c_zblas3.c c_z3chke.c auxiliary.c c_xerbla.c)



if(BUILD_SINGLE)
  add_executable(xscblat1 c_sblat1.f ${STESTL1O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xscblat2 c_sblat2.f ${STESTL2O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xscblat3 c_sblat3.f ${STESTL3O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)

  if(HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xscblat1 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xscblat2 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xscblat3 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
  endif()

  target_link_libraries(xscblat1 ${CBLASLIB})
  target_link_libraries(xscblat2 ${CBLASLIB})
  target_link_libraries(xscblat3 ${CBLASLIB})

  add_cblas_test(stest1.out ""   xscblat1)
  add_cblas_test(stest2.out sin2 xscblat2)
  add_cblas_test(stest3.out sin3 xscblat3)
endif()

if(BUILD_DOUBLE)
  add_executable(xdcblat1 c_dblat1.f ${DTESTL1O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xdcblat2 c_dblat2.f ${DTESTL2O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xdcblat3 c_dblat3.f ${DTESTL3O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)

  if(HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xdcblat1 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xdcblat2 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xdcblat3 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
  endif()

  target_link_libraries(xdcblat1 ${CBLASLIB})
  target_link_libraries(xdcblat2 ${CBLASLIB})
  target_link_libraries(xdcblat3 ${CBLASLIB})

  add_cblas_test(dtest1.out ""   xdcblat1)
  add_cblas_test(dtest2.out din2 xdcblat2)
  add_cblas_test(dtest3.out din3 xdcblat3)
endif()

if(BUILD_COMPLEX)
  add_executable(xccblat1 c_cblat1.f ${CTESTL1O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xccblat2 c_cblat2.f ${CTESTL2O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xccblat3 c_cblat3.f ${CTESTL3O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)

  if(HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xccblat1 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xccblat2 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xccblat3 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
  endif()

  target_link_libraries(xccblat1 ${CBLASLIB} ${BLAS_LIBRARIES})
  target_link_libraries(xccblat2 ${CBLASLIB})
  target_link_libraries(xccblat3 ${CBLASLIB})

  add_cblas_test(ctest1.out ""   xccblat1)
  add_cblas_test(ctest2.out cin2 xccblat2)
  add_cblas_test(ctest3.out cin3 xccblat3)
endif()

if(BUILD_COMPLEX16)
  add_executable(xzcblat1 c_zblat1.f ${ZTESTL1O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xzcblat2 c_zblat2.f ${ZTESTL2O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)
  add_executable(xzcblat3 c_zblat3.f ${ZTESTL3O} ${LAPACK_BINARY_DIR}/test_include/cblas_test.h)

  if(HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xzcblat1 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xzcblat2 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(xzcblat3 PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
  endif()

  target_link_libraries(xzcblat1 ${CBLASLIB})
  target_link_libraries(xzcblat2 ${CBLASLIB})
  target_link_libraries(xzcblat3 ${CBLASLIB})

  add_cblas_test(ztest1.out ""   xzcblat1)
  add_cblas_test(ztest2.out zin2 xzcblat2)
  add_cblas_test(ztest3.out zin3 xzcblat3)
endif()
