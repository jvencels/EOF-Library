INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

ADD_ELMERTEST_MODULE(interpolationTest OpenFOAM2Elmer OpenFOAM2Elmer.F90)
ADD_ELMERTEST_MODULE(interpolationTest Elmer2OpenFOAM Elmer2OpenFOAM.F90)
ADD_ELMERTEST_MODULE(interpolationTest InterpolationTester InterpolationTester.F90)

file(COPY 0 constant system travis_run.sh meshElmer elem.sif dg.sif ip.sif OpenFOAM2Elmer.F90 Elmer2OpenFOAM.F90 InterpolationTester.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(interpolationTest NPROCS 2 LABELS EOF)
