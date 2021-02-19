cmake_minimum_required(VERSION 3.18.0)


project(googletest-git NONE)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

################################
# Tests
################################
# Add test cpp file
#add_executable(blank blank.cpp)
add_executable(fb fermi_bose.cpp)
add_executable(getE get_energy.cpp)
add_executable(pole_equiv pole_equiv.cpp)


target_link_libraries(fb amibase gtest gtest_main)
target_link_libraries(getE amibase gtest gtest_main)
target_link_libraries(pole_equiv amibase gtest gtest_main)


