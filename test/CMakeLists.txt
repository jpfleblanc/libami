cmake_minimum_required(VERSION 3.18.0)
#include(GoogleTest)
#add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
#                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
#                 EXCLUDE_FROM_ALL)

## The gtest/gtest_main targets carry header search path
## dependencies automatically when using CMake 2.8.11 or
## later. Otherwise we have to add them here ourselves.
#if (CMAKE_VERSION VERSION_LESS 2.8.11)
#  include_directories("${gtest_SOURCE_DIR}/include")
#endif()

## Now simply link against gtest or gtest_main as needed. Eg
#add_executable(pole_equiv pole_equiv.cpp)
#target_link_libraries(pole_equiv gtest_main)
#add_test(NAME example_test COMMAND pole_equiv)

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
add_executable(${TESTNAME} blank.cpp)
add_executable(${TESTNAME} secondtest.cpp)
add_executable(${TESTNAME} pole_equiv.cpp)
# Link test executable against gtest & gtest_main
target_link_libraries(${TESTNAME} gtest gtest_main)







#include(GoogleTest) #defined in ./cmake
#set(test_src
#  blank
#  secondtest
#  pole_equiv
#)
#foreach(test ${test_src})
#  add_gtest(${test})
#endforeach(test)

