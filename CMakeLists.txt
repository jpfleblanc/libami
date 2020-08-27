cmake_minimum_required(VERSION 3.11)

project (ami CXX)
#find_package(ALPSCore REQUIRED)
#For efficiency, compile with optimization by default
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()



list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
link_directories("/usr/lib")

#find_package(Boost 1.36.0)
#if(Boost_FOUND)
#  include_directories(${Boost_INCLUDE_DIRS})
#endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11 -lstdc++fs")# -ffast-math -funsafe-math-optimizations -msse4.2 ")#-ffast-math")# -ffast-math")

#include_directories(/home/jpfleblanc/software/boost_1_55_0)

#include_directories(/home/jpfleblanc/software/boost_1_55_0)

set(LIB_FILES 
  ./src/ami_construct.cpp
  ./src/ami_evaluate.cpp
  ./src/ami_helper.cpp
  ./src/ami_example.cpp
  ./src/ami_load.cpp
  )

ADD_LIBRARY(ami SHARED ${LIB_FILES})


#target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES} m)
target_link_libraries(${PROJECT_NAME}  m stdc++fs)


INSTALL(TARGETS ami 
   ARCHIVE DESTINATION "lib/"
   LIBRARY DESTINATION "lib/"
   COMPONENT library)

#FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/." # source directory
        DESTINATION include # target directory
        FILES_MATCHING PATTERN "*.hpp" # select header files
)

add_subdirectory("doc")



#testing setup
#option(Testing "Enable testing" ON)

#if (Testing)
#  add_library(gtest ./test/gtest-all.cc ./test/gtest_main.cc)
#  list(APPEND  libami)
#  enable_testing(test)
#  add_subdirectory(test)
#endif (Testing)

###########################################################
# New testing code


add_subdirectory(test)
enable_testing()
add_test(initial test/initial)
add_test(pole_equiv_fun test/pole_equiv)



############################################################
#option(Documentation "Build documentation" ON)
#if (DocumentationOnly)
#  set(Documentation ON)
#endif (DocumentationOnly)

#set(_am_building_documentation ${Documentation})
#if (Documentation)
#   message(STATUS "Enable build of documentation")
#    set(Documentation OFF) # not to spawn make doc jobs for each module
#endif(Documentation)


#Doxygen building is a function to prevent namespace damage
#function(build_documentation_)
#    set(DOXYFILE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
#	set(DOXYFILE_EXTRA_SOURCES "${CMAKE_SOURCE_DIR}/doc")
#    message(STATUS "Building documentation for sources in ${DOXYFILE_EXTRA_SOURCES}")  
#    set(DOXYFILE_IN "${CMAKE_SOURCE_DIR}/doc/Doxyfile.in") 
#    set(CMAKE_PROJECT_NAME "LIBAMI") # to provide name for the documentation 
#    include(UseDoxygen)
#    unset(CMAKE_PROJECT_NAME)
#endfunction()
#if(_am_building_documentation)
#  build_documentation_()
#endif(_am_building_documentation)


#########################################################
#New documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  # exclude sqlite code
  set(DOXYGEN_EXCLUDE_PATTERNS
        */sqlite3/*
  )
  # doxygen settings can be set here, prefixed with "DOXYGEN_"
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/LibAMI-docs")
  set(DOXYGEN_GENERATE_LATEX YES)
  set(DOXYGEN_GENERATE_XML YES)
  # this target will only be built if specifically asked to.
  # run "make LibAMI-docs" to create the doxygen documentation
  doxygen_add_docs(
    LibAMI-docs
    ${PROJECT_SOURCE_DIR}/src
    COMMENT "Generate API-documents for LibAMI."
  )
endif(DOXYGEN_FOUND)






