cmake_minimum_required(VERSION 3.11)

project (amibase CXX)
#find_package(ALPSCore REQUIRED)
#For efficiency, compile with optimization by default
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

if (NOT BUILD_DOC)
	set(BUILD_DOC OFF)
endif()



list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
link_directories("/usr/lib")


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11 -lstdc++fs ") 


set(LIB_FILES 
   ./src/ami_base.cpp
  ./src/ami_base_evaluate.cpp
  ./src/ami_base_construct.cpp 
	./src/ami_base_optimize.cpp 
	./src/ami_base_terms.cpp
	./src/ami_base_terms_evaluate.cpp
		./src/ami_base_terms_optimize.cpp
  ./src/ami_calc_io.cpp
  ./src/ami_calc_evaluate.cpp
	./src/ami_calc_spectral.cpp
	./src/ami_calc_terms_spectral.cpp
	./src/ami_spec.cpp
	./src/ami_spec_io.cpp
  )

ADD_LIBRARY(amibase SHARED ${LIB_FILES})


#target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES} m)
target_link_libraries(${PROJECT_NAME}  m stdc++fs)

   
   INSTALL(TARGETS amibase 
   ARCHIVE DESTINATION "lib/"
   LIBRARY DESTINATION "lib/"
   COMPONENT library)

INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/." # source directory
        DESTINATION include # target directory
        FILES_MATCHING PATTERN "*.hpp" # select header files
)




add_subdirectory(test)
enable_testing()
add_test(fb test/fb)
add_test(getE test/getE)
add_test(pole_equiv_fun test/pole_equiv)


if(BUILD_DOC)
message("Adding doc subdirectory")
add_subdirectory("doc")

#########################################################
#Documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  # exclude sqlite code
  set(DOXYGEN_EXCLUDE_PATTERNS
        */sqlite3/*
  )
  # doxygen settings can be set here, prefixed with "DOXYGEN_"
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/DOCS")
  set(DOXYGEN_GENERATE_LATEX YES)
  set(DOXYGEN_GENERATE_XML YES)
  # this target will only be built if specifically asked to.
  # run "make LibAMI-docs" to create the doxygen documentation
  doxygen_add_docs(
    DOCS
    ${PROJECT_SOURCE_DIR}/src
    COMMENT "Generate API-documents for LibAMI."
  )
endif(DOXYGEN_FOUND)


endif(BUILD_DOC)





