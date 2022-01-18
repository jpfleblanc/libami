cmake_minimum_required(VERSION 3.18...3.22.1)

project (ami CXX)
#find_package(ALPSCore REQUIRED)
#For efficiency, compile with optimization by default
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

if (NOT BUILD_DOC)
	set(BUILD_DOC OFF)
endif()


# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
link_directories("/usr/lib")


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -stdlib=libc++ -lc++experimental")
else
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -lstdc++fs ") 
endif()



set(LIB_FILES 
   ./src/ami_base.cpp
  ./src/ami_base_evaluate.cpp
  ./src/ami_base_construct.cpp 
	./src/ami_base_optimize.cpp 
	./src/ami_base_terms.cpp
	./src/ami_base_terms_evaluate.cpp
		./src/ami_base_terms_optimize.cpp
  )

ADD_LIBRARY(ami SHARED ${LIB_FILES})


#target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES} m)
target_link_libraries(${PROJECT_NAME}  m stdc++fs)

   
   INSTALL(TARGETS ami 
   ARCHIVE DESTINATION "lib/"
   LIBRARY DESTINATION "lib/"
   COMPONENT library)

INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/." # source directory
        DESTINATION include # target directory
        FILES_MATCHING PATTERN "*.hpp" # select header files
)





if(CMAKE_BUILD_TYPE MATCHES Release OR TEST)

add_subdirectory(test)
enable_testing()
add_test(fb test/fb)
add_test(helper_functions test/helper_functions)
add_test(num_test test/num_test)
add_test(construct_tests test/construct_tests)

endif(CMAKE_BUILD_TYPE MATCHES Release OR TEST)

if(BUILD_DOC)


message("Adding doc subdirectory")
add_subdirectory("doc")

#########################################################
#Documentation
# find_package(Doxygen)
# if(DOXYGEN_FOUND)
  # # exclude sqlite code
  # set(DOXYGEN_EXCLUDE_PATTERNS
        # */sqlite3/*
  # )
  # message("run 'make DOCS' to generate documentation")
  # # doxygen settings can be set here, prefixed with "DOXYGEN_"
  # set(DOXYGEN_EXTRACT_PRIVATE YES)
  # set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/DOCS")
  # set(DOXYGEN_GENERATE_LATEX YES)
  # set(DOXYGEN_GENERATE_XML YES)
  # set(DOXYGEN_SORT_MEMBER_DOCS YES)
  # set(DOXYGEN_USE_MATHJAX YES)
  # #set(DOXYGEN_USE_MDFILE_AS_MAINPAGE main_page.dox)
  # # this target will only be built if specifically asked to.
  # # run "make LibAMI-docs" to create the doxygen documentation
  # # run "make DOCS" to create the doxygen documentation
  # doxygen_add_docs(
    # DOCS
    # ${PROJECT_SOURCE_DIR}/src
    # COMMENT "Generate API-documents for LibAMI."
  # )
# endif(DOXYGEN_FOUND)


endif(BUILD_DOC)





