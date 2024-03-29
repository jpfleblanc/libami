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

find_package(Boost 1.71.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  IF(BOOST_MP)
    ADD_DEFINITIONS(-DBOOST_MP)
	message("Using boost::multiprecision.")
  else()
	message("Multiprecision is likely available but not enabled")
  ENDIF(BOOST_MP)
  
else()

message("Boost not found, so using standard precision. To enable multiprecision use -DBOOST_MP=ON. ")

endif()


# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
link_directories("/usr/lib")

# Changed to compile with c++14 compiler - Maybe can stay as c++11 since these are compiled separately
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -stdlib=libc++")
else()
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



IF(BOOST_MP)
target_link_libraries(${PROJECT_NAME} m quadmath  )
else()
target_link_libraries(${PROJECT_NAME} m)
ENDIF(BOOST_MP)

# if(MACOSX)
# target_link_libraries(${PROJECT_NAME}  m)
# else()
# target_link_libraries(${PROJECT_NAME}  m)
# endif()
   
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

#########################################################
# Pyami - requiring pybind11

# See if user has pybind11 installed
find_package(pybind11)
if(pybind11_FOUND)
  
  if(MAKE_PYAMI)
    project(pyami)

    set(CMAKE_CXX_RELEASE "-O3") # not sure if needed. It was in the example
    set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ version selection")

    include_directories("${CMAKE_SOURCE_DIR}/src")
    include_directories("${CMAKE_SOURCE_DIR}/pyami_src")

    file (GLOB SOURCE_FILES "src/*.cpp")
    file (GLOB HEADER_FILES "src/*.hpp")
    file (GLOB PYTHON_FILES "pyami_src/*.cpp" "pyami_src/*.hpp")

    pybind11_add_module(pyami
      ${SOURCE_FILES}
      ${HEADER_FILES}
      ${PYTHON_FILES}
    )

    target_link_libraries(pyami PUBLIC)


    install(TARGETS pyami
      COMPONENT python
      LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
    )

  else()
    message("PYAMI is available. But is not being compiled.  To Enable set MAKE_PYAMI=ON and set PYTHON_LIBRARY_DIR path for installation.")
  ENDIF(MAKE_PYAMI)

else()
  message("pybind11 was not found.  Python library will not be compiled. Continuing without.")
endif()




