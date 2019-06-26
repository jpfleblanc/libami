cmake_minimum_required(VERSION 2.8)

project (ami CXX)
#find_package(ALPSCore REQUIRED)
#For efficiency, compile with optimization by default
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

link_directories("/usr/lib")

#find_package(Boost 1.36.0)
#if(Boost_FOUND)
#  include_directories(${Boost_INCLUDE_DIRS})
#endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11 -lstdc++fs")# -ffast-math")

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

#INSTALL(FILES ${files} DESTINATION include)




