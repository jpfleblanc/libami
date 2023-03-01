cmake_minimum_required(VERSION 3.1)
project(example CXX)

add_executable(${PROJECT_NAME} 
./src/examples_main.cpp
./src/ami_example.cpp 
)


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -stdlib=libc++")
else()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -lstdc++fs ") 
endif()


#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -DUSE_MPI -lstdc++fs")# -msse4.2")# -ffast-math -funsafe-math-optimizations -msse4.2")# -ffast-math")

if(LIBAMI_ROOT)
link_directories(${LIBAMI_ROOT}/lib)
include_directories(${LIBAMI_ROOT}/include)
target_link_libraries(${PROJECT_NAME}  ${LIBAMI_ROOT}/lib/libami.so
 m )
else()
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../libami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../libami/install/include)
target_link_libraries(${PROJECT_NAME}  ${CMAKE_CURRENT_SOURCE_DIR}/../../libami/install/lib/libami.so
 m )
ENDIF(LIBAMI_ROOT)




