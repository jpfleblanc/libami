cmake_minimum_required(VERSION 3.1)
project(example CXX)

add_executable(${PROJECT_NAME} 
./src/examples_main.cpp
./src/ami_example.cpp 
)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11 -DUSE_MPI -lstdc++fs")# -msse4.2")# -ffast-math -funsafe-math-optimizations -msse4.2")# -ffast-math")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../libami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../libami/install/include)

target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/../../libami/install/lib/libami.so
 m )


