cmake_minimum_required(VERSION 3.27)
project(MAKAsolve VERSION 0.0 LANGUAGES CXX CUDA)

find_package(MPI REQUIRED)
find_package(CUDAToolkit REQUIRED)
find_package(SCOREC 4 REQUIRED)
link_libraries(MPI::MPI_CXX)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(src)
if(BUILD_TESTING OR ${PROJECT_NAME}_BUILD_TESTING)
add_subdirectory(test)
endif()
