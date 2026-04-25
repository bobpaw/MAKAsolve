cmake_minimum_required(VERSION 3.12)
project(MAKAsolve VERSION 0.0 LANGUAGES CXX)

include(CTest)
include(GNUInstallDirs)
include(CheckLanguage)

check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
	enable_language(CUDA)
	find_package(CUDAToolkit REQUIRED)
endif()
find_package(MPI REQUIRED)
find_package(SCOREC 4 REQUIRED)

# Link everything with MPI
link_libraries(MPI::MPI_CXX)

# Top level installed include directory
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
install(DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/include/MAKAsolve" TYPE INCLUDE
)

add_subdirectory(src)
if(BUILD_TESTING OR ${PROJECT_NAME}_BUILD_TESTING)
  add_subdirectory(test)
endif()
