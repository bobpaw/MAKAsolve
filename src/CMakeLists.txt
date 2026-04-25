set(SOURCES
	Solver.cc
	Input.cc
)

add_library(MAKAsolve ${SOURCES})
target_link_libraries(MAKAsolve PUBLIC SCOREC::core)
if(CMAKE_CUDA_COMPILER)
  target_compile_definitions(MAKAsolve PUBLIC USE_CUDA)
  target_sources(MAKAsolve PRIVATE SparseSolver.cu)
  target_link_libraries(MAKAsolve PUBLIC CUDA::cusparse CUDA::cusolver)
endif()
install(TARGETS MAKAsolve)
