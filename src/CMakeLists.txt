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
if(USE_AMGX)
  target_compile_definitions(MAKAsolve PUBLIC USE_AMGX)
  target_sources(MAKAsolve PRIVATE ParallelSolver.cu)
  target_link_libraries(MAKAsolve PUBLIC amgx CUDA::cublas)
endif()
if(CMAKE_CUDA_COMPILER OR USE_AMGX)
  set_target_properties(MAKAsolve PROPERTIES CUDA_ARCHITECTURES OFF)
endif()
install(TARGETS MAKAsolve)
