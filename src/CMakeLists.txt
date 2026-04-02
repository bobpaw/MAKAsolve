set(SOURCES
	MAKAsolve.cc
	Solver.cu
)

add_library(MAKAsolve ${SOURCES})
set_target_properties(MAKAsolve PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(MAKAsolve
  CUDA::cusparse
  CUDA::cusolver
  MPI::MPI_CXX
)
install(TARGETS MAKAsolve)
