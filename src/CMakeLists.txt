set(SOURCES
	MAKAsolve.cc
	SparseSolver.cu
)

add_library(MAKAsolve ${SOURCES})
set_target_properties(MAKAsolve PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(MAKAsolve
  CUDA::cusparse
  CUDA::cusolver
)
install(TARGETS MAKAsolve)
