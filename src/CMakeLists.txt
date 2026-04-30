set(SOURCES
	Solver.cc
	Input.cc
  Timer.cc
)

add_library(MAKAsolve ${SOURCES})
target_link_libraries(MAKAsolve PUBLIC SCOREC::core)

if(USE_HYPRE)
  target_compile_definitions(MAKAsolve PUBLIC USE_HYPRE)
  target_sources(MAKAsolve PRIVATE HYPRESolver.cc)
  target_link_libraries(MAKAsolve PUBLIC HYPRE::HYPRE)
  if(CMAKE_CUDA_COMPILER)
  target_sources(MAKAsolve PRIVATE CUDALand.cu)
    target_link_libraries(MAKAsolve PUBLIC 
      -Wl,--start-group 
      HYPRE::HYPRE
      CUDA::cusparse
      -Wl,--end-group
    )
  else()
  endif()
endif()

install(TARGETS MAKAsolve)