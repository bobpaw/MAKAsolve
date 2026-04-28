set(SOURCES
	Solver.cc
	Input.cc
	Timer.cc
)

add_library(MAKAsolve ${SOURCES})
target_link_libraries(MAKAsolve PUBLIC SCOREC::core)
if(CMAKE_CUDA_COMPILER)
	target_link_libraries(MAKAsolve PUBLIC 
		-Wl,--start-group 
		HYPRE::HYPRE
		CUDA::cusparse
		-Wl,--end-group
	)
else()
	target_link_libraries(MAKAsolve PUBLIC HYPRE::HYPRE)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64le")
    target_compile_definitions(MAKAsolve PRIVATE POWER9)
endif()

install(TARGETS MAKAsolve)

