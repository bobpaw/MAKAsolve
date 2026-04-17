set(SOURCES
	MAKAsolve.cc
	Input.cc
	LinearSolver.cc
	CPULinearSolver.cc
)

add_library(MAKAsolve ${SOURCES})
target_link_libraries(MAKAsolve PUBLIC SCOREC::core)
install(TARGETS MAKAsolve)
