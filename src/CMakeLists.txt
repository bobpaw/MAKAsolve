set(SOURCES
	MAKAsolve.cc
	Input.cc
)

add_library(MAKAsolve ${SOURCES})
target_link_libraries(MAKAsolve PUBLIC SCOREC::core)
install(TARGETS MAKAsolve)
