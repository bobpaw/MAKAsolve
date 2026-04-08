set(${PROJECT_NAME}_TEST_LAUNCHER
	"" CACHE STRING
	"Test launcher, e.g. valgrind"
)

# Add a MAKAsolve test target executable.
macro(add_maka_test_exe EXE)
	add_executable(${EXE} ${ARGN})
	target_link_libraries(${EXE} MAKAsolve)
endmacro(add_maka_test_exe)

# Add a MAKAsolve test case.
macro(maka_test TESTNAME EXE NPROCS)
	if(${PROJECT_NAME}_TEST_LAUNCHER)
		add_test(NAME ${TESTNAME}
			COMMAND "${${PROJECT_NAME}_TEST_LAUNCHER}"
			"${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${NPROCS}
			"${EXE}" ${ARGN}
		)
	else()
		add_test(NAME ${TESTNAME}
			COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${NPROCS}
			"${EXE}" ${ARGN}
		)
	endif()
	set_tests_properties(${TESTNAME} PROPERTIES PROCESSORS ${NPROCS})
endmacro(maka_test)

# Add a simple MAKAsolve test case build from a single source file.
macro(maka_test_single TESTTARGET TESTSOURCE NPROCS)
	add_maka_test_exe(${TESTTARGET} ${TESTSOURCE})
	maka_test(${TESTTARGET} "$<TARGET_FILE:${TESTTARGET}>" ${NPROCS} ${ARGN})
endmacro(maka_test_single)

maka_test_single(SparseSolverTest SparseSolverTest.cc 1)
