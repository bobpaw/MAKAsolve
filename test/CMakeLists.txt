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
macro(maka_test TESTNAME NPROCS EXE)
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
	maka_test(${TESTTARGET} ${NPROCS} "$<TARGET_FILE:${TESTTARGET}>" ${ARGN})
endmacro(maka_test_single)

set(TESTDATA "${CMAKE_CURRENT_SOURCE_DIR}/data")

# maka_test_single(empty empty.cc 1)
maka_test_single(load-mesh load-mesh.cc 1
  "${TESTDATA}/square.dmg" "${TESTDATA}/square.smb"
)
maka_test_single(read-input read-input.cc 1
  "${TESTDATA}/square-stable.rc"
)
add_maka_test_exe(stiffness stiffness.cc)
maka_test(stiffness0 1
	./stiffness "${TESTDATA}/square.dmg" "${TESTDATA}/square.smb" 0
)

maka_test(stiffness1 1
	./stiffness "${TESTDATA}/square.dmg" "${TESTDATA}/square.smb" 1
)

maka_test(stiffness2 1
	./stiffness "${TESTDATA}/square.dmg" "${TESTDATA}/square.smb" 2
)

maka_test_single(solver-test solver-test.cc 1
  "${TESTDATA}/square.dmg"
  "${TESTDATA}/square.smb"
  0
  "${TESTDATA}/solve-data.rc"
)
