################################################################################
# General Informations
################################################################################

CMake_Minimum_Required (VERSION 3.0)
Project (MeshFix)

################################################################################

If (NOT CMAKE_BUILD_TYPE)
	Message (STATUS "No build type selected, default to Release")
	Set (CMAKE_BUILD_TYPE "Release")
Endif ()

################################################################################

Set (SOURCES
		src/TMesh/edge.cpp
		src/TMesh/vertex.cpp
		src/TMesh/io.cpp
		src/TMesh/triangle.cpp
		src/TMesh/tin.cpp
		src/Algorithms/detectIntersections.cpp
		src/Algorithms/holeFilling.cpp
		src/Algorithms/marchIntersections.cpp
		src/Algorithms/checkAndRepair.cpp
		src/Kernel/heap.cpp
		src/Kernel/matrix.cpp
		src/Kernel/orientation.c
		src/Kernel/list.cpp
		src/Kernel/coordinates.cpp
		src/Kernel/tmesh.cpp
		src/Kernel/graph.cpp
		src/Kernel/point.cpp
		src/Kernel/jqsort.cpp
		src/MeshFix/meshfix.cpp
)

################################################################################

Add_Executable (MeshFix ${SOURCES})

Target_Include_Directories(MeshFix PUBLIC
		include/TMesh
		include/Kernel
)

################################################################################

If (CMAKE_SIZEOF_VOID_P EQUAL 8)
	Message (STATUS "64 bits compiler detected")
	Target_Compile_Definitions (MeshFix PUBLIC -DIS64BITPLATFORM)
	Set (BIN_FOLDER "bin64")
Else ()
	Message (STATUS "32 bits compiler detected")
	Set (BIN_FOLDER "bin")
Endif ()

################################################################################

Set_Target_Properties (MeshFix
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${BIN_FOLDER}"
)

# Use C++11
Set_Target_Properties(MeshFix PROPERTIES CXX_STANDARD 11)
Set_Target_Properties(MeshFix PROPERTIES CXX_STANDARD_REQUIRED ON)
