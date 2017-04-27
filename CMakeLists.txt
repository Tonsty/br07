cmake_minimum_required(VERSION 2.8)

PROJECT (br07)

INCLUDE_DIRECTORIES(
	"D:/glut-3.7/include"
	"C:/Program Files/GMNR/include"
	${PROJECT_SOURCE_DIR}/include 
	${PROJECT_SOURCE_DIR}/third_party/trimesh2_2.11/include
	${PROJECT_SOURCE_DIR}/third_party/prebuilt-dll-2-9-1-release/include
	${PROJECT_SOURCE_DIR}/third_party/tnt_126
	${PROJECT_SOURCE_DIR}/third_party/jama125
	${PROJECT_SOURCE_DIR}/third_party/glui-2.35/src/include/GL)
LINK_DIRECTORIES(
	"D:/glut-3.7/lib/x64"
	${PROJECT_SOURCE_DIR}/third_party/trimesh2_2.11/lib.Win64
	${PROJECT_SOURCE_DIR}/third_party/prebuilt-dll-2-9-1-release/lib/x64
	${PROJECT_SOURCE_DIR}/third_party/glui-2.35/src/lib.Win64
	"C:/Program Files/freeglut/lib"
	"C:/Program Files/GMNR/lib")

ADD_DEFINITIONS(-DLIBMATH_DYNAMIC -DLIBTPS_DYNAMIC -DLIBMULTITPS_DYNAMIC)

find_package(PCL 1.8 REQUIRED COMPONENTS common io filters)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

ADD_EXECUTABLE(autorient
					./src/autorient.cc)
TARGET_LINK_LIBRARIES(autorient debug trimeshd optimized trimesh)

ADD_EXECUTABLE(correspond
					./include/ICP.h
					./include/EdgeMesh.h
					./include/rand48.h 
					./src/correspond.cc
					./src/ICP.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(correspond debug trimeshd optimized trimesh)

ADD_EXECUTABLE(distance
					./src/distance.cc)
TARGET_LINK_LIBRARIES(distance debug trimeshd optimized trimesh)

ADD_EXECUTABLE(distance2
					./src/distance2.cc)
TARGET_LINK_LIBRARIES(distance2 debug trimeshd optimized trimesh)

ADD_EXECUTABLE(ply2vtk
				   ./src/ply2vtk.cc)
TARGET_LINK_LIBRARIES(ply2vtk ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES})

ADD_EXECUTABLE(vtks2vtk
					./src/vtks2vtk.cc)
TARGET_LINK_LIBRARIES(vtks2vtk ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES})

ADD_EXECUTABLE(distort
				   ./src/distort.cc)
TARGET_LINK_LIBRARIES(distort debug trimeshd optimized trimesh)

ADD_EXECUTABLE(global_reg 
					./include/global_reg.h
					./include/tps.h
					./src/global_reg.cc
					./src/global_reg_util.cc
					./src/tps.cc)
TARGET_LINK_LIBRARIES(global_reg pthreadVC2 debug trimeshd optimized trimesh)

ADD_EXECUTABLE(global_reg2 
					./include/global_reg.h
					./include/tps.h
					./src/global_reg2.cc
					./src/global_reg_util.cc
					./src/tps.cc)
TARGET_LINK_LIBRARIES(global_reg2 pthreadVC2 math tps multitps debug trimeshd optimized trimesh)

ADD_EXECUTABLE(interactive 
					./include/ICP.h
					./include/EdgeMesh.h
					./include/rand48.h 
					./src/interactive.cc
					./src/ICP.cc
					./src/EdgeMesh.cc)
SET_TARGET_PROPERTIES(interactive 
							PROPERTIES 
								LINK_FLAGS_DEBUG /LIBPATH:"D:/glut-3.7/lib/x64/Debug"
								LINK_FLAGS_RELEASE /LIBPATH:"D:/glut-3.7/lib/x64/Release"
								LINK_FLAGS_RELWITHDEBINFO /LIBPATH:"D:/glut-3.7/lib/x64/Release")
TARGET_LINK_LIBRARIES(interactive debug glui32d optimized glui32 debug trimeshd optimized trimesh ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES})

ADD_EXECUTABLE(make_session
					./include/EdgeMesh.h
					./src/make_session.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(make_session debug trimeshd optimized trimesh)


ADD_EXECUTABLE(ply2vn ./src/ply2vn.cc)
TARGET_LINK_LIBRARIES(ply2vn debug trimeshd optimized trimesh)

ADD_EXECUTABLE(preprocess 
					./include/EdgeMesh.h
					./include/rand48.h 
					./include/ICP.h
					./src/preprocess.cc
					./src/ICP.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(preprocess debug trimeshd optimized trimesh ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${PCL_FILTERS_LIBRARIES})

ADD_EXECUTABLE(pre_to_ply
					./include/EdgeMesh.h
					./src/pre_to_ply.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(pre_to_ply debug trimeshd optimized trimesh)






