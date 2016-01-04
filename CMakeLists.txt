cmake_minimum_required(VERSION 2.8)

PROJECT (br07)

INCLUDE_DIRECTORIES(
	"D:/glut-3.7/include"
	${PROJECT_SOURCE_DIR}/include 
	${PROJECT_SOURCE_DIR}/third_party/trimesh2_2.11/include
	${PROJECT_SOURCE_DIR}/third_party/Pre-built.2/include
	${PROJECT_SOURCE_DIR}/third_party/tnt_126
	${PROJECT_SOURCE_DIR}/third_party/jama125
	${PROJECT_SOURCE_DIR}/third_party/glui-2.35/src/include/GL)
LINK_DIRECTORIES(
	${PROJECT_SOURCE_DIR}/third_party/trimesh2_2.11/lib.Win32
	${PROJECT_SOURCE_DIR}/third_party/Pre-built.2/lib
	${PROJECT_SOURCE_DIR}/third_party/glui-2.35/src/msvc/lib
	"C:/Program Files (x86)/freeglut/lib")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

ADD_EXECUTABLE(autorient
					./src/autorient.cc)
TARGET_LINK_LIBRARIES(autorient debug trimeshd optimized trimesh)

ADD_EXECUTABLE(correspond 
					./src/correspond.cc
					./src/ICP.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(correspond debug trimeshd optimized trimesh)

ADD_EXECUTABLE(distance
					./src/distance.cc)
TARGET_LINK_LIBRARIES(distance debug trimeshd optimized trimesh)

ADD_EXECUTABLE(global_reg 
					./src/global_reg.cc
					./src/global_reg_util.cc
					./src/tps.cc)
TARGET_LINK_LIBRARIES(global_reg pthreadVC2 debug trimeshd optimized trimesh)

ADD_EXECUTABLE(interactive 
					./src/interactive.cc
					./src/ICP.cc
					./src/EdgeMesh.cc)
SET_TARGET_PROPERTIES(interactive 
							PROPERTIES 
								LINK_FLAGS_DEBUG /LIBPATH:"D:/glut-3.7/lib/win32/Debug"
								LINK_FLAGS_RELEASE /LIBPATH:"D:/glut-3.7/lib/win32/Release"
								LINK_FLAGS_RELWITHDEBINFO /LIBPATH:"D:/glut-3.7/lib/win32/Release")
TARGET_LINK_LIBRARIES(interactive debug glui32d optimized glui32 debug trimeshd optimized trimesh)

ADD_EXECUTABLE(make_session
					./src/make_session.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(make_session debug trimeshd optimized trimesh)


ADD_EXECUTABLE(ply2vn ./src/ply2vn.cc)
TARGET_LINK_LIBRARIES(ply2vn debug trimeshd optimized trimesh)

ADD_EXECUTABLE(preprocess 
					./src/preprocess.cc
					./src/ICP.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(preprocess debug trimeshd optimized trimesh)

ADD_EXECUTABLE(pre_to_ply
					./src/pre_to_ply.cc
					./src/EdgeMesh.cc)
TARGET_LINK_LIBRARIES(pre_to_ply debug trimeshd optimized trimesh)






