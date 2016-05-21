cmake_minimum_required(VERSION 2.6)

#Configuration du projet
project(puffinVisuBall)
set(EXECUTABLE_OUTPUT_PATH binCmake/${CMAKE_BUILD_TYPE})

#Inclusion de Boost
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(. lib/include)
link_directories(. lib/lib)

#Configuration de l'exécutable
file(
	GLOB
	source_files
	*.h *.cpp
)
add_executable(
	puffinVisuBall
	${source_files}
)

#Configuration de l'édition de liens
if(VTK_LIBRARIES)
  target_link_libraries(
	puffinVisuBall
	fiberPuff
	${VTK_LIBRARIES})
else()
  target_link_libraries(
	puffinVisuBall
	fiberPuff
	 vtkHybrid)
endif()

