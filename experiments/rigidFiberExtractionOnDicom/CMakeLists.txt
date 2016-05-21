cmake_minimum_required(VERSION 2.6)

#Configuration du projet
project(rigidFiberExtractionOnDicom100pts)
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(. /usr/include/gdcm ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO  ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/tractography ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/synthetic ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils)
link_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO/lib  ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/tractography/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/synthetic/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils/lib)

#Configuration de l'exécutable
file(
	GLOB
	source_files
	main.cpp
)
add_executable(
	rigidFiberExtractionOnDicom100pts
	${source_files}
)

#Configuration de l'édition de liens
target_link_libraries(
	rigidFiberExtractionOnDicom100pts
	utilsPuff #some tools that can be of interest like image container rawData.h
	tractPuff #the actual tractigraphy library
	syntheticData #synthetic data generation and tensor computation
	fileIO #load data form files
	pthread #multithread library of POSIX
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()



