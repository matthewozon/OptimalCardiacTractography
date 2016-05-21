cmake_minimum_required(VERSION 2.6)

#Configuration du projet
project(dataCylinder)
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/synthetic ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils) # lib/include /usr/include/gdcm)
link_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/synthetic/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils/lib) #lib/lib /usr/lib64)

#Configuration de l'exécutable
file(
	GLOB
	source_files
	main.cpp
)
add_executable(
	dataCylinder
	${source_files}
)

#Configuration de l'édition de liens
target_link_libraries(
	dataCylinder
	utilsPuff
	#fiberPuff
	syntheticData
	fileIO
	pthread
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
