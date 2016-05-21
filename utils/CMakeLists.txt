cmake_minimum_required(VERSION 2.8)

#Configuration du projet
project(utilsPuff)
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils)
link_directories(.)

#Configuration de l'exécutable
#file(
#     	GLOB
#	source_files
#        *.h *.cpp
#)
add_library(
        utilsPuff
	SHARED
	C_process.cpp
	C_thread.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
        utilsPuff
        pthread
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
