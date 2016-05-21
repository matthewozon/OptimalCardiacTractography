cmake_minimum_required(VERSION 2.8)

#Configuration du projet
project(syntheticData)
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils)
link_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils/lib)

#Configuration de l'exécutable
#file(
#     	GLOB
#	source_files
#        *.h *.cpp
#)
add_library(
        syntheticData
	SHARED
	C_synthetic_DWI.cpp
	C_synthetic_DWI_cylinder.cpp
	C_synthetic_DWI_LV.cpp
	C_synthetic_DWI_patho.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
        syntheticData
	toolBoxPuff
        utilsPuff
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
