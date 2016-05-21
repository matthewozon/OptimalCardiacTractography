cmake_minimum_required(VERSION 2.8)

#Configuration du projet
project(fileIO)
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(. /usr/include/gdcm ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils )
link_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox/lib ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils/lib)

#Configuration de l'exécutable
#file(
#     	GLOB
#	source_files
#        *.h *.cpp
#)
add_library(
        fileIO
	SHARED
	C_load_fiber.cpp
	C_ReadAnalyze75.cpp
	C_ReadDicomDiffInfo.cpp
	C_readDicomDiffData.cpp
	C_getDicomFromDirectory.cpp
	C_sort_files.cpp
	C_maskMaker.cpp
	C_tensorMaker.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
        fileIO
        #pthread
        gdcmCommon
        gdcmMSFF
        gdcmIOD
        gdcmDSED
        gdcmDICT
        gdcmjpeg8
        gdcmjpeg16
        gdcmjpeg12
	toolBoxPuff
        utilsPuff
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
