cmake_minimum_required(VERSION 2.8)

#Configuration du projet
project(similarityRawTranslation)
#graphCount
#distributionCurvatureLength
#distributionU2meanLength
#distributionU2
#distributionCurvature
#meanPosition
#meanErrorFiberRawTranslation
#(meanPosition)
#fiberMeasureSpline
#meanErrorFiberLenRawTranslation
#similarityRawTranslation
#similaritySplineTranslation
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

#les path pour les header
include_directories(.
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/tractography 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/synthetic 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils
)

#les path pour les
link_directories(. 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox/lib 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO/lib  
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/tractography/lib 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/synthetic/lib 
	~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils/lib
)

#Configuration de l'exécutable
file(
	GLOB
	source_files
	*.h *.cpp
)
add_executable(
        similarityRawTranslation
	${source_files}
)

target_link_libraries(
        similarityRawTranslation
	utilsPuff #some tools that can be of interest like image container rawData.h
	tractPuff #the actual tractigraphy library
	syntheticData #synthetic data generation and tensor computation
	fileIO #load data form files
        toolBoxPuff #toolbox... well for quite many little thing to be all right
	pthread #multithread library of POSIX
)


if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
