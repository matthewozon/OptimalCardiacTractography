cmake_minimum_required(VERSION 2.8)

#Configuration du projet
project(toolBoxPuff)
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils)
link_directories(. ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils/lib)

#Configuration de l'exécutable
#file(
#     	GLOB
#	source_files
#        *.h *.cpp
#)
add_library(
        toolBoxPuff
	SHARED
	C_toolbox_SVD.cpp
	C_toolbox_rand2.cpp
	C_toolbox_integrate.cpp
	C_toolbox_eigen_sym.cpp
	C_toolbox.cpp
	C_thread_spline.cpp
	C_toolbox_spline.cpp
	C_toolbox_PCA_3D.cpp
	C_kmean2D.cpp
	C_kmean3D.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
        toolBoxPuff
        utilsPuff
        #pthread
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
