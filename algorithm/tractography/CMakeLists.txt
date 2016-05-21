cmake_minimum_required(VERSION 2.8) #en fait je n'ai absolument aucune idee de la version qu'il faut, mais c'etait la dans l'exemple que j'avais et c'est reste

#Configuration du projet
project(tractPuff)
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

#path des headers
include_directories(.
    ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox
    ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO
    ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/defs
    ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils
)
#path des lib
link_directories(.
    ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/toolbox/lib
    ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/algorithm/fileIO/lib
    ~/Documents/codeblocks_workspace/ozon_code_thesis_tracto/utils/lib
)

#les sources
add_library(
        tractPuff
	SHARED
	C_point_array_data.cpp
	C_cooling_schedule.cpp
	C_bundle.cpp
	C_toolbox_interpolate.cpp
	C_graph.cpp
	C_SA.cpp
	C_temperature.cpp
	C_SC.cpp
	C_move.cpp
	C_energy.cpp
	C_minimization.cpp
	C_ICM.cpp
	C_streamline.cpp
	C_thread_deltaU_move_vertex.cpp
	C_thread_init_energy.cpp
	C_measure.cpp
)

#les libraries a inclure
target_link_libraries(
        tractPuff
	toolBoxPuff
	fileIO
        utilsPuff
)

#quelques options pour gcc
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
