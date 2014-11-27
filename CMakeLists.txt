# default build-type (release)
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

# check cmake version
cmake_minimum_required(VERSION 2.8)

# project name
project (MPLP)

# set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# set source files
SET(MPLP_SRC_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/muldim_arr.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/read_model_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mplp_alg.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/matrix.cpp
)

# create library
add_library (mplp ${MPLP_SRC_FILES})
add_library (mplp-shared SHARED ${MPLP_SRC_FILES})

# point to include directory
target_include_directories (mplp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories (mplp-shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# build example
add_executable (mplpSolver ${CMAKE_CURRENT_SOURCE_DIR}/src/cycle_tighten_main.cpp)
# point to include directory
target_include_directories (mplpSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# link mplp library
target_link_libraries (mplpSolver LINK_PUBLIC mplp)

# install
if(NOT DEFINED MPLP_INSTALL_LIB_DESTINATION)
   set(MPLP_INSTALL_LIB_DESTINATION lib)
endif(NOT DEFINED MPLP_INSTALL_LIB_DESTINATION)
if(NOT DEFINED MPLP_INSTALL_HEADER_DESTINATION)
   set(MPLP_INSTALL_HEADER_DESTINATION include)
endif(NOT DEFINED MPLP_INSTALL_HEADER_DESTINATION)

install(TARGETS mplp DESTINATION ${MPLP_INSTALL_LIB_DESTINATION})
install(TARGETS mplp-shared DESTINATION ${MPLP_INSTALL_LIB_DESTINATION})
install(DIRECTORY include/MPLP DESTINATION ${MPLP_INSTALL_HEADER_DESTINATION})