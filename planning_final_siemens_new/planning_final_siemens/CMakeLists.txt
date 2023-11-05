# make a directory build, go there and run “cmake ..”.
# may need to specify generator in Windows:
# e.g.: for Visual Studio 2013 on x64, run cmake -G "Visual Studio 12 Win64" ..
# need separate build folders for Win32 and Win64
# The makefile will be created there.

CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

#PROJECT(CRICO_CPP)
project(QuadProg++ CXX)
set(LIBRARY_OUTPUT_PATH ./build)

#set(CMAKE_CXX_STANDARD 14)

INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(src/fmincon)
INCLUDE_DIRECTORIES(src/lib)
INCLUDE_DIRECTORIES(src/lib/QP_lib)
INCLUDE_DIRECTORIES(src/include)
#INCLUDE_DIRECTORIES(/usr/local/include/optim)
INCLUDE_DIRECTORIES(/usr/local/include/Eigen)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)

add_subdirectory(src)
# add_subdirectory(src/lib)

set(yaml_cxx_flags "-Wall ${GCC_EXTRA_OPTIONS} -pedantic -Wno-long-long -Wno-unused-variable -std=c++11 ${yaml_cxx_flags}")

add_executable(safe_polish src/safe_polish_main.cpp ${headers})
set_target_properties(safe_polish PROPERTIES
  COMPILE_FLAGS "${yaml_cxx_flags}"
)
target_link_directories(safe_polish PRIVATE /usr/local/lib)

target_link_libraries(safe_polish nlopt quadprog qhullcpp qhullstatic_r yaml-cpp armadillo openblas ${CMAKE_DL_LIBS} ${ARMADILLO_LIBRARIES})

