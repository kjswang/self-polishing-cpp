# make a directory build, go there and run “cmake ..”.
# may need to specify generator in Windows:
# e.g.: for Visual Studio 2013 on x64, run cmake -G "Visual Studio 12 Win64" ..
# need separate build folders for Win32 and Win64
# The makefile will be created there.

CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

#PROJECT(CRICO_CPP)
project(QuadProg++ CXX)

INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(src/lib)
INCLUDE_DIRECTORIES(src/lib/QP_lib)
INCLUDE_DIRECTORIES(src/include)



find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(src)
# add_subdirectory(src/lib)

set(yaml_cxx_flags "-Wall ${GCC_EXTRA_OPTIONS} -pedantic -Wno-long-long -Wno-unused-variable -std=c++11 ${yaml_cxx_flags}")

add_executable(safe_polish src/safe_polish_main.cpp ${headers})
set_target_properties(safe_polish PROPERTIES
  COMPILE_FLAGS "${yaml_cxx_flags}"
)
target_link_libraries(safe_polish quadprog Eigen3::Eigen yaml-cpp ${CMAKE_DL_LIBS})

