# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

file(GLOB sources lib/QP_lib/*.cc)
file(GLOB headers lib/QP_lib/*.hh)
# file(GLOB yaml_sources "lib/yaml/src/[a-zA-Z]*.cpp")
# file(GLOB_RECURSE public_headers "lib/yaml/include/yaml-cpp/[a-zA-Z]*.h")

INCLUDE_DIRECTORIES(lib/QP_lib)
# INCLUDE_DIRECTORIES(lib/yaml)
# INCLUDE_DIRECTORIES(lib/yaml/include)
# INCLUDE_DIRECTORIES(lib/yaml/src)

set(library_sources ${sources})
list(REMOVE_ITEM library_sources ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_library(quadprog STATIC ${library_sources} ${headers})
set_property(TARGET quadprog PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET quadprog PROPERTY PUBLIC_HEADER ${headers})


install(TARGETS quadprog
  EXPORT quadprog-targets
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/QuadProg++
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(EXPORT quadprog-targets DESTINATION cmake)
# add_subdirectory(lib)

set(yaml_cxx_flags "-Wall ${GCC_EXTRA_OPTIONS} -pedantic -Wno-long-long -std=c++11 ${yaml_cxx_flags}")

file(GLOB yaml_sources "src/[a-zA-Z]*.cpp")
file(GLOB_RECURSE public_headers "include/yaml-cpp/[a-zA-Z]*.h")


set(yaml_library_sources
  ${yaml_sources}
  ${public_headers}
)

add_library(yaml-cpp ${yaml_library_sources})
set_target_properties(yaml-cpp PROPERTIES
  COMPILE_FLAGS "${yaml_cxx_flags}"
)