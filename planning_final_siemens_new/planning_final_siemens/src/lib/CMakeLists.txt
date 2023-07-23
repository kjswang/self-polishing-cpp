# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

# INCLUDE_DIRECTORIES(yaml)
# INCLUDE_DIRECTORIES(yaml/include)
# INCLUDE_DIRECTORIES(yaml/src)

set(yaml_cxx_flags "-Wall ${GCC_EXTRA_OPTIONS} -pedantic -Wno-long-long -std=c++11 ${yaml_cxx_flags}")

file(GLOB yaml_sources "yaml/src/[a-zA-Z]*.cpp")
file(GLOB_RECURSE public_headers "yaml/include/yaml-cpp/[a-zA-Z]*.h")


set(yaml_library_sources
  ${yaml_sources}
  ${public_headers}
)

add_library(yaml-cpp ${yaml_library_sources})
set_target_properties(yaml-cpp PROPERTIES
  COMPILE_FLAGS "${yaml_cxx_flags}"
)
