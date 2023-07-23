#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "quadprog" for configuration "Debug"
set_property(TARGET quadprog APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(quadprog PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libquadprog.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS quadprog )
list(APPEND _IMPORT_CHECK_FILES_FOR_quadprog "${_IMPORT_PREFIX}/lib/libquadprog.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
