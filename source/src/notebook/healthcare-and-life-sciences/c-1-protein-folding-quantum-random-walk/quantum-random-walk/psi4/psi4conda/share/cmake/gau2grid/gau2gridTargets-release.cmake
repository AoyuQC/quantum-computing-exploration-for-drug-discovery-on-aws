#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "gau2grid::gg" for configuration "Release"
set_property(TARGET gau2grid::gg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(gau2grid::gg PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgg.so.2"
  IMPORTED_SONAME_RELEASE "libgg.so.2"
  )

list(APPEND _cmake_import_check_targets gau2grid::gg )
list(APPEND _cmake_import_check_files_for_gau2grid::gg "${_IMPORT_PREFIX}/lib/libgg.so.2" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
