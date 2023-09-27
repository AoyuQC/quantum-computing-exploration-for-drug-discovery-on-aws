#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Libint2::int2" for configuration "Release"
set_property(TARGET Libint2::int2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Libint2::int2 PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libint2.so.2"
  IMPORTED_SONAME_RELEASE "libint2.so.2"
  )

list(APPEND _cmake_import_check_targets Libint2::int2 )
list(APPEND _cmake_import_check_files_for_Libint2::int2 "${_IMPORT_PREFIX}/lib/libint2.so.2" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
