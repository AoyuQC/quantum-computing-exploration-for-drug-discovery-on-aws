#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "dkh::dkh" for configuration "Release"
set_property(TARGET dkh::dkh APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(dkh::dkh PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdkh.so"
  IMPORTED_SONAME_RELEASE "libdkh.so"
  )

list(APPEND _cmake_import_check_targets dkh::dkh )
list(APPEND _cmake_import_check_files_for_dkh::dkh "${_IMPORT_PREFIX}/lib/libdkh.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
