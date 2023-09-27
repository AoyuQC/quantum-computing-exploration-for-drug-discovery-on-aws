# libint2-config.cmake
# --------------------
#
# Libint2 cmake module.
# This module sets the following variables in your project:
#
# ::
#
#   Libint2_FOUND - true if Libint2 and all required components found on the system
#   Libint2_VERSION - Libint2 version in format Major.Minor.Release
#   Libint2_EXT_VERSION - Libint2 version including the (optional) buildid, such as beta.3
#   Libint2_MAX_AM_ERI - maximum angular momentum level of Libint2 libraries
#
#
# Available components:
#
# ::
#
#   shared - search for only shared library
#   static - search for only static library
#
#   onebody_dD_lL - search for library including 1-body integrals with derivative order D (D=0..4) and max angular momentum up to L (L=2..10)
#   eri_cC_dD_lL  - search for library including 2-body integrals with C (C=2,3,4) centers, derivative order D (D=0..4), and max angular momentum up to L (L=2..10)
#   g12_dD-lL - search for library including F12 integrals with Gaussian factors with derivative order D and max angular momentum up to L
#   g12dkh_dD-lL - search for library including F12 integrals with Gaussian factors and DKH with derivative order D and max angular momentum up to L (NYI)
#
#   impure_sh - search for library that doesn't assume 2- and 3-center integrals involve pure solid harmonics
#
#                    sph        cart       shell_set  used_by
#                    --------   --------   ---------  -------
#   sss - search for standard + standard + standard = mpqc4, cp2k
#   sso - search for                     + orca
#   sis - search for          + intv3    + standard = mpqc3
#   sio - search for                     + orca
#   sgs - search for          + gamess   + standard = gamess
#   sgo - search for                     + orca
#   sos - search for          + orca     + standard
#   soo - search for                     + orca     = orca
#   sbs - search for          + bagel    + standard = bagel
#   sbo - search for                     + orca
#   gss - search for gaussian + standard + standard = psi4 (since v1.4)
#   gso - search for                     + orca
#   gis - search for          + intv3    + standard
#   gio - search for                     + orca
#   ggs - search for          + gamess   + standard
#   ggo - search for                     + orca
#   gos - search for          + orca     + standard
#   goo - search for                     + orca
#   gbs - search for          + bagel    + standard
#   gbo - search for                     + orca
#
#   C - search for at least Libint2::int2 target
#   CXX_ho - search for at least Libint2::cxx target
#   CXX - search for at least Libint2::int2-cxx target
#   Fortran - search for at least libint_f target (NYI)
#
#
# Exported targets:
#
# ::
#
# If Libint2 is found and no language components are requested, this module
# defines at least the following :prop_tgt:`IMPORTED` target. ::
#
#   Libint2::int2 - library with C API
#
# If Libint2 is found, depending on components requested, available
# dependencies, and fullness of the installation, this module defines up to the
# following :prop_tgt:`IMPORTED` targets. ::
#
#   Libint2::int2 - library with C API (COMPONENT C)
#   Libint2::cxx - Libint2::int2 plus interface to header-only C++11 API (COMPONENT CXX_ho)
#   Libint2::int2-cxx - Libint2::int2 plus compiled C++11 API (COMPONENT CXX)
#   Libint2::fortran (NYI)
#
#
# Suggested usage:
#
# ::
#
#   find_package(Libint2)
#   find_package(Libint2 2.7.1 CONFIG REQUIRED COMPONENTS shared sss eri_c4_d0_l5 eri_c4_d1_l4)
#
#
# The following variables can be set to guide the search for this package:
#
# ::
#
#   Libint2_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package
#   CMAKE_DISABLE_FIND_PACKAGE_Libint2 - CMake variable, disables
#     find_package(Libint2) when not REQUIRED, perhaps to force internal build


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was libint2-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set(pnv libint2)  # projectnameversion
set(L2 Libint2)  # NameSpace

set(Libint2_EXT_VERSION "2.7.2-post1")

# make detectable the various cmake modules exported alongside
# * prepend to trump any pre-target FindEigen3.cmake modules lying around
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

# check library style component
if(1)  # BUILD_SHARED_LIBS
    set(${L2}_shared_FOUND 1)
endif()
if()  # BUILD_STATIC_LIBS
    set(${L2}_static_FOUND 1)
endif()
list(FIND ${L2}_FIND_COMPONENTS "shared" _seek_shared)
list(FIND ${L2}_FIND_COMPONENTS "static" _seek_static)

# check library language component
include(CMakeFindDependencyMacro)

set(${L2}_C_FOUND 1)
list(FIND ${L2}_FIND_COMPONENTS "C" _seek_C)
if(ON)  # LIBINT2_REQUIRE_CXX_API
    if(NOT TARGET Eigen3::Eigen)
        find_dependency(Eigen3 REQUIRED)
    endif()

    if (1)
        # Boost headers _not_ unpacked to within `include/libint2/`
        if (NOT TARGET Boost::headers)
            find_dependency(Boost 1.57 REQUIRED)
        endif()
    else()
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "Boost detected. satisfied by headers bundled with ${L2} distribution")
        endif()
    endif()

    set(${L2}_CXX_ho_FOUND 1)

    if(OFF)  # LIBINT2_REQUIRE_CXX_API_COMPILED
        set(${L2}_CXX_FOUND 1)
    endif()
endif()
list(FIND ${L2}_FIND_COMPONENTS "CXX_ho" _seek_CXX_ho)
list(FIND ${L2}_FIND_COMPONENTS "CXX" _seek_CXX)
if(OFF)  # LIBINT2_ENABLE_FORTRAN
    set(${L2}_Fortran_FOUND 1)
endif()
list(FIND ${L2}_FIND_COMPONENTS "Fortran" _seek_Fortran)

# check AM & derivative component
set(${L2}_MAX_AM_ERI 5)
foreach(_eri eri_c4_d0_l2;eri_c4_d0_l3;eri_c4_d0_l4;eri_c4_d0_l5;eri_c4_d1_l2;eri_c4_d1_l3;eri_c4_d1_l4;eri_c4_d2_l2;eri_c4_d2_l3;eri_c3_d0_l2;eri_c3_d0_l3;eri_c3_d0_l4;eri_c3_d0_l5;eri_c3_d0_l6;eri_c3_d1_l2;eri_c3_d1_l3;eri_c3_d1_l4;eri_c3_d1_l5;eri_c3_d2_l2;eri_c3_d2_l3;eri_c3_d2_l4;eri_c2_d0_l2;eri_c2_d0_l3;eri_c2_d0_l4;eri_c2_d0_l5;eri_c2_d0_l6;eri_c2_d1_l2;eri_c2_d1_l3;eri_c2_d1_l4;eri_c2_d1_l5;eri_c2_d2_l2;eri_c2_d2_l3;eri_c2_d2_l4;onebody_d0_l2;onebody_d0_l3;onebody_d0_l4;onebody_d0_l5;onebody_d0_l6;onebody_d1_l2;onebody_d1_l3;onebody_d1_l4;onebody_d1_l5;onebody_d2_l2;onebody_d2_l3;onebody_d2_l4;g12_d0_l2;g12_d0_l3;g12_d0_l4;g12_d1_l2;g12_d1_l3;g12_d1_l4)
    set(${L2}_${_eri}_FOUND 1)
endforeach()

# check pure restriction component
if((0 EQUAL 0) AND (0 EQUAL 0))  # ERI3/ERI2_PURE_SH
    set(${L2}_impure_sh_FOUND 1)
endif()
list(FIND ${L2}_FIND_COMPONENTS "impure_sh" _seek_impure_sh)

# check orderings component:
# LIBINT_SHGSHELL_ORDERING = 1
# LIBINT_CGSHELL_ORDERING = 1
# LIBINT_SHELL_SET = 1
set(${L2}_sss_FOUND 1)


# thanks, https://stackoverflow.com/a/9328525
function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    set(founds "")
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
            if (NOT ${${_variableName}})
                continue()
            endif()
        endif()
        list(APPEND found ${CMAKE_MATCH_1})
    endforeach()
    message(STATUS "${ARGV1}${found}")
endfunction()


if(NOT CMAKE_REQUIRED_QUIET)
    list(SORT ${L2}_FIND_COMPONENTS COMPARE STRING)
    message(STATUS "${L2}Config components requested: ${${L2}_FIND_COMPONENTS}")
    dump_cmake_variables("^Libint2_([A-Za-z0-9_]+)_FOUND$" "${L2}Config components found: ")
endif()

check_required_components(${L2})

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built this as a subproject
#-----------------------------------------------------------------------------
if(NOT TARGET ${L2}::int2)
    if(_seek_static GREATER -1)
        include("${CMAKE_CURRENT_LIST_DIR}/${pnv}-targets-static.cmake")
    elseif(_seek_shared GREATER -1)
        include("${CMAKE_CURRENT_LIST_DIR}/${pnv}-targets-shared.cmake")
    elseif(1)  # BUILD_SHARED_LIBS
        include("${CMAKE_CURRENT_LIST_DIR}/${pnv}-targets-shared.cmake")
    elseif()  # BUILD_STATIC_LIBS
        include("${CMAKE_CURRENT_LIST_DIR}/${pnv}-targets-static.cmake")
    endif()

    get_property(_loc TARGET ${L2}::int2 PROPERTY LOCATION)
    get_property(_ill TARGET ${L2}::int2 PROPERTY INTERFACE_LINK_LIBRARIES)
    get_property(_id TARGET ${L2}::int2 PROPERTY INCLUDE_DIRECTORIES)
    get_property(_iid TARGET ${L2}::int2 PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    message(DEBUG "Libint2::int2")
    message(DEBUG "loc ${_loc}")
    message(DEBUG "ill ${_ill}")
    message(DEBUG "id  ${_id}")
    message(DEBUG "iid ${_iid}")
endif()
