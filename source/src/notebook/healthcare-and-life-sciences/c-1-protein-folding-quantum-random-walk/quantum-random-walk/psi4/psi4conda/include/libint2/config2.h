/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* This file is automatically processed by configure script.
   It MUST NOT be changed manually after configuration, otherwise
   the library will likely fail to compile or produce erroneous results
 */

/* EXTRA DEFINES DETERMINED BY CONFIGURE OF THE EXPORTED LIBRARY */

#ifndef _libint2_include_libint2_config_h_1
#define _libint2_include_libint2_config_h_1

#undef LIBINT_ALIGN_SIZE
#undef LIBINT_HAS_MPFR
#undef HAVE_POSIX_MEMALIGN
#undef LIBINT_USER_DEFINED_REAL
/* retired Jan 2023 #undef LIBINT_SHGSHELL_ORDERING */
#undef LIBINT_HAS_EIGEN
#undef LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS

/* if can be controlled with posix_memalign, alignment size */
#define LIBINT2_ALIGN_SIZE 0

/* User-defined real type */
#define LIBINT2_REALTYPE double

/* Specifies the ordering of solid harmonics Gaussians in a shell. Retired Jan 2023 in favor of becoming a generator-time-only option */
/* #define LIBINT_SHGSHELL_ORDERING 1 */

/* define if Eigen library is available. */
#define LIBINT_HAS_EIGEN 1

/* define if system-wide Boost.Preprocessor is available */
#define LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 1

/* have posix_memalign ? */
#define HAVE_POSIX_MEMALIGN 1

#endif /* header guard */
