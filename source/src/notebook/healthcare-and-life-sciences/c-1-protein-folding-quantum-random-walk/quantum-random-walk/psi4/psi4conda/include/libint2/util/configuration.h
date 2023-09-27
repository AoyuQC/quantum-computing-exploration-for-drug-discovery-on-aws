/*
 *  Copyright (C) 2004-2023 Edward F. Valeev
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

#ifndef _libint2_include_libint2_util_configuration_h_
#define _libint2_include_libint2_util_configuration_h_

#include <string>

namespace libint2 {

  /// Runtime accessor for the library configuration:
  /// integral derivatives, AM, orderings, etc.
  /// @return the semicolon-separated strings from CMake components
  inline std::string configuration_accessor() {
    std::string components = "eri_c2_d0_l2;eri_c2_d0_l3;eri_c2_d0_l4;eri_c2_d0_l5;eri_c2_d0_l6;eri_c2_d1_l2;eri_c2_d1_l3;eri_c2_d1_l4;eri_c2_d1_l5;eri_c2_d2_l2;eri_c2_d2_l3;eri_c2_d2_l4;eri_c3_d0_l2;eri_c3_d0_l3;eri_c3_d0_l4;eri_c3_d0_l5;eri_c3_d0_l6;eri_c3_d1_l2;eri_c3_d1_l3;eri_c3_d1_l4;eri_c3_d1_l5;eri_c3_d2_l2;eri_c3_d2_l3;eri_c3_d2_l4;eri_c4_d0_l2;eri_c4_d0_l3;eri_c4_d0_l4;eri_c4_d0_l5;eri_c4_d1_l2;eri_c4_d1_l3;eri_c4_d1_l4;eri_c4_d2_l2;eri_c4_d2_l3;g12_d0_l2;g12_d0_l3;g12_d0_l4;g12_d1_l2;g12_d1_l3;g12_d1_l4;impure_sh;onebody_d0_l2;onebody_d0_l3;onebody_d0_l4;onebody_d0_l5;onebody_d0_l6;onebody_d1_l2;onebody_d1_l3;onebody_d1_l4;onebody_d1_l5;onebody_d2_l2;onebody_d2_l3;onebody_d2_l4;sss";
    return components;
  }

}

#endif /* header guard */
