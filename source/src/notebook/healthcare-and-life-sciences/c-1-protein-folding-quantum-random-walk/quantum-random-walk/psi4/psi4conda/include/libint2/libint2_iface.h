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

#ifndef _libint2_libint2iface_h_
#define _libint2_libint2iface_h_

#ifdef __cplusplus
# include <cstddef>
#else
# include <stddef.h>
#endif
#ifdef __cplusplus
LIBINT_PRAGMA_CLANG(diagnostic push)
LIBINT_PRAGMA_CLANG(diagnostic ignored "-Wunused-variable")
LIBINT_PRAGMA_GCC(diagnostic push)
LIBINT_PRAGMA_GCC(diagnostic ignored "-Wunused-variable")
extern "C" {
#endif
extern void (*libint2_build_default[7][7][7][7])(const Libint_t*);
extern void (*libint2_build_default1)(const Libint_t*);
extern void (*libint2_build_default2)(const Libint_t*);
extern void (*libint2_build_overlap[7][7])(const Libint_t*);
extern void (*libint2_build_kinetic[7][7])(const Libint_t*);
extern void (*libint2_build_elecpot[7][7])(const Libint_t*);
extern void (*libint2_build_1emultipole[7][7])(const Libint_t*);
extern void (*libint2_build_2emultipole[7][7])(const Libint_t*);
extern void (*libint2_build_3emultipole[7][7])(const Libint_t*);
extern void (*libint2_build_sphemultipole[7][7])(const Libint_t*);
extern void (*libint2_build_overlap1[7][7])(const Libint_t*);
extern void (*libint2_build_kinetic1[7][7])(const Libint_t*);
extern void (*libint2_build_elecpot1[7][7])(const Libint_t*);
extern void (*libint2_build_1emultipole1[7][7])(const Libint_t*);
extern void (*libint2_build_2emultipole1[7][7])(const Libint_t*);
extern void (*libint2_build_3emultipole1[7][7])(const Libint_t*);
extern void (*libint2_build_sphemultipole1[7][7])(const Libint_t*);
extern void (*libint2_build_overlap2[7][7])(const Libint_t*);
extern void (*libint2_build_kinetic2[7][7])(const Libint_t*);
extern void (*libint2_build_elecpot2[7][7])(const Libint_t*);
extern void (*libint2_build_1emultipole2[7][7])(const Libint_t*);
extern void (*libint2_build_2emultipole2[7][7])(const Libint_t*);
extern void (*libint2_build_3emultipole2[7][7])(const Libint_t*);
extern void (*libint2_build_sphemultipole2[7][7])(const Libint_t*);
extern void (*libint2_build_eri[6][6][6][6])(const Libint_t*);
extern void (*libint2_build_eri1[5][5][5][5])(const Libint_t*);
extern void (*libint2_build_eri2[4][4][4][4])(const Libint_t*);
extern void (*libint2_build_3eri[7][7][7])(const Libint_t*);
extern void (*libint2_build_3eri1[6][6][6])(const Libint_t*);
extern void (*libint2_build_3eri2[5][5][5])(const Libint_t*);
extern void (*libint2_build_2eri[7][7])(const Libint_t*);
extern void (*libint2_build_2eri1[6][6])(const Libint_t*);
extern void (*libint2_build_2eri2[5][5])(const Libint_t*);
extern void (*libint2_build_r12kg12[5][5][5][5])(const Libint_t*);
void libint2_static_init();
void libint2_static_cleanup();
void libint2_init_default(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_default(int max_am);
void libint2_cleanup_default(Libint_t* inteval);
void libint2_init_default1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_default1(int max_am);
void libint2_cleanup_default1(Libint_t* inteval);
void libint2_init_default2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_default2(int max_am);
void libint2_cleanup_default2(Libint_t* inteval);
void libint2_init_overlap(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_overlap(int max_am);
void libint2_cleanup_overlap(Libint_t* inteval);
void libint2_init_kinetic(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_kinetic(int max_am);
void libint2_cleanup_kinetic(Libint_t* inteval);
void libint2_init_elecpot(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_elecpot(int max_am);
void libint2_cleanup_elecpot(Libint_t* inteval);
void libint2_init_1emultipole(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_1emultipole(int max_am);
void libint2_cleanup_1emultipole(Libint_t* inteval);
void libint2_init_2emultipole(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_2emultipole(int max_am);
void libint2_cleanup_2emultipole(Libint_t* inteval);
void libint2_init_3emultipole(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_3emultipole(int max_am);
void libint2_cleanup_3emultipole(Libint_t* inteval);
void libint2_init_sphemultipole(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_sphemultipole(int max_am);
void libint2_cleanup_sphemultipole(Libint_t* inteval);
void libint2_init_overlap1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_overlap1(int max_am);
void libint2_cleanup_overlap1(Libint_t* inteval);
void libint2_init_kinetic1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_kinetic1(int max_am);
void libint2_cleanup_kinetic1(Libint_t* inteval);
void libint2_init_elecpot1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_elecpot1(int max_am);
void libint2_cleanup_elecpot1(Libint_t* inteval);
void libint2_init_1emultipole1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_1emultipole1(int max_am);
void libint2_cleanup_1emultipole1(Libint_t* inteval);
void libint2_init_2emultipole1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_2emultipole1(int max_am);
void libint2_cleanup_2emultipole1(Libint_t* inteval);
void libint2_init_3emultipole1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_3emultipole1(int max_am);
void libint2_cleanup_3emultipole1(Libint_t* inteval);
void libint2_init_sphemultipole1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_sphemultipole1(int max_am);
void libint2_cleanup_sphemultipole1(Libint_t* inteval);
void libint2_init_overlap2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_overlap2(int max_am);
void libint2_cleanup_overlap2(Libint_t* inteval);
void libint2_init_kinetic2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_kinetic2(int max_am);
void libint2_cleanup_kinetic2(Libint_t* inteval);
void libint2_init_elecpot2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_elecpot2(int max_am);
void libint2_cleanup_elecpot2(Libint_t* inteval);
void libint2_init_1emultipole2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_1emultipole2(int max_am);
void libint2_cleanup_1emultipole2(Libint_t* inteval);
void libint2_init_2emultipole2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_2emultipole2(int max_am);
void libint2_cleanup_2emultipole2(Libint_t* inteval);
void libint2_init_3emultipole2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_3emultipole2(int max_am);
void libint2_cleanup_3emultipole2(Libint_t* inteval);
void libint2_init_sphemultipole2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_sphemultipole2(int max_am);
void libint2_cleanup_sphemultipole2(Libint_t* inteval);
void libint2_init_eri(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_eri(int max_am);
void libint2_cleanup_eri(Libint_t* inteval);
void libint2_init_eri1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_eri1(int max_am);
void libint2_cleanup_eri1(Libint_t* inteval);
void libint2_init_eri2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_eri2(int max_am);
void libint2_cleanup_eri2(Libint_t* inteval);
void libint2_init_3eri(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_3eri(int max_am);
void libint2_cleanup_3eri(Libint_t* inteval);
void libint2_init_3eri1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_3eri1(int max_am);
void libint2_cleanup_3eri1(Libint_t* inteval);
void libint2_init_3eri2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_3eri2(int max_am);
void libint2_cleanup_3eri2(Libint_t* inteval);
void libint2_init_2eri(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_2eri(int max_am);
void libint2_cleanup_2eri(Libint_t* inteval);
void libint2_init_2eri1(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_2eri1(int max_am);
void libint2_cleanup_2eri1(Libint_t* inteval);
void libint2_init_2eri2(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_2eri2(int max_am);
void libint2_cleanup_2eri2(Libint_t* inteval);
void libint2_init_r12kg12(Libint_t* inteval, int max_am, void* buf);
size_t libint2_need_memory_r12kg12(int max_am);
void libint2_cleanup_r12kg12(Libint_t* inteval);
#ifdef __cplusplus
};
LIBINT_PRAGMA_CLANG(diagnostic pop)
LIBINT_PRAGMA_GCC(diagnostic pop)
#endif

/** Use LIBINT2_PREFIXED_NAME(fncname) to form properly prefixed function name from LIBINT2 API */
#define LIBINT2_PREFIXED_NAME(name) __libint2_prefixed_name__(LIBINT2_API_PREFIX,name)
#define __libint2_prefixed_name__(prefix,name) __prescanned_prefixed_name__(prefix,name)
#define __prescanned_prefixed_name__(prefix,name) prefix##name
/** Use LIBINT2_PREFIXED_NAME(fncname) to form properly prefixed function name from LIBINT2 API */
#define LIBINT2_DEFINED(taskname,symbol) __prescanned_libint2_defined__(taskname,symbol)
#define __prescanned_libint2_defined__(taskname,symbol) LIBINT2_DEFINED_##symbol

#endif

