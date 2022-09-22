/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_MOEAD_DE_H
#define MOGEN_MGN_MOEAD_DE_H

#include "mgn_types.h"
#include "mgn_moead_types.h"

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "decomposition/mgn_scalarization.h"


mgnMoa* mgn_moead_de_init(gsl_matrix *W,size_t nobj, size_t T
                       ,mgn_popl *epop
                       ,mgnMop *mop
                       ,void (*apply)(void*, void*), void* params
                       ,bool external);

//mgnMoa* mgn_moead_init(gsl_matrix *W,size_t nobj, size_t T
//                       ,mgn_popl *epop
//                       ,mgnMop *mop
//                       ,void (*apply)(void*, void*), void* params
//                       ,bool external);

void mgn_moead_free(mgnMoa *moead);

void mgn_moead_set_scalarization(mgnMoa* moead
                                 ,mgnf_decomp_scalar f);

mgn_pop* mgn_moead_getpop(mgnMoa* moead);

gsl_matrix* mgn_moead_get_w(mgnMoa* moead);

#endif //MOGEN_MGN_MOEAD_DE_H
