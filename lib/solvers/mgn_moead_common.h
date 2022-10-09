/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_MOEAD_COMMON_H
#define MOGEN_MGN_MOEAD_COMMON_H

#include "mgn_moead_types.h"

#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "mgn_random.h"
#include "mgn_moa.h"
#include "mgn_mop.h"
#include "population.h"
#include "mgn_poplist.h"
#include "mgn_pareto.h"
#include "mgn_initializer.h"

#include "decomposition/mgn_scalarization.h"
#include "operators/mgn_vector_distance.h"
#include "operators/mgn_gen_operator.h"


struct moead_features {
    bool isalloc;
    bool ismopset;
    bool isprobset;
    mgn_ga_sets *ga_set;
    size_t size_nei;
    double *z;
    mgn_pop *pop;
    mgn_popl *epop;
    mgnMop *mop;
    mgnMop *_mop; // private
    gsl_matrix *wei;
    gsl_matrix *dist;
    gsl_matrix_int *dindex;
    mgnf_decomp_scalar scalarize;
//    double (*scalarize)(gsl_vector*, gsl_vector*, gsl_vector*);
};


moeadf* mgn_moead_getfeatures(mgnMoa* moead);


void moead_run(mgnMoa* moead);

// instead of checking dominance on every loop
// do this after every run, else just add to EP
// use sort by dominance
void moead_update_ep(moeadf *set, mgn_pop *lpop);


bool moead_stop(mgnMoa *moa);


double* mgnp_moead_alloc_z(size_t size);


// xfg param
void mgnp_moead_update_z(gsl_vector* x, gsl_vector* f, gsl_vector* g, moeadf* param);

// rpop = external reference pop
moeadf* mgn_moead_alloc_features(gsl_matrix *W
                                 , size_t nobj
                                 , size_t T
                                 , mgn_popl *rpop
                                 , mgnMop *mop
                                 , bool external);


void mgn_moead_free(mgnMoa *moead);


void moead_set_prob(mgnMoa*moead, mgn_ga_sets *gasets);


void moead_pop_evaluate(mgn_pop *pop, moeadf* set);


mgnMoa* mgn_moead_common_init(gsl_matrix *W,size_t nobj, size_t T
                       ,mgn_popl *epop
                       ,mgnMop *mop
                       ,void (*apply)(void*, void*), void* params
                       ,bool external);


void mgn_moead_pop_init_eval(mgnMoa *moa, mgn_initializer *init);


void mgn_moead_set_scalarization(mgnMoa* moead
                                 ,mgnf_decomp_scalar f);


gsl_matrix* mgn_moead_get_w(mgnMoa* moead);


mgn_pop* mgn_moead_getpop(mgnMoa* moead);


void mgn_moead_set_mop(mgnMoa* moead, mgnMop *mop, int type);

#endif //MOGEN_MGN_MOEAD_COMMON_H
