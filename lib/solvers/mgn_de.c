/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_de.h"

#include <string.h>
#include <stdbool.h>

//#include <gsl/gsl_permutation.h>
//#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_randist.h>

#include "population.h"
#include "individual.h"
#include "mgn_mop.h"

#include "mgn_random.h"
#include "mgn_moa.h"
#include "mgn_selection.h"
//#include "gsl_vector_additional.h"
//#include "mgn_sphere.h"

/*
 * Need
 *  P_x Current population : Np
 *  P_v intermediary population (mutants) : Np
 *
 *  each vector of P_x is recombined with P_v to produce
 *  P_u <- trial population : Np
 *
 *  During recomb, trial overwrite mutant, so P_v can hold both
 */

typedef struct mgn_de_param de_param;

// TODO add setter and getter for F, Cr, etc
struct mgn_de_param {
    size_t size;
    size_t D;
    double F;
    double cr;
    bool popeval;
    bool mopset;
    mgn_pop* pop;
    mgn_pop* trial;
    mgnMop *mop;
    mgn_de_ef(ef);
    mgn_de_ef_param* ef_param;
};

typedef struct mgnp_de_cross_param {
    double factor;
    gsl_vector *r0;
    gsl_vector *r1;
    gsl_vector *r2;
    gsl_vector *tu;
} mgnp_crossparam;

//void mgn_de_run(mgnMoa* de);
//void mgnp_de_mut(mgnp_crossparam *param, size_t pos);
//void mgn_de_crossover(mgn_pop *pop, size_t *rsel, size_t jrand, de_param *dparam);
//void mgn_de_selection(mgn_pop *pop, mgn_pop *trial, mgn_de_ef(ef));


de_param* mgn_de_getfeatures(mgnMoa* moa)
{
    return (de_param*)moa->features;
}

mgn_pop_proto* mgn_de_pop_get(mgnMoa* moa)
{
    return (mgn_pop_proto*)mgn_de_getfeatures(moa)->pop;
}

/*
 * Selection
 * if u <= x, replace from target, compare with parent
 *      -- u    if f(u) <= f(x)
 * x = -|
 *      -- x    otherwise
 */
void mgn_de_selection(mgn_pop *pop, mgn_pop *trial, mgn_de_ef(ef), mgn_de_ef_param* ef_p)
{
//    mgn_mop_eval_pop(mop, trial, mop->params);
    for (size_t i = 0; i < pop->size; ++i) {
        gsl_vector *x = pop_get_iparam(pop,i).f;
        gsl_vector *v = pop_get_iparam(trial,i).f;
        ef_p->pos = i;
        if (ef(x,v,ef_p) == 1){
//            printf("true x > v\n");
            mgn_pop_copy(pop,trial,i,i,1);
        }
    }
}

void mgnp_de_rndseq_except(gsl_vector_ulong *seq, size_t except)
{
    // generate seq
    size_t s_i = 0;
    for (size_t i = 0; i <= seq->size; ++i) {
        if(i == except) continue; // skip except element
        gsl_vector_ulong_set(seq,s_i,i);
        s_i++;
    }
    gsl_ran_shuffle(rnd_get_generator(), seq->data, seq->size, sizeof(seq->data[0]));
}

/*
 * mutation
 * P_X -> P_v + P_x == P_u
 * Differential mutation adds a scaled, randomly sampled,
 * vector difference to third vector
 *
 * v_i = x_r0 + F(x_r1 -xr2) -> F positive factor > 1
 *
 */
void mgnp_de_mut(mgnp_crossparam *param, size_t pos)
{
    double *r0 = param->r0->data;
    double *r1 = param->r1->data;
    double *r2 = param->r2->data;
    double *cur = param->tu->data;

    cur[pos] = r0[pos] + (param->factor * (r1[pos] - r2[pos]));
}

/*
 * Crossover
 * Uniform crossover (discrete recombination)
 * Build _u from 2 vectors
 *
 * u = -- v     if (rand <= Cr or j = j_rand)
 *     |
 *     -- x     otherwise
 *
 *     j == trial parameter choosen rand from mutant
 */
void mgn_de_crossover(mgn_pop *pop, size_t *rsel, size_t jrand, de_param *dparam)
{
//    printf("tu %zu, %zu, %zu, %zu\n", rsel[0], rsel[1], rsel[2], rsel[3]);
    gsl_vector *r0 = pop->ops->get_iparams(mgn_pop_get(pop,rsel[0])).x;
    gsl_vector *r1 = pop->ops->get_iparams(mgn_pop_get(pop,rsel[1])).x;
    gsl_vector *r2 = pop->ops->get_iparams(mgn_pop_get(pop,rsel[2])).x;
    gsl_vector *xi = pop->ops->get_iparams(mgn_pop_get(pop,rsel[3])).x;
    gsl_vector *tu = dparam->trial->ops->get_iparams(mgn_pop_get(dparam->trial,rsel[3])).x;

    struct mgnp_de_cross_param fparam = {dparam->F, r0, r1, r2, tu};

//    gsl_vector_map(tu, mgnp_de_mut, &fparam);
    for (size_t j = 0; j < tu->size; ++j) {
        if (rnd_getUniform() <= dparam->cr || j == jrand) {
            mgnp_de_mut(&fparam, j);
        } else {
            gsl_vector_set(tu,j,gsl_vector_get(xi,j));
        }
    }
}

void mgn_de_run(mgnMoa* de)
{
    de_param* m_p = mgn_de_getfeatures(de);
    if(m_p->mopset == false || m_p->popeval == false) return;

    // generate trial pop
    gsl_vector_ulong *rindex = gsl_vector_ulong_alloc(m_p->size - 1);

    for (size_t i = 0; i < m_p->size; ++i) {
        // select 4 random position and set last to current one.
        mgnp_de_rndseq_except(rindex,i);
        gsl_vector_ulong_view rsel = gsl_vector_ulong_subvector(rindex,0,4);
        gsl_vector_ulong_set(&rsel.vector,3,i);

        // seleccionamos param x_i,j aleatorio
        unsigned long jrand = rnd_getUniform_int((int)m_p->D);
        mgn_de_crossover(m_p->pop, rsel.vector.data, jrand, m_p);
        mgn_mop_eval_pop_index(m_p->mop, m_p->pop, m_p->mop->params,i,1);
        de->tot_exec++;
    }

    // original Selection
    mgn_mop_eval_pop(m_p->mop, m_p->trial, m_p->mop->params);
    mgn_de_selection(m_p->pop,m_p->trial,m_p->ef, m_p->ef_param);
//    mgn_pop_copy(m_p->trial, m_p->pop, 0, 0, m_p->trial->size);

    // A sel
//    mgn_sel_sub(m_p->pop, m_p->trial);

// B sel
//    mgn_pop_prank_sort(m_p->pop);
//    mgn_pop_prank_sort(m_p->trial);
//    mgn_sel_lambda_mu(m_p->pop, m_p->trial,5);
//    pop_sort_1d(m_p->pop);
//    pop_sort_1d(m_p->trial);

// C sel
    mgn_sel_lambda_mu_plus(m_p->pop, m_p->trial, mgn_pop_prank_sort);

//    for (size_t i = 0; i < m_p->trial->size; ++i) {
//        gsl_vector *f = pop_get_iparam(m_p->trial,i).f;
//        printf("%.6f ",f->data[0]
//        );
//    }
//    puts("");


    gsl_vector_ulong_free(rindex);
}


// ============== public functions ==
/*
 * Initialization need
 *  bounds for x_
 *  used for initialization, random uniform
 *
 *
 */
mgnMoa* mgn_moa_de_alloc(size_t Np
                         ,void* iops
                         ,void* iparams
                         ,double factor
                         ,double cr
)
{
    mgnMoa *de = malloc(sizeof(*de));
    de_param* params = malloc(sizeof(*params));

    params->size = Np;
    params->pop = mgn_pop_alloc(Np,iops, iparams);
    params->popeval = false;
    params->trial = mgn_pop_alloc(Np,iops,iparams);
    params->D = params->pop->ops->get_iparams(params->pop->I).x->size;
    params->F = factor;
    params->cr = cr;
    params->mopset = false;

    strncpy(de->name, "Differential Evolution", MOA_NAME_LEN);
    de->tot_exec = 0;
    de->run = mgn_de_run;
    de->pop_get = mgn_de_pop_get;
    de->features = params;

    return de;
}
//
//de_param* mgn_de_alloc(size_t Np, void* iops, void* iparams)
//{
//    de_param* params = malloc(sizeof(*params));
//    params->size = Np;
//    params->pop = mgn_pop_alloc(Np,iops, iparams);
//    params->popeval = false;
//    params->trial = mgn_pop_alloc(Np,iops,iparams);
//    params->D = params->pop->ops->get_iparams(params->pop->I).x->size;
//    params->F = 1.10;
//
//    params->mopset = false;
//
//    return params;
//}
//
void mgn_de_init(mgnMoa* moa
                 ,void (*apply)(void*, void*)
                 ,void* apply_params
)
{
    de_param *de_p = mgn_de_getfeatures(moa);
    mgn_pop_init(de_p->pop, apply, apply_params);
}

// TODO move to moa;
//      how to allow inheritance (no this)
void mgn_de_setmop(mgnMoa *de, mgnMop *mop, mgn_de_ef(ef), mgn_de_ef_param* ef_p)
{
    de_param *prm = mgn_de_getfeatures(de);
    prm->mop = mop;
    prm->ef = ef;
    prm->ef_param = ef_p;
    prm->mopset = true;
}

void mgn_de_eval(mgnMoa *de)
{
    de_param *prm = mgn_de_getfeatures(de);
    mgn_mop_eval_pop(prm->mop, prm->pop, prm->mop->params);
    de->tot_exec += prm->pop->size;
    prm->popeval = true;
}

void mgn_moa_de_free(mgnMoa *de)
{
    de_param *prm = mgn_de_getfeatures(de);
    mgn_pop_free(prm->pop);
    mgn_pop_free(prm->trial);
    free(prm);
    free(de);
}
