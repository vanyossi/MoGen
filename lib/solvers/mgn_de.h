//
// Created by Iván Yossi on 10/05/22.
//

#ifndef MOGEN_MGN_DE_H
#define MOGEN_MGN_DE_H

#include <string.h>
#include <stdbool.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_randist.h>

#include "population.h"
#include "mgn_mop.h"

#include "mgn_random.h"
#include "mgn_moa.h"
#include "mgn_selection.h"
#include "gsl_vector_additional.h"
#include "mgn_sphere.h"

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

// tru refactor to make de_param private
struct mgn_de_param {
    size_t size;
    size_t D;
    double F;
    bool popeval;
    bool mopset;
    mgn_pop* pop;
    mgn_pop* trial;
    mgnMop *mop;
};

typedef struct mgnp_de_cross_param {
    double factor;
    gsl_vector *r0;
    gsl_vector *r1;
    gsl_vector *r2;
    gsl_vector *tcur;
} mgnp_crossparam;

void mgn_de_run(mgnMoa* de);
void mgnp_de_mut(mgnp_crossparam *param, size_t pos);
void mgn_de_crossover(mgn_pop *pop, size_t *rsel, size_t jrand, de_param *dparam);
void mgn_de_selection(mgn_pop *pop, mgn_pop *trial, mgnMop *mop);


/*
 * Initialization need
 *  bounds for x_
 *  used for initialization, random uniform
 *
 *
 */
de_param* mgn_de_alloc(size_t Np, void* iops, void* iparams)
{
    de_param* params = malloc(sizeof(*params));
    params->size = Np;
    params->pop = mgn_pop_alloc(Np,iops, iparams);
    params->popeval = false;
    params->trial = mgn_pop_alloc(Np,iops,iparams);
    params->D = params->pop->ops->get_iparams(params->pop->I).x->size;
    params->F = 1.15;

    params->mopset = false;

    return params;
}

mgnMoa* mgn_de_init(de_param* dparam
    ,void (*apply)(void*, void*), void* params)
{
    mgnMoa *de = malloc(sizeof(*de));
    de->tot_exec = 0;
    de->run = mgn_de_run;
    strncpy(de->name, "Differential Evolution", MOA_NAME_LEN);
    de->features = dparam;

    mgn_pop_init(dparam->pop, apply, params);

    return de;
}

de_param* mgn_de_getfeatures(mgnMoa* moa)
{
    return (de_param*)moa->features;
}

void mgn_de_setmop(mgnMoa *de, mgnMop *mop)
{
    de_param *prm = mgn_de_getfeatures(de);
    prm->mop = mop;
    prm->mopset = true;
}

void mgn_de_eval(mgnMoa *de)
{
    de_param *prm = mgn_de_getfeatures(de);
    mgn_mop_eval_pop(prm->mop, prm->pop, prm->mop->params);
    prm->popeval = true;
}

void mgn_de_free(mgnMoa *de)
{
    de_param *prm = mgn_de_getfeatures(de);
    mgn_pop_free(prm->pop);
    mgn_pop_free(prm->trial);
    free(prm);
    free(de);
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

void mgn_de_run(mgnMoa* de)
{
    de_param* m_p = mgn_de_getfeatures(de);
    if(m_p->mopset == false || m_p->popeval == false) return;

    // generate trial pop
    gsl_vector_ulong *rindex = gsl_vector_ulong_alloc(m_p->size - 1);

    for (size_t i = 0; i < m_p->size; ++i) {
        mgnp_de_rndseq_except(rindex,i);
        gsl_vector_ulong_view rsel = gsl_vector_ulong_subvector(rindex,0,4);
        gsl_vector_ulong_set(&rsel.vector,3,i);

        unsigned long jrand = rnd_getUniform_int((int)m_p->D);
        mgn_de_crossover(m_p->pop, rsel.vector.data, jrand, m_p);
        mgn_mop_eval_pop_index(m_p->mop, m_p->pop, m_p->mop->params,i,1);
    }
//    mgn_de_selection(m_p->pop, m_p->trial, m_p->mop);
    mgn_mop_eval_pop(m_p->mop, m_p->trial, m_p->mop->params);

    // A sel
//    mgn_sel_sub(m_p->pop, m_p->trial);

// B sel
//    pop_sort_1d(m_p->pop);
//    pop_sort_1d(m_p->trial);
//    mgn_sel_lambda_mu(m_p->pop, m_p->trial,5);

// C sel
    mgn_sel_lambda_mu_plus(m_p->pop, m_p->trial, pop_sort_1d);

//    for (size_t i = 0; i < m_p->trial->size; ++i) {
//        gsl_vector *f = pop_get_iparam(m_p->trial,i).f;
//        printf("%.6f ",f->data[0]
//        );
//    }
//    puts("");


    gsl_vector_ulong_free(rindex);
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
    double *tcur = param->tcur->data;

    tcur[pos] = r0[pos] + param->factor * (r1[pos] - r2[pos]);
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
//    printf("cur %zu, %zu, %zu, %zu\n", rsel[0], rsel[1], rsel[2], rsel[3]);
    gsl_vector *r0 = pop->ops->get_iparams(mgn_pop_get(pop,rsel[0])).x;
    gsl_vector *r1 = pop->ops->get_iparams(mgn_pop_get(pop,rsel[1])).x;
    gsl_vector *r2 = pop->ops->get_iparams(mgn_pop_get(pop,rsel[2])).x;
    gsl_vector *cur = pop->ops->get_iparams(mgn_pop_get(pop,rsel[3])).x;

    gsl_vector *tcur = dparam->trial->ops->get_iparams(mgn_pop_get(dparam->trial,rsel[3])).x;
    double Cr = 0.9;
    struct mgnp_de_cross_param fparam = {dparam->F,r0,r1,r2,tcur};

//    gsl_vector_map(tcur, mgnp_de_mut, &fparam);
    for (size_t j = 0; j < tcur->size; ++j) {
        if (rnd_getUniform() <= Cr || j == jrand) {
            mgnp_de_mut(&fparam, j);
        } else {
            gsl_vector_set(tcur,j,gsl_vector_get(cur,j));
        }
    }
}

/*
 * Selection
 * if u <= x, replace from target, compare with parent
 *      -- u    if f(u) <= f(x)
 * x = -|
 *      -- x    otherwise
 */
void mgn_de_selection(mgn_pop *pop, mgn_pop *trial, mgnMop *mop)
{
    mgn_mop_eval_pop(mop, trial, mop->params);
    for (size_t i = 0; i < pop->size; ++i) {
        gsl_vector *x = pop_get_iparam(pop,i).f;
        gsl_vector *v = pop_get_iparam(trial,i).f;
        if (mgn_mop_sphere_min(v->data,x->data,v->size) == -1){
            mgn_pop_copy(pop,trial,i,i,1);
        }
    }
}

#endif //MOGEN_MGN_DE_H
