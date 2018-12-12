/*
 *  Copyright (c) 2018 Iván Yossi <ghevan@gmail.com>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "mgf_moead.h"

#include <float.h>
#include <stdio.h>

#include "mogen_mop.h"
#include "moea_math.h"
#include "mgf_operators.h"
#include "rand.h"
#include "mgf_pareto.h"

//#include <stdlib.h>

static void moa_moead_alloc(Moa *moa){
    int nobj = moa->mop->set.nobj;
    struct moa_moead_t *moead = mgf_moa_moead_data(moa);
    moead->cross_eta = 5.0; // big index for children close to parent
    moead->cross_prob = 1.0;
    moead->mut_eta = 5.0;
    moead->mut_prob = 1.0 / (moa->mop->set.xsize + moa->mop->set.bsize + moa->mop->set.isize);
    moead->st_id = 0;
    moead->st_p = NULL;
    moead->st_nadir = calloc(nobj, sizeof(double));
    /* Setting the utopian vector */
    moead->z_ref = calloc(nobj, sizeof(double));
    for (int i = 0; i < nobj; ++i) {
        moead->z_ref[i] = DBL_MAX;
    }
    moead->t_size = 20;   //!< Neibourghood size
    moead->st_indv = mgf_indv_new(mgf_indvtype_std(moa));
    moead->weights.transformation = 0;
    moead->weights.method = (nobj < 20) ? WM_WEIGHT_SLD: WM_WEIGHT_FILE;
    mgf_moa_set_scalarization(&moead->s_m, SCLM_TCH_NORM);
}

static void moa_moead_free(Moa *moa){
    struct moa_moead_t *moead = mgf_moa_moead_data(moa);
    if(moead->st_nadir) free(moead->st_nadir);
    if(moead->z_ref) free(moead->z_ref);
    mgf_moa_free_std(moa);
}

MoaType* mgf_moatype_moead(){
    struct moa_type_t *moat = mgf_moatype_new(sizeof(struct moa_moead_t),
        moa_moead_alloc,
        0,
        moa_moead_first_run,
        moa_moead_free);
    moat->get_crossover_vals = moead_get_crossvals;
    moat->get_mutation_vals = moead_get_mutvals;
    return moat;
}

struct moa_moead_t* mgf_moa_moead_data(Moa* moa){
    return mgf_moa_buffer(moa);
}

struct moa_cm_container moead_get_crossvals(Moa *self){
    struct moa_moead_t* data = mgf_moa_moead_data(self);
    return (struct moa_cm_container){data->cross_prob, data->cross_eta};
}

struct moa_cm_container moead_get_mutvals(Moa *self){
    struct moa_moead_t* data = mgf_moa_moead_data(self);
    return (struct moa_cm_container){data->mut_prob, data->mut_eta};
}


// move to sort_methods.h
static int compare_inc(const void **_a, const void **_b)
{
    const double *a, *b;
    a = (const double *) *_a;
    b = (const double *) *_b;
    return (a[1] < b[1]) ? -1 : (a[1] > b[1]) ? 1 : 0;
}

static void moead_neighborhood_set(int nobj, struct mgf_neighborhood_set *ns, struct mgf_weight_set *ws){
    unsigned int i, j;
    double **dist, **sort;

    dist = (double**) malloc(sizeof(double*) * ws->size);
    sort = (double**) malloc(sizeof(double*) * ws->size);

    for (i = 0; i < ws->size; i++)
    {
        dist[i] = (double*) malloc(sizeof(double) * ws->size);
        sort[i] = (double*) malloc(sizeof(double) * 2);
    }

    // O(kxN^2)
    for (int i = 0; i < ws->size; i++) {
        dist[i][i] = 0.0;

        for (int j = i + 1; j < ws->size; j++) {
            dist[i][j] = dist[j][i] = pDistance(ws->weight[i].w, ws->weight[j].w, nobj, 2.0);
        }
    }


    for (i = 0; i < ws->size; i++) {
        for (j = 0; j < ws->size; j++)
        {
            sort[j][0] = j;
            sort[j][1] = dist[i][j];
        }
        qsort(sort, (size_t) ws->size, sizeof(sort[0]), (int (*)(const void *, const void *)) &compare_inc);
        for (j = 0; j < ns->neighborhood[i].size; j++)
        {
            ns->neighborhood[i].members[j] = (int) sort[j][0];
        }
    }

    for (i = 0; i < ws->size; i++)
    {
        free(dist[i]);
        free(sort[i]);
    }

    free(dist);
    free(sort);
//    vprint(&moea_report, "\tget weight neighborhood...OK\n");

    return;
}

/**
 * @brief Update neighborhood and define the problem with best convergence
 * @param pop 		        The population
 * @param W 		        The weight set
 * @param nb 		        The neighborhood
 * @param y	    	        The new individual
 * @param nadir             Nadir point
 * @param dup               Duplication flag
 * @return                  1 if updated
 */
static unsigned int moead_neighborhood_update(struct moa_moead_t *moead_param,
    MoeazPop *pop, unsigned int nb_id, unsigned int *dup)
{
    struct mgf_neighborhood *nb = &(moead_param->neighbor_set.neighborhood[nb_id]);
    struct mgf_weight_set *w = moead_param->weights.wset;
    double *z_ref = moead_param->z_ref;
    double *st_nadir = moead_param->st_nadir;
    Individual *indv = moead_param->st_indv;

    int p_id;
    double g1, g2;
    unsigned int flag = 0;

    for (int i = 0; i < nb->size; i++) {
        p_id = nb->members[i];

        // scalar function
        g1 = moead_param->s_m.func(indv->type->fsize, indv->f, w->weight[p_id].w, z_ref, st_nadir);
        g2 = moead_param->s_m.func(indv->type->fsize, mgf_pop_get_indv(pop, p_id)->f, w->weight[p_id].w, z_ref, st_nadir);

        // On high unfeasible populations including the best one always
        // will cause a big loss in diversity and no parete will be found
        // in a resaonable time. Avoiding adding always the best allows
        // the crossover to mantain effectiveness.
        if (g1 < g2) {
            flag = 1;
        }

        UNUSED(dup);

        if (flag) {
            mgf_indv_copy(mgf_pop_get_indv(pop, p_id), indv);
            flag = 0;
        }
//        if (ind_realInPop(y, pop, mop.nreal) > -1) {
//            (*dup)++;
//        } else {
//            *dup = 0;
//        }
//
//        if(flag && *dup <= mop.nreal) {
//            cpy_ind(&(pop->ind[p_id]), y, &mop, NULL);
//            flag = 0;
//        }
    }
    return flag;
}

/**
 * @brief Recombination MOEA/D
 * @param pop 		Population
 * @param child 	The new individual
 * @param nb		The neighborhood
 * @param eta_c	    The crossover index
 * @param eta_m	    The mutation index
 * @param Pc		The crossover probability
 * @param Pm		The mutation probability
 */
static void mgf_moead_recombination(Mop *mop, MoeazPop *pop, unsigned int nb_id)
{
    Individual *tmp_child;
    Individual *child1, *child2;
    Individual *parent1, *parent2;

    struct moa_moead_t *moead_data = mgf_moa_moead_data(mop->solver);
    struct mgf_neighborhood *nb = &(moead_data->neighbor_set.neighborhood[nb_id]);

    Individual *child = moead_data->st_indv;

    int l, k;

    do {
        l = nb->members[rnd_int(0, (int)nb->size - 1)];
        k = nb->members[rnd_int(0, (int)nb->size - 1)];
    } while (l == k);

    parent1 = mgf_pop_get_indv(pop, l);
    parent2 = mgf_pop_get_indv(pop, k);

    tmp_child = mgf_indv_new(mgf_indvtype_std(mop->solver));

    child1 = child;
    child2 = tmp_child;

    if (rnd_perc() >= 0.5)
    {
        child1 = tmp_child;
        child2 = child;
    }

    mgf_operator()->crossover(mop,
        parent1, parent2,
        child1, child2);

    mgf_operator()->mutation(child, mop->solver->type->get_mutation_vals(mop->solver), mop->limits);

    mgf_indv_free(tmp_child);
//    free_ind(&tmp_child, NULL);
//    return;
}

Moa *moa_moead(Mop *mop, WeightFileOption foption){
    Moa *moa = mgf_moa_new(mop, "MOEAD", mgf_moatype_moead());
    WeightParams *weiParam = &mgf_moa_moead_data(moa)->weights;

    unsigned int psize = weights_settings(mop->set.nobj, weiParam, foption);

    mgf_moa_new_pop(moa, psize, mgf_indvtype_std(moa));
    mgf_pop_init(moa); // init population as it was created inside

    mgf_moa_moead_data(moa)->st_p = calloc(psize, sizeof(int));
    unsigned int *st_p = mgf_moa_moead_data(moa)->st_p;
    for (unsigned int i = 0; i < psize; ++i) {
        st_p[i] = i;
    }

//    alloc_weight_set(&st_W, moead.pop_size);
    weiParam->wset = calloc(1, sizeof(struct mgf_weight_set));
    weiParam->wset->size = psize;
    weiParam->wset->weight = calloc(weiParam->wset->size, sizeof(struct mgf_weight));
    for (int i = 0; i < weiParam->wset->size; i++) {
        weiParam->wset->weight[i].id = i;
        weiParam->wset->weight[i].w = calloc(mop->set.nobj, sizeof(double));
    }
//    alloc_neighborhood_set(&st_B, moead.pop_size, moead.T);
    struct mgf_neighborhood_set *nb = &mgf_moa_moead_data(moa)->neighbor_set;
    nb->size = psize;
    nb->neighborhood = calloc(nb->size, sizeof(struct mgf_neighborhood));
    for (int i = 0; i < nb->size; i++) {
        nb->neighborhood[i].size = mgf_moa_moead_data(moa)->t_size;
        nb->neighborhood[i].members = calloc(mgf_moa_moead_data(moa)->t_size, sizeof(int));
    }

    /* Initialize weight vectors ******************************************** */
    weights_initialize_set(weiParam, mop->set.nobj);
    moead_neighborhood_set(mop->set.nobj, nb, weiParam->wset);

    return moa;
}

// this is more like an initializer function.
mbool moa_moead_first_run(Mop *mop){
    struct moa_moead_t *moead_data = mgf_moa_moead_data(mop->solver);

    mgf_pop_evaluate(mop->pop, mop);

    for (unsigned int i = 0; i < mop->pop->size; ++i) {
        // only for elite
//        update_ind_list_continuous(&st_elite, &st_pop.ind[i], NULL, NULL, NULL);
        mgf_pareto_update_z( mgf_pop_get_indv(mop->pop,i), moead_data->z_ref, mop->set.nobj);
    }

    mop->solver->type->run = moa_moead_run;
    return mtrue;
}

mbool moa_moead_run(Mop *mop){
    struct moa_moead_t *moead_data = mgf_moa_moead_data(mop->solver);

    rnd_shuffle_vector(moead_data->st_p, mop->pop->size, sizeof(int));
    mgf_pop_get_nadir(mop->pop, moead_data->st_nadir, mop->set.nobj);

    // evaluate pop
    MoeazPop *pop = mop->pop;
    Individual *indv;

    for (int i = 0; i < pop->size; ++i){
        moead_data->st_id = moead_data->st_p[i];
        mgf_moead_recombination(mop, pop, moead_data->st_id);

//        indv = mgf_pop_get_indv(pop, i);
        indv = moead_data->st_indv;
        if (!mop_evaluate(mop, indv)) return mfalse; // stop after N evals

        mgf_pareto_update_z(indv, moead_data->z_ref, mop->set.nobj);
        moead_neighborhood_update(moead_data, pop, moead_data->st_id, NULL);
//        MOEAD_neighborhood_update(&st_pop, &st_W, &(st_B.neighborhood[st_id]), &st_child, st_nadir, &dup);
    }
    return mtrue;
}
