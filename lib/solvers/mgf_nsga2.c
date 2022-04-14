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


#include "mgf_nsga2.h"

#include <assert.h>

#include "rand.h"

#include "mogen_mop.h"
#include "mgf_operators.h"
#include "mgf_rank.h"

void indv_nsga_alloc(Mop *mop, struct indv_t *indv) {
    UNUSED(mop);
    mgf_indv_nsga2_data(indv)->rank = 0;
    mgf_indv_nsga2_data(indv)->crow_dist = 0.0;
}

void mgf_indv_nsga2_copy(struct indv_t *to, struct indv_t *from){
    indv_nsga2_set_crowdist(to, indv_nsga2_crowdist(from));
    indv_nsga2_set_rank(to, indv_nsga2_rank(from));
}

struct indv_type_t* mgf_indvtype_nsga2(Moa *moa){
    struct indv_type_t* new_indv = mgf_indvtype_new(
        moa, sizeof(struct indv_t_nsga2_type), indv_nsga_alloc, mgf_indv_nsga2_copy, mgf_indv_free_std);
    new_indv->get_rank = indv_nsga2_rank;
    new_indv->set_rank = indv_nsga2_set_rank;
    new_indv->get_crowdist = indv_nsga2_crowdist;
    new_indv->set_crowdist = indv_nsga2_set_crowdist;
    return new_indv;
}

struct indv_t_nsga2_type* mgf_indv_nsga2_data(struct indv_t *self){
    return mgf_indv_buffer(self);
}

int indv_nsga2_rank(struct indv_t *indv){
    return mgf_indv_nsga2_data(indv)->rank;
}

void indv_nsga2_set_rank(struct indv_t *indv, int rank){
    mgf_indv_nsga2_data(indv)->rank = rank;
}

double indv_nsga2_crowdist(struct indv_t *indv){
    return mgf_indv_nsga2_data(indv)->crow_dist;
}

void indv_nsga2_set_crowdist(struct indv_t *indv, double crowdist){
    mgf_indv_nsga2_data(indv)->crow_dist = crowdist;
}

static void moa_nsga2_alloc(Moa *moa){
    struct moa_nsga2_t *nsga2 = mgf_moa_nsga2_data(moa);
    nsga2->cross_eta = 20;
    nsga2->cross_prob = 0.9;
    nsga2->mut_eta = 20;
    nsga2->mut_prob = 1.0 / moa->mop->set.xsize;
}

static void moa_nsga2_free(Moa *moa){
    struct moa_nsga2_t *nsga2 = mgf_moa_nsga2_data(moa);
    if(nsga2->Q) mgf_pop_free(nsga2->Q);
    if(nsga2->mix) mgf_pop_free(nsga2->mix);
    if(nsga2->pop_ranks) free(nsga2->pop_ranks);
    mgf_moa_free_std(moa);
}

MoaType* mgf_moatype_nsga2(){
    struct moa_type_t *moat = mgf_moatype_new(sizeof(struct moa_nsga2_t),
        moa_nsga2_alloc,
        0,
        moa_nsga2_first_run,
        moa_nsga2_free);
    moat->get_crossover_vals = ngsa2_get_crossvals;
    moat->get_mutation_vals = ngsa2_get_mutvals;
    return moat;
}

struct moa_nsga2_t* mgf_moa_nsga2_data(Moa* moa){
    return mgf_moa_buffer(moa);
}

void moa_nsga2_setparams(struct moa_nsga2_t *nsga2, double cx_eta, double cxp, double m_eta, double m_p) {
    moa_nsga2_cross_eta(nsga2,cx_eta);
    moa_nsga2_cross_prob(nsga2,cxp);
    moa_nsga2_mut_eta(nsga2,m_eta);
    moa_nsga2_mut_prob(nsga2,m_p);
}

struct moa_cm_container ngsa2_get_crossvals(Moa *self){
    struct moa_nsga2_t* data = mgf_moa_nsga2_data(self);
    return (struct moa_cm_container){data->cross_prob, data->cross_eta};
}

struct moa_cm_container ngsa2_get_mutvals(Moa *self){
    struct moa_nsga2_t* data = mgf_moa_nsga2_data(self);
    return (struct moa_cm_container){data->mut_prob, data->mut_eta};
}

double moa_nsga2_cross_eta(struct moa_nsga2_t *nsga2, double value) {
    if(value > 0) nsga2->cross_eta = value;
    return nsga2->cross_eta;
}

double moa_nsga2_cross_prob(struct moa_nsga2_t *nsga2, double value) {
    if(value > 0) nsga2->cross_prob = value;
    return nsga2->cross_prob;
}

double moa_nsga2_mut_eta(struct moa_nsga2_t *nsga2, double value) {
    if(value > 0) nsga2->mut_eta = value;
    return nsga2->mut_eta;
}

double moa_nsga2_mut_prob(struct moa_nsga2_t *nsga2, double value) {
    if(value > 0) nsga2->mut_prob = value;
    return nsga2->mut_prob;
}

/** nsga2 methods
 */

Moa *moa_nsga2(Mop *mop){
    Moa *moa = mgf_moa_new(mop, "NSGA2", mgf_moatype_nsga2());
//    moa_nsga2_setparams(mgf_moa_nsga2_data(moa), 20, 0.5, 30, 0.5);
    return moa;
}
//void mgf_nsga2_settings(char *moea_name, int max_gen, double seed){}
//
//void NSGA2_report_headers(MOEA_REPORT *report){}

int mgf_nsga2_dominance(Individual *a, Individual *b) {

    int nobjs = a->type->fsize;
    int flag1 = 0;
    int flag2 = 0;

    /* Handling constraint mechanism */
    if (a->CV < 0 && b->CV < 0) {
        if (a->CV > b->CV) {
            return (1);
        } else {
            if (a->CV < b->CV) {
                return (-1);
            } else {
                return (0);
            }
        }
    } else {
        if (a->CV < 0 && b->CV == 0) {
            return (-1);
        } else {
            if (a->CV == 0 && b->CV < 0) {
                return (1);
            } else {
                for (int i = 0; i < nobjs; i++) {
                    if (a->f[i] < b->f[i]) {
                        flag1 = 1;
                    } else {
                        if (a->f[i] > b->f[i]) {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1 == 1 && flag2 == 0) {
                    return (1);
                } else {
                    if (flag1 == 0 && flag2 == 1) {
                        return (-1);
                    } else {
                        return (0);
                    }
                }
            }
        }
    }
}

unsigned int mgf_nsga2_binary_tournament(MoeazPop *P, unsigned int p1, unsigned int p2){
    Individual *ind1 = mgf_pop_get_indv(P, p1);
    Individual *ind2 = mgf_pop_get_indv(P, p2);
    int flag = mgf_nsga2_dominance(ind1, ind2);

    if (flag == 1) {
        return p1;

    } else if (flag == -1) {
        return p2;

    } else if (mgf_indv_nsga2_data(ind1)->crow_dist > mgf_indv_nsga2_data(ind2)->crow_dist) {
        return p1;

    } else if (mgf_indv_nsga2_data(ind2)->crow_dist > mgf_indv_nsga2_data(ind1)->crow_dist) {
        return p2;

    } else {
        return (rnd_perc() < 0.5) ? p1 : p2;
    }
}

static void set_integer_sequence(int *dest, int from, int to) {
    int seq_size = to - from;
    int j = from;
    for (int i = 0; i < seq_size; i++, j++) {
        dest[i] = j;
    }
}


void mgf_nsga2_parent_selection(MoeazPop *P, unsigned int *parent, unsigned int n_parent)
{
    unsigned int *p;
    unsigned int p1, p2;
    int j = 0;
    int N = (P->size % 2 == 0) ? P->size : P->size - 1;

    UNUSED(n_parent);
    // @TODO check if only even pop is still necessary
    // only: even population
    assert(n_parent == (int )(P->size / 2));

    p = calloc(P->size, sizeof(int));
    set_integer_sequence((int*)p, 0, P->size);
    rnd_shuffle_vector(p, P->size, sizeof(int));

    for (unsigned int i = 0; i < N; i += 2)
    {
        p1 = p[i];
        p2 = p[i + 1];

        parent[j] = mgf_nsga2_binary_tournament(P, p1, p2);
        j++;
    }
    free(p);
}

// #include <stdio.h>
int mgf_nsga2_offspring(Mop *mop, struct moa_nsga2_t *moa_data) {
    Individual *p1 = 0, *p2 = 0, tmp_ind;
    MoeazPop *P = mop->pop;
    MoeazPop *Q = moa_data->Q;

    MutationSettings mutset = mop->solver->type->get_mutation_vals(mop->solver);
//     for (int k = 0; k < P->size; ++k) {
//         if ( mgf_pop_get_indv(P,k)->type == NULL) {
//             printf("type %d, %d, %p\n", P->size, k, mgf_pop_get_indv(P,k)->type);
//         }
//     }
// //    mutset.prob = moa_data->mut_prob;
// //    mutset.eta = moa_data->mut_eta;

    int j = 0;
    unsigned int *parent1, *parent2;
    unsigned int n_parent;
    unsigned int N;

    /* Even population */
    n_parent = P->size / 2;
    parent1 = malloc(sizeof(int) * n_parent);
    parent2 = malloc(sizeof(int) * n_parent);
    N = ((P->size % 2) == 0) ? (P->size) : (P->size - 1);

    mgf_nsga2_parent_selection(P, parent1, n_parent);
    mgf_nsga2_parent_selection(P, parent2, n_parent);

    for (unsigned int i = 0; i < N; i += 2, j++)
    {
        p1 = mgf_pop_get_indv(P, parent1[j]);
        p2 = mgf_pop_get_indv(P, parent2[j]);

        mgf_operator()->crossover(mop,
            p1, p2,
            mgf_pop_get_indv(Q, i),
            mgf_pop_get_indv(Q, i+1));

        mgf_operator()->mutation(mgf_pop_get_indv(Q,i), mutset, mop->limits);
        mgf_operator()->mutation(mgf_pop_get_indv(Q,i +1), mutset, mop->limits);

        mop_evaluate(mop, mgf_pop_get_indv(Q, i));
        mop_evaluate(mop, mgf_pop_get_indv(Q, i+1));

//        ELITE_update_ind_list_continuous(&elite, &Q->ind[i], NULL, NULL, NULL);
//        ELITE_update_ind_list_continuous(&elite, &Q->ind[i + 1], NULL, NULL, NULL);
    }
    /* Odd population */
    if (N != Q->size)
    {
        // @TODO add error message asserts
        assert(N == (Q->size - 1));
        tmp_ind = *mgf_indv_new(P->indv_type);
        p1 = mgf_pop_get_indv(P, (unsigned int)rnd_int(0, P->size - 1)); // no tournament
        p2 = mgf_pop_get_indv(P, (unsigned int)rnd_int(0, P->size - 1)); // no tournament

        mgf_operator()->crossover(mop,
            p1, p2,
            mgf_pop_get_indv(Q, N),
            &tmp_ind);

        mgf_operator()->mutation(mgf_pop_get_indv(Q,N), mutset, mop->limits);
        mop_evaluate(mop, mgf_pop_get_indv(Q, N));

//        ELITE_update_ind_list_continuous(&elite, &Q->ind[N], NULL, NULL, NULL);
        mgf_indv_free(&tmp_ind);
    }
    free(parent1);
    free(parent2);

    return 0;
}

void mgf_nsga2_offspring_generation(Mop *mop, struct moa_nsga2_t *moa)
{
    // change offsprint generation depending on type
    mgf_nsga2_offspring(mop, moa);
}



void mgf_nsga2_reduce(MoeazPop *P, Moa *moa)
{
    MoeazPop *mixed_pop = mgf_moa_nsga2_data(moa)->mix;
    struct mgf_pop_ranks_t *mix_ranks = mgf_new_pop_ranks(mixed_pop);

    int rank;
    int j, selection_count;
    int delete;
    int sol, id_cur;

    pop_fast_non_dominated_sort(mix_ranks);

    selection_count = 0;
    for (rank = 0; rank < mix_ranks->ranks; rank++)
    {
        pop_crowding_assignment(mix_ranks, rank);
        selection_count += mix_ranks->front[rank].size;

        if (selection_count >= (P->size)) {
            break;
        }
    }
    delete = selection_count - P->size;
    assert(delete >= 0);
    UNUSED(delete);

    sol = 0;
    /* Copying the first ranks */
    for (int i = 0; i < rank; i++)
    {
        for (j = 0; j < mix_ranks->front[i].size; j++)
        {
            assert(sol <= P->size);
            id_cur = mix_ranks->front[i].idx[j];
            mgf_indv_copy(mgf_pop_get_indv(P,sol), mgf_pop_get_indv(mixed_pop, id_cur));
            sol++;
        }
    }

    struct mgf_indv_crowdist_a *crowdists = mgf_front_to_crowdist_array(&(mix_ranks->front[rank]), mixed_pop);
    crowding_sort(crowdists, mix_ranks->front[rank].size);
    mgf_crowdist_array_to_front(crowdists, &(mix_ranks->front[rank]));

    /* Copying the best solutions according to crowding distance */
    for (j = 0; j < mix_ranks->front[rank].size && sol < P->size; j++)
    {
        id_cur = mix_ranks->front[rank].idx[j];
        mgf_indv_copy(mgf_pop_get_indv(P,sol), mgf_pop_get_indv(mixed_pop, id_cur));
        sol++;
    }
    assert(sol == P->size);


    free(crowdists);
    mgf_free_pop_ranks(mix_ranks);
}

// run step
mbool moa_nsga2_first_run(Mop *mop) {
    MoeazPop *ori_pop = mop->pop;

    // initialize Offsrping and Mix pops
    struct moa_nsga2_t *nsga2 = mgf_moa_nsga2_data(mop->solver);
    nsga2->Q = mgf_pop_alloc(ori_pop->size, ori_pop->indv_type);
    nsga2->mix = mgf_pop_alloc(ori_pop->size * 2, ori_pop->indv_type);
    nsga2->pop_ranks = mgf_new_pop_ranks(ori_pop);

//    struct mgf_pop_ranks_t *p_rank = mgf_moa_nsga2_data(mop->solver)->pop_ranks;

    mgf_pop_evaluate(mop->pop, mop);
    // we need ranks and Front
    // rank is integer number
    // Front is int array.
    pop_fast_non_dominated_sort(nsga2->pop_ranks);
    for (int i = 0; i < nsga2->pop_ranks->ranks; ++i) {
        pop_crowding_assignment(nsga2->pop_ranks, i);
    }

    mop->solver->type->run = moa_nsga2_run;
    return mtrue;
}

mbool moa_nsga2_run(Mop *mop){

    // We assume moa is of type nsga2
    struct moa_nsga2_t *nsga2_data =  mgf_moa_nsga2_data(mop->solver);

    mgf_nsga2_offspring_generation(mop, nsga2_data);
    mgf_pop_merge(nsga2_data->mix, mop->pop, nsga2_data->Q);
    mgf_nsga2_reduce(mop->pop, mop->solver);

    return mtrue;
}
