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


#include "mgf_rank.h"

#include <stdlib.h>
#include "lib/mgf_population.h"
#include "mgf_dominance.h"

struct mgf_pop_ranks_t* mgf_new_pop_ranks(MoeazPop *pop){
    struct mgf_pop_ranks_t* new_ranks = calloc(1, sizeof(struct mgf_pop_ranks_t));
    new_ranks->pop = pop;
    return new_ranks;
};


void mgf_free_pop_ranks(struct mgf_pop_ranks_t* pop_ranks){
    free(pop_ranks);
}


struct mgf_indv_crowdist_a *mgf_front_to_crowdist_array(struct mgf_pop_front_t *front, MoeazPop *pop) {
    IndvidualType *indv_type = pop->indv_type;
    int size = front->size;
    struct mgf_indv_crowdist_a* new_crowd = calloc(size, sizeof(struct mgf_indv_crowdist_a));

    for (int i = 0; i < front->size; ++i) {
        new_crowd[i].index = front->idx[i];
        new_crowd[i].crowdist = indv_type->get_crowdist(mgf_pop_get_indv(pop, front->idx[i]));
    }
    return new_crowd;
}


void mgf_crowdist_array_to_front(struct mgf_indv_crowdist_a *cd, struct mgf_pop_front_t *front){
    int size = front->size;
    for (int i = 0; i < size; ++i) {
        front->idx[i] = cd[i].index;
    }
}


struct mgf_indv_crowdist_a *mgf_front_to_objidx_array(struct mgf_pop_front_t *front, MoeazPop *pop) {
    int size = front->size;
    struct mgf_indv_crowdist_a* new_crowd = calloc((size_t)size, sizeof(struct mgf_indv_crowdist_a));

    for (int i = 0; i < front->size; ++i) {
        new_crowd[i].index = front->idx[i];
        new_crowd[i].objs = mgf_indv_get_realdatapointer(mgf_pop_get_indv(pop, front->idx[i]));
    }
    return new_crowd;
}


static void mgf_crowdist_a_set_obj(struct mgf_indv_crowdist_a *cd, int size, int obj){
    for (int i = 0; i < size; ++i) {
        cd[i].obj = obj;
    }
}


/**
* Fast non.dominated sort. Following the description in [Deb et. al, 2002]
* @param pop:	The population to be ranked
*/
void pop_fast_non_dominated_sort(struct mgf_pop_ranks_t *pop_ranks)
{
    MoeazPop *pop = pop_ranks->pop;
    IndvidualType *indv_type = mgf_pop_get_indv(pop,0)->type;
    int qq, pp;
    int p, q, j;
    int dominance_flag;
    int rank;
    int **S, *nS;
    int *Q, nQ;
    int *n;

    int more_fronts_flag;

    /* Getting the ranks of the population (fast non-dominated sort) ******** */
    //printf("size: %d\n", pop->ranks);
    if (pop_ranks->front != NULL) {
        //vprint("\n\n Clean memory in fast nondominated sort\n");
//        assert(pop_ranks->front != NULL && pop_ranks->ranks > 0);
        for (j = 0; j < pop_ranks->ranks; ++j) {
            free(&pop_ranks->front[j]);
        }
        free(pop_ranks->front);
        pop_ranks->front = NULL;
        pop_ranks->ranks = 0;
    }

    rank = 1;
    S = (int**) malloc(sizeof(int*) * pop->size);
    nS = (int*) malloc(sizeof(int) * pop->size);
    n = (int*) malloc(sizeof(int) * pop->size);

    pop_ranks->ranks = rank;
    pop_ranks->front = realloc(pop_ranks->front, sizeof(struct mgf_pop_front_t) * pop_ranks->ranks);
    pop_ranks->front[0].size = 0;
    pop_ranks->front[0].idx = NULL;

    for (p = 0; p < pop->size; p++) {
        S[p] = NULL;
        nS[p] = 0;
        n[p] = 0;
        for (q = 0; q < pop->size; ++q)
        {
            if (q != p)
            {   // @TODO make dominance for bin a mixed
                dominance_flag = dominance(
                    mgf_indv_get_realdatapointer(mgf_pop_get_indv(pop,p)),
                    mgf_indv_get_realdatapointer(mgf_pop_get_indv(pop,q)),
                    indv_type->xsize);
                if (dominance_flag == 1) /* p dominates q */
                {
                    nS[p]++;
                    S[p] = (int*) realloc(S[p], sizeof(int) * nS[p]);
                    S[p][nS[p] - 1] = q;
                }
                else if (dominance_flag == -1) /* q dominates p */
                {
                    n[p]++;
                }
            }
        }
        if (n[p] == 0) {
            indv_type->set_rank(mgf_pop_get_indv(pop,p),rank);
            pop_ranks->front[rank - 1].size++;
            pop_ranks->front[rank - 1].idx = realloc(pop_ranks->front[rank - 1].idx,
                sizeof(int) * pop_ranks->front[rank - 1].size);
            pop_ranks->front[rank - 1].idx[pop_ranks->front[rank - 1].size - 1] = (unsigned int) p;

            /*
             pop->last_rank_size++;
             pop->last_rank = (int*) realloc(pop->last_rank,
             sizeof(int) * pop->last_rank_size);
             pop->last_rank[pop->last_rank_size - 1] = p;
             */
        }
    }
    pop_ranks->ranks = rank;

    //more_fronts_flag = (pop->last_rank_size == pop->size) ? 0 : 1;
    more_fronts_flag = (pop_ranks->front[rank - 1].size == pop->size) ? 0 : 1;
    while (more_fronts_flag == 1) {
        more_fronts_flag = 0;
        Q = NULL;
        nQ = 0;
        for (pp = 0; pp < pop_ranks->front[rank - 1].size; ++pp) {
            p = pop_ranks->front[rank - 1].idx[pp];
            for (qq = 0; qq < nS[p]; ++qq) {
                q = S[p][qq];
                n[q]--;
                if (n[q] == 0) {
                    more_fronts_flag = 1;
                    indv_type->set_rank(mgf_pop_get_indv(pop,q),rank + 1);
                    nQ++;
                    Q = (int*) realloc(Q, sizeof(int) * nQ);
                    Q[nQ - 1] = q;
                }
            }
        }
        if (more_fronts_flag == 1) {
            rank++;
            //free(pop->last_rank);
            pop_ranks->ranks = rank;
            pop_ranks->front = realloc(pop_ranks->front, sizeof(struct mgf_pop_front_t) * pop_ranks->ranks);
            pop_ranks->front[rank - 1].idx = (unsigned int*)Q;
            pop_ranks->front[rank - 1].size = nQ;

        }
    }
    free(n);
    free(nS);
    for (p = 0; p < pop->size; ++p) {
        if (S[p] != NULL) {
            free(S[p]);
        }
    }
    free(S);

    return;
}

int qsort_compare_obj(const void *a_, const void *b_)
{
    struct mgf_indv_crowdist_a *a, *b;
    a = (struct mgf_indv_crowdist_a*) a_;
    b = (struct mgf_indv_crowdist_a*) b_;

    int obj = a->obj;

    if (a->objs[obj] < b->objs[obj])
    {
        return -1;
    }
    else if (a->objs[obj] > b->objs[obj])
    {
        return 1;
    }
    return 0;
}

int qsort_compare_crowd(const void *a_, const void *b_)
{
    struct mgf_indv_crowdist_a *a, *b;
    a = (struct mgf_indv_crowdist_a*) a_;
    b = (struct mgf_indv_crowdist_a*) b_;

    if (a->crowdist < b->crowdist) {
        return -1;

    } else if (a->crowdist > b->crowdist) {
        return 1;
    }
    return 0;
}

void pop_crowding_assignment(struct mgf_pop_ranks_t *pop_ranks){
    MoeazPop *pop = pop_ranks->pop;
    IndvidualType *indv_type = pop->indv_type;
    struct mgf_pop_front_t *F = pop_ranks->front;

    double f_max, f_min;
    unsigned int id_max, id_min;
    unsigned int id_prev, id_next, id_cur;

//    crw_struct.pop = pop;

    //printf("F.size: %d\n", F->size);

    /* One element in the front*/
    if (F->size == 1) {
        id_cur = F->idx[0];
        indv_type->set_crowdist(mgf_pop_get_indv(pop,id_cur), 1.0e14);
        return;
    }
    /* Two element in the front*/
    if (F->size == 2) {
        id_cur = F->idx[0];
        indv_type->set_crowdist(mgf_pop_get_indv(pop,id_cur), 1.0e14);
        id_cur = F->idx[1];
        indv_type->set_crowdist(mgf_pop_get_indv(pop,id_cur), 1.0e14);
        return;
    }

    /* Cleaning crowding distance */
    for (int i = 0; i < F->size; ++i) {
        id_cur = F->idx[i];
        indv_type->set_crowdist(mgf_pop_get_indv(pop,id_cur), 0.0);
    }
    if (F->size > 2) {
        struct mgf_indv_crowdist_a *crowdists = mgf_front_to_objidx_array(F, pop);
        for (int obj = 0; obj < indv_type->xsize; obj++) {
            mgf_crowdist_a_set_obj(crowdists, F->size, obj);

            qsort(crowdists, F->size, sizeof(struct mgf_indv_crowdist_a), qsort_compare_obj);
//            qsort(F->idx, F->size, sizeof(F->idx[0]), _compare_obj);
//            obj = crw_struct.obj;

            mgf_crowdist_array_to_front(crowdists, F);
            id_min = F->idx[0];
            id_max = F->idx[F->size - 1];

            // @TODO convert to general, bin mixed int op.
            f_min = mgf_pop_get_indv(pop,id_min)->real[obj];
            f_max = mgf_pop_get_indv(pop,id_max)->real[obj];

            indv_type->set_crowdist(mgf_pop_get_indv(pop,id_min), 1.0e14);
            indv_type->set_crowdist(mgf_pop_get_indv(pop,id_max), 1.0e14);
            for (int i = 1; i < F->size - 1; ++i) {
                id_cur = F->idx[i];
                if (indv_type->get_crowdist(mgf_pop_get_indv(pop,id_cur)) != 1.0e14) {
                    id_prev = F->idx[i - 1];
                    id_next = F->idx[i + 1];
                    indv_type->set_crowdist(mgf_pop_get_indv(pop,id_cur),
                        indv_type->get_crowdist(mgf_pop_get_indv(pop,id_cur)) + (
                        (mgf_pop_get_indv(pop,id_next)->real[obj] -
                        mgf_pop_get_indv(pop,id_prev)->real[obj]) / (f_max - f_min))
                    );
                }
            }

        }
        /* Following lines are like the original implementation of NSGA-II */
        for (int i = 1; i < F->size - 1; ++i) {
            id_cur = F->idx[i];
            if (indv_type->get_crowdist(mgf_pop_get_indv(pop,id_cur)) != 1.0e14) {
                indv_type->set_crowdist(mgf_pop_get_indv(pop, id_cur),
                    indv_type->get_crowdist(mgf_pop_get_indv(pop,id_cur)) / indv_type->xsize
                );
            }
        }
    }
}

void crowding_sort(struct mgf_indv_crowdist_a *pop_crowdist, int size) {
    qsort(pop_crowdist, (size_t) size, sizeof(struct mgf_indv_crowdist_a), qsort_compare_crowd);
}
