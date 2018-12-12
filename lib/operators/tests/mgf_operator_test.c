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

#include <stdio.h>

#include "rand.h"

#include "mogen_mop.h"
#include "mgf_nsga2.h"

void print_test_pop(MoeazPop*);
void print_front(struct mgf_pop_front_t *front);

void mgf_crowdist_a_set_obj(struct mgf_indv_crowdist_a *cd, int size, int obj);
int qsort_compare_obj(const void *a_, const void *b_);
int qsort_compare_crowd(const void *a_, const void *b_);

int main(int argc, char const *argv[]) {
    set_random(0.141516);

    Mop *mop = mogen_mop("Op_test", MOP_REAL, 0);
    mop_set_params(mop, 5, 0, 0, 2, 0);

    double xmin = 0.0;
    double xmax = 1.0;

    mop_set_limits_ndec(mop, &xmin, &xmax, 1, NULL, NULL, 0);

    // initialize pop for individual type.
    Moa *moa = moa_nsga2(mop);
    mgf_moa_new_pop(moa, 5, mgf_indvtype_nsga2(moa));
    mgf_pop_init(moa);

    Individual *indv = NULL;
    for (int i = 0; i < mop->pop->size; ++i) {
        indv = mgf_pop_get_indv(mop->pop, i);
        indv_nsga2_set_crowdist(indv, 0); // crowdist initial is 0
        for (int j = 0; j < 5; ++j) {
            mgf_indv_set_double(indv, j, 5 - i);
        }
        for (int j = 0; j < mop->set.xsize; ++j) {
            indv->f[j] = (2 - j) * (5 - i);
        }
    }

    MoeazPop *obj_sorted_pop = mgf_pop_alloc(mop->pop->size, mgf_indvtype_nsga2(moa));
    for (int i = 0; i < mop->pop->size; ++i) {
        mgf_indv_copy(
            mgf_pop_get_indv(obj_sorted_pop, i),
            mgf_pop_get_indv(mop->pop, i)
        );
    }

    print_test_pop(obj_sorted_pop);

    // testing ranks
    struct mgf_pop_ranks_t* prank = mgf_new_pop_ranks(mop->pop);

    prank->ranks = 1;
    prank->front = realloc(prank->front, sizeof(struct mgf_pop_front_t) * prank->ranks);
    struct mgf_pop_front_t *F = prank->front;

    F->size = 5;
    F->idx = calloc(5, sizeof(int));
    for (unsigned int j = 0; j < 5; ++j) {
        F->idx[j] = j;
    }

    print_front(F);

    struct mgf_indv_crowdist_a *cdist = mgf_front_to_objidx_array(F, obj_sorted_pop);
    mgf_crowdist_a_set_obj(cdist, F->size, 1);
    qsort(cdist, F->size, sizeof(struct mgf_indv_crowdist_a), qsort_compare_obj);
    mgf_crowdist_array_to_front(cdist, F);

    print_front(F);

    unsigned int id_min = F->idx[0];
    unsigned int id_max = F->idx[F->size - 1];
    // @TODO convert to general, bin mixed int op.
    double f_min = mgf_pop_get_indv(obj_sorted_pop,id_min)->f[1];
    double f_max = mgf_pop_get_indv(obj_sorted_pop,id_max)->f[1];

    printf("id_min %d, min: %.5f\n", id_min, f_min);
    printf("id_max %d, max: %.5f", id_max, f_max);

    mop->pop->indv_type->set_crowdist(mgf_pop_get_indv(obj_sorted_pop,id_min), 1.0e14);
    mop->pop->indv_type->set_crowdist(mgf_pop_get_indv(obj_sorted_pop,id_max), 1.0e14);

    int id_cur, id_next, id_prev;
    for (int i = 1; i < F->size - 1; ++i) {
        id_cur = F->idx[i];
        if (mop->pop->indv_type->get_crowdist(mgf_pop_get_indv(obj_sorted_pop,id_cur)) != 1.0e14) {
            id_prev = F->idx[i - 1];
            id_next = F->idx[i + 1];
            mop->pop->indv_type->set_crowdist(mgf_pop_get_indv(obj_sorted_pop,id_cur),
                mop->pop->indv_type->get_crowdist(mgf_pop_get_indv(obj_sorted_pop,id_cur)) + (
                    (mgf_pop_get_indv(obj_sorted_pop,id_next)->f[1] -
                     mgf_pop_get_indv(obj_sorted_pop,id_prev)->f[1]) / (f_max - f_min))
            );
        }
    }

    print_test_pop(obj_sorted_pop);




//    MoeazPop *crowd_sorted_pop = mgf_pop_alloc(mop->pop->size, mgf_indvtype_nsga2(moa));
////    for (int i = 0; i < mop->pop->size; ++i) {
////        mgf_indv_copy(
////            mgf_pop_get_indv(crowd_sorted_pop, i),
////            mgf_pop_get_indv(mop->pop, i)
////        );
////    }

    print_test_pop(obj_sorted_pop);

    // testing ranks
    prank->ranks = 1;
    prank->front = realloc(prank->front, sizeof(struct mgf_pop_front_t) * prank->ranks);
    F = prank->front;

    F->size = 5;
    F->idx = calloc(5, sizeof(int));
    for (unsigned int j = 0; j < 5; ++j) {
        F->idx[j] = j;
    }

    print_front(F);

    struct mgf_indv_crowdist_a *crow_dist = mgf_front_to_crowdist_array(F, obj_sorted_pop);
    qsort(crow_dist, F->size, sizeof(struct mgf_indv_crowdist_a), qsort_compare_crowd);
    mgf_crowdist_array_to_front(crow_dist, F);

    print_front(F);

    return 0;
}

void print_test_pop(MoeazPop *pop){
    Individual *indv = NULL;
    for (int i = 0; i < pop->size; ++i) {
        indv = mgf_pop_get_indv(pop, i);
        for (int j = 0; j < 5; ++j) {
            printf("%.6f, ", mgf_indv_get_double(indv,j));
        }
        puts("");
        for (int j = 0; j < pop->indv_type->fsize; ++j) {
            printf("%.6f ", indv->f[j]);
        }
        printf("\ncrowd: %.6f", indv_nsga2_crowdist(indv));
        puts("");
    }
}

void print_front(struct mgf_pop_front_t *front){
    for (int i = 0; i < front->size; ++i) {
        printf("%d", front->idx[i]);
    }
    puts("");
}
