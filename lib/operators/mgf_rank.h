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


#ifndef MOGEN_RANK_H
#define MOGEN_RANK_H

#include "mgf_global_types.h"

struct mgf_pop_front_t {
    int size;
    unsigned int *idx;
};

struct mgf_pop_ranks_t {
    MoeazPop *pop;
    int ranks;
    struct mgf_pop_front_t *front;

};

struct mgf_indv_crowdist_a {
    unsigned int index;
    double *objs; //!> NULL by default
    double crowdist;
    int obj;
};

struct mgf_pop_ranks_t* mgf_new_pop_ranks(MoeazPop *pop);

void mgf_free_pop_ranks(struct mgf_pop_ranks_t* pop_ranks);

struct mgf_indv_crowdist_a *mgf_front_to_crowdist_array(struct mgf_pop_front_t *front, MoeazPop *pop);

void mgf_crowdist_array_to_front(struct mgf_indv_crowdist_a *cd, struct mgf_pop_front_t *front);

struct mgf_indv_crowdist_a *mgf_front_to_objidx_array(struct mgf_pop_front_t *front, MoeazPop *pop);

/**
 * Fast non.dominated sort as in [Deb et. al, 2002]
 * @param pop:	The population to be ranked
 */
void pop_fast_non_dominated_sort(struct mgf_pop_ranks_t *pop_ranks);

void pop_crowding_assignment(struct mgf_pop_ranks_t *pop_ranks, int front_idx);

void crowding_sort(struct mgf_indv_crowdist_a *pop_crowdist, int size);

void nondominated_sorting_and_crowding_assignment(struct mgf_pop_ranks_t *pop_ranks);


#endif //MOGEN_RANK_H
