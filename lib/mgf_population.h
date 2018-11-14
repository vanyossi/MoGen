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


#ifndef MOGEN_MGF_POPULATION_H
#define MOGEN_MGF_POPULATION_H

#include <stdlib.h>

#include "global_types.h"
#include "mgf_individual.h"

struct mgf_pop_t {
    unsigned int size;
    unsigned int current;
    size_t indv_size;
    IndvidualType *indv_type;
    Individual* indv;
};

MoeazPop *mgf_pop_alloc(unsigned int size, IndvidualType *indvtype);

void mgf_moa_new_pop(Moa *moa, unsigned int size, IndvidualType *indvtype);

void mgf_pop_init(Moa *moa);

void mgf_pop_free(MoeazPop *pop);

mbool mgf_pop_evaluate(MoeazPop *pop, Mop *mop);

Individual* mgf_pop_get_indv(MoeazPop *pop, unsigned int pos);

Individual* mgf_pop_current_indv(MoeazPop *pop);

// void until NULL
void * mgf_pop_next(MoeazPop *pop);

void mgf_pop_reset_cursor(MoeazPop *pop);

void mgf_pop_merge(MoeazPop *from1, MoeazPop *from2, MoeazPop *to);


#endif //MOGEN_MGF_POPULATION_H
