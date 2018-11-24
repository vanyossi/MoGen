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


#ifndef MOGEN_POP_H
#define MOGEN_POP_H

#include "mgf_global_types.h"
#include "multi_array.h"

#define indv_real_data double
#define indv_bin_data unsigned short
#define get_data_real(var) (indv_real_data*)(var->x.data)
#define get_data_bin(var) (indv_bin_data*)(var.x.data)

struct mop_t;

typedef void(*indv_init_extra)(Mop*, Individual*);

union multi_data_t {
    void *data;             //!> real and binary
    Multiarray *mix;  // que sea continuo. alguna cod para cambiar
};

struct moeaz_indv_t {
    int type:16;
    int feasible:16;
    double CV;
    unsigned int xsize;
    unsigned int fsize;
    unsigned int gsize;
    MultiData x; // real, bin, or mixed
    double* f;
    double* g;
};

void mogen_indv_alloc(Individual *indv, struct mop_t *mop);

void mogen_indv_init(Individual *indv, struct mop_t *mop);

void mogen_indv_init_extra(Mop *mop, Individual *indv, indv_init_extra init_f);

void mogen_indv_free(Individual *indv);


struct moeaz_pop_t {
    int size;
    size_t indv_size;
    Individual* indv;
};

MoeazPop *mogen_pop_alloc(Moa *moa, unsigned int size, size_t indv_size);

void mogen_pop_init(Moa *moa);

void mogen_pop_init_extra(Moa *moa, indv_init_extra init_f);

void mogen_pop_free(MoeazPop *pop);

#endif //MOGEN_POP_H
