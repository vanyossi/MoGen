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

#include "multi_array.h"

struct mop_t;

typedef union multi_data_t {
    double* real;
    unsigned short* bin;
    Multiarray mix;
} MultiData;

typedef struct moeaz_indv_t {
    int type;
    unsigned int xsize;
    unsigned int fsize;
    unsigned int gsize;
    MultiData x; // real, bin, or mixed
    double* f;
    double* g;
} MoeazIndv;

void moeaz_indv_init(MoeazIndv *indv, struct mop_t *mop);

void moeaz_indv_free(MoeazIndv *indv);


typedef struct moeaz_pop_t {
    int size;
    MoeazIndv* indv;
    unsigned int *front;
} MoeazPop;

MoeazPop *moeaz_pop_init(struct mop_t *mop, unsigned int size);

void moeaz_pop_free(MoeazPop* pop);

#endif //MOGEN_POP_H
