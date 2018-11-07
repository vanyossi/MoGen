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


#include "mogen_population.h"
#include <stdlib.h>
#include <memory.h>

#include "mogen_mop.h"
#include "rand.h"

void mogen_indv_alloc(Individual *indv, struct mop_t *mop){
//    Individual* indv = calloc(1,sizeof(Individual));
    indv->type = mop->set.type;
    indv->feasible = 1;
    memcpy(&indv->xsize, &mop->set.ndec, sizeof(int) * 3); // copy sizes

    if (mop->set.type == MOP_REAL) {
        indv->x.data = malloc(sizeof(indv_real_data) * mop->set.ndec);
    }
    else if (mop->set.type == MOP_BIN) {
        indv->x.data = malloc(sizeof(indv_bin_data) * mop->set.ndec);
    }
    else {
        Multiarray *mda = malloc(sizeof(Multiarray));
        mda->size = mop->set.ndec;
        mua_multiarray(mda, MOP_MIX);
        indv->x.mix = mda;
    }

    indv->f = calloc(sizeof(double),mop->set.nobj);
    indv->g = calloc(sizeof(double),mop->set.ncons);

};


void mogen_indv_init(Individual *indv, struct mop_t *mop){

    if (indv->type == MOP_REAL){
        indv_real_data *data = (indv_real_data*) indv->x.data;
        for (int i = 0; i < mop->set.ndec; ++i) {
            data[i] = rnd_real(mop->limits.xmin[i], mop->limits.xmax[i]);
        }
    } else if (indv->type == MOP_BIN) {
        indv_bin_data *data = (indv_bin_data*) indv->x.data;
        for (int i = 0; i < mop->set.ndec; ++i) {
            data[i] = (indv_bin_data)((rnd_perc() < 0.5) ? 0 : 1);
        }
    } else {
        double value;
        for (unsigned int i = 0; i < mop->set.ndec; ++i) {
            // bin or real
            if(rnd_perc() < 0.5) {
                value = rnd_real(mop->limits.xmin[i], mop->limits.xmax[i]);
                mua_set_double(indv->x.mix, i, value);
            } else {
                mua_set_int(indv->x.mix, i, (rnd_perc() < 0.5) ? 0 : 1);
            }
        }
    }

    if (mop->solver->extra_indv_alloc != NULL) {
        mogen_indv_init_extra(mop, indv, mop->solver->extra_indv_alloc);
    }
}

void mogen_indv_init_extra(Mop *mop, Individual *indv, indv_init_extra init_f) {
    init_f(mop, indv);
}

void mogen_indv_free(Individual *indv){
    if (indv->type != MOP_MIX) {
        free(indv->x.data);
    }
    else {
        mua_multiarray_del(indv->x.mix);
    }
};

#include <stdio.h>
MoeazPop* mogen_pop_alloc(Moa *moa, unsigned int size, size_t indv_size) {
    MoeazPop* new_pop = calloc(1,sizeof(MoeazPop));

    size_t alloc_size = sizeof(Individual);
    if (indv_size > alloc_size) {
        alloc_size = indv_size;
    }

    new_pop->size = size;
    new_pop->indv_size = alloc_size;
    new_pop->indv = calloc(size, alloc_size);

    char *indv_adress = (char*)new_pop->indv;
    for (int i = 0; i < size; ++i) {
        mogen_indv_alloc((Individual*)indv_adress, moa->mop);
        indv_adress += alloc_size;
    }
    moa->mop->pop = new_pop;

    for (int i = 1; i < 2; ++i) {
        printf("%p, %d", &moa->mop->pop->indv[i], moa->mop->pop->indv[i].type);
    }
    return new_pop;
}


#define mogen_pop_preploop(pop) \
    char *indv = (char*)(pop->indv); \
    size_t indvsize = pop->indv_size; \
    int size = pop->size;

void mogen_pop_init(Moa *moa){
    mogen_pop_preploop(moa->mop->pop);

    for (int i = 0; i < size; ++i) {
        mogen_indv_init( (Individual*)indv, moa->mop);
        indv += indvsize;
    }
}

void mogen_pop_init_extra(Moa *moa, indv_init_extra init_f){
    mogen_pop_preploop(moa->mop->pop);

    for (int i = 0; i < size; ++i) {
        mogen_indv_init_extra(moa->mop, (Individual *) indv, init_f);
        indv += indvsize;
    }
}

void mogen_pop_free(MoeazPop *pop){
    mogen_pop_preploop(pop);

    for (int i = 0; i < size; ++i) {
        mogen_indv_free((Individual*)indv);
        indv += indvsize;
    }
    free(pop->indv);
    free(pop);
}
