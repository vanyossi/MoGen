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

void moeaz_indv_alloc(MoeazIndv *indv, struct mop_t *mop){
//    MoeazIndv* indv = calloc(1,sizeof(MoeazIndv));
    indv->type = mop->set.type;
    memcpy(&indv->xsize, &mop->set.ndec, sizeof(int) * 3);

    if (mop->set.type == MOP_REAL) {
        indv->x.real = malloc(sizeof(double) * mop->set.ndec);
    }
    else if (mop->set.type == MOP_BIN) {
        indv->x.bin = malloc(sizeof(int) * mop->set.ndec);
    }
    else {
        indv->x.mix = (Multiarray){mop->set.ndec};
        mua_multiarray(&indv->x.mix, MOP_MIX);
    }

    indv->f = malloc(sizeof(double) * mop->set.nobj);
    indv->g = malloc(sizeof(double) * mop->set.ncons);

};


void moeaz_indv_init(MoeazIndv *indv, struct mop_t *mop){

    if (indv->type == MOP_REAL){
        for (int i = 0; i < mop->set.ndec; ++i) {
            indv->x.real[i] = rnd_real(mop->limits.xmin[i], mop->limits.xmax[i]);
        }
    } else if (indv->type == MOP_BIN) {
        for (int i = 0; i < mop->set.ndec; ++i) {
            indv->x.bin[i] = (short)((rnd_perc() < 0.5) ? 0 : 1);
        }
    } else {
        double value;
        for (int i = 0; i < mop->set.ndec; ++i) {
            // bin or real
            if(rnd_perc() < 0.5) {
                value = rnd_real(mop->limits.xmin[i], mop->limits.xmax[i]);
                mua_set_double(&indv->x.mix, i, value);
            } else {
                mua_set_int(&indv->x.mix, i, (rnd_perc() < 0.5) ? 0 : 1);
            }
        }
    }
}


void moeaz_indv_free(MoeazIndv *indv){
    if (indv->type == MOP_REAL) {
      free(indv->x.real);
    }
    else if (indv->type == MOP_BIN) {
        free(indv->x.bin);
    }
    else {
        mua_multiarray_del(&indv->x.mix);
    }
};


MoeazPop* moeaz_pop_alloc(struct mop_t *mop, unsigned int size) {
    MoeazPop* new_pop = calloc(1,sizeof(MoeazPop));

    new_pop->size = size;
    new_pop->indv = calloc(size, sizeof(MoeazIndv));
    for (int i = 0; i < size; ++i) {
        moeaz_indv_alloc(&new_pop->indv[i], mop);
    }
    mop->pop = new_pop;
    return new_pop;
}


void moeaz_pop_init(struct mop_t *mop){
    for (int i = 0; i < mop->pop->size; ++i) {
        moeaz_indv_init(&mop->pop->indv[i], mop);
    }
}

void moeaz_pop_free(MoeazPop* pop){

    for (int i = 0; i < pop->size; ++i) {
        moeaz_indv_free(&pop->indv[i]);
    }
    free(pop->indv);
    free(pop);
}
