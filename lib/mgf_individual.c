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


#include "mgf_individual.h"

#include <stdlib.h>
#include <memory.h>
#include <limits.h>

#include "mgf_moa.h"
#include "mogen_mop.h"
#include "rand.h"

static const unsigned int INT_BITSIZE = WORD_BIT;
static int INDV_TRUE = 1;
static int INDV_FALSE = 0;

//#define indv_standard mgf_indvtype_new(0, mgf_indv_free_std);

struct indv_type_t *mgf_indvtype_new(
    Moa *moa, // @TODO remove need for moa
    int data_size,
    void (*typealloc)(Mop *mop, struct indv_t *),
    void (*free)(struct indv_t *))
{
    struct indv_type_t *indv = calloc(1, sizeof(struct indv_type_t));
    memcpy(&indv->xsize, &moa->mop->set.ndec, sizeof(int) * 3);
    indv->mop_type = moa->mop->set.type;
    indv->data_size = data_size;
    indv->typealloc = typealloc;
    indv->free = (!free)? free : mgf_indv_free_std;

    return indv;
};

struct indv_type_t* mgf_indvtype_std(Moa *moa){
    return mgf_indvtype_new(moa, 0, 0, mgf_indv_free_std);
}

// alloc
struct indv_t* mgf_indv_new(struct indv_type_t *type){
    int size = sizeof(struct indv_t) + type->data_size;
    struct indv_t* indv = calloc(1,size);
    indv->type = type;
    indv->xtype = indv->type->mop_type;
    // assign 1 bit for each element in variables
    indv->type_idx = calloc(1, indv->type->xsize / INT_BITSIZE);
    if ((indv->type->mop_type & MOP_REAL) == MOP_REAL) {
        indv->real = malloc(sizeof (double) * indv->type->xsize);
    }
    if ((indv->type->mop_type & MOP_BIN) == MOP_BIN) {
        indv->bin = calloc(1,indv->type->xsize / INT_BITSIZE);
    }

    indv->f = calloc(sizeof(double), type->fsize);
    indv->g = calloc(sizeof(double), type->gsize);

    return indv;
};


void mgf_indv_init(struct indv_t *indv, Mop *mop){

    if (indv->xtype == MOP_REAL){
        double *data = mgf_indv_get_realdatapointer(indv);
        for (int i = 0; i < mop->set.ndec; ++i) {
            data[i] = rnd_real(mop->limits.xmin[i], mop->limits.xmax[i]);
        }
    } else if (indv->xtype == MOP_BIN) {
        for (int i = 0; i < mop->set.ndec; ++i) {
            mgf_indv_set_bin(indv, i, ((rnd_perc() < 0.5) ? 0 : 1));
        }
    } else {
        double value;
        for (unsigned int i = 0; i < mop->set.ndec; ++i) {
            // bin or real
            if(rnd_perc() < 0.5) {
                value = rnd_real(mop->limits.xmin[i], mop->limits.xmax[i]);
                mgf_indv_set_double(indv, i, value);
            } else {
                mgf_indv_set_bin(indv, i, (rnd_perc() < 0.5) ? 0 : 1);
            }
        }
    }

    if (indv->type->typealloc != NULL) {
        indv->type->typealloc(mop, indv);
    }
}

// get extended data pointer
struct indv_type_t* mgf_indv_type(struct indv_t* self){
    return self->type;
}

void* mgf_indv_buffer(struct indv_t* self){
    return (void*) &(self->buffer_start);
}

void mgf_indv_free(struct indv_t* indv){
    if (indv) {
        indv->type->free(indv);
    }
}

/** MixedData API **/
#define DOUBLE_NAN __builtin_nan("0xFFF0000000000000")

union mixed {
    double val;
    int data[sizeof(double)/sizeof(int)];
    unsigned long long i;
};

double mgf_indv_get_double(struct indv_t *indv, unsigned int pos){
    if (pos >= indv->type->xsize) {
        union mixed ret = {DOUBLE_NAN};
        ret.data[0] = 1;
        return ret.val;
    }
    return indv->real[pos];
}

int mgf_indv_get_bin(struct indv_t *indv, unsigned int pos){
    int bitpos = pos / INT_BITSIZE;
    int shift = pos % INT_BITSIZE;

    return (indv->bin[bitpos] >> shift & 1)? INDV_TRUE : INDV_FALSE;
}

void mgf_indv_set_double(struct indv_t *indv, unsigned int pos, double value){
    if (pos < indv->type->xsize) {
        indv->real[pos] = value;
    }
}

void mgf_indv_set_bin(struct indv_t *indv, unsigned int pos, int value){
    int bitpos = pos / INT_BITSIZE;
    int shift = pos % INT_BITSIZE;

    indv->type_idx[bitpos] = indv->type_idx[bitpos] | (1 << shift);

    if (value) {
        indv->bin[bitpos] = indv->bin[bitpos] | (1 << shift);
    } else {
        indv->bin[bitpos] = indv->bin[bitpos] & ~(1 << shift);
    }
}



double mgf_indv_value_at(struct indv_t *indv, unsigned int pos){
    union mixed value = {DOUBLE_NAN};
    value.data[0] = 1; // regular NaN this data is 0

    int bitpos = pos / INT_BITSIZE;
    int shift = pos % INT_BITSIZE;

    if (pos < indv->type->xsize) {
        if((indv->type_idx[bitpos] >> shift) & 1) {
            value.val = ((indv->bin[bitpos] >> shift & 1)? INDV_TRUE: INDV_FALSE);
        } else {
            value.val = indv->real[pos];
        }
    }
    return value.val;
}

int mgf_indv_value_isbin(struct indv_t *indv, unsigned int pos){
    int bitpos = pos / INT_BITSIZE;
    int shift = pos % INT_BITSIZE;
    return indv->type_idx[bitpos] >> shift & 1;
}

double* mgf_indv_get_realdatapointer(struct indv_t *indv){
    return indv->real;
}

/** Standard individual functions **/
void mgf_indv_free_std(struct indv_t* indv){
    if (indv->real) free(indv->real);
    if (indv->bin)  free(indv->bin);
    if (indv->type_idx)  free(indv->type_idx);
}

