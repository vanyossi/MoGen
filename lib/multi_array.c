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

#include "multi_array.h"

#define MULTYPE_REAL 1
#define MULTYPE_BIN 2
#define MULTYPE_MIXED 3

#define MUA_BYTESIZE 8

static const unsigned int int_bytesize = (MUA_BYTESIZE * sizeof(int));

static int mua_true = 1;
static int mua_false = 0;

void mua_multiarray(Multiarray* ma, int type){
    ma->muatype = type;
    ma->type_idx = calloc(ma->size / int_bytesize,int_bytesize); // assign 1 bit per position
    if ((type & MULTYPE_REAL) == MULTYPE_REAL) {
        ma->real = malloc(sizeof (double) * ma->size);
    }
    if ((type & MULTYPE_BIN) == MULTYPE_BIN) {
        ma->bin = calloc(ma->size / int_bytesize, int_bytesize);
    }
};


void mua_set_double(Multiarray* ma, unsigned int pos, double value){
    if (pos < ma->size) {
        ma->real[pos] = value;
    }
};


void mua_set_int(Multiarray* ma, unsigned int pos, int value){
    int bitpos = pos / int_bytesize;
    int shift = pos % int_bytesize;

    ma->type_idx[bitpos] = ma->type_idx[bitpos] | (1 << shift);

    if (value) {
        ma->bin[bitpos] = ma->bin[bitpos] | (1 << shift);
    } else {
        ma->bin[bitpos] = ma->bin[bitpos] & ~(1 << shift);
    }
}


void* mua_value_at(Multiarray* ma, unsigned int pos){

    void* value = NULL;

    int bitpos = pos / int_bytesize;
    int shift = pos % int_bytesize;

    if (pos < ma->size) {
        if((ma->type_idx[bitpos] >> shift) & 1) {
            value = ((ma->bin[bitpos] >> shift & 1)? (void*) &mua_true: (void*) &mua_false);
        } else {
            value = (void*) &ma->real[pos];
        }
    }
    return value;
};


int mua_value_isbin(Multiarray* ma, unsigned int pos){
    int bitpos = pos / int_bytesize;
    int shift = pos % int_bytesize;
    return ma->type_idx[bitpos] >> shift & 1;
};


void mua_multiarray_del(Multiarray *ma){
    free(ma->type_idx);
    if (ma->muatype == MULTYPE_REAL) {
        free(ma->real);
    } else if (ma->muatype == MULTYPE_BIN) {
        free(ma->bin);
    } else {
        free(ma->real);
        free(ma->bin);
    }
};

