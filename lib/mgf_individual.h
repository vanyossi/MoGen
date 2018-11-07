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


#ifndef MOGEN_MGF_INDIVIDUAL_H
#define MOGEN_MGF_INDIVIDUAL_H

#include "global_types.h"

struct indv_t;
struct indv_type_t;

struct indv_type_t {
    int data_size;
    unsigned int mop_type;
    unsigned int xsize;
    unsigned int fsize;
    unsigned int gsize;
    // interface
    void (*typealloc)(Mop *mop, struct indv_t *);
    void (*free)(struct indv_t*);
    // traits: initialized null
    void (*get_decomp)(struct indv_t*);
    void (*get_crossindex)(struct indv_t*);
};

struct indv_t {
    struct indv_type_t *type;
    int xtype:16;
    int feasible:16;
    int        *type_idx;   //!< Type index identify which position are binary
    int        *bin;        //!< Int array, each bit represents one position
    double     *real;       //!< Double array for real values
    double      *f;
    double      *g;
    double CV;
    char buffer_start;
};


struct indv_type_t *mgf_indvtype_new(
    Moa *moa,
    int data_size,
    void (*typealloc)(Mop *mop, struct indv_t *),
    void (*free)(struct indv_t *));

struct indv_type_t* mgf_indvtype_std();

void mgf_indv_set_double(struct indv_t *indv, unsigned int pos, double value);

void mgf_indv_set_bin(struct indv_t *indv, unsigned int pos, int value);

void* mgf_indv_value_at(struct indv_t *indv, unsigned int pos);

int mgf_indv_value_isbin(struct indv_t *indv, unsigned int pos);

double* mgf_indv_get_realdatapointer(struct indv_t *indv);


struct indv_t* mgf_indv_new(struct indv_type_t *type);

void mgf_indv_init(struct indv_t *indv, Mop *mop);

struct indv_type_t* mgf_indv_type(struct indv_t* self);

void* mgf_indv_buffer(struct indv_t* self);

void mgf_indv_free_std(struct indv_t* indv);

void mgf_indv_free(struct indv_t* indv);


#endif //MOGEN_MGF_INDIVIDUAL_H