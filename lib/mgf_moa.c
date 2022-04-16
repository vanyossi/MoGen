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


#include "mgf_moa.h"

#include <string.h>

#include "mogen_mop.h"

struct moa_type_t *mgf_moatype_new(
    int data_size,
    void (*typealloc)(struct moa_t *),
    mbool (*stop)(struct moa_t* moa, MoaStopCriterion criterion),
    mbool (*run)(struct mop_t* mop),
    void (*free)(struct moa_t *))
{
    struct moa_type_t *moat = calloc(1, sizeof(struct moa_type_t));

    moat->data_size = data_size;
    moat->typealloc = (typealloc)? typealloc : 0;
    moat->stop = (stop)? stop : moa_stop;
    moat->run = (run)? run : moa_run;
    moat->free = (free)? free : mgf_moa_free;

    return moat;
}

struct moa_type_t* mgf_moa_std(){
    return mgf_moatype_new(0, 0, moa_stop, moa_run, mgf_moa_free_std);
}

struct moa_t *mgf_moa_new(Mop *mop, char *name, struct moa_type_t *type) {
    int size = sizeof(struct moa_t) + type->data_size;
    struct moa_t* moa = calloc(1,size);

    moa->type = type;

    strcpy(moa->name, name);
    moa->mop = mop;
    mop->solver = moa;

    if (moa->type->typealloc) {
        moa->type->typealloc(moa);
    }

    return moa;
}

//void mgf_moa_init(struct moa_t *moa, Mop *mop);

struct moa_type_t* mgf_moa_type(struct moa_t* self){
    return self->type;
}

mbool mgf_moa_run(Moa* moa){
    return moa->type->run(moa->mop);
}

void* mgf_moa_buffer(struct moa_t* self){
    return (void*) &(self->buffer_start);
}

scalarization_f mgf_moa_scalarization(struct moa_t *self){
    return mgf_moa_type(self)->scalarize;
}

void mgf_moa_free(struct moa_t* moa){
    if (moa) {
        moa->type->free(moa);
    }
}

void mgf_moa_free_std(struct moa_t* moa){
    UNUSED(moa);
    return;
}


mbool moa_run(struct mop_t *mop){
    for (size_t i = 0; i < mop->pop->size; ++i){
        Individual* indv = mgf_pop_get_indv(mop->pop,i);
        if (!mop_evaluate(mop, indv)){
            return mfalse;
        }
    }
    return mtrue;
}

void moa_stopat_gen(Moa *moa, unsigned int gen){
    moa->stops.gens = gen;
}

void moa_stopat_eval(Moa *moa, unsigned int eval){
    moa->stops.evals = eval;
}

mbool moa_stop(Moa *moa, MoaStopCriterion criterion){
    mbool stop = mfalse;
    switch (criterion) {
        case MGN_STOPIF_EXEC:
            if (moa->stops.evals > 0)
                stop = (moa->mop->report.total.evals >= moa->stops.evals) ? mtrue : mfalse;
            break;
        case MGN_STOPIF_GEN:
            if (moa->stops.gens > 0)
                stop = (moa->mop->report.total.gens >= moa->stops.gens) ? mtrue : mfalse;
            break;
        default:
            break;
    }
    return stop;
}
