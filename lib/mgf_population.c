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


#include "mgf_population.h"

#include <memory.h>

#include "mgf_moa.h"
#include "mogen_mop.h"

#include <stdio.h> // @TODO remove

MoeazPop *mgf_pop_alloc(unsigned int size, IndvidualType *indvtype) {
    MoeazPop* new_pop = calloc(1,sizeof(MoeazPop));

    size_t alloc_size = sizeof(Individual);
    if (indvtype->data_size > 0) {
        alloc_size += indvtype->data_size;
    }

    new_pop->size = size;
    new_pop->indv_size = alloc_size;
    new_pop->indv_type = indvtype;
    new_pop->indv = calloc(size, alloc_size);

    char *indv_adress = (char*)new_pop->indv;
    for (int i = 0; i < size; ++i) {
        // @TODO assigning array value missing
        memcpy(indv_adress, mgf_indv_new(indvtype), alloc_size);
        indv_adress += alloc_size;
    }
    return new_pop;
}

void mgf_moa_new_pop(Moa *moa, unsigned int size, IndvidualType *indvtype){
    MoeazPop *new_pop = mgf_pop_alloc(size, indvtype);

    if(moa->mop->pop) {
        mgf_pop_free(moa->mop->pop);
    }

    moa->mop->pop = new_pop;

    Individual* new_indv;
    for (int i = 1; i < 3; ++i) {
        new_indv = mgf_pop_get_indv(new_pop, i);
        printf("%p, %d\n", new_indv, new_indv->xtype);
    }
}


#define mgf_pop_preploop(pop) \
    char *indv = (char*)(pop->indv); \
    size_t indvsize = pop->indv_size; \
    int size = pop->size;


void mgf_pop_init(Moa *moa){
    mgf_pop_preploop(moa->mop->pop);

    for (int i = 0; i < size; ++i) {
        mgf_indv_init( (Individual*)indv, moa->mop);
        indv += indvsize;
    }
}

void mgf_pop_free(MoeazPop *pop){
    mgf_pop_preploop(pop);

    for (int i = 0; i < size; ++i) {
        mgf_indv_free((Individual*)indv);
        indv += indvsize;
    }
    free(pop->indv);
    free(pop);
}

mbool mgf_pop_evaluate(MoeazPop *pop, Mop *mop) {
    Individual *indv;
    for (int i = 0; i < pop->size; ++i){
        indv = mgf_pop_get_indv(pop, i);
            if (!mop_evaluate(mop, indv)) return mfalse;
    }
    return mtrue;
}

Individual* mgf_pop_get_indv(MoeazPop *pop, unsigned int pos){
    mgf_pop_preploop(pop);
    UNUSED(size);
    return (Individual*)(indv + (indvsize * pos));
}

Individual* mgf_pop_current_indv(MoeazPop *pop){
    mgf_pop_preploop(pop);
    UNUSED(size);
    return (Individual*)(indv + (indvsize * pop->current));
}

// void until NULL
void * mgf_pop_next(MoeazPop *pop){
    mgf_pop_preploop(pop);
    UNUSED(size);
    pop->current++;
    if (pop->current <= pop->size) {
        return (Individual*)(indv + (indvsize * pop->current));
    } else {
        return 0;
    }
}

void mgf_pop_reset_cursor(MoeazPop *pop)
{
    pop->current = 0;
}

void mgf_pop_merge(MoeazPop *to, MoeazPop *from1, MoeazPop *from2)
{
    unsigned int i, k;
//    assert((pop1->size + pop2->size) == pop3->size);
    for (i = 0; i < from1->size; i++) {
        mgf_indv_copy(mgf_pop_get_indv(to,i),mgf_pop_get_indv(from1,i));
    }
    for (k = 0; k < from2->size; i++, k++) {
        mgf_indv_copy(mgf_pop_get_indv(to,i),mgf_pop_get_indv(from2,k));
    }
}

