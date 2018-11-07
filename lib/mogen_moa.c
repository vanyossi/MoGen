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

#include "mogen_moa.h"
#include "mogen_mop.h"

#include <stdlib.h>
#include <string.h>

Moa *moa_init(struct mop_t *mop, char *name, MoaTypes type, size_t mem_size) {
    size_t alloc_size = sizeof(Moa);
    if(mem_size > alloc_size){
        alloc_size = mem_size;
    }
    Moa* new_moa = calloc(1, alloc_size);

    new_moa->type = type;
    strcpy(new_moa->name, name);
    new_moa->mop = mop;
    mop->solver = new_moa;

    mop->solver->stop = moa_stop;

    return new_moa;
}

mbool moa_run(struct mop_t *mop){
    for (int i = 0; i < mop->pop->size; ++i){
        Individual* indv = mgf_pop_get_indv(mop->pop,i);
        if (!mop_evaluate(mop, indv)){
            return mfalse;
        }
    }
    return mtrue;
}

void moa_stopat_gen(Moa *moa, unsigned int gen) {
    moa->stops.gens = gen;
}

void moa_stopat_eval(Moa *moa, unsigned int eval){
    moa->stops.evals = eval;
}

mbool moa_stop(struct mogen_moa_t* moa, MoaStopCriterion criterion){
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
