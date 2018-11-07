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


#include "mgf_mono_solver.h"


void indv_mono_alloc(Mop *mop, struct indv_t *indv) {
    ((struct indv_t_mono_type*)mgf_indv_buffer(indv))->error = MAXFLOAT;
}

struct indv_type_t* mgf_indvtype_mono(Moa *moa){
    return mgf_indvtype_new(moa, sizeof(struct indv_t_mono_type), indv_mono_alloc, mgf_indv_free_std);
}

// to mop_monoobj assign fx
void mop_mono_assign_fx(Mop *mop, mono_fx f){
    ((MopMono*)mop)->fx = f;
}

Moa *moa_mono(Mop *mop, char *name, double epsilon, void (*evaluate)(Mop*, Individual*)) {
    MoaMono* moamono = (MoaMono*) moa_init(mop, name, MOA_MONO, sizeof(MoaMono));

    moamono->moa.run = moa_mono_run;

//    moamono->error = malloc(sizeof(double) * mop->set.nobj);
    moamono->epsilon = epsilon;

    mop->evaluate = evaluate;
    mop->solver = (Moa*)moamono;

    return (Moa*)moamono;
}

// to mop mono_run
mbool moa_mono_run(Mop *mop) {
    MoaMono* secant = (MoaMono*)(mop->solver);

    struct indv_t_mono_type* indv_data;
    Individual* indv;
    int ev = 0;
    for (int i = 0; i < mop->pop->size; ++i){
        indv = mgf_pop_get_indv(mop->pop, i);
        indv_data = mgf_indv_buffer(indv);
        if ( indv_data->error > secant->epsilon) {
            ev++;
            if (!mop_evaluate(mop, indv)){
                return mfalse;
            }
        }
    }
    if (!ev) {
        return mfalse;
    }
    return mtrue;
}
