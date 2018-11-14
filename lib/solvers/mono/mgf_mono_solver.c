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
    mgf_indv_get_mono_buffer(indv)->error = MAXFLOAT;
}

void moa_mono_alloc(Moa *moa){
    mgf_moa_get_mono_buffer(moa)->epsilon = 1e-5;
}

// @TODO add missisng moa to moa_mono and back
// @TODO data access to epsilon should take moa_mono as arg

struct moa_mono_type* mgf_moa_get_mono_buffer(struct moa_t *moa){
    return mgf_moa_buffer(moa);
}

struct moa_type_t* mgf_moatype_mono(){
    return mgf_moatype_new(sizeof(struct moa_mono_type), moa_mono_alloc, 0, moa_mono_run, 0);
}

struct indv_t_mono_type* mgf_indv_get_mono_buffer(struct indv_t *indv){
    return mgf_indv_buffer(indv);
}

struct indv_type_t* mgf_indvtype_mono(Moa *moa){
    return mgf_indvtype_new(moa, sizeof(struct indv_t_mono_type), indv_mono_alloc, NULL, mgf_indv_free_std);
}

// to mop_monoobj assign fx
void mop_mono_assign_fx(Mop *mop, mono_fx f){
    ((MopMono*)mop)->fx = f;
}

void mgf_moa_mono_set_solver(Moa* moa, void (*evaluate)(Mop*, Individual*)){
    moa->mop->evaluate = evaluate;
}

void mgf_moa_mono_set_epsilon(Moa *moa, double epsilon){
    mgf_moa_get_mono_buffer(moa)->epsilon = epsilon;
}

// to mop mono_run
mbool moa_mono_run(Mop *mop) {
    MoaMono* secant = mgf_moa_get_mono_buffer(mop->solver);

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
