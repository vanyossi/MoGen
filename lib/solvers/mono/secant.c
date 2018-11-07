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

#include "secant.h"

#include <stdlib.h>
#include <math.h>


#include <stdio.h>

void indv_mono_alloc(Mop *mop, struct indv_t *indv) {
    ((struct indv_t_mono_type*)mgf_indv_buffer(indv))->error = MAXFLOAT;
}

struct indv_type_t* mgf_indvtype_mono(Moa *moa){
    return mgf_indvtype_new(moa, sizeof(struct indv_t_mono_type), indv_mono_alloc, mgf_indv_free_std);
}

Moa *moa_secant(Mop *mop, double epsilon) {
    MoaSecant* secant = (MoaSecant*) moa_init(mop, "Secant", MOA_MONO, sizeof(MoaSecant));

    secant->moa.run = moa_secant_run;
//    secant->moa.extra_indv_alloc = moa_mono_indv_alloc;

//    secant->error = malloc(sizeof(double) * mop->set.nobj);
    secant->epsilon = epsilon;

    mop->evaluate = moa_secant_solver;
    mop->solver = (Moa*)secant;

    return (Moa*)secant;
}

//#include <stdio.h>
//Moa *moa_secant(Mop *mop, double epsilon) {
//    MoaSecant* secant = malloc(sizeof(MoaSecant));
//    secant->moa = moa_init(mop, "Secant", MOA_MONO);
//
//    secant->moa->run = moa_secant_run;
//
//    secant->error = malloc(sizeof(double) * mop->set.nobj);
//    secant->epsilon = epsilon;
//
//    mop->evaluate = moa_secant_solver;
//    mop->solver = (Moa*)secant;
//
//    return (Moa*)secant;
//}


// to mop_monoobj assign fx
void mop_secant_assign_fx(Mop *mop, mono_fx f){
    ((MopMono*)mop)->fx = f;
}

// to mop mono_run
mbool moa_secant_run(Mop *mop) {
    MoaSecant* secant = (MoaSecant*)(mop->solver);

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

//#include <stdio.h>
void moa_secant_solver(Mop* mop, Individual *indv){
    MopMono *sec_mop = (MopMono*)mop;
    double *x = mgf_indv_get_realdatapointer(indv);

    mono_fx fx = sec_mop->fx;
    double xn = solve(fx, x);


    if( !isnan(xn) || isfinite(xn)) {
        indv->f[0] = xn;
        x[0] = x[1];
        x[1] = xn;
    }
    // calculate error
    ((struct indv_t_mono_type*)mgf_indv_buffer(indv))->error = fabs(fx(xn));
}


//mbool moa_mono_stop(Moa *moa, MoaStopCriterion criterion){
//    if (criterion == MGN_STOPIF_EPSILON){
//        return ( ((MoaSecant*)moa)->epsilon >= (unsigned int) criterion) ? mtrue : mfalse;
//    }
//    return moa_stop(moa, criterion);
//}
