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

#include <stdio.h>

Moa *moa_secant(Mop *mop, double epsilon) {
    MoaMono* secant = (MoaMono*) moa_init(mop, "Secant", MOA_MONO, sizeof(MoaMono));

    secant->moa.run = moa_mono_run;
//    secant->moa.extra_indv_alloc = moa_mono_indv_alloc;

//    secant->error = malloc(sizeof(double) * mop->set.nobj);
    secant->epsilon = epsilon;

    mop->evaluate = moa_secant_solver;
    mop->solver = (Moa*)secant;

    return (Moa*)secant;
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
//        return ( ((MoaMono*)moa)->epsilon >= (unsigned int) criterion) ? mtrue : mfalse;
//    }
//    return moa_stop(moa, criterion);
//}
