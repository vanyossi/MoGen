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
    Moa *moa = mgf_moa_new(mop, "Secant", mgf_moatype_mono());
    mgf_moa_mono_set_solver(moa, moa_secant_solver);
    mgf_moa_mono_set_epsilon(moa, epsilon);
    return moa;
}

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
    mgf_indv_get_mono_buffer(indv)->error = fabs(fx(xn));
}

//mbool moa_mono_stop(Moa *moa, MoaStopCriterion criterion){
//    if (criterion == MGN_STOPIF_EPSILON){
//        return ( (mgf_moa_get_mono_buffer(moa)->epsilon >= (unsigned int) criterion) ? mtrue : mfalse;
//    }
//    return moa_stop(moa, criterion);
//}
