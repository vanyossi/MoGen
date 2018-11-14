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

#include "mgf_mutation.h"

#include <math.h>
#include "rand.h"

#include "mgf_operators.h"

mgf_op_mutation mut_PBM[3] = {
    PBM_real, mut_bitwise, PBM_mut
};

static double PBM_real_on_x(double y, MutationSettings mutation, double xmin, double xmax);
static inline unsigned int mut_bitwise_on_x(unsigned int x);

void moa_mutation_setup(Moa *moa, MutType mut_type) {
    switch(mut_type){
        case MUT_PBM:
        default:
            mgf_opset_mutation(mut_PBM[moa->mop->set.type - 1]);
            break;
    }
}

void PBM_mut(Individual *indv, MutationSettings mutation, Mop_limit limits) {
    IndvidualType *indv_type = mgf_indv_type(indv);

    for (unsigned int j = 0; j < indv_type->xsize; j++) {
        if (rnd_perc() <= mutation.prob) {
            if (mgf_indv_value_isbin(indv, j)) {
                mgf_indv_get_bin( indv, mut_bitwise_on_x(mgf_indv_get_bin(indv, j)) );

            } else {
                mgf_indv_set_double( indv, j, PBM_real_on_x(
                    mgf_indv_get_double(indv, j), mutation, limits.xmin[j], limits.xmax[j]));
            }
        }
    }
}

static double PBM_real_on_x(double y, MutationSettings mutation, double xmin, double xmax) {

    double rnd, delta1, delta2, mut_pow, deltaq;
    double val, xy;

    delta1 = (y - xmin) / (xmax - xmin);
    delta2 = (xmax - y) / (xmax - xmin);
    rnd = rnd_perc();
    mut_pow = 1.0 / (mutation.eta + 1.0);
    if (rnd <= 0.5) {
        xy = 1.0 - delta1;
        val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (mutation.eta + 1.0)));
        deltaq = pow(val, mut_pow) - 1.0;

    } else {
        xy = 1.0 - delta2;
        val = 2.0 * (1.0 - rnd)
              + 2.0 * (rnd - 0.5) * (pow(xy, (mutation.eta + 1.0)));
        deltaq = 1.0 - (pow(val, mut_pow));
    }
    y = y + deltaq * (xmax - xmin);
    if (y < xmin) {
        y = xmin;
    }
    if (y > xmax) {
        y = xmax;
    }
    return y;
}

void PBM_real(Individual *indv, MutationSettings mutation, Mop_limit limits) {
    IndvidualType *indv_type = mgf_indv_type(indv);

    for (unsigned int j = 0; j < indv_type->xsize; j++) {
        if (rnd_perc() <= mutation.prob) {
            mgf_indv_set_double( indv, j, PBM_real_on_x(
                mgf_indv_get_double(indv, j), mutation, limits.xmin[j], limits.xmax[j]));
        }
    }
}

void mut_bitwise(Individual *indv, MutationSettings mutation, Mop_limit limits){
    UNUSED(limits);
    IndvidualType *indv_type = mgf_indv_type(indv);

    for (unsigned int j = 0; j < indv_type->xsize; j++) {
        if (rnd_perc() <= mutation.prob) {
            mgf_indv_get_bin( indv, mut_bitwise_on_x(mgf_indv_get_bin(indv, j)) );
        }
    }

}

static inline unsigned int mut_bitwise_on_x(unsigned int x) {
    return x ^ 1;
}
