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


#include "zdt_test.h"

#include <stdio.h>
#include <math.h>

#include "../zdt.h"
#include "lib/mogen_mop.h"

#define POP_SIZE 3

extern void mgf_zdt1(Mop *mop, Individual *indv);

int double_isequal(double a, double b, double delta);

struct indv_results {
    double res[2];
};

int main(int argc, char const *argv[]) {

    int fail = 0;

    Mop *mop = mogen_mop("test_zdt", MOP_REAL | MOP_CONTIGUOUS, 0);
    mop_set_params(mop, 4, 0, 0, 2, 0);
    Moa *moa = mgf_moa_new(mop, "zdt_testmoa", mgf_moa_std());

    MoeazPop *pop = mgf_pop_alloc(POP_SIZE, mgf_indvtype_std(moa));

    struct indv_results indv_res[POP_SIZE] = {
        {1, 22.708497},
        {2, 28.397675},
        {3, 34.252659}
    };
    // manual init
    double *F;
    for (int i = 0; i < pop->size; ++i) {
        for (int j = 0; j < mop->set.xsize; ++j) {
            mgf_indv_set_double(mgf_pop_get_indv(pop, i), j, (i+j+1));
//            printf("val %d is %.7f\n", j, mgf_indv_get_double(mgf_pop_get_indv(pop, i),j));
        }
        mgf_zdt1(mop, mgf_pop_get_indv(pop, i));

        for (int j = 0; j < mop->set.nobj; ++j) {
            F = mgf_indv_get_solution_pointer(mgf_pop_get_indv(pop, i));
            if (!double_isequal(F[j], indv_res[i].res[j], 0.00001) ){
                fail = 1;
            } else {
            //    printf("f %d,%d are equal\n", i, j);
            }
        }
    }

    if (!fail){
        printf("Zdt1 test passed!\n");
    }
    return fail;
}

int double_isequal(double a, double b, double delta){
    return (fabs(a - b) < delta);
}
