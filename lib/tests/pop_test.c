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

#include <stdio.h>

#include "rand.h"  // @TODO integrar main api

#include "mgf_moa.h"
#include "mogen_mop.h"
#include "mgf_population.h"

#include "solvers/mgf_nsga2.h"

int main(int argc, char const *argv[]) {
    set_random(0.141559);

    Mop* mop = mogen_mop("test_pop", MOP_REAL | MOP_INT,0);
    mop_set_params(mop, 5, 0, 5, 1, 0);

    double min = 0.2;
    double max = 0.8;
    int imin = 0;
    int imax = 1;
    mop_set_limits_ndec(mop, &min, &max, 1, &imin, &imax, 1);
    Moa* moa = mgf_moa_new(mop, "tests_moa_pop", mgf_moatype_nsga2());

    mgf_moa_new_pop(moa, 20, mgf_indvtype_nsga2(moa));
    mgf_pop_init(moa);

    MoeazPop *copy_pop = mgf_pop_alloc(20, mgf_indvtype_nsga2(moa));
    for (int i = 0; i < mop->pop->size; ++i) {
        mgf_indv_copy(
            mgf_pop_get_indv(copy_pop, i),
            mgf_pop_get_indv(mop->pop, i)
        );
    }

    printf("xsize %d, indv type %d\n", mop->pop->indv->type->xsize, mop->pop->indv->xtype);

    Individual *indv;
    Individual *indv2;
    IndvidualType *type;
    int non_equal = 0;
    for (int j = 0; j < copy_pop->size - 1; ++j) {
        indv = mgf_pop_get_indv(copy_pop, j);
        indv2 = mgf_pop_get_indv(mop->pop, j);
        type = indv->type;

        non_equal += (indv->real[type->xsize -1] != indv2->real[type->xsize -1])? 1 : 0;
        non_equal += (indv->type != indv2->type)? 10: 0;
        non_equal += (indv->buffer_start != indv2->buffer_start)? 100: 0;
        non_equal += (indv_nsga2_crowdist(indv) != indv_nsga2_crowdist(indv2))? 1000: 0;
        non_equal += (indv->xtype != indv2->xtype)? 10000: 0;
        non_equal += (indv->integer[type->isize - 1] != indv2->integer[type->isize - 1])? 100000 : 0;

        if(non_equal){
            printf("Copy %d, %d is not perfect %0.15f %0.15f\n", j, non_equal, indv->real[1], indv2->real[1]);
        }

        if(type == NULL){
            printf("indv %d is null crowd = %p %g\n", j, type, type->get_crowdist(indv));
        }
    }
    mgf_pop_free(copy_pop);
    mgf_pop_free(mop->pop);

    return 0;
}
