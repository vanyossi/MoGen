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

int main(int argc, char const *argv[]) {
    set_random(0.141559);

    Mop* mop = mogen_mop("test_pop", MOP_REAL,0);
    mop_set_params(mop, 5,1,0);
    double min = 0.2;
    double max = 0.8;
    mop_set_limits_ndec(mop, &min, &max, 1);
    Moa* moa = mgf_moa_new(NULL, mop, "tests_moa_pop");

    mgf_pop_alloc(moa, 5, mgf_indvtype_std(moa));
    mgf_pop_init(moa);

    printf("xsize %d, indv type %d\n", mop->pop->indv->type->xsize, mop->pop->indv->xtype);

    for (int j = 0; j < mop->pop->size; ++j) {
        printf("indv %d:\n", j);
        for (int i = 0; i < mop->set.ndec; ++i) {
            printf("x%d: %0.5f, ", i, mgf_pop_get_indv(mop->pop, j)->real[i]);
        }
        puts("\n");
    }

    mgf_pop_free(mop->pop);

    return 0;
}
