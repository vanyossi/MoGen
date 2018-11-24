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

#include "mops/zdt.h"
#include "solvers/mgf_nsga2.h"
#include "operators/mgf_operators.h"

#include "rand.h"

#include <stdio.h>

int main(int argc, char const *argv[]) {
    set_random(0.123452);

    Mop *zdt = mop_zdt(ZDT1, MOP_REAL | MOP_CONTIGUOUS);
    Moa *nsga2 = moa_nsga2(zdt);

    mgf_moa_new_pop(nsga2, 100, mgf_indvtype_nsga2(nsga2));
    mgf_pop_init(nsga2);

    moa_cross_setup(nsga2, CX_PNX);
    moa_mutation_setup(nsga2, MUT_PBM);
    moa_nsga2_cross_eta(mgf_moa_nsga2_data(nsga2), 10);

    mop_solve(zdt, 300);

    double res[2];
    for (int i = 0; i < zdt->pop->size; ++i) {
        res[0] = mogen_mop_getindv(zdt,i)->f[0];
        res[1] = mogen_mop_getindv(zdt,i)->f[1];
        printf("%.10f, %.10f\n", res[0], res[1]);
    }

    printf("total eval/runs: %d, %d \n", zdt->report.total.evals, zdt->report.total.gens);
    printf("time taken %s: %.8f, %ld\n",
        zdt->solver->name,
        mogen_time_ms2sec(zdt->report.total.t_elapsed),
        zdt->report.total.t_elapsed);
    return 0;
}
