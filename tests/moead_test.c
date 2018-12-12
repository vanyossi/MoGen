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
#include "solvers/mgf_moead.h"
#include "operators/mgf_operators.h"

#include "rand.h"

#include <stdio.h>

int main(int argc, char const *argv[]) {
    set_random(0.123452);

    Mop *zdt = mop_zdt(ZDTM1);
    // moead creates its own population and initializes it.
    Moa *moead = moa_moead(zdt, W_RES_210K);
    mgf_moa_set_scalarization(&mgf_moa_moead_data(moead)->s_m, SCLM_TCH_NORM); // @TODO need API
//    mop_set_params(zdt, 10, 30, 3, 2, 0);
//    double min = 0;
//    double max = 1;
//    int imin = 0;
//    int imax = 60;
//    mop_set_limits_ndec(zdt, &min, &max, 1, &imin, &imax, 1);

    moa_cross_setup(moead, CX_PNX);
    moa_mutation_setup(moead, MUT_PBM);
//    moa_moead_cross_eta(mgf_moa_moead_data(moead), 10)

    mop_solve(zdt, 1000);

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

    mgf_pop_free(zdt->pop);

    return 0;
}
