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

#include "mogen_mop.h"
#include "solvers/mgf_moead.h"
#include "solvers/mgf_nsga2.h"
#include "operators/mgf_operators.h"
#include "mops/zdt.h"

#include "rand.h"


Mop *mop_houseoptim();

#define HOUSE_CARAC_SIZE 96
static double carac[] = {6.1101, 5.5277, 8.5186, 7.0032, 5.8598, 8.3829, 7.4764, 8.5781, 6.4862, 5.0546, 5.7107, 14.164, 5.734, 8.4084, 5.6407, 5.3794, 6.3654, 5.1301, 6.4296, 7.0708, 6.1891, 20.27, 5.4901, 6.3261, 5.5649, 18.945, 12.828, 10.957, 13.176, 22.203, 5.2524, 6.5894, 9.2482, 5.8918, 8.2111, 7.9334, 8.0959, 5.6063, 12.836, 6.3534, 5.4069, 6.8825, 11.708, 5.7737, 7.8247, 7.0931, 5.0702, 5.8014, 11.7, 5.5416, 7.5402, 5.3077, 7.4239, 7.6031, 6.3328, 6.3589, 6.2742, 5.6397, 9.3102, 9.4536, 8.8254, 5.1793, 21.279, 14.908, 18.959, 7.2182, 8.2951, 10.236, 5.4994, 20.341, 10.136, 7.3345, 6.0062, 7.2259, 5.0269, 6.5479, 7.5386, 5.0365, 10.274, 5.1077, 5.7292, 5.1884, 6.3557, 9.7687, 6.5159, 8.5172, 9.1802, 6.002, 5.5204, 5.0594, 5.7077, 7.6366, 5.8707, 5.3054, 8.2934, 13.394, 5.4369};
static double result[] = {17.592, 9.1302, 13.662, 11.854, 6.8233, 11.886, 4.3483, 12.0, 6.5987, 3.8166, 3.2522, 15.505, 3.1551, 7.2258, 0.71618, 3.5129, 5.3048, 0.56077, 3.6518, 5.3893, 3.1386, 21.767, 4.263, 5.1875, 3.0825, 22.638, 13.501, 7.0467, 14.692, 24.147, -1.22, 5.9966, 12.134, 1.8495, 6.5426, 4.5623, 4.1164, 3.3928, 10.117, 5.4974, 0.55657, 3.9115, 5.3854, 2.4406, 6.7318, 1.0463, 5.1337, 1.844, 8.0043, 1.0179, 6.7504, 1.8396, 4.2885, 4.9981, 1.4233, -1.4211, 2.4756, 4.6042, 3.9624, 5.4141, 5.1694, -0.74279, 17.929, 12.054, 17.054, 4.8852, 5.7442, 7.7754, 1.0173, 20.992, 6.6799, 4.0259, 1.2784, 3.3411, -2.6807, 0.29678, 3.8845, 5.7014, 6.7526, 2.0576, 0.47953, 0.20421, 0.67861, 7.5435, 5.3436, 4.2415, 6.7981, 0.92695, 0.152, 2.8214, 1.8451, 4.2959, 7.2029, 1.9869, 0.14454, 9.0551, 0.61705};

/**
 * ZDT1 Mop
 * @param mop Multiobjective problem
 * @param indv Individual to be evaluated
 */
void ex_house_cost_optim(Mop *mop, Individual *indv)
{
//    IndvidualType *indv_type = mgf_indv_type(indv);
    double *F = mgf_indv_get_solution_pointer(indv);

    double h, sum;

    sum = 0.0;
    for (int k = 0; k < HOUSE_CARAC_SIZE; k++) {
        h = (mgf_indv_get_double(indv, 0) + mgf_indv_get_double(indv, 1) * carac[k]) - result[k];
        sum += h*h;
    }
    sum *= (0.5 * 1/(double)HOUSE_CARAC_SIZE);

    F[0] = sum;
    F[1] = 0;

    return;
}

/**
 * The parameter settings for the unconstrained DTLZ test problems
 * @param mop Multiobjectibe problem reference
 */
Mop *mop_houseoptim()
{
    Mop *mop = mogen_mop("House optim", (MopSpecs)(MOP_REAL | MOP_CONTIGUOUS), 0);
    mop_set_params(mop, 2, 0, 0, 2, 0); // default params

    double xmin = -10;
    double xmax = 10;
    mop_set_limits_ndec(mop, &xmin, &xmax, 1, 0, 0, 0);

    mop->evaluate = ex_house_cost_optim;

    return mop;
}

int main(int argc, char const *argv[]) {
    set_random(0.123452);

    Mop *problem = mop_houseoptim();
    // moead creates its own population and initializes it.
    Moa *moead = moa_moead(problem, W_RES_210K);
    mgf_moa_set_scalarization(&mgf_moa_moead_data(moead)->s_m, SCLM_TCH); // @TODO need API
    moa_stopat_eval(moead,1500);
    // mop_set_params(problem, 2, 0, 0, 2, 0);
//    double min = 0;
//    double max = 1;
//    int imin = 0;
//    int imax = 60;
//    mop_set_limits_ndec(zdt, &min, &max, 1, &imin, &imax, 1);

    moa_cross_setup(moead, CX_PNX);
    moa_mutation_setup(moead, MUT_PBM);
//    moa_moead_cross_eta(mgf_moa_moead_data(moead), 10)

    mop_solve(problem, 100);

    double res[3];
    for (int i = 0; i < problem->pop->size; ++i) {
        res[0] = mogen_mop_getindv(problem,i)->real[0];
        res[1] = mogen_mop_getindv(problem,i)->real[1];
        res[2] = mogen_mop_getindv(problem,i)->f[0];
        printf("%.10f, %.10f :: %.10f\n", res[0], res[1], res[2]);
    }

    printf("total eval/runs: %d, %d \n", problem->report.total.evals, problem->report.total.gens);
    printf("time taken %s: %.8f, %ld\n",
        problem->solver->name,
        mogen_time_ms2sec(problem->report.total.t_elapsed),
        problem->report.total.t_elapsed);

    mgf_pop_free(problem->pop);

    return 0;
}
