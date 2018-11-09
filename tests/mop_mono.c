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

#include "rand.h"
#include "secant.h"

double cuadratic(double x) {
    return (x*x) - 6;
} // +/- 2,44944897428

double polinom(double x){
    return (x*x*x) + (4*x*x) - 10;
} // res = 1.365230013

double polinom3(double x){
    return 3*(x +1)*(x -.5)*(x -1);
} // res = -1, 0.5, 1

int main(int argc, char const *argv[]) {

    // initialize seed
    set_random(0.141559);

    Mop* mypolin = mogen_mop("Polinom_3x", MOP_REAL | MOP_CONTIGUOUS, sizeof(MopMono));
    mop_set_params(mypolin, 2, 1, 0);

    double xmin = -6;
    double xmax = 4;

    mop_set_limits_ndec(mypolin, &xmin, &xmax, 1);
    moa_secant(mypolin, .000001);

//    printf("epsilon = %.6f", mgf_moa_get_mono_buffer(mypolin->solver)->epsilon );

    mgf_pop_alloc(mypolin->solver, 10, mgf_indvtype_mono(mypolin->solver));

    mono_fx eval_funcs[3] =  {cuadratic, polinom, polinom3};

    for (int j = 0; j < 3; ++j) {
        mono_fx eval_func = eval_funcs[j];

        mgf_pop_init(mypolin->solver);
        mop_mono_assign_fx(mypolin, eval_func);
        moa_stopat_eval(mypolin->solver, 900);

        mop_solve(mypolin,50);

        double res;
        for (int i = 0; i < mypolin->pop->size; ++i) {
            res = mogen_mop_getindv(mypolin,i)->f[0];
            printf("result %.10f, eval = %.10f\n", res, eval_func(res));
        }
        printf("total eval/runs: %d, %d \n", mypolin->report.total.evals, mypolin->report.total.gens);
        printf("end. %d\n\n", j);
    }
    printf("time taken %s: %.8f, %ld\n", mypolin->solver->name, mogen_time_ms2sec(mypolin->report.total.t_elapsed), mypolin->report.total.t_elapsed);
    return 0;
}
