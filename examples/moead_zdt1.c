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
#include <stdlib.h>
#include <assert.h>

#include <math.h>
#include "rand.h"

#include "mogen_mop.h"
#include "lib/operators/crossover.h"

int moead_run(Mop* mop){
    for (unsigned int i = 0; i < mop->pop->size; ++i) {
        mop_evaluate(mop, mogen_mop_getindv(mop, i));
        mop->solver->cross(mop,
            mogen_mop_getindv(mop, (unsigned) rnd_int(0,mop->pop->size)),
            mogen_mop_getindv(mop, (unsigned) rnd_int(0,mop->pop->size)),
            mogen_mop_getindv(mop, (unsigned) rnd_int(0,mop->pop->size)),
            mogen_mop_getindv(mop, (unsigned) rnd_int(0,mop->pop->size))
        );
    }
    return mtrue;
}

Moa* mymoead_init(Mop *mop, double cp, double ci){

    Moa* new_moa = moa_init(mop, "MOEAD", MOA_DECOMP, 0);
    new_moa->bias.cxprob = cp;
    new_moa->run = moead_run;

    return new_moa;
}

void mopeval_zdt1(Mop* mop, MoeazIndv* cind){

    indv_real_data *x = cind->x.data;

    double f1, f2, g, h, sum;
    unsigned int i;

    f1 = x[0];
    sum = 0.0;
    for (i = 1; i < cind->xsize; i++)
    {
        sum += x[i];
    }
    g = 1.0 + 9.0 * sum / (cind->xsize - 1.0);
    h = 1.0 - sqrt(f1 / g);
    f2 = g * h;
    cind->f[0] = f1;
    cind->f[1] = f2;

    return;
}

Mop* init_zdt(unsigned int nreal, unsigned int nobjs){

    assert(nreal > nobjs);
    Mop* mop = mogen_mop("ZDT", MOP_REAL | MOP_CONTIGUOUS, 0);

    mop_set_params(mop, nreal, nobjs, 0);

    double xmin_limits[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    double xmax_limits[5] = {1.0, 2.0, 1.0, 0.5, 1.5};

    mop_set_limits_ndec(mop, xmin_limits, xmax_limits, 5);

//    int n_var[5]; /* number of variables */
//
//    enum { zdt1, zdt2, zdt3, zdt4, zdt6 };
//    n_var[zdt1] = n_var[zdt2] = n_var[zdt3] = 30;
//    n_var[zdt4] = n_var[zdt6] = 10;

    mop->evaluate = mopeval_zdt1;

    return mop;
}


int main(int argc, char const *argv[]) {

    // initialize seed
    set_random(0.56);

    // declare mop problem
    Mop* my_mop = init_zdt(10,2);
    mop_print(my_mop);

    // declare algorithm
    Moa *my_moa = mymoead_init(my_mop, 1.0,1.0);

    // initialize population
    MoeazPop* m_pop = mogen_pop_alloc(my_mop, 20);
    mogen_pop_init(my_mop);

    moa_cross_setup(my_moa, CX_PNX);

    printf("size of ind %d and type %d\n", m_pop->indv[0].xsize, m_pop->indv[0].type);
    printf("size of ind %d and type %d\n", m_pop->indv[1].xsize, m_pop->indv[1].type);

    double *vals;

    for (int i = 0; i < 3; ++i) {
        printf("ind %d: ", i);
        vals = get_data_real((&my_mop->pop->indv[i]));
        for (int j = 0; j < my_mop->pop->indv[i].xsize; ++j) {
            printf("%d:%.3f, ", j, vals[j]);
        }
        puts("");
    }
    puts("");
    // evaluate runs
    my_moa->run(my_mop);

    // analize pop
    for (int i = 0; i < 3; ++i) {
        printf("ind %d: ", i);
        vals = get_data_real((&my_mop->pop->indv[i]));
        for (int j = 0; j < my_mop->pop->indv[i].xsize; ++j) {
            printf("%d:%.3f, ", j, vals[j]);
        }
        puts("");
    }

    for (int i = 0; i < 2; ++i) {
        printf("zdt ind %d: ", i);
        vals = my_mop->pop->indv[i].f;
        for (int j = 0; j < my_mop->pop->indv[i].fsize; ++j) {
            printf("%d:%.3f, ",j, vals[j]);
        }
        puts("");
    }

    // free resources.
    mogen_pop_free(m_pop);

    return 0;
}

