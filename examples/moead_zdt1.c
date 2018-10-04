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

#include "mogen_mop.h"

void mopeval_zdt1(Mop* mop, unsigned int idx){

    MoeazIndv* cind = &mop->pop->indv[idx];
    double f1, f2, g, h, sum;
    unsigned int i;

    f1 = cind->x.real[0];
    sum = 0.0;
    for (i = 1; i < cind->xsize; i++)
    {
        sum += cind->x.real[i];
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
    Mop* mop = mogen_mop("ZDT", MOP_REAL | MOP_CONTIGUOUS);
    mop->set.ndec = nreal;
    mop->set.nobj = nobjs;

//    int n_var[5]; /* number of variables */
//
//    enum { zdt1, zdt2, zdt3, zdt4, zdt6 };
//    n_var[zdt1] = n_var[zdt2] = n_var[zdt3] = 30;
//    n_var[zdt4] = n_var[zdt6] = 10;

    mop->set.ncons = 0;
    mop->limits.xmin = (double*) malloc(sizeof(double) * mop->set.ndec);
    mop->limits.xmax = (double*) malloc(sizeof(double) * mop->set.ndec);

    for (int j = 0; j < mop->set.ndec; j++)
    {
        mop->limits.xmin[j] = 0.0;
        mop->limits.xmax[j] = 1.0;
    }

    mop->evaluate = mopeval_zdt1;

    return mop;
}



int main(int argc, char const *argv[]) {

    // declare mop problem
    Mop* my_mop = init_zdt(14,2);
    mop_print(my_mop);

    // declare algorithm

    // initialize population
    MoeazPop* m_pop = moeaz_pop_init(my_mop, 3);

    printf("size of ind %d and type %d", m_pop->indv[0].xsize, m_pop->indv[0].type);
    printf("size of ind %d and type %d", m_pop->indv[1].xsize, m_pop->indv[1].type);


    // evaluate runs


    // analize pop


    // free resources.
    moeaz_pop_free(m_pop);

    return 0;
}

