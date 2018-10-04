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


#include "crossover.h"

#include <memory.h>
#include <math.h>

#include "rand.h"

void crossover(Mop *mop, unsigned int p1, unsigned int p2, unsigned int c1, unsigned int c2){
    mop->cross(mop, mogen_mop_getindv(mop, p1), mogen_mop_getindv(mop, p2), mogen_mop_getindv(mop, c1), mogen_mop_getindv(mop, c2));
}

void PNX(Mop *mop, MoeazIndv *p1, MoeazIndv *p2, MoeazIndv *c1, MoeazIndv *c2){
    if (mop->set.type != 1) {
        // mop is not real, crossover might give unexpected results.
    }
    unsigned int j;
    double s;
    double eta = 0.5;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    double *par1 = p1->x.real;
    double *par2 = p2->x.real;
    double *chi1 = c1->x.real;
    double *chi2 = c2->x.real;

    if (rnd_perc() > mop->solver->algorithm.moea.cross_prob)
    {
        memcpy(chi1, par1, sizeof(double) * mop->set.ndec);
        memcpy(chi2, par2, sizeof(double) * mop->set.ndec);
        return;
    }

    for (j = 0; j < mop->set.ndec; j++)
    {
        /* crossover */
        s = fabs(par2[j] - par1[j]) / eta;
        chi1[j] = (rnd_perc() < .5) ? rnd_normal(par1[j], s) : rnd_normal(par2[j], s);
        chi2[j] = (rnd_perc() < .5) ? rnd_normal(par1[j], s) : rnd_normal(par2[j], s);

        /* child1 bounds */
        chi1[j] = (chi1[j] < low_bound[j]) ? low_bound[j] : chi1[j];
        chi1[j] = (chi1[j] > up_bound[j]) ? up_bound[j] : chi1[j];

        /* child2 bounds */
        chi2[j] = (chi2[j] < low_bound[j]) ? low_bound[j] : chi2[j];
        chi2[j] = (chi2[j] > up_bound[j]) ? up_bound[j] : chi2[j];
    }
    return;
}

