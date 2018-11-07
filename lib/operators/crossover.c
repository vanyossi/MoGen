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

#include "mogen_mop.h"
#include "mogen_moa.h"

#include <memory.h>
#include <limits.h>
#include <math.h>

#include "rand.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

static void (*cross_PNX[3])(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2) = {
    PNX_real, PNX_bin, PNX_mixed
};

void moa_cross_setup(Moa *moa, CXType cx_type) {
    switch(cx_type){
        case CX_PNX:
            moa->cross = cross_PNX[moa->mop->set.type - 1];
            break;
        default:
            break;
    }
}


void PNX_real(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2){
//    assert(mop->set.type != 1);
    unsigned int j;
    double s;
    double eta = 0.5;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    double *par1 = mgf_indv_get_realdatapointer(parent1);
    double *par2 = mgf_indv_get_realdatapointer(parent2);
    double *chi1 = mgf_indv_get_realdatapointer(c1);
    double *chi2 = mgf_indv_get_realdatapointer(c2);

    if (rnd_perc() > mop->solver->bias.cxprob)
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
}

void PNX_bin(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2) {
    //    assert(mop->set.type != 2);
    unsigned int j;
    double s;
    double eta = 0.5;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    if (rnd_perc() > mop->solver->bias.cxprob)
    {
        memcpy(c1->bin, parent1->bin, mop->set.ndec / WORD_BIT);
        memcpy(c2->bin, parent2->bin, mop->set.ndec / WORD_BIT);
        return;
    }

    for (j = 0; j < mop->set.ndec; j++)
    {
        /* crossover */
        s = abs(mgf_indv_get_bin(parent2, j) - mgf_indv_get_bin(parent1, j)) / eta;
        mgf_indv_set_bin(c1,j, (rnd_perc() < .5) ? (int)rnd_normal(mgf_indv_get_bin(parent1, j), s) : (int)rnd_normal(mgf_indv_get_bin(parent2, j), s) );
        mgf_indv_set_bin(c2,j, (rnd_perc() < .5) ? (int)rnd_normal(mgf_indv_get_bin(parent1, j), s) : (int)rnd_normal(mgf_indv_get_bin(parent2, j), s) );

        /* child1 bounds */
        mgf_indv_set_bin(c1,j, (mgf_indv_get_bin(c1, j) < low_bound[j]) ? (int)low_bound[j] : mgf_indv_get_bin(c1, j) );
        mgf_indv_set_bin(c1,j, (mgf_indv_get_bin(c1, j) > up_bound[j]) ? (int)up_bound[j] : mgf_indv_get_bin(c1, j) );

        /* child2 bounds */
        mgf_indv_set_bin(c2,j, (mgf_indv_get_bin(c2, j) < low_bound[j]) ? (int)low_bound[j] : mgf_indv_get_bin(c2, j) );
        mgf_indv_set_bin(c2,j, (mgf_indv_get_bin(c2, j) > up_bound[j]) ? (int)up_bound[j] : mgf_indv_get_bin(c2, j) );
    }
}

void PNX_mixed(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2) {
    //    assert(mop->set.type != 3);
    unsigned int j;
    double s;
    double eta = 0.5;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    if (rnd_perc() > mop->solver->bias.cxprob)
    {
        memcpy(c1->type_idx, parent1->type_idx, (2 * mop->set.ndec / WORD_BIT) + sizeof(double) * mop->set.ndec);
        memcpy(c2->type_idx, parent2->type_idx, (2 * mop->set.ndec / WORD_BIT) + sizeof(double) * mop->set.ndec);
        return;
    }
    // @TODO after crossing all variables are real!
    for (j = 0; j < mop->set.ndec; j++)
    {
        /* crossover */
        s = fabs(mgf_indv_value_at(parent2, j) - mgf_indv_value_at(parent1, j)) / eta;
        mgf_indv_set_double(c1,j, (rnd_perc() < .5) ? rnd_normal(mgf_indv_value_at(parent1, j), s) : rnd_normal(mgf_indv_value_at(parent2, j), s) );
        mgf_indv_set_double(c2,j, (rnd_perc() < .5) ? rnd_normal(mgf_indv_value_at(parent1, j), s) : rnd_normal(mgf_indv_value_at(parent2, j), s) );

        /* child1 bounds */
        mgf_indv_set_double(c1,j, (mgf_indv_get_bin(c1, j) < low_bound[j]) ? low_bound[j] : mgf_indv_value_at(c1, j) );
        mgf_indv_set_double(c1,j, (mgf_indv_get_bin(c1, j) > up_bound[j]) ? up_bound[j] : mgf_indv_value_at(c1, j) );

        /* child2 bounds */
        mgf_indv_set_double(c2,j, (mgf_indv_get_bin(c2, j) < low_bound[j]) ? low_bound[j] : mgf_indv_value_at(c2, j) );
        mgf_indv_set_double(c2,j, (mgf_indv_get_bin(c2, j) > up_bound[j]) ? up_bound[j] : mgf_indv_value_at(c2, j) );
    }
}
