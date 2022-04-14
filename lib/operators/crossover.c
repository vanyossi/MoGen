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
#include <limits.h>
#include <math.h>

#include "rand.h"

#include "mogen_mop.h"
#include "mgf_moa.h"
#include "mgf_operators.h"

//static void (*cross_PNX[3])(Mop*, Individual*, Individual*, Individual*, Individual*) = {
//    PNX_real, PNX_bin, PNX_mixed
//};

void moa_cross_setup(Moa *moa, CXType cx_type) {
    UNUSED(moa);
    switch(cx_type){
        case CX_PNX:
            // @TODO this wont work for other types than real.
            mgf_opset_crossover(PNX_mixed);
            break;
        default:
            break;
    }
}

void PNX_real(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2){
    struct moa_cm_container cross = mop->solver->type->get_crossover_vals(mop->solver);
//    assert(mop->set.type != 1);
    unsigned int j;
    double s;
    double eta = (cross.eta == 0)? 0.5 : cross.eta;
    double cprob = (cross.prob == 0)? 0.5 : cross.prob;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    double *par1 = mgf_indv_get_realdatapointer(parent1);
    double *par2 = mgf_indv_get_realdatapointer(parent2);
    double *chi1 = mgf_indv_get_realdatapointer(c1);
    double *chi2 = mgf_indv_get_realdatapointer(c2);

    if (rnd_perc() > cprob) {
        memcpy(chi1, par1, sizeof(double) * mop->set.xsize);
        memcpy(chi2, par2, sizeof(double) * mop->set.xsize);
        return;
    }

    for (j = 0; j < mop->set.xsize; j++)
    {
        /* crossover */
        s = fabs(par2[j] - par1[j] + 1e-14) / eta;
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
    struct moa_cm_container cross = mop->solver->type->get_crossover_vals(mop->solver);
    //    assert(mop->set.type != 2);
    unsigned int j;
    double s;
    double eta = (cross.eta == 0)? 0.5 : cross.eta;
    double cprob = (cross.prob == 0)? 0.5 : cross.prob;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    if (rnd_perc() > cprob) {
        memcpy(c1->bin, parent1->bin, mop->set.xsize / WORD_BIT);
        memcpy(c2->bin, parent2->bin, mop->set.xsize / WORD_BIT);
        return;
    }

    for (j = 0; j < mop->set.bsize; j++)
    {
        /* crossover */
        s = mgf_indv_get_bin(parent2, j) - mgf_indv_get_bin(parent1, j) / eta;
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
    struct moa_cm_container cross = mop->solver->type->get_crossover_vals(mop->solver);
    //    assert(mop->set.type != 3);
    unsigned int j;
    double eta = (cross.eta == 0)? 0.5 : cross.eta;
    double cprob = (cross.prob == 0)? 0.5 : cross.prob;

    if (rnd_perc() > cprob) {
        memcpy(c1->real, parent1->real, sizeof(double) * mop->set.xsize);
        memcpy(c1->integer, parent1->integer, sizeof(int) * mop->set.isize);
        memcpy(c1->bin, parent1->bin, mop->set.bsize / WORD_BIT);

        memcpy(c2->real, parent2->real, sizeof(double) * mop->set.xsize);
        memcpy(c2->integer, parent2->integer, sizeof(int) * mop->set.isize);
        memcpy(c2->bin, parent2->bin, mop->set.bsize / WORD_BIT);
        return;
    }

    double s;
    double *par1 = mgf_indv_get_realdatapointer(parent1);
    double *par2 = mgf_indv_get_realdatapointer(parent2);
    double *chi1 = mgf_indv_get_realdatapointer(c1);
    double *chi2 = mgf_indv_get_realdatapointer(c2);

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    for (j = 0; j < mop->set.xsize; j++)
    {
        /* crossover */
        s = fabs(par2[j] - par1[j] + 1e-14) / eta;
        chi1[j] = (rnd_perc() < .5) ? rnd_normal(par1[j], s) : rnd_normal(par2[j], s);
        chi2[j] = (rnd_perc() < .5) ? rnd_normal(par1[j], s) : rnd_normal(par2[j], s);

        /* child1 bounds */
        chi1[j] = (chi1[j] < low_bound[j]) ? low_bound[j] : chi1[j];
        chi1[j] = (chi1[j] > up_bound[j]) ? up_bound[j] : chi1[j];

        /* child2 bounds */
        chi2[j] = (chi2[j] < low_bound[j]) ? low_bound[j] : chi2[j];
        chi2[j] = (chi2[j] > up_bound[j]) ? up_bound[j] : chi2[j];
    }

    //missing integer cross
    int s_i;
    int *par1_i = mgf_indv_get_integerdatapointer(parent1);
    int *par2_i = mgf_indv_get_integerdatapointer(parent2);
    int *chi1_i = mgf_indv_get_integerdatapointer(c1);
    int *chi2_i = mgf_indv_get_integerdatapointer(c2);

    int *low_bound_i = mop->limits.imin;
    int *up_bound_i = mop->limits.imax;

    for (j = 0; j < mop->set.isize; j++)
    {
        /* crossover */
        s_i = (abs(par2_i[j] - par1_i[j]) / (int) round(eta));
        chi1_i[j] = (rnd_perc() < .5) ? rnd_int(par1_i[j], s_i) : rnd_int(par2_i[j], s_i);
        chi2_i[j] = (rnd_perc() < .5) ? rnd_int(par1_i[j], s_i) : rnd_int(par2_i[j], s_i);

        /* child1 bounds */
        chi1_i[j] = (chi1_i[j] < low_bound_i[j]) ? low_bound_i[j] : chi1_i[j];
        chi1_i[j] = (chi1_i[j] > up_bound_i[j]) ? up_bound_i[j] : chi1_i[j];

        /* child2 bounds */
        chi2_i[j] = (chi2_i[j] < low_bound_i[j]) ? low_bound_i[j] : chi2_i[j];
        chi2_i[j] = (chi2_i[j] > up_bound_i[j]) ? up_bound_i[j] : chi2_i[j];

    }

    for (j = 0; j < mop->set.bsize; j++)
    {
        /* crossover */
        mgf_indv_set_bin(c1,j, (rnd_perc() < .5) ? mgf_indv_get_bin(parent1, j) : mgf_indv_get_bin(parent2, j) );
        mgf_indv_set_bin(c2,j, (rnd_perc() < .5) ? mgf_indv_get_bin(parent1, j) : mgf_indv_get_bin(parent2, j) );
    }
}
