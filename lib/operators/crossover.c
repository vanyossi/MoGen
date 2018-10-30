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
#include <math.h>

#include "rand.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

static void (*cross_PNX[3])(Mop* mop, MoeazIndv *parent1, MoeazIndv* parent2, MoeazIndv* c1, MoeazIndv* c2) = {
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


void PNX_real(Mop* mop, MoeazIndv *parent1, MoeazIndv* parent2, MoeazIndv* c1, MoeazIndv* c2){
//    assert(mop->set.type != 1);
    unsigned int j;
    double s;
    double eta = 0.5;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    indv_real_data *par1 = parent1->x.data;
    indv_real_data *par2 = parent2->x.data;
    indv_real_data *chi1 = c1->x.data;
    indv_real_data *chi2 = c2->x.data;

    if (rnd_perc() > mop->solver->bias.cxprob)
    {
        memcpy(chi1, par1, sizeof(indv_real_data) * mop->set.ndec);
        memcpy(chi2, par2, sizeof(indv_real_data) * mop->set.ndec);
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

void PNX_bin(Mop* mop, MoeazIndv *parent1, MoeazIndv* parent2, MoeazIndv* c1, MoeazIndv* c2) {
    //    assert(mop->set.type != 2);
    unsigned int j;
    double s;
    double eta = 0.5;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    indv_bin_data *par1 = parent1->x.data;
    indv_bin_data *par2 = parent2->x.data;
    indv_bin_data *chi1 = c1->x.data;
    indv_bin_data *chi2 = c2->x.data;

    if (rnd_perc() > mop->solver->bias.cxprob)
    {
        memcpy(chi1, par1, sizeof(indv_bin_data) * mop->set.ndec);
        memcpy(chi2, par2, sizeof(indv_bin_data) * mop->set.ndec);
        return;
    }

    for (j = 0; j < mop->set.ndec; j++)
    {
        /* crossover */
        s = (indv_bin_data) abs(par2[j] - par1[j]) / eta;
        chi1[j] = (indv_bin_data)((rnd_perc() < .5) ? rnd_normal(par1[j], s) : rnd_normal(par2[j], s));
        chi2[j] = (indv_bin_data)((rnd_perc() < .5) ? rnd_normal(par1[j], s) : rnd_normal(par2[j], s));

        /* child1 bounds */
        chi1[j] = (indv_bin_data)((chi1[j] < low_bound[j]) ? low_bound[j] : chi1[j]);
        chi1[j] = (indv_bin_data)((chi1[j] > up_bound[j]) ? up_bound[j] : chi1[j]);

        /* child2 bounds */
        chi2[j] = (indv_bin_data)((chi2[j] < low_bound[j]) ? low_bound[j] : chi2[j]);
        chi2[j] = (indv_bin_data)((chi2[j] > up_bound[j]) ? up_bound[j] : chi2[j]);
    }
}

void PNX_mixed(Mop* mop, MoeazIndv *parent1, MoeazIndv* parent2, MoeazIndv* c1, MoeazIndv* c2) {
    //    assert(mop->set.type != 3);
    double s;
    double eta = 0.5;

    double *low_bound = mop->limits.xmin;
    double *up_bound = mop->limits.xmax;

    Multiarray *par1 = parent1->x.mix;
    Multiarray *par2 = parent2->x.mix;
    Multiarray *chi1 = c1->x.mix;
    Multiarray *chi2 = c2->x.mix;

    if (rnd_perc() > mop->solver->bias.cxprob)
    {
        memcpy(chi1, par1, sizeof(Multiarray) * mop->set.ndec);
        memcpy(chi2, par2, sizeof(Multiarray) * mop->set.ndec);
        return;
    }

    double m_p1;
    double m_p2;
    double m_c1;
    double m_c2;
    for (unsigned int j = 0; j < mop->set.ndec; j++)
    {
        // Generalize to double to work with GSL funcs.
        m_p1 = *(double*)mua_value_at(par1, j);
        m_p2 = *(double*)mua_value_at(par2, j);
        /* crossover */
        s = fabs(m_p2 - m_p1) / eta;
        m_c1 = (rnd_perc() < .5) ? rnd_normal(m_p1, s) : rnd_normal(m_p2, s);
        m_c2 = (rnd_perc() < .5) ? rnd_normal(m_p1, s) : rnd_normal(m_p2, s);

        /* child1 bounds */
        m_c1 = (m_c1 < low_bound[j]) ? low_bound[j] : m_c1;
        m_c1 = (m_c1 > up_bound[j]) ? up_bound[j] : m_c1;

        /* child2 bounds */
        m_c2 = (m_c2 < low_bound[j]) ? low_bound[j] : m_c2;
        m_c2 = (m_c2 > up_bound[j]) ? up_bound[j] : m_c2;

        // After cross all are double, as no Bin/Real cross is defined yet
        // possibly depending on cross op, a flag must be set.
        mua_set_double(par1, j, m_p1);
        mua_set_double(par2, j, m_p2);
        mua_set_double(chi1, j, m_c1);
        mua_set_double(chi1, j, m_c2);
    }
}
