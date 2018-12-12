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


#ifndef MOGEN_MGF_MOEAD_H
#define MOGEN_MGF_MOEAD_H

#include "mgf_moa.h"
#include "decomposition/mgf_scalarization.h"
#include "decomposition/mgf_weights.h"

struct moa_moead_t {
    double cross_eta;
    double cross_prob;
    double mut_eta;
    double mut_prob;
    unsigned int t_size;   //!< Neibourghood size
    unsigned int st_id;
    unsigned int *st_p;
    double *st_nadir;
    double *z_ref;
    Individual *st_indv;
    WeightParams weights;
    struct mgf_neighborhood_set neighbor_set;
    struct mgf_scalar_method_t s_m;
};

MoaType* mgf_moatype_moead();

struct moa_moead_t* mgf_moa_moead_data(Moa* moa);

struct moa_cm_container moead_get_crossvals(Moa *self);

struct moa_cm_container moead_get_mutvals(Moa *self);

double moa_moead_cross_eta(struct moa_moead_t *moead, double value);

double moa_moead_cross_prob(struct moa_moead_t *moead, double value);

double moa_moead_mut_eta(struct moa_moead_t *moead, double value);

double moa_moead_mut_prob(struct moa_moead_t *moead, double value);

// moead methods
Moa *moa_moead(Mop *mop, WeightFileOption foption);

mbool moa_moead_first_run(Mop *mop);

mbool moa_moead_run(Mop *mop);

#endif //MOGEN_MGF_MOEAD_H
