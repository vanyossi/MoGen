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


#ifndef MOGEN_MGF_NSGA2_H
#define MOGEN_MGF_NSGA2_H

#include "global_types.h"

#include "mgf_moa.h"
#include "operators/mgf_rank.h"

struct indv_t_nsga2_type {
    int rank;
    double crow_dist;
};

struct indv_type_t* mgf_indvtype_nsga2(Moa *moa);

struct indv_t_nsga2_type* mgf_indv_nsga2_data(struct indv_t *indv);

int indv_nsga2_rank(struct indv_t *indv);
void indv_nsga2_set_rank(struct indv_t *indv, int rank);
double indv_nsga2_crowdist(struct indv_t *indv);
void indv_nsga2_set_crowdist(struct indv_t *indv, double crowdist);


struct moa_nsga2_t {
    double cross_eta;
    double cross_prob;
    double mut_eta;
    double mut_prob;
    MoeazPop *mix;
    MoeazPop *Q;
    struct mgf_pop_ranks_t *pop_ranks;
};


MoaType* mgf_moatype_nsga2();

struct moa_nsga2_t* mgf_moa_nsga2_data(Moa* moa);


// set to -1 to ignore
void moa_nsga2_setparams(struct moa_nsga2_t *nsga2, double cx_eta, double cxp, double m_eta, double m_p);

struct moa_cm_container ngsa2_get_crossvals(Moa* self);

struct moa_cm_container ngsa2_get_mutvals(Moa* self);

double moa_nsga2_cross_eta(struct moa_nsga2_t *nsga2, double value);

double moa_nsga2_cross_prob(struct moa_nsga2_t *nsga2, double value);

double moa_nsga2_mut_eta(struct moa_nsga2_t *nsga2, double value);

double moa_nsga2_mut_prob(struct moa_nsga2_t *nsga2, double value);

// nsga2 methods
Moa *moa_nsga2(Mop *mop);
//
//void NSGA2_report_headers(MOEA_REPORT *report);

int mgf_nsga2_dominance(Individual *a, Individual *b);

unsigned int mgf_nsga2_binary_tournament(MoeazPop *P, unsigned int p1, unsigned int p2);

void mgf_nsga2_parent_selection(MoeazPop *P, unsigned int *parent, unsigned int n_parent);

int mgf_nsga2_offspring(Mop *mop, struct moa_nsga2_t *moa);

void mgf_nsga2_offspring_generation(Mop *mop, struct moa_nsga2_t *moa);

void mgf_nsga2_reduce(MoeazPop *P, Moa *moa);

// run step
mbool moa_nsga2_first_run(Mop *mop);
mbool moa_nsga2_run(Mop *mop);

#endif //MOGEN_MGF_NSGA2_H
