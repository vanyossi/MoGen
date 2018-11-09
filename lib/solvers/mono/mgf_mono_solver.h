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


#ifndef MOGEN_MGF_MONO_SOLVER_H
#define MOGEN_MGF_MONO_SOLVER_H

#include "mgf_moa.h"
#include "mogen_mop.h"
#include <math.h>

typedef double(*mono_fx)(double);

typedef enum moa_stop_criterion_mono_t {
    MGN_STOPIF_EPSILON = MGN_STOP_SIZE_MARKER,
} MoaStopCriterionMono;

typedef struct mop_mono_f {
    Mop mop;
    mono_fx fx;
} MopMono;

typedef struct moa_mono_type {
    double *error;
    double epsilon;
} MoaMono;

struct indv_t_mono_type {
    double error;
};

struct moa_mono_type* mgf_moa_get_mono_buffer(struct moa_t *moa);

struct moa_type_t* mgf_moatype_mono();

struct indv_t_mono_type* mgf_indv_get_mono_buffer(struct indv_t *indv);

struct indv_type_t* mgf_indvtype_mono(Moa *moa);

void mop_mono_assign_fx(Mop *mop, mono_fx f);

void mgf_moa_mono_set_solver(Moa* moa, void (*evaluate)(Mop*, Individual*));

void mgf_moa_mono_set_epsilon(Moa *moa, double epsilon);

mbool moa_mono_run(Mop *mop);

//mbool moa_mono_stop(Moa *moa, MoaStopCriterion criterion);

#endif //MOGEN_MGF_MONO_SOLVER_H
