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


#ifndef MOGEN_MOGEN_MOA_H
#define MOGEN_MOGEN_MOA_H

#include "global_types.h"
#include "mogen_population.h"
#include "mop_report.h"

#define mbool int
#define mfalse 0
#define mtrue 1

struct mop_t;
struct mogen_moa_t;
struct moeaz_indv_t;

enum moa_stop_criterion_t {
    MGN_STOPIF_GEN,
    MGN_STOPIF_EXEC,
    MGN_STOP_SIZE_MARKER            // not used, only for extending enum
};

enum mogen_moa_types_e {
    MOA_MONO,
    MOA_DECOMP,
    MOA_INDIC,
    MOA_REND,
    MOA_PARETO
};

struct mogen_moa_bias_t {
    double cxprob;
    double cxidx;
};

struct mogen_moa_t {
    MoaTypes type;
    char name[64];
    struct mop_t* mop;
    struct mogen_moa_bias_t bias;
    MopReportStats stops;
    indv_init_extra extra_indv_alloc;
    void* cross;

    mbool (*stop)(struct mogen_moa_t* moa, MoaStopCriterion criterion);
    mbool (*run)(struct mop_t* mop);
};


Moa *moa_init(struct mop_t *mop, char *name, MoaTypes type, size_t mem_size);

mbool moa_run(struct mop_t *mop);

void moa_stopat_gen(Moa *moa, unsigned int gen);

void moa_stopat_eval(Moa *moa, unsigned int eval);

mbool moa_stop(struct mogen_moa_t* moa, MoaStopCriterion criterion);

//void moa_crossover(Moa* moa, enum moa_cx_types cx_func);

#endif //MOGEN_MOGEN_MOA_H
