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

/**
 * @file mogen_mop.h
 * @brief
*/

#ifndef MOGEN_MOP_H
#define MOGEN_MOP_H

#include <stdbool.h>

#include "global_types.h"

#include "mgf_population.h"
#include "mogen_moa.h"
#include "mop_report.h"

enum mop_specs_e {
    MOP_REAL =          1 << 0,     // 0b000000001
    MOP_BIN =           1 << 1,     // 0b000000010
    MOP_MIX =           11 << 0,    // 0b000000011
    MOP_CONTIGUOUS =    1 << 6,     // 0b001000000
    MOP_RESTRICTED =    1 << 7,     // 0b010000000
    MOP_DYNAMIC =       1 << 8      // 0b100000000
};

struct mop_t;
struct mogen_moa_t;

struct mop_base_t {
    char name[128];
    unsigned int type:8;    //!< real, bin or mixed
    unsigned int continuous:8;
    unsigned int restricted:8;
    unsigned int dynamic:8;

    unsigned int ndec;     //!< Number of real decision variables
    unsigned int nobj;    //!< Number of objectives
    unsigned int ncons;   //!< Constrain number
};


struct mop_limit_t {
    double *xmin;
    double *xmax;
};


struct mop_extra_t {
    int size;
};

struct mop_t {
    Mop_base set;
    Mop_limit limits;
    MoeazPop *pop;
    MoeazPop *ext_pop;
    struct mogen_moa_t* solver;
    void (*evaluate)(struct mop_t *mop, Individual *indv);
    MopReport report;
};

// void (*cross)(struct mop_t *mop, Individual* p1, Individual* p2, Individual* c1, Individual* c2);

Mop *mogen_mop(char *name, MopSpecs mop_specs, size_t mem_size);

void mop_set_params(Mop *mop, unsigned int nreal, unsigned int nobjs, unsigned int ncons);

void mop_set_limits_ndec(Mop *mop, double *min, double *max, unsigned int size);

Individual *mogen_mop_getindv(Mop *mop, unsigned int pos);

mbool mop_evaluate(Mop *mop, Individual *indv);

void mop_solve(Mop *mop, int steps);

void mop_print(Mop *mop);

#endif
