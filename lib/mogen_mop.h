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

#include "mgf_global_types.h"

#include "mgf_population.h"
#include "mgf_moa.h"
#include "mop_report.h"

struct mop_t;
struct moa_t;

struct mop_base_t {
    char name[128];
    unsigned int type:8;    //!< real, bin or mixed
    unsigned int continuous:8;
    unsigned int restricted:8;
    unsigned int dynamic:8;

    unsigned int xsize;     //!< Number of real decision variables
    unsigned int isize;
    unsigned int bsize;
    unsigned int nobj;    //!< Number of objectives
    unsigned int ncons;   //!< Constrain number
};

struct mop_limit_t {
    double *xmin;
    double *xmax;
    int *imin;
    int *imax;
};


struct mop_extra_t {
    int size;
};

struct mop_t {
    Mop_base set;
    Mop_limit limits;
    MoeazPop *pop;
    MoeazPop *ext_pop;
    Moa *solver;
    void (*evaluate)(struct mop_t *mop, Individual *indv);
    MopReport report;
};

// void (*cross)(struct mop_t *mop, Individual* p1, Individual* p2, Individual* c1, Individual* c2);

Mop *mogen_mop(char *name, MopSpecs mop_specs, size_t mem_size);

void mop_set_params(
    Mop *mop, unsigned int real, unsigned int binary, unsigned int integer, unsigned int nobjs, unsigned int ncons);

void mop_set_limits_ndec(
    Mop *mop, double *xmin, double *xmax, unsigned int lim_xsize, int *imin, int *imax, unsigned int lim_isize);

Individual *mogen_mop_getindv(Mop *mop, unsigned int pos);

mbool mop_evaluate(Mop *mop, Individual *indv);

void mop_solve(Mop *mop, int steps);

void mop_print(Mop *mop);

#endif
