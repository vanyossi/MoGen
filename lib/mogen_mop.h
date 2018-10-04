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

#include "mogen_population.h"
#include "mogen_moa.h"

typedef enum mop_specs_e {
    MOP_REAL =          1 << 0,     // 0b00001
    MOP_BIN =           1 << 1,     // 0b00010
    MOP_MIX =           11 << 0,    // 0b00011
    MOP_CONTIGUOUS =    1 << 2,     // 0b00100
    MOP_RESTRICTED =    1 << 3,     // 0b01000
    MOP_DYNAMIC =       1 << 4      // 0b10000
} MopSpecs;

struct mop_t;
struct mogen_moa_t;

typedef struct mop_base_t {
    char name[128];
    unsigned int type:8;    //!< real, bin or mixed
    unsigned int continuous:8;
    unsigned int restricted:8;
    unsigned int dynamic:8;

    unsigned int ndec;     //!< Number of real decision variables
    unsigned int nobj;    //!< Number of objectives
    unsigned int ncons;   //!< Constrain number

    double prob_cross;
    void (*x_func)(struct mop_t mop, int p1, int p2, int i1, int i2);

} Mop_base;


typedef struct mop_limit_t {
    double *xmin;
    double *xmax;
} Mop_limit;


typedef struct mop_t {
    Mop_base set;
    Mop_limit limits;
    MoeazPop* pop;
    struct mogen_moa_t* solver;
    void (*evaluate)(struct mop_t* mop, unsigned int idx);
    void (*cross)(struct mop_t *mop, MoeazIndv* p1, MoeazIndv* p2, MoeazIndv* c1, MoeazIndv* c2);
} Mop;


Mop *mogen_mop(char *name, MopSpecs mop_specs);

MoeazIndv *mogen_mop_getindv(Mop *mop, unsigned int pos);

void mop_print(Mop *mop);

#endif
