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


#ifndef MOGEN_CROSSOVER_H
#define MOGEN_CROSSOVER_H

#include "mogen_population.h"

struct mop_t;
struct moeaz_indv_t;

enum moa_cx_types {
    CX_ONE_POINT,
    CX_TWO_POINT,
    CX_UNIFORM,
    CX_SBX,
    CX_PNX
};

void crossover(struct mop_t *mop, unsigned int p1, unsigned int p2, unsigned int c1, unsigned int c2);

void PNX(struct mop_t *mop, struct moeaz_indv_t* p1, struct moeaz_indv_t* p2,
    struct moeaz_indv_t* c1, struct moeaz_indv_t* c2);


#endif //MOGEN_CROSSOVER_H
