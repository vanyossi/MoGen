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

#include "mgf_global_types.h"

typedef enum moa_cx_types {
    CX_ONE_POINT,
    CX_TWO_POINT,
    CX_UNIFORM,
    CX_SBX,
    CX_PNX
} CXType;

void moa_cross_setup(Moa *moa, CXType cx_type);

/**
 * @brief PNX CrossOver
 * @param mop
 * @param parent1
 * @param parent2
 * @param c1
 * @param c2
 */
void PNX_real(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2);
void PNX_bin(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2);
void PNX_mixed(Mop* mop, Individual *parent1, Individual* parent2, Individual* c1, Individual* c2) ;

#endif //MOGEN_CROSSOVER_H
