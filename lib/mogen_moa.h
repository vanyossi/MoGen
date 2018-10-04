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

#include "mogen_moea.h"
#include "crossover.h"

struct mop_t;

typedef union mogen_moa_type_t {
    Moea moea;
} MoaType;


typedef struct mogen_moa_t {
    char name[64];
    struct mop_t* mop;
    void (*run)(struct mop_t* mop, int steps);          //!< Moa evaluate step function pointer
    MoaType algorithm;
} Moa;

Moa* moa_init(struct mop_t *mop, char* name);

void moa_crossover(Moa* moa, enum moa_cx_types cx_func);

#endif //MOGEN_MOGEN_MOA_H
