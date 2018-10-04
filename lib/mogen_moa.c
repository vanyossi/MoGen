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

#include "mogen_moa.h"
#include "mogen_mop.h"

#include <stdlib.h>
#include <string.h>

Moa* moa_init(struct mop_t *mop, char* name){
    Moa* new_moa = calloc(1, sizeof(Moa));

    strcpy(new_moa->name, name);
    new_moa->mop = mop;
    mop->solver = new_moa;

    return new_moa;
}

void moa_crossover(Moa* moa, enum moa_cx_types cx_func){
    moa->mop->cross = PNX;
}
